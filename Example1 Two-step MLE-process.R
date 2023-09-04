library(Sim.DiffProc)
library(ggplot2)
library(deSolve)
library(tictoc)

set.seed(1234) #used to create reproducible results
f <- 3
a <- 0.25
sigma <- 100
b <- 100
x0 <- 0 #initial value of X
y0 <- 100 #initial value of Y
theta <- 1 #real value of theta
realizations <- 100 #we will simulate 100 processes each

eps_vec <- c(0.2)
#eps_vec <- c(0.9, 0.5, 0.3, 0.25, 0.2, 0.15, 0.125, 0.1, 0.05)
#eps_vec <- seq(from = 0.995, to = 0.005, by = -0.005) 

tic("total")
#----------
# Simulate processes X_t and Y_t and calculate its median values
#----------
print("Simulate processes X_t and Y_t")
tic("simulation of processes X_t and Y_t")
x_for_epsilon_list <- c()
y_for_epsilon_list <- c()
for(epsilon in eps_vec){
  sim_x_y_drift <- expression(theta*f*y,a*y)
  sim_x_y_diffusion <- expression(epsilon*sigma, epsilon*b)
  #The following function simulated a 2-dimensional SDE of Itô-type; more details can be found in the bibliography.
  sim_x_y <- snssde2d(drift=sim_x_y_drift, diffusion=sim_x_y_diffusion, M=realizations, x0=c(x0,y0), t0=0, T=3, Dt=0.003, method="euler")
  x_for_epsilon_list[[as.character(round(epsilon,3))]] <- sim_x_y$X 
  y_for_epsilon_list[[as.character(round(epsilon,3))]] <- sim_x_y$Y
}

#We need the values of X_t and Y_t at the corresponding times 0, 0.003, ..., 3.
get_median_values_of_process <- function(t, epsilon, process, realizations, time_steps){
  vec_at_t <- c()
  index <- (t*time_steps)+1 #we need +1 since we start at time 0, which has index 1
  for(i in 0:(realizations-1)){
    vec_at_t <- c(vec_at_t, process[[as.character(round(epsilon,3))]][[i*(time_steps+1) + index]])
  }
  return(median(vec_at_t))
}

Xt_values <- list()
Yt_values <- list()
for(epsilon in eps_vec){
  for(i in 0:1000){
    Xt_values[[as.character(round(epsilon,3))]][[as.character(round(i/1000*3,3))]] <- get_median_values_of_process(round(i/1000,3), epsilon, x_for_epsilon_list, realizations, time_steps=1000)
    Yt_values[[as.character(round(epsilon,3))]][[as.character(round(i/1000*3,3))]] <- get_median_values_of_process(round(i/1000,3), epsilon, y_for_epsilon_list, realizations, time_steps=1000)
  }
  print(epsilon)
}

#Let us delete results that we no longer need in the further course to save memory
rm(x_for_epsilon_list)
rm(y_for_epsilon_list)
toc()
#----------
# Calculate preliminary estimators
#----------
tic("calculation of preliminary estimators")
prelim_est_list <- c()
for(epsilon in round(eps_vec,3)){
  tau_eps <- epsilon^(1.25)
  #We have to check whether round(tau_eps,3) is divisible by 3;
  #if not, we have to round to the next number which is divisible by 3.
  calculated_median <- Xt_values[[as.character(epsilon)]][[as.character(ceiling(round(tau_eps,3)*1000/3)*3*0.001)]]
  #The following calculation is based on the results of the second example.
  prelim_est_list[[as.character(round(epsilon,3))]] <- calculated_median/(f*y0*((exp(a*tau_eps)/a)-(1/a)))
}
toc()

#----------
# For the One-step MLE-process:
#----------

#----------
# Solve Riccati equation for different preliminary estimators
# to solve: 
# \partial \gamma(\vartheta,t)/\partial t = 2*a*\gamma(\vartheta,t) 
#       - (\gamma(\vartheta,t)^2*\vartheta^2*f^2)/(\sigma^2)+b^2
# initial value: \gamma(\vartheta,0)=0
# Note that we have to use the prelim est at tau_eps instead of epsilon!
#----------
print("Calculate the solutions of the Riccati equations")
tic("calculation of solutions of Riccati equation")
#The following function is represents the ODE we have to solve:
riccati_eq <- function(t,y, parms){
  list(2*a*y-((y^2*v^2*f^2)/(sigma^2))+b^2)
}

#We will save the results in the following list
sol_riccati_eq_for_diff_prelim_est_list <- list()

times <- seq(from = 0, to = 3, by = 0.003)

for(epsilon in eps_vec){
  epsilon_rounded <- round(epsilon,3)
  #We have to insert the preliminary estimator at time tau_eps
  v <- prelim_est_list[[as.character(epsilon_rounded)]]
  #The following vector includes all parameters needed in the riccati_eq function
  parameter <- c(a=a, f=f, sigma=sigma, b=b, v = v)
  sol_riccati_eq <- ode(y = 0, times = times, func = riccati_eq, parms = parameter)
  sol_riccati_eq_for_diff_prelim_est_list[[as.character(epsilon_rounded)]] <- sol_riccati_eq
}
toc()
#We need the values of the solution of the Riccati equation for corresponding t 
#and preliminary estimator for the calculation of the conditional expectation
get_values_of_sol_of_riccati_eq_at_t_and_epsilon <- function(t, epsilon){
  index <- (t/3)*1000 + 1002
  value <- sol_riccati_eq_for_diff_prelim_est_list[[as.character(round(epsilon,3))]][[index]]
  return(value)
}

#----------
# Calculate conditional expectation
# For each preliminary estimator, we have to solve the equation just once!
#----------
print("Calculate value of conditional expectation")
tic("calculation of conditional expectation")
#The following function will simulate the one-dimensional SDE m(.,t)
calculate_cond_exp_at_eps <- function(f, a, sigma, b, y0, realizations, epsilon){
  epsilon_rounded <- round(epsilon,3)
  sim_cond_exp_drift <- expression((a-(get_values_of_sol_of_riccati_eq_at_t_and_epsilon(t,epsilon)*prelim_est_list[[as.character(epsilon_rounded)]]*f^2)/sigma^2)*x+(get_values_of_sol_of_riccati_eq_at_t_and_epsilon(t,epsilon)*f^2*prelim_est_list[[as.character(epsilon_rounded)]])/(sigma^2)*(Yt_values[[as.character(epsilon_rounded)]][[as.character(t)]]))
  sim_cond_exp_diffusion <- expression(epsilon_rounded*get_values_of_sol_of_riccati_eq_at_t_and_epsilon(t,epsilon)*f/sigma)
  sim_cond_exp <- snssde1d(drift=sim_cond_exp_drift, diffusion=sim_cond_exp_diffusion, M=realizations, x0=y0, t0=0, T=3, Dt=0.003, method="euler")
  return(sim_cond_exp)
}

#We will save our results in the list "cond_exp_list". The other list exists for technical reasons only.
cond_exp_list <- list()
cond_exp_list_before_calculating_median <- list()
for(epsilon in eps_vec){
  epsilon_rounded <- round(epsilon,3)
  cond_exp_list_before_calculating_median[[as.character(epsilon_rounded)]] <- calculate_cond_exp_at_eps(f, a, sigma, b, y0, realizations, epsilon)$X
  for(i in 0:1000){
    cond_exp_list[[as.character(epsilon_rounded)]][[as.character(round(i/1000*3,3))]] <- get_median_values_of_process(round(i/1000,3), epsilon, cond_exp_list_before_calculating_median, realizations, 1000)
  }
  print(epsilon)
}
rm(cond_exp_list_before_calculating_median)
toc()
#----------
# Calculate dot(y)(theta,t) as we will need it several times
#----------
tic("calculation of dot(y)")
D_at_theta_and_t <- list()
for(epsilon in eps_vec){
  for(i in 0:1000){
    D_at_theta_and_t[[as.character(round(epsilon,3))]][[as.character(round(i/1000*3,3))]] <- get_values_of_sol_of_riccati_eq_at_t_and_epsilon(round(i/1000*3,3), round(epsilon,3))*prelim_est_list[[as.character(round(epsilon,3))]]*f/sigma^2
  }
}
dot_y_exp_integral <- list()
for(epsilon in eps_vec){
  epsilon_rounded <- round(epsilon,3)
  for(i in 0:1000){
    dot_y_exp_integral[[as.character(epsilon_rounded)]][[as.character(round(i/1000*3,3))]][[as.character(round(i/1000*3,3))]] <- 0.001*(a-D_at_theta_and_t[[as.character(epsilon_rounded)]][[as.character(round(i/1000*3,3))]]*prelim_est_list[[as.character(epsilon_rounded)]]*f)
  }
  for(i in 0:999){
    accumulated_value <- dot_y_exp_integral[[as.character(epsilon_rounded)]][[as.character(round(i/1000*3, 3))]][[as.character(round(i/1000*3, 3))]]
    for(j in (i+1):1000){
      if(j<=1000){
        accumulated_value <- accumulated_value + 0.001 * 0.001*(a-D_at_theta_and_t[[as.character(epsilon_rounded)]][[as.character(round(j/1000*3,3))]]*prelim_est_list[[as.character(epsilon_rounded)]]*f)
        dot_y_exp_integral[[as.character(epsilon_rounded)]][[as.character(round(i/1000*3,3))]][[as.character(round(j/1000*3,3))]] <- accumulated_value
      }
    }
  }
}

dot_y <- list()
for(epsilon in eps_vec){
  epsilon_rounded <- round(epsilon,3)
  for(upper_bound in 1:1000){
    value_integral <- 0
    for(i in 0:(upper_bound-1)){
      value_integral <- value_integral + 0.001*(-exp(dot_y_exp_integral[[as.character(epsilon_rounded)]][[as.character(round(i/1000*3,3))]][[as.character(round(upper_bound/1000*3,3))]])*D_at_theta_and_t[[as.character(round(epsilon,3))]][[as.character(round(i/1000*3,3))]]*f*y0*exp(a*round(i/1000*3)))
    }
    dot_y[[as.character(epsilon_rounded)]][[as.character(round(upper_bound/1000*3,3))]] <- value_integral
  }
}
toc()
#----------
# Calculate dot(M)
#----------
tic("calculation of dot(M)")
dot_M <- list()
for(epsilon in eps_vec){
  epsilon_rounded <- round(epsilon,3)
  for(i in 0:1000){
    i_rounded <- round(i/1000*3,3)
    dot_M[[as.character(epsilon_rounded)]][[as.character(i_rounded)]] <- f*y0*exp(a*i_rounded)+prelim_est_list[[as.character(epsilon_rounded)]]*f*dot_y[[as.character(epsilon_rounded)]][[as.character(i_rounded)]]
  }
}
toc()
#----------
# Calculate Fisher information
#----------
tic("calculation of Fisher information")
fisher_info <- list()
for(epsilon in eps_vec){
  epsilon_rounded <- round(epsilon,3)
  tau_eps_rounded <- round(epsilon^(1.25),3)
  for(upper_bound in (ceiling((tau_eps_rounded*1000/3)/3)*3):1000){
    value_integral <- 0
    for(i in 1:(upper_bound-round((epsilon^(1.25))/3*1000))){
      value_integral <- value_integral + 0.001*(((f*y0*exp(a*round(i/1000*3,3))+f*dot_y[[as.character(epsilon_rounded)]][[as.character(round(i/1000*3,3))]])/sigma)^2)
    }
    fisher_info[[as.character(epsilon_rounded)]][[as.character(round(upper_bound/1000*3,3))]] <- value_integral
  }
}
toc()
#----------
# Calculate stochastic integral
#----------
print("Calculate values of stochastic integral")
tic("calculation of stochastic integral")
stochastic_integral <- list()
for(epsilon in eps_vec){
  tau_eps <- epsilon^(1.25)
  tau_eps_rounded <- round(tau_eps,3)
  epsilon_rounded <- round(epsilon,3)
  value_integral <- 0
  if(round(tau_eps,3)*1000/3+2<1000){
    number_of_simulation_steps <- 1000-ceiling(round(epsilon^(1.25),3)*1000/3)
    #number_of_simulation_steps <- length(unlist(dot_M[[as.character(round(epsilon,3))]], use.names= FALSE))-1
    ito_integrand <- expression(dot_M[[as.character(round(epsilon,3))]][[as.character(t)]]*epsilon/sigma)
    sim_of_stochastic_integral <- Sim.DiffProc::st.int(expr=ito_integrand,type="ito",M=realizations,lower=(ceiling(round(epsilon^(1.25),3)*1000/3)*3/1000),upper=3, subdivisions=number_of_simulation_steps)
    #We now, once again, have to calculate the median of all realizations. 
    #This is more difficult here as we also have to consider the varying lower bound depending on tau_eps
    for(i in 1:(number_of_simulation_steps+1)){
      stochastic_integral_at_specific_sim_step <- c()
      for(j in 0:99){
        stochastic_integral_at_specific_sim_step <- c(stochastic_integral_at_specific_sim_step, sim_of_stochastic_integral$X[[i+(number_of_simulation_steps+1)*j]])
      }
      stochastic_integral[[as.character(epsilon_rounded)]][[as.character(ceiling(tau_eps_rounded*1000/3)*3/1000+(i-1)*0.003)]] <- median(stochastic_integral_at_specific_sim_step)
    }
  }
}
toc()
#----------
# Calculate One-step MLE-process
#----------
print("Calculate One-step MLE-process")
tic("calculation of One-step MLE-process")
ds_integral <- list()
one_step_mle_process <- list()
for(epsilon in eps_vec){
  epsilon_rounded <- round(epsilon,3)
  tau_epsilon_rounded <- round(epsilon^(1.25),3)
  #Let us first calculate the value of the integral involving Y_s and the conditional expectation
  for(upper_bound in (ceiling((round((epsilon^(1.25)),3)*1000/3))+2):1000){
    upper_div_by_1000_times_3_rounded <- round(upper_bound/1000*3,3)
    value_integral <- 0
    for(i in (tau_epsilon_rounded*1000+2):upper_bound){
      value_integral <- value_integral + 0.001*((dot_M[[as.character(epsilon)]][[as.character(upper_div_by_1000_times_3_rounded)]])/sigma^2)*f*prelim_est_list[[as.character(epsilon_rounded)]]*(Yt_values[[as.character(epsilon_rounded)]][[as.character(upper_div_by_1000_times_3_rounded)]]-cond_exp_list[[as.character(epsilon_rounded)]][[as.character(upper_div_by_1000_times_3_rounded)]])
    }
    ds_integral[[as.character(epsilon_rounded)]][[as.character(upper_div_by_1000_times_3_rounded)]] <- value_integral
  }
  for(upper_bound in ceiling((round((epsilon^(1.25)),3)*1000/3+2)):1000){
    upper_div_by_1000_times_3_rounded <- round(upper_bound/1000*3,3)
    one_step_mle_process[[as.character(epsilon,3)]][[as.character(upper_div_by_1000_times_3_rounded)]] <- prelim_est_list[[as.character(epsilon_rounded)]] + 1/fisher_info[[as.character(epsilon_rounded)]][[as.character(upper_div_by_1000_times_3_rounded)]]*ds_integral[[as.character(epsilon_rounded)]][[as.character(upper_div_by_1000_times_3_rounded)]]+1/fisher_info[[as.character(epsilon_rounded)]][[as.character(upper_div_by_1000_times_3_rounded)]]*stochastic_integral[[as.character(epsilon_rounded)]][[as.character(upper_div_by_1000_times_3_rounded)]]
  }
}
toc()

#----------
# Calculate Two-step MLE-process
#----------

#----------
# For the first preliminary estimator process:
#----------

#----------
# Solve Riccati equation for different preliminary estimators
# The preliminary estimators now are the values of the One-step MLE-process!
# to solve: 
# \partial \gamma(\vartheta,t)/\partial t = 2*a*\gamma(\vartheta,t) 
#       - (\gamma(\vartheta,t)^2*\vartheta^2*f^2)/(\sigma^2)+b^2
# initial value: \gamma(\vartheta,0)=0
# Note that we have to use the prelim est at tau_eps instead of epsilon!
#----------
print("Calculate the solutions of the Riccati equations")
tic("Two-step MLE: calculation of solutions of Riccati equations")
#The following function is represents the ODE we have to solve:
riccati_eq <- function(t,y, parms){
  list(2*a*y-((y^2*v^2*f^2)/(sigma^2))+b^2)
}

#We will save the results in the following list
sol_riccati_eq_for_diff_prelim_est_list_two_step <- list()

times <- seq(from = 0, to = 3, by = 0.003)

for(epsilon in eps_vec){
  epsilon_rounded <- round(epsilon,3)
  time_steps <- length(one_step_mle_process[[as.character(epsilon_rounded)]])
  for(time in (1000-time_steps+1):1000){
    v <- one_step_mle_process[[as.character(epsilon_rounded)]][[as.character(round(time/1000*3,3))]]
    parameter <- c(a=a, f=f, sigma=sigma, b=b, v = v)
    sol_riccati_eq <- ode(y = 0, times = times, func = riccati_eq, parms = parameter)
    sol_riccati_eq_for_diff_prelim_est_list_two_step[[as.character(epsilon_rounded)]][[as.character(round(time/1000*3,3))]] <- sol_riccati_eq
  }
}
toc()
get_values_of_sol_of_riccati_eq_at_t_and_epsilon_two_step <- function(t, epsilon, time_point){
  index <- (t/3)*1000 + 1002
  value <- sol_riccati_eq_for_diff_prelim_est_list_two_step[[as.character(round(epsilon,3))]][[as.character(round(time_point/1000*3,3))]][[index]]
  return(value)
}

#----------
# Calculate conditional expectation
# For each preliminary estimator, we have to solve the equation just once!
# We need a different function to calculate the median of the realizations as we will not have 1000 steps anymore
#----------
tic("Two-step MLE: calculation of conditional expectation")
get_median_values_of_process_two_step <- function(t, epsilon, process, realizations, time_steps, time_point){
  vec_at_t <- c()
  index <- ceiling(t*time_steps)+1 #we need +1 since we start at time 0, which has index 1
  for(i in 0:(realizations-1)){
    vec_at_t <- c(vec_at_t, process[[as.character(round(epsilon,3))]][[as.character(round(time_point,3))]][[i*(time_steps) + index]])
  }
  return(median(vec_at_t))
}

print("Calculate value of conditional expectation")
calculate_cond_exp_at_eps_two_step <- function(f, a, sigma, b, y0, realizations, epsilon, time_point){
  epsilon_rounded <- round(epsilon,3)
  sim_cond_exp_drift <- expression((a-(get_values_of_sol_of_riccati_eq_at_t_and_epsilon_two_step(t,epsilon,time_point)*one_step_mle_process[[as.character(epsilon_rounded)]][[as.character(t)]]*f^2)/sigma^2)*x+(get_values_of_sol_of_riccati_eq_at_t_and_epsilon_two_step(t,epsilon,time_point)*f^2*one_step_mle_process[[as.character(epsilon_rounded)]][[as.character(t)]])/(sigma^2)*(Yt_values[[as.character(epsilon_rounded)]][[as.character(t)]]))
  sim_cond_exp_diffusion <- expression(epsilon_rounded*get_values_of_sol_of_riccati_eq_at_t_and_epsilon_two_step(t,epsilon,time_point)*f/sigma)
  sim_cond_exp <- snssde1d(drift=sim_cond_exp_drift, diffusion=sim_cond_exp_diffusion, N=953, M=realizations, x0=y0, t0=0.141, T=3, Dt=0.003, method="euler") #we have to start at 0.141 as we do not have values for the One-step MLE-process before that time
  return(sim_cond_exp)
}

#We will save our results in the list "cond_exp_list". The other list exists for technical reasons only.
cond_exp_list_two_step <- list()
cond_exp_list_before_calculating_median_two_step <- list()
for(epsilon in eps_vec){
  epsilon_rounded <- round(epsilon,3)
  time_steps <- length(one_step_mle_process[[as.character(epsilon_rounded)]])
  for(time_point in (1000-time_steps+1):1000){
    cond_exp_list_before_calculating_median_two_step[[as.character(epsilon_rounded)]][[as.character(round(time_point,3))]] <- calculate_cond_exp_at_eps_two_step(f, a, sigma, b, y0, realizations, epsilon, time_point)$X
    for(i in 0:998){
      cond_exp_list_two_step[[as.character(epsilon_rounded)]][[as.character(round(time_point/1000*3,3))]][[as.character(round(i/1000*3,3))]] <- get_median_values_of_process_two_step(round(i/1000,3), epsilon, cond_exp_list_before_calculating_median_two_step, realizations, 954, time_point)
    }
  }
}
rm(cond_exp_list_before_calculating_median_two_step)
toc()
#----------
# Calculate dot(y)(theta,t) as we will need it several times
#----------
tic("calculation of dot(y)")
D_at_theta_and_t_two_step <- list()
for(epsilon in eps_vec){
  epsilon_rounded <- round(epsilon,3)
  time_steps <- length(one_step_mle_process[[as.character(epsilon_rounded)]])
  for(time_point in ((1000-time_steps)+1):1000){
    time_point_rounded <- round(time_point/1000*3,3)
    for(i in 0:time_point){
      i_rounded <- round(i/1000*3,3)
      D_at_theta_and_t_two_step[[as.character(epsilon_rounded)]][[as.character(time_point_rounded)]][[as.character(i_rounded)]] <- get_values_of_sol_of_riccati_eq_at_t_and_epsilon(i_rounded, epsilon_rounded)*one_step_mle_process[[as.character(epsilon_rounded)]][[as.character(time_point_rounded)]]*f/sigma^2
    }
  }
}

dot_y_exp_integral_two_step <- list()
for(epsilon in eps_vec){
  epsilon_rounded <- round(epsilon,3)
  time_steps <- length(one_step_mle_process[[as.character(epsilon_rounded)]])
  for(time_point in ((1000-time_steps)+1):1000){
    time_point_rounded <- round(time_point/1000*3,3)
    for(lower_bound in 0:time_point){
      value_integral <-0
      for(i in lower_bound:time_point){
        value_integral <- value_integral + 0.001*(a-D_at_theta_and_t_two_step[[as.character(epsilon_rounded)]][[as.character(time_point_rounded)]][[as.character(round(i/1000*3,3))]]*one_step_mle_process[[as.character(epsilon_rounded)]][[as.character(time_point_rounded)]]*f)
      }
      dot_y_exp_integral_two_step[[as.character(epsilon_rounded)]][[as.character(time_point_rounded)]][[as.character(round(lower_bound/1000*3,3))]] <- value_integral
    }
  }
}

dot_y_two_step <- list()
for(epsilon in eps_vec){
  epsilon_rounded <- round(epsilon,3)
  time_steps <- length(one_step_mle_process[[as.character(epsilon_rounded)]])
  for(time_point in ((1000-time_steps)+1):1000){
    time_step_rounded <- round(time_point/1000*3,3)
    value_integral <- 0
    for(i in 0:time_point){
      i_rounded <- round(i/1000*3,3)
      value_integral <- value_integral - 0.001*exp(dot_y_exp_integral_two_step[[as.character(epsilon_rounded)]][[as.character(time_step_rounded)]][[as.character(i_rounded)]])*D_at_theta_and_t_two_step[[as.character(epsilon_rounded)]][[as.character(time_point_rounded)]][[as.character(i_rounded)]]*f*y0*exp(a*i_rounded)
      dot_y_two_step[[as.character(epsilon_rounded)]][[as.character(time_step_rounded)]][[as.character(i_rounded)]] <- value_integral
    }
  }
}
toc()
#----------
# Calculate Fisher information
#----------
tic("Two-step MLE: calculation of Fisher information")
fisher_info_two_step <- list()
for(epsilon in eps_vec){
  epsilon_rounded <- round(epsilon,3)
  tau_eps_rounded <- round(epsilon^(1.25),3)
  for(upper_bound in (ceiling((tau_eps_rounded*1000/3)/3)*3+2):1000){
    upper_bound_rounded <- round(upper_bound/1000*3,3)
    value_integral <- 0
    for(i in (ceiling((tau_eps_rounded*1000/3)/3)*3+2):upper_bound){
      i_rounded <- round(i/1000*3,3)
      value_integral <- value_integral + 0.001*(((f*y0*exp(a*i_rounded)+one_step_mle_process[[as.character(epsilon_rounded)]][[as.character(i_rounded)]]*f*dot_y_two_step[[as.character(epsilon_rounded)]][[as.character(upper_bound_rounded)]][[as.character(as.character(i_rounded))]])/sigma)^2)
    }
    fisher_info_two_step[[as.character(epsilon_rounded)]][[as.character(upper_bound_rounded)]] <- value_integral
  }
}
toc()
#----------
# Calculate ds-integral
#----------
tic("Two-step MLE: calculation of ds-integral")
print("Calculate ds-integral")
ds_integral_two_step <- list()
for(epsilon in eps_vec){
  epsilon_rounded <- round(epsilon,3)
  tau_epsilon_rounded <- round(epsilon^(1.25),3)
  time_steps <- length(one_step_mle_process[[as.character(epsilon_rounded)]])
  for(time_point in ((1000-time_steps)+2):1000){
    time_rounded <- round(time_point/1000*3,3)
    value_integral <- 0
    for(i in (ceiling((round((epsilon^(1.25)),3)*1000/3))+2):time_point){
      i_rounded <- round(i/1000*3,3)
      value_integral <- value_integral + 0.001*(dot_M[[as.character(epsilon_rounded)]][[as.character(i_rounded)]]/sigma^2*f*one_step_mle_process[[as.character(epsilon_rounded)]][[as.character(i_rounded)]]*(Yt_values[[as.character(epsilon_rounded)]][[as.character(i_rounded)]]-cond_exp_list_two_step[[as.character(epsilon_rounded)]][[as.character(time_point_rounded)]][[as.character(i_rounded)]]))
    }
    ds_integral_two_step[[as.character(epsilon_rounded)]][[as.character(time_rounded)]] <- value_integral
  }
}
#----------
# Calculate Two-step MLE-process
#----------
#Note that we have to plug in dot(M)(preliminary estimator), which we have already calculated in order to calculate the (first) preliminary estimator. Hence, we can reuse these results now.
two_step_mle_process <- list()
for(epsilon in eps_vec){
  epsilon_rounded <- round(epsilon,3)
  tau_epsilon_rounded <- round(epsilon^(1.25),3)
  time_steps <- length(one_step_mle_process[[as.character(epsilon_rounded)]])
  for(at_t in ((1000-time_steps)+2):1000){ 
    at_t_rounded <- round(at_t/1000*3,3)
    two_step_mle_process[[as.character(epsilon_rounded)]][[as.character(at_t_rounded)]] <- one_step_mle_process[[as.character(epsilon_rounded)]][[as.character(at_t_rounded)]] + (1/fisher_info_two_step[[as.character(epsilon_rounded)]][[as.character(at_t_rounded)]])*ds_integral_two_step[[as.character(epsilon_rounded)]][[as.character(at_t_rounded)]]+(1/fisher_info_two_step[[as.character(epsilon_rounded)]][[as.character(at_t_rounded)]])*stochastic_integral[[as.character(epsilon_rounded)]][[as.character(at_t_rounded)]]
  }
}
toc()

#----------
# Create a plot regarding the consistency with respect to t for constant coefficients, vartheta_0=0, epsilon = 0.2
#----------
#The following definitions are technical necessities to layer the lines. 
#We define all unnecessary times to have values 0 so that these values do not appear in the plot.
values_one_step <- c()
values_two_step <- c()
values_prelim <- c()
df_consistency_wrt_t_list_two_step <- list()
values_real <- rep(1, 1000)

for(epsilon in eps_vec){
  one_step_temp <- c(rep(-100, round(epsilon^(1.25)/3+0.001,3)*1000),unlist(one_step_mle_process[[as.character(epsilon)]], use.names=FALSE))
  if(length(one_step_temp)<1000){
    one_step_temp <- c(one_step_temp, rep(-100, 1000-length(one_step_temp)))
  }
  values_one_step[[as.character(epsilon)]] <- one_step_temp
  two_step_temp <- c(rep(-100, round(epsilon^(1.25)/3+0.001,3)*1000),unlist(two_step_mle_process[[as.character(epsilon)]], use.names = FALSE))
  if(length(two_step_temp)<1000){
    two_step_temp <- c(two_step_temp, rep(-100, 1000-length(two_step_temp)))
  }
  values_two_step[[as.character(epsilon)]] <- two_step_temp
  values_prelim[[as.character(epsilon)]] <- c(rep(prelim_est_list[[as.character(epsilon)]], ceiling(round(epsilon^(1.25)+0.001,3)*1000/3)), rep(-100, 1000-ceiling(round(epsilon^(1.25)+0.001,3)*1000/3)))
  df_consistency_wrt_t_list_two_step[[as.character(epsilon)]] <- data.frame(times=seq(0.003,3,0.003), values_prelim = values_prelim[[as.character(epsilon)]], values_one_step = values_one_step[[as.character(epsilon)]], values_two_step = values_two_step[[as.character(epsilon)]], values_real = values_real)
}

plot_consistency_wrt_t_full <- ggplot(df_consistency_wrt_t_list_two_step[[as.character(0.2)]], aes(x=times)) +
  geom_line(aes(y=values_prelim,linetype="solid"))+ #It should be the other way around, but this works
  geom_line(aes(y=values_one_step,linetype="dotted"))+
  geom_line(aes(y=values_two_step,linetype="longdash"))+
  geom_line(aes(y=values_real,linetype="dashed"))+
  labs(x ="t", y = "One-step and Two-step MLE-processes") +
  ylim(0.5,1.15)+
  theme(panel.background = element_rect(fill = "white", colour = "black",
                                        size = 0.25, linetype = "solid"))+
  labs(linetype='Estimators')+
  scale_linetype_manual(values = c(1, 2, 3, 4),
                        labels = c("Real value", "One-step MLE-process", "Two-step MLE-process", "Preliminary Estimator"))+
  theme(legend.position = c(0.8, 0.2))+
  theme(legend.background = element_rect(colour ="black"))

#----------
# Create a plot of the conditional expectation 
#----------
cond_exp_two_step <- list()
for(i in 0:46){
  cond_exp_two_step[[as.character(0.2)]][[as.character(0+i*0.003)]]<- -100
}
for(i in 0:998){
  cond_exp_two_step[[as.character(0.2)]][[as.character(0.141+i*0.003)]] <- cond_exp_list_two_step[[as.character(0.2)]][[as.character(0.141+i*0.003)]][[as.character(0.141+i*0.003)]]
}
cond_exp_two_step[[as.character(0.2)]][[as.character(2.997)]] <- -100
cond_exp_two_step[[as.character(0.2)]][[as.character(3)]] <- -100

df_cond_exp_wrt_t <- data.frame(time=seq(0,3,0.003),values_yt=unlist(Yt_values[[as.character(0.2)]], use.names = FALSE), values_with_two_step=unlist(cond_exp_two_step[[as.character(0.2)]], use.names = FALSE))
plot_cond_exp_wrt_t <- ggplot(df_cond_exp_wrt_t, aes(x=time))+
  geom_line(aes(y=values_yt, linetype="solid"))+
  geom_line(aes(y=values_with_two_step), linetype="dotted")+
  ylim(100, 225)

#We will now simulate the conditional expectation for our original construction with delta^(1/6)
#in order to compare the results and decide whether the new approximation really performs better.
prelim_est_list_plot <- c()
for(epsilon in c(0.2)){
  tau_eps <- epsilon^(1/6)
  calculated_median <- Xt_values[[as.character(epsilon)]][[as.character(ceiling(round(tau_eps,3)*1000/3)*3*0.001)]]
  prelim_est_list_plot[[as.character(round(epsilon,3))]] <- calculated_median/(f*y0*((exp(a*tau_eps)/a)-(1/a)))
}
riccati_eq <- function(t,y, parms){
  list(2*a*y-((y^2*v^2*f^2)/(sigma^2))+b^2)
}
sol_riccati_eq_for_diff_prelim_est_list_plot <- list()
times <- seq(from = 0, to = 3, by = 0.003)
for(epsilon in c(0.2)){
  epsilon_rounded <- round(epsilon,3)
  #We have to insert the preliminary estimator at time tau_eps
  v <- prelim_est_list_plot[[as.character(epsilon_rounded)]]
  #The following vector includes all parameters needed in the riccati_eq function
  parameter <- c(a=a, f=f, sigma=sigma, b=b, v = v)
  sol_riccati_eq <- ode(y = 0, times = times, func = riccati_eq, parms = parameter)
  sol_riccati_eq_for_diff_prelim_est_list_plot[[as.character(epsilon_rounded)]] <- sol_riccati_eq
}
get_values_of_sol_of_riccati_eq_at_t_and_epsilon <- function(t, epsilon){
  index <- (t/3)*1000 + 1002
  value <- sol_riccati_eq_for_diff_prelim_est_list_plot[[as.character(round(epsilon,3))]][[index]]
  return(value)
}
calculate_cond_exp_at_eps <- function(f, a, sigma, b, y0, realizations, epsilon){
  epsilon_rounded <- round(epsilon,3)
  sim_cond_exp_drift <- expression((a-(get_values_of_sol_of_riccati_eq_at_t_and_epsilon(t,epsilon)*prelim_est_list[[as.character(epsilon_rounded)]]*f^2)/sigma^2)*x+(get_values_of_sol_of_riccati_eq_at_t_and_epsilon(t,epsilon)*f^2*prelim_est_list[[as.character(epsilon_rounded)]])/(sigma^2)*(Yt_values[[as.character(epsilon_rounded)]][[as.character(t)]]))
  sim_cond_exp_diffusion <- expression(epsilon_rounded*get_values_of_sol_of_riccati_eq_at_t_and_epsilon(t,epsilon)*f/sigma)
  sim_cond_exp <- snssde1d(drift=sim_cond_exp_drift, diffusion=sim_cond_exp_diffusion, M=realizations, x0=y0, t0=0, T=3, Dt=0.003, method="euler")
  return(sim_cond_exp)
}
cond_exp_list_plot <- list()
cond_exp_list_before_calculating_median_plot <- list()
for(epsilon in eps_vec){
  epsilon_rounded <- round(epsilon,3)
  cond_exp_list_before_calculating_median_plot[[as.character(epsilon_rounded)]] <- calculate_cond_exp_at_eps(f, a, sigma, b, y0, realizations, epsilon)$X
  for(i in 0:1000){
    cond_exp_list_plot[[as.character(epsilon_rounded)]][[as.character(round(i/1000*3,3))]] <- get_median_values_of_process(round(i/1000,3), epsilon, cond_exp_list_before_calculating_median_plot, realizations, 1000)
  }
}
cond_exp_original_delta <- list()
for(i in 0:254){
  cond_exp_original_delta[[as.character(0.2)]][[as.character(0+i*0.003)]] <- -100
}
for(i in 0:745){
  cond_exp_original_delta[[as.character(0.2)]][[as.character(0.765+i*0.003)]] <- cond_exp_list_plot[[as.character(0.2)]][[as.character(0.765+i*0.003)]]
}

df_cond_exp_wrt_t_with_delta_comparison <- data.frame(time=seq(0,3,0.003),values_yt=unlist(Yt_values[[as.character(0.2)]], use.names = FALSE), values_with_two_step=unlist(cond_exp_two_step[[as.character(0.2)]], use.names = FALSE), values_with_original_delta=unlist(cond_exp_original_delta, use.names = FALSE))
plot_cond_exp_wrt_t_with_delta_comparison <- ggplot(df_cond_exp_wrt_t_with_delta_comparison, aes(x = time)) +
  geom_line(aes(y = values_yt, linetype = "Real Value of Process Yt")) +
  geom_line(aes(y = values_with_two_step, linetype = "Two-step MLE-process")) +
  geom_line(aes(y = values_with_original_delta, linetype = "One-step MLE-process")) +
  ylim(100, 225) +
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.25, linetype = "solid")) +
  labs(linetype = 'Estimators') +
  scale_linetype_manual(values = c("Real Value of Process Yt" = "solid",
                                   "One-step MLE-process" = "dotted",
                                   "Two-step MLE-process" = "dashed")) +
  theme(legend.position = c(0.8, 0.2)) +
  theme(legend.background = element_rect(colour = "black"))+
  ylab("values")

