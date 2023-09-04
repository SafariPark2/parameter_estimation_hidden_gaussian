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

#eps_vec <- c(0.9, 0.5, 0.3, 0.25, 0.2, 0.15, 0.125, 0.1, 0.05)
eps_vec <- c(0.2)
#eps_vec <- seq(from = 0.3, to = 0.25, by = -0.005) 

tic("total")
#----------
# Simulate processes X_t and Y_t and calculate its median values
#----------
tic("simulation of processes X_t and Y_t")
print("Simulate processes X_t and Y_t")
x_for_epsilon_list <- c()
y_for_epsilon_list <- c()
for(epsilon in eps_vec){
  sim_x_y_drift <- expression(f*y,theta*a*y)
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
  tau_eps <- epsilon^(1/6)
  #We have to check whether round(tau_eps,3) is divisible by 3;
  #if not, we have to round to the next number which is divisible by 3.
  calculated_median <- Xt_values[[as.character(epsilon)]][[as.character(ceiling(round(tau_eps,3)*1000/3)*3*0.001)]]
  #The following calculation is based on the results of the second example.
  prelim_est_list[[as.character(round(epsilon,3))]] <- 2/(f*a*y0*(tau_eps^2))*(calculated_median-y0*f*tau_eps)
}
toc()
#-------------
# For the One-step MLE-process:
#-------------
#----------
# solve Riccati equation
# to solve: 
# \partial \gamma(\vartheta,t)/\partial t 
#    = 2\vartheta*a*\gamma(\vartheta,t) 
#        - (\gamma(\vartheta,t)^2*f^2)/(\varepsilon^2\sigma^2)+\varepsilon^2b^2
# initial value: \gamma(\vartheta,0)=0
# Note that we have to use the prelim est at tau_eps instead of epsilon!
#----------
print("Calculate the solutions of the Riccati equations")
tic("calculation of solutions of Riccati equations")
#The following function is represents the ODE we have to solve:
riccati_eq <- function(t,y, parms){
  list(2*v*a*y-((y^2*f^2)/(epsilon^2*sigma^2))+epsilon^2*b^2)
}

#We will save the results in the following list
sol_riccati_eq_for_diff_prelim_est_list <- list()

times <- seq(from = 0, to = 3, by = 0.003)

for(epsilon in eps_vec){
  epsilon_rounded <- round(epsilon,3)
  #We have to insert the preliminary estimator at time tau_eps
  v <- prelim_est_list[[as.character(epsilon_rounded)]]
  #The following vector includes all parameters needed in the riccati_eq function
  parameter <- c(a=a, f=f, sigma=sigma, b=b, v = v, epsilon=epsilon)
  sol_riccati_eq <- ode(y = 0, times = times, func = riccati_eq, parms = parameter)
  sol_riccati_eq_for_diff_prelim_est_list[[as.character(epsilon_rounded)]] <- sol_riccati_eq
}

#We need the values of the solution of the Riccati equation for corresponding t 
#and preliminary estimator for the calculation of the conditional expectation
get_values_of_sol_of_riccati_eq_at_t_and_epsilon <- function(t, epsilon){
  index <- (t/3)*1000 + 1002
  value <- sol_riccati_eq_for_diff_prelim_est_list[[as.character(round(epsilon,3))]][[index]]
  return(value)
}
toc()
#----------
# Calculate conditional expectation
# For each preliminary estimator, we have to solve the equation just once!
#----------
print("Calculate value of conditional expectation")
tic("calculation of conditonal expectation")
#The following function will simulate the one-dimensional SDE m(.,t)
calculate_cond_exp_at_eps <- function(f, a, sigma, b, y0, realizations, epsilon){
  epsilon_rounded <- round(epsilon,3)
  sim_cond_exp_drift <- expression(prelim_est_list[[as.character(epsilon_rounded)]]*a*x+(get_values_of_sol_of_riccati_eq_at_t_and_epsilon(t,epsilon)*f^2)/(epsilon^2*sigma^2)*(Yt_values[[as.character(epsilon_rounded)]][[as.character(t)]])-(get_values_of_sol_of_riccati_eq_at_t_and_epsilon(t,epsilon)*f^2)/(epsilon^2*sigma^2)*x)
  sim_cond_exp_diffusion <- expression((1/epsilon_rounded)*get_values_of_sol_of_riccati_eq_at_t_and_epsilon(t,epsilon)*f/sigma)
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
}
rm(cond_exp_list_before_calculating_median)
toc()
#----------
# Calculate dot(M)
#----------
print("Calculate integral of gamma star")
tic("calculation of integral of gamma star")
integral_gamma_star <- list()
for(epsilon in eps_vec){
  epsilon_rounded <- round(epsilon,3)
  for(i in 0:1000){
    integral_gamma_star[[as.character(epsilon_rounded)]][[as.character(round(i/1000*3,3))]][[as.character(round(i/1000*3,3))]] <- 0.001*sol_riccati_eq_for_diff_prelim_est_list[[as.character(epsilon_rounded)]][[1002+i]]
  }
  for(i in 0:999){
    accumulated_value <- integral_gamma_star[[as.character(epsilon_rounded)]][[as.character(round(i/1000*3, 3))]][[as.character(round(i/1000*3, 3))]]
    for(j in (i+1):1000){
      if(j<=1000){
        accumulated_value <- accumulated_value + 0.001 * sol_riccati_eq_for_diff_prelim_est_list[[as.character(epsilon_rounded)]][[1001 + j]]
        integral_gamma_star[[as.character(epsilon_rounded)]][[as.character(round(i/1000*3,3))]][[as.character(round(j/1000*3,3))]] <- accumulated_value
      }
    }
  }
}
toc()
#We can now calculate dot(M). Note that we only need its values for times 
#[tau_eps, 1], hence we can neglect the calculcation of all values of dot(M) 
#for times < tau_eps, which will save a lot of running time.
#Hence, only calculate dot(M) for UPPER bounds >= tau_eps.
tic("calculation of dot(M)")
print("Calculate values of dot(M)")
dot_M <- list()
for(epsilon in eps_vec){
  epsilon_rounded <- round(epsilon,3)
  for(upper_bound in ceiling((round(epsilon^(1/6),3)*1000/3)):1000){
    value_integral <- 0
    for(i in 0:(upper_bound-1)){
      value_integral <- value_integral + 0.001*exp(-f^2/sigma^2*integral_gamma_star[[as.character(epsilon_rounded)]][[as.character(round(0.001*i*3,3))]][[as.character(round(upper_bound/1000*3,3))]])
    }
    dot_M[[as.character(epsilon_rounded)]][[as.character(round(upper_bound/1000*3,3))]] <- f*a*y0*exp(prelim_est_list[[as.character(epsilon_rounded)]]*a*(upper_bound/1000))*value_integral
  }
}
toc()
#----------
# Calculate Fisher information
#----------
print("Calculate Fisher information")
tic("calculation of Fisher information")
fisher_info <- list()
for(epsilon in eps_vec){
  epsilon_rounded <- round(epsilon,3)
  tau_eps_rounded <- round(epsilon^(1/6),3)
  for(upper_bound in (ceiling((tau_eps_rounded*1000/3)/3)*3):1000){
    value_integral <- 0
    for(i in 0:(upper_bound-round((epsilon^(1/6))/3*1000))){
      value_integral <- value_integral + 0.001*((dot_M[[as.character(epsilon_rounded)]][[as.character(round(upper_bound/1000*3,3))]])^2/sigma^2)
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
  tau_eps <- epsilon^(1/6)
  tau_eps_rounded <- round(tau_eps,3)
  epsilon_rounded <- round(epsilon,3)
  value_integral <- 0
  if(round(tau_eps,3)*1000/3+2<1000){
    number_of_simulation_steps <- length(unlist(dot_M[[as.character(round(epsilon,3))]], use.names= FALSE))-1
    ito_integrand <- expression(dot_M[[as.character(round(epsilon,3))]][[as.character(t)]]*epsilon/sigma)
    sim_of_stochastic_integral <- Sim.DiffProc::st.int(expr=ito_integrand,type="ito",M=realizations,lower=(ceiling(round(epsilon^(1/6),3)*1000/3)*3/1000),upper=3, subdivisions=number_of_simulation_steps)
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
  tau_epsilon_rounded <- round(epsilon^(1/6),3)
  #Let us first calculate the value of the integral involving Y_s and the conditional expectation
  for(upper_bound in (ceiling((round((epsilon^(1/6)),3)*1000/3))+2):1000){
    upper_div_by_1000_times_3_rounded <- round(upper_bound/1000*3,3)
    value_integral <- 0
    for(i in (tau_epsilon_rounded*1000+2):upper_bound){
      value_integral <- value_integral + 0.001*((dot_M[[as.character(epsilon)]][[as.character(upper_div_by_1000_times_3_rounded)]])/sigma^2)*f*(Yt_values[[as.character(epsilon_rounded)]][[as.character(upper_div_by_1000_times_3_rounded)]]-cond_exp_list[[as.character(epsilon_rounded)]][[as.character(upper_div_by_1000_times_3_rounded)]])
    }
    ds_integral[[as.character(epsilon_rounded)]][[as.character(upper_div_by_1000_times_3_rounded)]] <- value_integral
  }
  for(upper_bound in ceiling((round((epsilon^(1/6)),3)*1000/3+2)):1000){
    upper_div_by_1000_times_3_rounded <- round(upper_bound/1000*3,3)
    one_step_mle_process[[as.character(epsilon,3)]][[as.character(upper_div_by_1000_times_3_rounded)]] <- prelim_est_list[[as.character(epsilon_rounded)]] + 1/fisher_info[[as.character(epsilon_rounded)]][[as.character(upper_div_by_1000_times_3_rounded)]]*ds_integral[[as.character(epsilon_rounded)]][[as.character(upper_div_by_1000_times_3_rounded)]]+1/fisher_info[[as.character(epsilon_rounded)]][[as.character(upper_div_by_1000_times_3_rounded)]]*stochastic_integral[[as.character(epsilon_rounded)]][[as.character(upper_div_by_1000_times_3_rounded)]]
  }
}
toc()
toc()
#----------
# Create a plot regarding the consistency of the preliminary and One-step MLE-process with respect to epsilon
#----------
one_step_mle_process_at_time_3 <- list()
for(epsilon in eps_vec){
  one_step_mle_process_at_time_3[[as.character(round(epsilon,3))]] <- one_step_mle_process[[as.character(epsilon,3)]][[as.character(3)]]
}
prelim_values <- unlist(prelim_est_list, use.names=FALSE)
one_step_values <- unlist(one_step_mle_process_at_time_3, use.names=FALSE)
df_consistency_wrt_eps <- data.frame(eps_vec, prelim_values, one_step_values)

plot_consistency_wrt_eps <- ggplot(df_consistency_wrt_eps, aes(x=eps_vec)) +
  geom_line(aes(y=prelim_values,linetype="dotted"))+
  geom_line(aes(y=one_step_values,linetype="solid"))+
  xlim(1,0) +
  #ylim(1,1.12)+
  geom_hline(yintercept=1, size=1) +
  labs(x ="ε", y = "Estimators") +
  theme(panel.background = element_rect(fill = "white", colour = "black",
                                        size = 0.25, linetype = "solid"))+
  labs(linetype='Estimators')+
  scale_linetype_manual(values = c(1, 3),
                        labels = c("Preliminary Estimator", "One-step MLE-process"))+
  theme(legend.position = c(0.8, 0.8))+
  theme(legend.background = element_rect(colour ="black"))

#----------
# Create a plot regarding the consistency with respect to t for constant coefficients, vartheta_0=0, epsilon = 0.5
#----------

#The following definitions are technical necessities to layer the lines. 
#We define all unnecessary times to have values 0 so that these values do not appear in the plot.
values_one_step <- c()
values_prelim <- c()
df_consistency_wrt_t_list <- list()
values_real <- rep(1, 1000)

for(epsilon in eps_vec){
  one_step_temp <- c(rep(-100, round(epsilon^(1/6)/3+0.001,3)*1000),unlist(one_step_mle_process[[as.character(epsilon)]], use.names=FALSE))
  if(length(one_step_temp)<1000){
    one_step_temp <- c(one_step_temp, rep(-100, 1000-length(one_step_temp)))
  }
  values_one_step[[as.character(epsilon)]] <- one_step_temp
  values_prelim[[as.character(epsilon)]] <- c(rep(prelim_est_list[[as.character(epsilon)]], ceiling(round(epsilon^(1/6)+0.001,3)*1000/3)), rep(-100, 1000-ceiling(round(epsilon^(1/6)+0.001,3)*1000/3)))
  df_consistency_wrt_t_list[[as.character(epsilon)]] <- data.frame(times=seq(0.003,3,0.003), values_prelim = values_prelim[[as.character(epsilon)]], values_one_step = values_one_step[[as.character(epsilon)]], values_real = values_real)
}
for(epsilon in eps_vec){
  print(df_consistency_wrt_t_list[[as.character(epsilon)]]$values_one_step[[999]])
}

plot_consistency_wrt_t_full <- ggplot(df_consistency_wrt_t_list[[as.character(0.285)]], aes(x=times)) +
  geom_line(aes(y=values_prelim,linetype="solid"))+ #It should be the other way around, but this works
  geom_line(aes(y=values_one_step,linetype="dotted"))+
  geom_line(aes(y=values_real,linetype="dashed"))+
  labs(x ="t", y = "One-step MLE-process") +
  ylim(-30,2)+
  theme(panel.background = element_rect(fill = "white", colour = "black",
                                        size = 0.25, linetype = "solid"))+
  labs(linetype='Estimators')+
  scale_linetype_manual(values = c(1, 2, 3),
                        labels = c("Real value", "One-step MLE-process", "Preliminary Estimator"))+
  theme(legend.position = c(0.8, 0.2))+
  theme(legend.background = element_rect(colour ="black"))

plot_consistency_wrt_t_partially <- ggplot(df_consistency_wrt_t_list[[as.character(0.285)]], aes(x=times)) +
  geom_line(aes(y=values_prelim,linetype="solid"))+ #It should be the other way around, but this works
  geom_line(aes(y=values_one_step,linetype="dotted"))+
  geom_line(aes(y=values_real,linetype="dashed"))+
  ylim(-4,2)+
  labs(x ="t", y = "One-step MLE-process") +
  theme(panel.background = element_rect(fill = "white", colour = "black",
                                        size = 0.25, linetype = "solid"))+
  labs(linetype='Estimators')+
  scale_linetype_manual(values = c(1, 2, 3),
                        labels = c("Real value", "One-step MLE-process", "Preliminary Estimator"))+
  theme(legend.position = c(0.8, 0.2))+
  theme(legend.background = element_rect(colour ="black"))

#----------
# Create a histogramm to display the asymptotic normality
#----------
histogram_list <- list()
values_for_histogram <- list()
trimmed_values_for_histogram <- list()
for(epsilon in eps_vec){
  epsilon_rounded <- round(epsilon,3)
  values_for_histogram[[as.character(epsilon_rounded)]] <- (unlist(one_step_mle_process[[as.character(epsilon_rounded)]], use.names = FALSE)-theta)/epsilon
  #Let us trim these values a bit in order to remove outliers
  trim_q <- function(x, lb, ub){
    x[(x > quantile(x, lb)) & (x < quantile(x, ub))]
  }
  trimmed_values_for_histogram[[as.character(epsilon_rounded)]] <- data.frame(values=trim_q(values_for_histogram[[as.character(epsilon_rounded)]], 0.0, 1))
  
  histogram_list[[as.character(epsilon)]] <- ggplot(trimmed_values_for_histogram[[as.character(epsilon_rounded)]], aes(x = values)) +
    geom_histogram(aes(y = ..density..), bins = 60, fill = "white", color = "black") +
    stat_function(
      fun = dnorm,
      args = list(mean = mean(trimmed_values_for_histogram[[as.character(epsilon_rounded)]]$values), sd = sd(trimmed_values_for_histogram[[as.character(epsilon_rounded)]]$values)),
      color = "black"
    ) +
    xlab("x") +
    ylab("Density Histogram")+
    theme(panel.background = element_rect(fill = "white", colour = "black",
                                          size = 0.25, linetype = "solid"))
}

histogram_list[[as.character(0.295)]] 


#set up QQ-Plot to (visually) check for normality
ggplot(trimmed_values_for_histogram[[as.character(0.295)]], aes(sample=values)) +
  stat_qq() + stat_qq_line() + theme(panel.background = element_rect(fill = "white", colour = "black",
                                                                     size = 0.25, linetype = "solid"))
#set up Kolmogorov-Smirnov test to check for normality
ks.test(unlist(trimmed_values_for_histogram[[as.character(epsilon_rounded)]],use.names = FALSE), "pnorm", mean=mean(unlist(trimmed_values_for_histogram[[as.character(epsilon_rounded)]],use.names = FALSE)), sd=sd(unlist(trimmed_values_for_histogram[[as.character(epsilon_rounded)]],use.names = FALSE)))
#set up Shapiro-Wilk test to check for normality
shapiro.test(unlist(trimmed_values_for_histogram[[as.character(epsilon_rounded)]],use.names = FALSE))


