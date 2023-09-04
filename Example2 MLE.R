library(Sim.DiffProc)
library(ggplot2)
library(deSolve)
library(dplyr)
library(tictoc)

set.seed(1234) #used to create reproducible results
f <- 3
a <- 0.25
sigma <- 100
b <- 100
x0 <- 0 #initial value of X
y0 <- 100 #initial value of Y
theta <- 1 #real value of theta; we will use this to simulate the processes X and Y
realizations <- 100 #we will simulate 100 processes each

#Note that there is one major difference: instead of different epsilons,
#we now simulate the processes for fixed epsilon, but different theta in Theta.
#Let us therefore fix epsilon = 0.2
epsilon <- 0.2

#Let us check the following theta:
theta_vec <- seq(1,100,0.1)

tic("total")
#----------
# Simulate processes X_t and Y_t for the real theta=1
# Since we do not have different epsilon, we only need a (median) list
# for the different times
#----------
tic("simulation of processes X_t and Y_t")
print("Simulate processes X_t and Y_t")
sim_x_y_drift <- expression(f*y,theta*a*y)
sim_x_y_diffusion <- expression(epsilon*sigma, epsilon*b)
#The following function simulated a 2-dimensional SDE of Itô-type; more details can be found in the bibliography.
sim_x_y <- snssde2d(drift=sim_x_y_drift, diffusion=sim_x_y_diffusion, M=realizations, x0=c(x0,y0), t0=0, T=3, Dt=0.003, method="euler")

#We need the values of X_t and Y_t at the corresponding times 0, 0.003, ..., 3.
get_median_values_of_process <- function(t, process, realizations, time_steps){
  vec_at_t <- c()
  index <- (t*time_steps)+1 #we need +1 since we start at time 0, which has index 1
  for(i in 0:(realizations-1)){
    vec_at_t <- c(vec_at_t, process[[i*(time_steps+1) + index]])
  }
  return(median(vec_at_t))
}

Xt_values <- list()
Yt_values <- list()
for(i in 0:1000){
  Xt_values[[as.character(round(i/1000*3,3))]] <- get_median_values_of_process(round(i/1000,3), sim_x_y$X, realizations, time_steps=1000)
  Yt_values[[as.character(round(i/1000*3,3))]] <- get_median_values_of_process(round(i/1000,3), sim_x_y$Y, realizations, time_steps=1000)
}
toc()
#----------
# Calculate solutions of the Riccati equations
#----------
tic("calculation of solutions of Riccati equations")
print("Calculate the solutions of the Riccati equations")
#The following function is represents the ODE we have to solve:
riccati_eq <- function(t,y, parms){
  list(2*v*a*y-((y^2*f^2)/(epsilon^2*sigma^2))+epsilon^2*b^2)
}

#We will save the results in the following list
sol_riccati_eq_for_diff_theta <- list()

times <- seq(from = 0, to = 3, by = 0.003)

for(theta in theta_vec){
  v <- theta
  #The following vector includes all parameters needed in the riccati_eq function
  parameter <- c(a=a, f=f, sigma=sigma, b=b, v = v, epsilon=epsilon)
  sol_riccati_eq <- ode(y = 0, times = times, func = riccati_eq, parms = parameter)
  sol_riccati_eq_for_diff_theta[[as.character(round(theta,3))]] <- sol_riccati_eq
}

#We need the values of the solution of the Riccati equation for corresponding t 
#and preliminary estimator for the calculation of the conditional expectation
get_values_of_sol_of_riccati_eq_at_t_and_theta <- function(t, theta){
  index <- (t/3)*1000 + 1002
  value <- sol_riccati_eq_for_diff_theta[[as.character(round(theta,3))]][[index]]
  return(value)
}
toc()
#----------
# Calculate values of the conditional expectation
#----------
tic("calculation of conditional expectations")
print("Calculate value of conditional expectation")
#The following function will simulate the one-dimensional SDE m(.,t)
calculate_cond_exp_at_theta <- function(f, a, sigma, b, y0, realizations, theta){
  sim_cond_exp_drift <- expression(theta*a*x+(get_values_of_sol_of_riccati_eq_at_t_and_theta(t,theta)*f^2)/(epsilon^2*sigma^2)*(Yt_values[[as.character(t)]]-x))
  sim_cond_exp_diffusion <- expression(get_values_of_sol_of_riccati_eq_at_t_and_theta(t,theta)*f/(sigma*epsilon))
  sim_cond_exp <- snssde1d(drift=sim_cond_exp_drift, diffusion=sim_cond_exp_diffusion, M=realizations, x0=y0, t0=0, T=3, Dt=0.003, method="euler")
  return(sim_cond_exp)
}

#We will save our results in the list "cond_exp_list". The other list exists for technical reasons only.
cond_exp_list <- list()
cond_exp_list_before_calculating_median <- list()
for(theta in theta_vec){
  cond_exp_list_before_calculating_median[[as.character(round(theta,3))]] <- calculate_cond_exp_at_theta(f, a, sigma, b, y0, realizations, theta)$X
  for(i in 0:1000){
    cond_exp_list[[as.character(theta)]][[as.character(round(i/1000*3,3))]] <- get_median_values_of_process(round(i/1000,3), cond_exp_list_before_calculating_median[[as.character(round(theta,3))]], realizations, 1000)
  }
}
rm(cond_exp_list_before_calculating_median)
toc()
#-----------
# Calculate (log-)likelihood function
# For that, we need to calculate the ds- and the dWs-integral
#-----------
tic("calculation of (log-)likelihood function and of the MLE")
print("Calculate ds- and dWs-integral and log-likelihood")
ds_integral <- list()
for(theta in theta_vec){
  value_integral <- 0
  for(i in 0:1000){
    value_integral <- value_integral + 0.001*((f^2*cond_exp_list[[as.character(round(theta,3))]][[as.character(round((i/1000)*3,3))]]*Yt_values[[as.character(round((i/1000)*3,3))]])/(epsilon^2*sigma^2)-(f^2*(cond_exp_list[[as.character(round(theta,3))]][[as.character(round((i/1000)*3,3))]])^2)/(2*epsilon^2*sigma^2))
  }
  ds_integral[[as.character(round(theta,3))]] <- value_integral
}

dWs_integral <- list()
for(theta in theta_vec){
  ito_integrand <- expression((f*cond_exp_list[[as.character(round(theta,3))]][[as.character(t)]])/(epsilon*sigma))
  sim_of_stochastic_integral <- Sim.DiffProc::st.int(expr=ito_integrand,type="ito",M=realizations,lower=0,upper=3, subdivisions=1000)
  median_stoch_integral <- list()
  for(i in 0:1000){
    median_stoch_integral[[as.character(theta)]][[as.character(round(i/1000*3,3))]] <- get_median_values_of_process(round(i/1000,3), sim_of_stochastic_integral$X, realizations, 1000)
  }
  value_integral <- 0
  for(i in 0:1000){
    value_integral <- value_integral + 0.001*(median_stoch_integral[[as.character(theta)]][[as.character(round(i/1000*3,3))]])
  }
  dWs_integral[[as.character(round(theta,3))]] <- value_integral
}

#Note that we can maximize the log-likelihood instead of the likelihood itself
log_likelihood_list <- list()
for(theta in theta_vec){
  log_likelihood_list[[as.character(round(theta,3))]] <- ds_integral[[as.character(round(theta,3))]]+dWs_integral[[as.character(round(theta,3))]]
}

mle_estimator_with_X <- data.frame(log_likelihood_list) %>%
  mutate(Theta = names(.)[max.col(.)])
mle_estimator <- gsub("X","",mle_estimator_with_X[["Theta"]])%>%
  as.numeric()
rm(mle_estimator_with_X)

print(paste0("The MLE is ", mle_estimator,"."))
toc()
toc()

log_likelihood_df <- data.frame(value=unlist(log_likelihood_list, use.names = TRUE), theta=seq(1,100,0.1))
log_likelihood_plot <- ggplot(log_likelihood_df, aes(x=theta))+
  geom_line(aes(y=value))+
  labs(x ="\U03D1", y = "value of log-likelihood function")+
  theme(panel.background = element_rect(fill = "white", colour = "black",
                                        size = 0.25, linetype = "solid"))


