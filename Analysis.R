#Data visualization.
library(ggplot2)
#Graph Legend
library(cowplot)
#Graphs
library(gridExtra)
library(grid)
#Volatility.
library(TTR)
#Spline functions.
library(splines)
#Merge ggplot2 plots.
library(patchwork)
#Bootstrap 
library(boot)
#Maximum likelihood estimation.
library(stats4)
#Nonlinear optimization.
library(nloptr)
#Anderson Darling tests.
library(nortest)
#Time series analysis.
library(tseries)
#MCMC diagnostics.
library(coda)
#Solve nonlinear equations.
library(nleqslv)
#Hawkes process modeling.
library(hawkes)
#Kolmogorov–Smirnov test
library(MASS)
#Testing for Hawkes
library(stats)
#Truncated Distribution
library(truncdist)
#QQ plot
library(QRM)
#Progress Bar
library(progress)

# Seed
set.seed(13793632)

# Template for output of graphs
plot_comparison <- function(data, title = "")
{
  ggp <- ggplot(data, aes(x = Day, y = Price, color = group)) +
    geom_line(linetype = "solid") + 
    labs(
      title = title,
      x = "Day",
      y = "Price"
    ) +
    theme_minimal(base_size = 15) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
  return(ggp)
}

# Load Data
raw <- read.csv("^IXIC.csv", sep = ",")

# Data Refine
data <- raw[head(which(raw$Volume > 0), 1):nrow(raw), ]
data <- data[!apply(data, 1, function(x) any(is.na(x) | is.nan(x) | is.infinite(x))), ]
rownames(data) <- NULL

# Data Process
date <- as.Date(data$Date, format="%Y-%m-%d")
open <- data$Open
low <- data$Low
high <- data$High
close <- data$Close
logreturn <- diff(log(close))
N <- length(date)
Close_data <- data.frame(Day = seq(1, length(close)), Price = close, group = "Close Data")

# 1 Group Student's t test 
t.test(logreturn, mu = 0)

# Calculation of Volatility
volatilities <- volatility(cbind(open, high, low, close), n = 30, calc = "rogers.satchell", N=1)
plot(seq(30:N), na.omit(volatilities), type="l", main="Volatility over time", xlab = "Days", ylab="Volatility")
vol <- mean(na.omit(volatilities))

# Drift
drift <- mean(logreturn) + 0.5 * vol^2

#Anderson-Darling Test
ad.test(logreturn)

#Histogram
returnhist <- hist(logreturn,breaks=2500,xlim=c(-0.15,0.15))
returnnorm <- curve(dnorm(x, drift - 0.5 * vol^2, vol), -0.15, 0.15)
plot(returnhist)
lines(returnnorm, col="red")


# GBM Simulation
GBM_path <- function(S0, mu, sigma, N) 
{
  GBMsim <- vector()
  GBMsim[1] <- S0
  t <- seq(1, N)
  for (i in 2:N)
  {
    GBMsim[i] <- GBMsim[i-1] * exp((mu - 0.5 * sigma^2) + sigma * rnorm(1))
  }
  return(data.frame(Day = t, Price = GBMsim, group = "GBM Simulation"))
}

# Comparison of GBM and real
GBMgraph <- list()
for(i in 1:6)
{
  GBMsim <- GBM_path(open[1], drift, vol, N)
  GBMgraph[[i]] <- plot_comparison(rbind(GBMsim, Close_data))
}
legend <- get_legend(GBMgraph[[1]] + theme(legend.position = "bottom"))
grid.arrange(
  do.call(arrangeGrob, c(GBMgraph, nrow = 3)),
  legend,
  ncol = 1,
  heights = c(10, 1),
  top = textGrob("GBM Simulation Comparisons", gp = gpar(fontsize = 16, fontface = "bold"))
)

#Expected value and Variance
GBMexp <- exp(drift*N)*open[1]
GBMvar <- exp(2*drift*N)*open[1]^2*(exp(vol^2*N)-1)

#Merton Jump Diffusion Model
#Jump Frequency
jump_detect <- function(alpha, returns) 
{
  return(sum(abs(returns) >= alpha))
}
alpha <- seq(0.00, max(logreturn), length.out = 100)
jumpcounts <- sapply(alpha, jump_detect, returns = logreturn)
cv.error <- rep(NA,80)
for (i in 1:length(cv.error)) 
{
  fit <- glm(jumpcounts ~ ns(alpha, df = i))
  cv.error[i] <- cv.glm(data.frame(alpha, jumpcounts), fit, K = 10)$delta[1]
}
optimaldf <- which.min(cv.error)

#cv.glm randomly separates the data and therefore optimal df varies
optimfit <- lm(jumpcounts ~ ns(alpha, df = optimaldf))
jumpplot <- predict(optimfit, newdata = data.frame(alpha))
curvature_func <- function(func, alpha) 
{
  d1 <- predict(func, data.frame(alpha), type = "terms")[, 1]
  d2 <- diff(d1) / diff(alpha)
  d2 <- c(d2, d2[length(d2)])
  curvature <- abs(d2) / (1 + d1^2)^(3/2)
  return(curvature)
}
curvature <- curvature_func(optimfit, alpha)
alpha_max <- alpha[which.max(curvature)]
detected_jumps <- which(abs(logreturn) > alpha_max)
plot(alpha, jumpplot, type = "l", col = "blue", 
     xlab = "alpha", ylab = "Jump Counts", main = "Jump Counts and Curvature")
points(alpha_max, jump_detect(alpha_max, logreturn), col = "red",pch = 19)
abline(v = alpha_max, col = "red", lty = 2)

#Jump frequency
frequency <- length(detected_jumps)/N

#Moments of jump size
jump <- vector()
for(i in seq_along(detected_jumps))
{
  if(logreturn[detected_jumps][i] > alpha_max)
  {
    jump[i] <- logreturn[detected_jumps][i] - alpha_max
  }
  else if(logreturn[detected_jumps][i] < alpha_max)
  {
    jump[i] <- logreturn[detected_jumps][i] + alpha_max
  }
}
jumpexp <- mean(jump)
jumpsd <- sd(jump)
jumphist <- hist(jump,breaks=350)
jumpnorm <- curve(dnorm(x,jumpexp,jumpsd),-0.10,0.10)
plot(jumphist)
lines(jumpnorm)

#Testing for autocorrelation of event frequency
whenjump <- numeric(N)
for(i in detected_jumps)
{
  whenjump[i] <- 1
}
acf(whenjump, main = "ACF plot for jumps")
finalp <- vector()
for(j in 1:50)
{
  ljung_box_test <- Box.test(jump, lag=j, type="Ljung-Box")
  finalp[j] <- ljung_box_test$p.value
}
plot(finalp,type="l", main = "p-value of Ljung Box test with different lags", xlab = "lags", ylab = "p-value")
abline(h = 0.05, col = "red", lty = 2)

#Testing for normality of non-jumps
ad.test(logreturn[-detected_jumps])

#Testing for normality of jumps
shapiro.test(jump)

#MJD simulation
MJD_path <- function(S0, drift, vol, N, lambda, jumpexp, jumpsd) 
{
  MJDsim <- vector()
  MJDsim[1] <- S0
  Poisson <- rpois(N, lambda)
  Wiener <- rnorm(N)
  
  diffusion <- drift -0.5*vol^2 - lambda * (exp(jumpexp + 0.5 * jumpsd^2) - 1)
  for (t in 2:N) 
  {
    dS <- diffusion + vol * Wiener[t] + Poisson[t]*rnorm(1, mean = jumpexp, sd = jumpsd)
    MJDsim[t] <- MJDsim[t - 1] * exp(dS)
  }
  return(data.frame(Day = seq(1:N), Price = MJDsim, group = "MJD Simulation"))
}

MJDgraph <- list()
for(i in 1:6)
{
  MJDsim <- MJD_path(open[1], drift, vol, N, frequency, jumpexp, jumpsd)
  MJDgraph[[i]] <- plot_comparison(rbind(MJDsim, Close_data))
}
legend <- get_legend(MJDgraph[[1]] + theme(legend.position = "bottom"))
grid.arrange(
  do.call(arrangeGrob, c(MJDgraph, nrow = 3)),
  legend,
  ncol = 1,
  heights = c(10, 1),
  top = textGrob("MJD Simulation Comparisons", gp = gpar(fontsize = 16, fontface = "bold"))
)
#MLE for EALD Parameters
Loglik_EALD <- function(params, data)
{
  p <- params[1]
  a <- params[2]
  b <- params[3]
  L <- 0
  
  for (i in data)
  {
    if (i>=0)
    {
      L <- L+log(p*a*exp(-a*i))
    }
    else
    {
      L <- L+log(-(1-p)*b*exp(-b*i))
    }
  }
  return(-L)
}
Initial_EALD <- c(0.5,0.5,-0.5)
Lower_EALD <- c(0, 0, -Inf)  
Upper_EALD <- c(1, Inf, 0)

Optim_EALD <- optim(par = Initial_EALD, fn = Loglik_EALD, data=jump, method = "L-BFGS-B", lower = Lower_EALD, upper = Upper_EALD)
Star_EALD <- Optim_EALD$par
Star_EALD

#Bootstrap for the confidence interval of p MLE
Loglik_EALD_p <- function(p, data)
{
  L <- 0
  a <- Star_EALD[2]
  b <- Star_EALD[3]
  for (i in data)
  {
    if (i>=0)
    {
      L <- L+log(p*a*exp(-a*i))
    }
    else
    {
      L <- L+log(-(1-p)*b*exp(-b*i))
    }
  }
  return(-L)
}
estimate_p_MLE <- function(data) {
  fit <- optim(par = 0.5, fn = Loglik_EALD_p, data = data, method = "Brent", lower = 0, upper = 1)
  return(fit$par)
}
bootstrap_MLE <- function(data, indices) {
  sample <- data[indices]
  return(estimate_p_MLE(sample))
}
results_p <- boot(data = jump, statistic = function(data, indices) bootstrap_MLE(data, indices), R = 1000)
conf_MLE <- boot.ci(results_p, type = "perc")
print(conf_MLE)
sd_p_MLE <- sd(results_p$t)

#Function A
A <- function(beta,event)
{
  T <- length(event)
  A <- numeric(T)
  A[1] <- 0
  for (i in 2:T) {
    A[i] <- exp(-beta * (event[i] - event[i-1])) * (1 + A[i-1])
  }
  return(A)
}

#MLE for Jump Intensity Parameters
Loglik_JIP <- function(params, T, NT, event) {
  lambda_inf <- params[1]
  alpha <- params[2]
  beta <- params[3]
  A <- A(beta,event)
  
  if (any(params < Lower_JIP)) {
    return(1e10)
  }
  term1 <- T - T * lambda_inf
  term2 <- - sum(alpha / beta * (1 - exp(-beta * (T - event))))
  term3 <- sum(log(lambda_inf + alpha * A))
  L <- term1 + term2 + term3
  
  return(-L)  
}
Initial_JIP <- c(frequency,runif(1),runif(1))
Lower_JIP <- c(.Machine$double.eps, .Machine$double.eps, .Machine$double.eps)
Optim_JIP <- optim(Initial_JIP, Loglik_JIP, T = N, NT = length(detected_jumps), event = detected_jumps, method = "Nelder-Mead")
Star_JIP <- Optim_JIP$par
Star_JIP


#Variance of beta
beta_variance <- function(lambda_inf,alpha,beta,event)
{
  T <- length(event)
  A <- A(beta,event)
  c <- numeric(T)
  c[1] <- 0
  B <- numeric(T)
  B[1] <- 0
  Last <- event[T]
  term1 <- sum(1/beta*(Last-event)^2*exp(-beta*(Last-event)))
  term2 <- sum(2/(beta^2)*(Last-event)*(exp(-beta*(Last-event))-1))
  term3 <- sum(2/(beta^3)*(exp(-beta*(Last-event))-1))
  for(i in 2:T)
  {
    valid <- event[1:i]
    c[i] <- sum((event[i]-valid)^2*exp(-beta*(event[i]-valid)))
  }
  term4 <- sum(alpha*c/(lambda_inf+alpha*A))
  for(i in 2:T)
  {
    valid <- event[event<event[i]]
    B[i] <- sum((event[i]-valid)*exp(-beta*(event[i]-valid)))
  }
  term5 <- sum((alpha*B/(lambda_inf+alpha*A))^2)
  return(alpha*(term1+term2+term3)+term4-term5)
}
beta_var <- 1/-beta_variance(Star_JIP[1],Star_JIP[2],Star_JIP[3],detected_jumps)
beta_sd <- sqrt(beta_var)

#MCMC for parameter beta
lambda_calculate <- function(mu,alpha,beta,data)
{
  sapply(1:length(data), function(i) {
    suming <- mu
    if (i > 1) {
      suming <- suming + sum(alpha * exp(-beta * (data[i] - data[1:(i-1)])))
    }
    if (suming == 0) {
      suming <- .Machine$double.eps
    }
    return(suming)
  })
}

B_Integral <- function(beta,data,T)
{
  sapply(1:length(data), function(i) {
    suming <- 0
    if (i > 1) {
      suming <- suming + sum(1/beta*exp(-beta * (T - data[1:(i-1)])))
    }
    return(suming)
  })
}

muMLE <- Star_JIP[3]
alphaMLE <- Star_JIP[1]
priorbeta <- function(beta)
{
  Star_JIP[2]^2/gamma(2)*beta^{-2-1}*exp(-Star_JIP[2]/2)
}
Excecute_MCMC <- function(initial_beta, iterations, sd)
{
  beta <- initial_beta
  beta_samples <- numeric(iterations)
  pb <- progress_bar$new(total = n, format = "[:bar] :percent Finishes in: :eta", clear = TRUE)
  for(i in 1:iterations)
  {
    pb$tick()
    beta_proposal <- rnorm(1, beta, sd)
    while (beta_proposal < 0)
    {
      beta_proposal <- rnorm(1, beta, sd)
    }
    original_lambda <- lambda_calculate(muMLE, alphaMLE, beta, detected_jumps)
    proposed_lambda <- lambda_calculate(muMLE, alphaMLE, beta_proposal, detected_jumps)
    Hastings <- priorbeta(beta_proposal) / priorbeta(beta) * prod(proposed_lambda / original_lambda) * exp(sum(alphaMLE * (B_Integral(beta, detected_jumps, N) - B_Integral(beta_proposal, detected_jumps, N))))
    if (runif(1) < min(1,Hastings)) {
      beta <- beta_proposal
    }
    beta_samples[i] <- beta
    Sys.sleep(1 / n)
  }
  return(beta_samples)
}
burnin <- 2000
Chain1_raw <- Excecute_MCMC(initial_beta = Star_JIP[2], iterations = 15000, sd = beta_sd)
Chain2_raw <- Excecute_MCMC(initial_beta = Star_JIP[2], iterations = 15000, sd = beta_sd)
Chain3_raw <- Excecute_MCMC(initial_beta = Star_JIP[2], iterations = 15000, sd = beta_sd)
Chain4_raw <- Excecute_MCMC(initial_beta = Star_JIP[2], iterations = 15000, sd = beta_sd)
Chain1 <- Chain1_raw[(burnin+1):length(Chain1_raw)]
Chain2 <- Chain2_raw[(burnin+1):length(Chain2_raw)]
Chain3 <- Chain3_raw[(burnin+1):length(Chain3_raw)]
Chain4 <- Chain4_raw[(burnin+1):length(Chain4_raw)]
Chain1_mcmc <- mcmc(Chain1)
Chain2_mcmc <- mcmc(Chain2)
Chain3_mcmc <- mcmc(Chain3)
Chain4_mcmc <- mcmc(Chain4)
Samples <- mcmc.list(Chain1_mcmc, Chain2_mcmc, Chain3_mcmc, Chain4_mcmc)

Gelman_Rubin <- gelman.diag(Samples)
print(Gelman_Rubin)

plot(Chain1_raw, type='l', col='blue', main='Trace Plot for MCMC Chains with burn-in', ylab='Value', xlab='Iterations')
lines(Chain2_raw, type='l', col='red')
lines(Chain3_raw, type='l', col='green')
lines(Chain4_raw, type='l', col='purple')
legend("topright", legend=c("Chain 1", "Chain 2", "Chain 3", "Chain 4"), col=c("blue", "red", "green", "purple"), lty=1)

plot(Chain1, type='l', col='blue', main='Trace Plot for MCMC Chains without burn-in', ylab='Value', xlab='Iterations')
lines(Chain2, type='l', col='red')
lines(Chain3, type='l', col='green')
lines(Chain4, type='l', col='purple')
legend("topright", legend=c("Chain 1", "Chain 2", "Chain 3", "Chain 4"), col=c("blue", "red", "green", "purple"), lty=1)

acf_colors <- c("blue", "red", "green", "purple")
acf_list <- list(acf(Chain1, plot=FALSE), acf(Chain2, plot=FALSE), acf(Chain3, plot=FALSE), acf(Chain4, plot=FALSE))
plot(acf_list[[1]], type="l", main="ACF Plot for MCMC Chains", col=acf_colors[1], xlab="Lag", ylab="ACF",ylim=c(0,1))
for (i in 2:length(acf_list)) {
  lines(acf_list[[i]]$lag, acf_list[[i]]$acf, type="l", col=acf_colors[i])
}
legend("bottomleft", legend=c("Chain 1", "Chain 2", "Chain 3", "Chain 4"), col=acf_colors, lty=1)

#Calculation of ESS and new estimated 
ESS <- effectiveSize(Samples)
print(ESS)
betaMCMC <- mean(cbind(Chain1,Chain2,Chain3,Chain4))

#Recalculation of MLE for alpha and lambda_infty
Score <- function(param, beta, event)
{
  lambda_inf <- param[1]
  alpha <- param[2]
  A <- A(beta,event)
  
  eq1 <- -max(event) + sum(1 / (lambda_inf + alpha * A))
  eq2 <- sum(1 / beta * (exp(-beta * (max(event) - event)) - 1)) + sum(A / (lambda_inf + alpha * A))
  
  return(c(eq1, eq2))
}
Initial_Score <- c(runif(1),runif(1))
Solution_Score <- nleqslv(Initial_Score, Score,beta=betaMCMC,event=detected_jumps)
Solution_lambda_inf <- Solution_Score$x[1]
Solution_alpha <- Solution_Score$x[2]

#HJD MLEMCMC Simulation
EALD_sample <- function(n,p,a,b)
{
  prob_pos <- rbinom(n, 1, p) == 1
  sample <- numeric(n)
  sample[prob_pos] <- rexp(sum(prob_pos), rate = a)
  sample[!prob_pos] <- -rexp(sum(!prob_pos), rate = -b)
  return(sample)
}

HJD_path_MLEMCMC <- function(S0, drift, vol, N, lambda_infty, alpha, beta, p, a, b) {
  HJDsim_MLEMCMC <- numeric(N)
  HJDsim_MLEMCMC[1] <- S0
  jumpindex <- unique(as.integer(simulateHawkes(lambda_infty, alpha, beta, 10000)[[1]]))
  jumplogic <- rep(0, N)
  jumplogic[jumpindex] <- 1
  lambda <- numeric(N)
  jumpsize <- EALD_sample(N, p, a, b)
  Wiener <- rnorm(N)
  
  for (i in 1:N) {
    valid <- jumpindex[jumpindex < i]
    if (length(valid) == 0) {
      lambda[i] <- lambda_infty
    } else {
      lambda[i] <- lambda_infty + sum(sapply(valid, function(ti) alpha * exp(-beta * (i - ti))))
    }
  }
  
  k <- p * a / (a - 1) + (1 - p) * b / (b - 1) - 1
  
  for (t in 2:N) {
    dS <- drift - 0.5 * vol^2 - lambda[t] * k + vol * Wiener[t] + jumplogic[t] * jumpsize[t]
    HJDsim_MLEMCMC[t] <- HJDsim_MLEMCMC[t - 1] * exp(dS)
  }
  
  return(data.frame(Day = seq(1:N), Price = HJDsim_MLEMCMC, group = "HJD Simulation with MLEMCMC"))
}
HJDgraph_MLEMCMC <- list()
for(i in 1:6)
{
  HJDsim_MLEMCMC <- HJD_path_MLEMCMC(open[1], drift, vol, N, Solution_lambda_inf, Solution_alpha, betaMCMC, Star_EALD[1], Star_EALD[2], Star_EALD[3])
  HJDgraph_MLEMCMC[[i]] <- plot_comparison(rbind(HJDsim_MLEMCMC, Close_data))
}
legend <- get_legend(HJDgraph_MLEMCMC[[1]] + theme(legend.position = "bottom"))
grid.arrange(
  do.call(arrangeGrob, c(HJDgraph_MLEMCMC, nrow = 3)),
  legend,
  ncol = 1,
  heights = c(10, 1),
  top = textGrob("HJD Simulation Comparisons with MLE & MCMC", gp = gpar(fontsize = 16, fontface = "bold"))
)

#Parameters for Distribution of Jumps
Moment_Distribution <- function(p,gam1, gam2,k)
{
  return((-1)^k*factorial(k)*p/(gam1^k)+factorial(k)*(1-p)/(gam2)^k)
}
GMM_Distribution <- function(par, data, k)
{
  p <- par[1]
  gam1 <- par[2]
  gam2 <- par[3]
  sample <- numeric(k)
  theoretical <- numeric(k)
  for (i in 1:k) {
    sample[i] <- mean(data^i)
    theoretical[i] <- Moment_Distribution(p,gam1,gam2,i)
  }
  diff <- sample - theoretical
  return(sum(diff^2))
}
Initial_Distribution <- c(Star_EALD[1],Star_EALD[2],-Star_EALD[3])
Lower_Distribution <- c(0, 0, 0)
Upper_Distribution <- c(1, Inf, Inf)

Optim_Distribution <- optim(par = Initial_Distribution, fn = GMM_Distribution, data=jump, k=3, method = "L-BFGS-B", lower = Lower_Distribution, upper = Upper_Distribution)
Star_Distribution <- Optim_Distribution$par
print(Star_Distribution)

#Bootstrap for the confidence interval of p GMM
GMM_Distribution_p <- function(p,data)
{
  a <- Star_Distribution[2]
  b <- Star_Distribution[3]
  
  m1 <- mean(data)-Moment_Distribution(p,a,b,1)
  return(abs(m1))
}
estimate_p <- function(data) {
  opt_res <- optim(par = 0.5, fn = GMM_Distribution_p, data = data, method = "Brent",lower = 0,upper = 1)
  return(opt_res$par)
}
bootstrap <- function(data, indices) {
  sample <- data[indices]
  return(estimate_p(sample))
}
results <- boot(jump, statistic = bootstrap, R = 100)
conf_GMM <- boot.ci(results, type = "perc")
print(conf_GMM)
sd_p_GMM <- sd(results$t)

#Parameters by GMM Estimations
Theoretical_GMM <- function(parameters) {
  mu <- parameters[1]
  sigma <- parameters[2]
  lambda_infty <- parameters[3]
  alpha <- parameters[4]
  beta <- parameters[5]
  
  lambda <- alpha*lambda_infty/(alpha-beta)
  delta <- 1
  
  first <- (mu + lambda*Moment_Distribution(Star_Distribution[1],Star_Distribution[2],Star_Distribution[3],1))*delta
  second <- (sigma+ lambda*Moment_Distribution(Star_Distribution[1],Star_Distribution[2],Star_Distribution[3],2))*delta+delta^2*beta*lambda*(2*alpha-beta)*Moment_Distribution(Star_Distribution[1],Star_Distribution[2],Star_Distribution[3],1)^2/(2*(alpha-beta))
  third <- lambda*Moment_Distribution(Star_Distribution[1],Star_Distribution[2],Star_Distribution[3],3)*delta+delta^2*3/2*(2*alpha-beta)*beta*lambda*Moment_Distribution(Star_Distribution[1],Star_Distribution[2],Star_Distribution[3],1)*Moment_Distribution(Star_Distribution[1],Star_Distribution[2],Star_Distribution[3],2)/(alpha-beta)
  fourth <- delta*lambda*Moment_Distribution(Star_Distribution[1],Star_Distribution[2],Star_Distribution[3],4)+3*delta^2*sigma^2+6*sigma*lambda*Moment_Distribution(Star_Distribution[1],Star_Distribution[2],Star_Distribution[3],2)+3*lambda*(lambda+(2*alpha-beta)*beta/(2*(alpha-beta)))*Moment_Distribution(Star_Distribution[1],Star_Distribution[2],Star_Distribution[3],2)^2+(2*(2*alpha-beta)*beta*lambda*Moment_Distribution(Star_Distribution[1],Star_Distribution[2],Star_Distribution[3],1)*Moment_Distribution(Star_Distribution[1],Star_Distribution[2],Star_Distribution[3],3))/(alpha-beta)
  autocov <- beta*lambda*(2*alpha-beta)/(2*(alpha-beta))*exp(-(alpha-beta)*10)*Moment_Distribution(Star_Distribution[1],Star_Distribution[2],Star_Distribution[3],1)^2*delta^2
  autocov2 <- beta*lambda*(2*alpha-beta)/(2*(alpha-beta))*exp(-(alpha-beta)*10)*Moment_Distribution(Star_Distribution[1],Star_Distribution[2],Star_Distribution[3],2)^2*delta^2
  return(c(first, second, third, fourth,autocov, autocov2))
}

Covariance <- function(delta_X, mean_delta_X, tau, r) {
  n <- length(delta_X)
  covariance <- 0
  tauincluded <- delta_X[tau:n]
  for (t in 1:(n-tau)) 
  {
    term1 <- (delta_X[t])^r - (mean_delta_X)^r
    term2 <- (delta_X[t + tau])^r  - (mean(tauincluded))^r
    covariance <- covariance + term1 * term2
    covariance <- covariance / (n - tau)
    return(covariance)
  }
}

GMM_Function <- function(parameters, data) {
  Sample_Moments <- c(mean(data),mean((data-mean(data))^2),mean((data-mean(data))^3),mean((data-mean(data))^4),Covariance(data,mean(data),10,1),Covariance(data,mean(data),5,2))
  Theoretical_Moments <- Theoretical_GMM(parameters)
  Diff <- Sample_Moments - Theoretical_Moments
  return(sum(Diff^2))
}

Initial_GMM <- c(drift, vol, Solution_lambda_inf, Solution_alpha, betaMCMC)
Optim_GMM <- optim(Initial_GMM, GMM_Function, data = logreturn, method = "BFGS")
Star_GMM <- Optim_GMM$par
print(Star_GMM)
Solution_alpha/betaMCMC
Star_GMM[4]/Star_GMM[5]

#Check for over-identification
Jstat <- N * GMM_Function(Star_GMM, logreturn)
#Follows a chi-squared distribution with df=1
p_val <- 1 - pchisq(Jstat, 1)
p_val
if (p_val < 0.05) 
{
  cat("Moment conditions may not be statistically compatible \n")
} else
{
  cat("Moment conditions are statistically compatible.\n")
}

#HJD GMM Simulation
Distribution_sample <- function(n,p,a,b)
{
  prob_pos <- rbinom(n, 1, p) == 1
  sample <- numeric(n)
  sample[prob_pos] <- rexp(sum(prob_pos), rate = a)
  sample[!prob_pos] <- -rexp(sum(!prob_pos), rate = b)
  return(sample)
}

HJD_path_GMM  <- function(S0, drift, vol, N, lambda_infty, alpha, beta, p, a, b) {
  HJDsim_GMM <- numeric(N)
  HJDsim_GMM[1] <- S0
  jumpindex <- unique(as.integer(simulateHawkes(lambda_infty, alpha, beta, 10000)[[1]]))
  jumplogic <- rep(0, N)
  jumplogic[jumpindex] <- 1
  jumpsize <- Distribution_sample(N, p, a, b)
  Wiener <- rnorm(N)
  
  for (t in 2:N) {
    dS <- drift+ sqrt(vol) * Wiener[t] + jumplogic[t] * jumpsize[t]
    HJDsim_GMM[t] <- HJDsim_GMM[t - 1] * exp(dS)
  }
  
  return(data.frame(Day = seq(1:N), Price = HJDsim_GMM, group = "HJD Simulation with GMM"))
}
HJDgraph_GMM <- list()
for(i in 1:6)
{
  HJDsim_GMM <- HJD_path_GMM(open[1], Star_GMM[1], Star_GMM[2], N, Star_GMM[3], Star_GMM[4], Star_GMM[5], Star_Distribution[1], Star_Distribution[2], Star_Distribution[3])
  HJDgraph_GMM[[i]] <- plot_comparison(rbind(HJDsim_GMM, Close_data))
}
legend <- get_legend(HJDgraph_GMM[[1]] + theme(legend.position = "bottom"))
grid.arrange(
  do.call(arrangeGrob, c(HJDgraph_GMM, nrow = 3)),
  legend,
  ncol = 1,
  heights = c(10, 1),
  top = textGrob("HJD Simulation Comparisons with GMM", gp = gpar(fontsize = 16, fontface = "bold"))
)

#Two sample Wald test
W <- (Star_EALD[1]-Star_Distribution[1])^2/(sd_p_MLE^2+sd_p_GMM^2)
crit_value <- qchisq(0.95, 1)
if (W > crit_value) {
  cat("The estimates are significantly different.\n")
} else {
  cat("The estimates are not significantly different.\n")
}

#Goodness of fit
EALD <- function(x, p, a, b) 
{
  ifelse(x < 0, (1 - p) * (1 - exp(b * x)), p * (1 - exp(-a * x)) + (1 - p))
}
KS_MLE <- ks.test(jump, function(x) EALD(x, Star_Distribution[1], Star_Distribution[2], Star_Distribution[3]))
KS_GMM <- ks.test(jump, function(x) EALD(x, Star_EALD[1], Star_EALD[2], Star_EALD[3]))
print(KS_MLE)
print(KS_GMM)

# Testing for Hawkes
Integration_Hawkes <- function(times, mu, alpha, beta)
{
  cum_intensity <- numeric(N)
  lambda <- numeric(N)
  for (i in 1:N) {
    valid <- detected_jumps[detected_jumps < i]
    if (length(valid) == 0) {
      lambda[i] <- mu
    } else {
      lambda[i] <- mu + sum(sapply(valid, function(ti) alpha * exp(-beta * (i - ti))))
    }
  }
  cum_intensity[1] <- lambda[1]
  for (t in 2:N)
  {
    cum_intensity[t] <- cum_intensity[t-1]+lambda[t]
  }
  return(cum_intensity)
}
#MLE
transformed_times_MLE <- Integration_Hawkes(detected_jumps, Solution_lambda_inf, Solution_alpha, betaMCMC)
Hawkes_Poisson_MLE <- transformed_times_MLE[detected_jumps]
intervals_MLE <- diff(Hawkes_Poisson_MLE)

qqplot(qexp(ppoints(length(intervals_MLE)), rate = 1), intervals_MLE, main = "Q-Q Plot of the interval vs Exponential with MLE")
abline(0, 1, col = "red")

ks.test(intervals_MLE, "pexp", rate = 1)

#GMM
transformed_times_GMM <- Integration_Hawkes(detected_jumps, Star_GMM[3], Star_GMM[4], Star_GMM[5])
Hawkes_Poisson_GMM <- transformed_times_GMM[detected_jumps]
intervals_GMM <- diff(Hawkes_Poisson_GMM)

qqplot(qexp(ppoints(length(intervals_GMM)), rate = 1), intervals_GMM, main = "Q-Q Plot of the interval vs Exponential with GMM")
abline(0, 1, col = "red")

ks.test(intervals_GMM, "pexp", rate = 1)

#Truncated MLE
truncated <- function(x, rate, q) {
  qexp_cdf <- pexp(q, rate)
  return(pexp(x, rate) / qexp_cdf)
}
results_MLE <- data.frame(Quantile = numeric(0), Size = numeric(0), Statistic = numeric(0), p_Value = numeric(0))
subintervals_MLE <- list()
options(repr.plot.width = 15, repr.plot.height = 6)
par(mfrow = c(2, 5), mar = c(4, 4, 2, 1))
for (i in 1:10)
{
  cutoff_MLE <- quantile(intervals_MLE, i / 10)
  subintervals_MLE[[i]] <- intervals_MLE[intervals_MLE <= cutoff_MLE]
  qqtrunc(subintervals_MLE[[i]], "exp", rate = 1, a = 0, b = cutoff_MLE, title = paste("MLE: D",i), xlabel = "Theoretical", ylabel = "Sample")
  abline(0, 1, col = "red")
  ks_result_MLE <- ks.test(subintervals_MLE[[i]], truncated, rate = 1, q = cutoff_MLE)
  results_MLE <- rbind(results_MLE, data.frame(Quantile = i / 10, Size = length(subintervals_MLE[[i]]), Statistic = ks_result_MLE$statistic, p_Value = ks_result_MLE$p.value))
}
par(mfrow = c(1, 1))
print(results_MLE)

#Truncated GMM
results_GMM <- data.frame(Quantile = numeric(0), Size = numeric(0), Statistic = numeric(0), p_Value = numeric(0))
subintervals_GMM <- list()
options(repr.plot.width = 15, repr.plot.height = 6)
par(mfrow = c(2, 5), mar = c(4, 4, 2, 1))
for (i in 1:10)
{
  cutoff_GMM <- quantile(intervals_GMM, i / 10)
  subintervals_GMM[[i]] <- intervals_GMM[intervals_GMM <= cutoff_GMM]
  qqtrunc(subintervals_GMM[[i]], "exp", rate = 1, a = 0, b = cutoff_GMM, title = paste("GMM: D",i), xlabel = "Theoretical", ylabel = "Sample")
  abline(0, 1, col = "red")
  ks_result_GMM <- ks.test(subintervals_GMM[[i]], truncated, rate = 1, q = cutoff_GMM)
  results_GMM <- rbind(results_GMM, data.frame(Quantile = i / 10, Size = length(subintervals_GMM[[i]]), Statistic = ks_result_GMM$statistic, p_Value = ks_result_GMM$p.value))
}
print(results_GMM)
par(mfrow = c(1, 1))

#MAPE
MAPE <- function(actual, simulation)
{
  return(sum(abs((actual-simulation)/actual)))
}
MAPEGBM <- vector()
MAPEMJD <- vector()
MAPEHJD_MLEMCMC <- vector()
MAPEHJD_GMM <- vector()
pb <- progress_bar$new(total = n, format = "[:bar] :percent Finishes in: :eta", clear = TRUE)
for (i in 1:10000)
{
  pb$tick()
  MAPEGBM[i] <- MAPE(close, GBM_path(open[1], drift, vol, N)$Price)
  MAPEMJD[i] <- MAPE(close,MJD_path(open[1], drift, vol, N, frequency, jumpexp, jumpsd)$Price)
  MAPEHJD_MLEMCMC[i] <- MAPE(close, HJD_path_MLEMCMC(open[1], drift, vol, N, Solution_lambda_inf, Solution_alpha, betaMCMC, Star_EALD[1], Star_EALD[2], Star_EALD[3])$Price)
  MAPEHJD_GMM[i] <- MAPE(close, HJD_path_GMM(open[1], Star_GMM[1], Star_GMM[2], N, Star_GMM[3], Star_GMM[4], Star_GMM[5], Star_Distribution[1], Star_Distribution[2], Star_Distribution[3])$Price)
  Sys.sleep(1 / 10000)
}
MCMAPEGBM <- mean(MAPEGBM)
MCMAPEMJD <- mean(MAPEMJD)
MCMAPEHJD_MLEMCMC <- mean(MAPEHJD_MLEMCMC)
MCMAPEHJD_GMM <- mean(MAPEHJD_GMM)
MCMAPEGBM
MCMAPEMJD
MCMAPEHJD_MLEMCMC
MCMAPEHJD_GMM
min(MAPEGBM)
min(MAPEMJD)
min(MAPEHJD_MLEMCMC)
min(MAPEHJD_GMM)
