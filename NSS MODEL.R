# Define the maturities and yields
maturities <- c(0.12, 0.45, 0.66, 0.71, 0.96, 1.65, 2.09, 2.84, 3.61, 3.99, 4.73, 6, 7, 8, 8.62, 8.99, 9.99, 14.08, 17.62, 29.99)
yields <- c(3.364, 3.233, 3.327, 3.217, 3.259, 2.844, 2.775, 2.882, 2.972, 3.109, 3.157, 3.278, 3.334, 3.410, 3.507, 3.504, 3.636, 3.895, 4.003, 4.185)

# Define the NSS model residuals function
nss_residuals <- function(parameters, maturities, yields) {
  beta0 <- parameters[1]
  beta1 <- parameters[2]
  beta2 <- parameters[3]
  beta3 <- parameters[4]
  tau1 <- parameters[5]
  tau2 <- parameters[6]
  
  fitted_yields <- beta0 +
    beta1 * (1 - exp(-maturities/tau1)) / (maturities/tau1) +
    beta2 * ((1 - exp(-maturities/tau1)) / (maturities/tau1) - exp(-maturities/tau1)) +
    beta3 * ((1 - exp(-maturities/tau2)) / (maturities/tau2) - exp(-maturities/tau2))
  
  residuals <- fitted_yields - yields
  return(sum(residuals^2))
}

# Initial guess for parameters
initial_parameters <- c(0.02, -0.02, -0.02, 0.02, 1, 5)

# Optimize parameters
nss_fit <- optim(initial_parameters, nss_residuals, maturities = maturities, yields = yields)

# Extract fitted parameters
fitted_parameters <- nss_fit$par
beta0 <- fitted_parameters[1]
beta1 <- fitted_parameters[2]
beta2 <- fitted_parameters[3]
beta3 <- fitted_parameters[4]
tau1 <- fitted_parameters[5]
tau2 <- fitted_parameters[6]

# Define function to calculate NSS yields
nss_yield <- function(maturities, beta0, beta1, beta2, beta3, tau1, tau2) {
  term1 <- beta0
  term2 <- beta1 * (1 - exp(-maturities/tau1)) / (maturities/tau1)
  term3 <- beta2 * ((1 - exp(-maturities/tau1)) / (maturities/tau1) - exp(-maturities/tau1))
  term4 <- beta3 * ((1 - exp(-maturities/tau2)) / (maturities/tau2) - exp(-maturities/tau2))
  return(term1 + term2 + term3 + term4)
}

# Calculate fitted yields for plotting
fitted_yields <- nss_yield(maturities, beta0, beta1, beta2, beta3, tau1, tau2)

# Plotting observed yields and fitted NSS curve
plot(maturities, yields, main = "Nelson-Siegel-Svensson Yield Curve", xlab = "Maturity", ylab = "Yield", col = "blue", pch = 16)
lines(maturities, fitted_yields, col = "red")
legend("topright", legend = c("Observed Yields", "Fitted Nelson-Siegel-Svensson Curve"), col = c("blue", "red"), pch = c(16, NA), lty = c(NA, 1))
