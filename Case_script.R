library(bootUR)
library(vars)
library(lpirfs)
library(readr) 
library(dplyr) 
library(lubridate)
fred_md <- read.csv("data/data.csv")
                    
case_data <- fred_md[, c("sasdate", "INDPRO", "CPIAUCSL", "FEDFUNDS")]

case_data <- case_data[-1, ]

# Convert date to Date format
case_data$sasdate <- as.Date(case_data$sasdate, format = "%m/%d/%Y")

# End sample in December 2019 to avoid the COVID-19 period
case_data <- subset(case_data, sasdate <= as.Date("2019-12-01"))

range(case_data$sasdate)

# Inspect the data

head(case_data)
str(case_data)
summary(case_data)

# Time series

INDPRO_ts <- ts(case_data$INDPRO, start = c(1959, 1), frequency = 12)
CPI_ts <- ts(case_data$CPIAUCSL, start = c(1959, 1), frequency = 12)
FEDFUNDS_ts <- ts(case_data$FEDFUNDS, start = c(1959, 1), frequency = 12)

# Plotting time series without transformation
plot(
  INDPRO_ts,
  col = "blue",
  main = "US Industrial Production",
  xlab = "Year",
)

plot(
  CPI_ts,
  col = "green",
  main = "US Consumer Price Index",
  xlab = "Year",
)

plot(
  FEDFUNDS_ts,
  col = "red",
  main = "Effective Federal Funds Rate",
  xlab = "Year",
)

INDPRO_trans <- diff(log(INDPRO_ts))
CPI_trans <- diff(log(CPI_ts), differences = 2)
FEDFUNDS_trans <- diff(FEDFUNDS_ts)


#Plotting transformed time series:
plot(
  INDPRO_trans,
  col = "blue",
  main = "US Industrial Production",
  xlab = "Year",
)

plot(
  CPI_trans,
  col = "green",
  main = "US Consumer Price Index",
  xlab = "Year",
)

plot(
  FEDFUNDS_trans,
  col = "red",
  main = "Effective Federal Funds Rate",
  xlab = "Year",
)

#log

log_INDPRO <- log(INDPRO_ts)
log_CPI <- log(CPI_ts)

# Unit root tests:

#INDPRO
d_log_INDPRO <- diff(log_INDPRO)
d2_log_INDPRO <- diff(log_INDPRO, differences = 2)

print(adf(d2_log_INDPRO, deterministics = "intercept"))
print(adf(d_log_INDPRO, deterministics = "intercept"))
print(adf(log_INDPRO, deterministics = "trend"))

#CPIAUCSL
d_log_CPI <- diff(log_CPI)
d2_log_CPI <- diff(log_CPI, differences = 2)

print(adf(d2_log_CPI, deterministics = "intercept"))
print(adf(d_log_CPI, deterministics = "intercept"))
print(adf(log_CPI, deterministics = "trend"))

#FEDFUNDS

d_FEDFUNDS  <- diff(FEDFUNDS_ts)
d2_FEDFUNDS <- diff(FEDFUNDS_ts, differences = 2)

print(adf(d2_FEDFUNDS, deterministics = "intercept"))
print(adf(d_FEDFUNDS,  deterministics = "intercept"))
print(adf(FEDFUNDS_ts, deterministics = "intercept"))

# Create the final stationary variables for the VAR

INDPRO_VAR <- diff(log(INDPRO_ts))
CPI_VAR <- diff(log(CPI_ts), differences = 2)
FEDFUNDS_VAR <- diff(FEDFUNDS_ts)

# Put the three stationary variables together

VAR_data <- cbind(
  INDPRO = INDPRO_VAR,
  CPI = CPI_VAR,
  FEDFUNDS = FEDFUNDS_VAR
)

head(VAR_data)

dim(VAR_data)


# Check the time-series structure.
VAR_data
VAR_data <- na.omit(VAR_data)

# SELECT THE VAR LAG ORDER

lag_selection <- VARselect(
  VAR_data,
  lag.max = 12,
  type = "const"
)

# Lag order selected by each information criterion:
lag_selection$selection
# Result: AIC = 10, HQ = 4, SC/BIC = 3 -> the criteria disagree.

lag_selection$criteria

# SEQUENTIAL TEST FOR THE 10TH LAG

# Unrestricted model: VAR(10)
VAR_10 <- VAR(
  VAR_data,
  p = 10,
  type = "const"
)

VAR_9 <- VAR(
  VAR_data[-1, ],
  p = 9,
  type = "const"
)

K <- ncol(VAR_data)
T <- VAR_10$obs

Sigma10 <- crossprod(residuals(VAR_10)) / T
Sigma9  <- crossprod(residuals(VAR_9)) / T

# Likelihood-ratio statistic:

LR_10 <- T * (
  log(det(Sigma9)) - log(det(Sigma10))
)

q <- K^2

p_value_10 <- 1 - pchisq(
  LR_10,
  df = q
)

LR_10        # 19.58
p_value_10   # 0.021 < 0.05 -> reject H0: the 10th lag is significant

# So the sequential test also points to VAR(10).

# FINAL CHOICE OF THE LAG ORDER: VAR(3)

K * (K * c(3, 4, 10) + 1)
# VAR(3) = 30, VAR(4) = 39, VAR(10) = 93 parameters.

p_final <- 3   # selected lag order, used in all sections below

VAR_3 <- VAR(
  VAR_data,
  p = p_final,
  type = "const"
)

summary(VAR_3)

# RESIDUAL DIAGNOSTICS

serial.test(
  VAR_3,
  lags.pt = 12,
  type = "PT.asymptotic"
)

# Test for multivariate normality of the residuals.
# H0: residuals are multivariate normally distributed.

normality.test(VAR_3)

# Test for ARCH effects.

arch.test(VAR_3)

# CHECK STABILITY OF THE VAR(3)


roots(VAR_3)   # largest root = 0.665 < 1 -> stable

# GRANGER CAUSALITY TESTS

# INDPRO AND CPI 

VAR_INDPRO_CPI <- VAR(
  VAR_data[, c("INDPRO", "CPI")],
  p = p_final,
  type = "const"
)

GC_INDPRO_CPI <- causality(
  VAR_INDPRO_CPI,
  cause = "INDPRO"
)

GC_INDPRO_CPI$Granger

GC_CPI_INDPRO <- causality(
  VAR_INDPRO_CPI,
  cause = "CPI"
)

GC_CPI_INDPRO$Granger

# INDPRO AND FEDFUNDS

VAR_INDPRO_FED <- VAR(
  VAR_data[, c("INDPRO", "FEDFUNDS")],
  p = p_final,
  type = "const"
)

GC_INDPRO_FED <- causality(
  VAR_INDPRO_FED,
  cause = "INDPRO"
)

GC_INDPRO_FED$Granger

GC_FED_INDPRO <- causality(
  VAR_INDPRO_FED,
  cause = "FEDFUNDS"
)

GC_FED_INDPRO$Granger

# CPI AND FEDFUNDS

VAR_CPI_FED <- VAR(
  VAR_data[, c("CPI", "FEDFUNDS")],
  p = p_final,
  type = "const"
)

GC_CPI_FED <- causality(
  VAR_CPI_FED,
  cause = "CPI"
)

GC_CPI_FED$Granger

GC_FED_CPI <- causality(
  VAR_CPI_FED,
  cause = "FEDFUNDS"
)

GC_FED_CPI$Granger

#Impulse Response Functions for reduced-form VARs:

irf_points <- irf(
  VAR_3,
  n.ahead = 24,
  ortho = FALSE,
  boot = FALSE
)

par(mfrow = c(3, 3), mar = c(4, 4, 2.5, 1))

#SHIOCK INDPRO -> INDPRO
y_val <- irf_points$irf$INDPRO[, "INDPRO"]
plot(0:24, y_val, type = "l", col = "blue", lwd = 2,
     ylim = range(c(0, y_val)),
     main = "INDPRO → INDPRO", xlab = "Months", ylab = "Response")
abline(h = 0, col = "gray40", lty = 2)

#SHOCK INDPRO -> CPI
y_val <- irf_points$irf$INDPRO[, "CPI"]
plot(0:24, y_val, type = "l", col = "blue", lwd = 2,
     ylim = range(c(0, y_val)),
     main = "INDPRO → CPI", xlab = "Months", ylab = "Response")
abline(h = 0, col = "gray40", lty = 2)

#SHOCK INDPRO -> FEDFUNDS
y_val <- irf_points$irf$INDPRO[, "FEDFUNDS"]
plot(0:24, y_val, type = "l", col = "blue", lwd = 2,
     ylim = range(c(0, y_val)),
     main = "INDPRO → FEDFUNDS", xlab = "Months", ylab = "Response")
abline(h = 0, col = "gray40", lty = 2)

#SHOCK CPI -> INDPRO
y_val <- irf_points$irf$CPI[, "INDPRO"]
plot(0:24, y_val, type = "l", col = "blue", lwd = 2,
     ylim = range(c(0, y_val)),
     main = "CPI → INDPRO", xlab = "Months", ylab = "Response")
abline(h = 0, col = "gray40", lty = 2)

#SHOCK CPI -> CPI
y_val <- irf_points$irf$CPI[, "CPI"]
plot(0:24, y_val, type = "l", col = "blue", lwd = 2,
     ylim = range(c(0, y_val)),
     main = "CPI → CPI", xlab = "Months", ylab = "Response")
abline(h = 0, col = "gray40", lty = 2)

#SHOCK CPI -> FEDFUNDS
y_val <- irf_points$irf$CPI[, "FEDFUNDS"]
plot(0:24, y_val, type = "l", col = "blue", lwd = 2,
     ylim = range(c(0, y_val)),
     main = "CPI → FEDFUNDS", xlab = "Months", ylab = "Response")
abline(h = 0, col = "gray40", lty = 2)

#SHOCK FEDFUNDS -> INDPRO
y_val <- irf_points$irf$FEDFUNDS[, "INDPRO"]
plot(0:24, y_val, type = "l", col = "blue", lwd = 2,
     ylim = range(c(0, y_val)),
     main = "FEDFUNDS → INDPRO", xlab = "Months", ylab = "Response")
abline(h = 0, col = "gray40", lty = 2)

#SHOCK FEDFUNDS -> CPI
y_val <- irf_points$irf$FEDFUNDS[, "CPI"]
plot(0:24, y_val, type = "l", col = "blue", lwd = 2,
     ylim = range(c(0, y_val)),
     main = "FEDFUNDS → CPI", xlab = "Months", ylab = "Response")
abline(h = 0, col = "gray40", lty = 2)

#SHOCK FEDFUNDS -> FEDFUNDS
y_val <- irf_points$irf$FEDFUNDS[, "FEDFUNDS"]
plot(0:24, y_val, type = "l", col = "blue", lwd = 2,
     ylim = range(c(0, y_val)),
     main = "FEDFUNDS → FEDFUNDS", xlab = "Months", ylab = "Response")
abline(h = 0, col = "gray40", lty = 2)

#Price puzzle, transforming the ts to level:

cpi_response <- irf_points$irf$FEDFUNDS[, "CPI"]

cpi_level_response <- 100 * cumsum(cumsum(cpi_response))

plot(
  0:24, cpi_level_response,
  type = "l",
  xlab = "Months after FEDFUNDS innovation",
  ylab = "Approximate CPI level response (%)",
  main = "CPI level response to a FEDFUNDS innovation"
)
abline(h = 0, col = "red")

#confidence intervals:
irf_ci <- irf(
  VAR_3,
  n.ahead = 24,
  ortho = FALSE,
  boot = TRUE,
  ci = 0.95,
  runs = 1000,
  seed = 123
)

#Plotting ts with confidence intervals this time:
par(mfrow = c(3, 3), mar = c(4, 4, 2.5, 1))

#Shock INDPRO -> response INDPRO
plot(0:24, irf_ci$irf$INDPRO[, "INDPRO"], type = "l",
     ylim = range(irf_ci$Lower$INDPRO[, "INDPRO"],
                  irf_ci$Upper$INDPRO[, "INDPRO"]),
     main = "INDPRO → INDPRO", xlab = "Months", ylab = "Response")
lines(0:24, irf_ci$Lower$INDPRO[, "INDPRO"], col = "red", lty = 2)
lines(0:24, irf_ci$Upper$INDPRO[, "INDPRO"], col = "red", lty = 2)
abline(h = 0)

#Shock INDPRO -> response CPI
plot(0:24, irf_ci$irf$INDPRO[, "CPI"], type = "l",
     ylim = range(irf_ci$Lower$INDPRO[, "CPI"],
                  irf_ci$Upper$INDPRO[, "CPI"]),
     main = "INDPRO → CPI", xlab = "Months", ylab = "Response")
lines(0:24, irf_ci$Lower$INDPRO[, "CPI"], col = "red", lty = 2)
lines(0:24, irf_ci$Upper$INDPRO[, "CPI"], col = "red", lty = 2)
abline(h = 0)

#Shock INDPRO -> response FEDFUNDS
plot(0:24, irf_ci$irf$INDPRO[, "FEDFUNDS"], type = "l",
     ylim = range(irf_ci$Lower$INDPRO[, "FEDFUNDS"],
                  irf_ci$Upper$INDPRO[, "FEDFUNDS"]),
     main = "INDPRO → FEDFUNDS", xlab = "Months", ylab = "Response")
lines(0:24, irf_ci$Lower$INDPRO[, "FEDFUNDS"], col = "red", lty = 2)
lines(0:24, irf_ci$Upper$INDPRO[, "FEDFUNDS"], col = "red", lty = 2)
abline(h = 0)

#Shock CPI -> response INDPRO
plot(0:24, irf_ci$irf$CPI[, "INDPRO"], type = "l",
     ylim = range(irf_ci$Lower$CPI[, "INDPRO"],
                  irf_ci$Upper$CPI[, "INDPRO"]),
     main = "CPI → INDPRO", xlab = "Months", ylab = "Response")
lines(0:24, irf_ci$Lower$CPI[, "INDPRO"], col = "red", lty = 2)
lines(0:24, irf_ci$Upper$CPI[, "INDPRO"], col = "red", lty = 2)
abline(h = 0)

#Shock CPI -> response CPI
plot(0:24, irf_ci$irf$CPI[, "CPI"], type = "l",
     ylim = range(irf_ci$Lower$CPI[, "CPI"],
                  irf_ci$Upper$CPI[, "CPI"]),
     main = "CPI → CPI", xlab = "Months", ylab = "Response")
lines(0:24, irf_ci$Lower$CPI[, "CPI"], col = "red", lty = 2)
lines(0:24, irf_ci$Upper$CPI[, "CPI"], col = "red", lty = 2)
abline(h = 0)

#Shock CPI -> response FEDFUNDS
plot(0:24, irf_ci$irf$CPI[, "FEDFUNDS"], type = "l",
     ylim = range(irf_ci$Lower$CPI[, "FEDFUNDS"],
                  irf_ci$Upper$CPI[, "FEDFUNDS"]),
     main = "CPI → FEDFUNDS", xlab = "Months", ylab = "Response")
lines(0:24, irf_ci$Lower$CPI[, "FEDFUNDS"], col = "red", lty = 2)
lines(0:24, irf_ci$Upper$CPI[, "FEDFUNDS"], col = "red", lty = 2)
abline(h = 0)

#Shock FEDFUNDS -> response INDPRO
plot(0:24, irf_ci$irf$FEDFUNDS[, "INDPRO"], type = "l",
     ylim = range(irf_ci$Lower$FEDFUNDS[, "INDPRO"],
                  irf_ci$Upper$FEDFUNDS[, "INDPRO"]),
     main = "FEDFUNDS → INDPRO", xlab = "Months", ylab = "Response")
lines(0:24, irf_ci$Lower$FEDFUNDS[, "INDPRO"], col = "red", lty = 2)
lines(0:24, irf_ci$Upper$FEDFUNDS[, "INDPRO"], col = "red", lty = 2)
abline(h = 0)

#Shock FEDFUNDS -> response CPI
plot(0:24, irf_ci$irf$FEDFUNDS[, "CPI"], type = "l",
     ylim = range(irf_ci$Lower$FEDFUNDS[, "CPI"],
                  irf_ci$Upper$FEDFUNDS[, "CPI"]),
     main = "FEDFUNDS → CPI", xlab = "Months", ylab = "Response")
lines(0:24, irf_ci$Lower$FEDFUNDS[, "CPI"], col = "red", lty = 2)
lines(0:24, irf_ci$Upper$FEDFUNDS[, "CPI"], col = "red", lty = 2)
abline(h = 0)

#Shock FEDFUNDS -> response FEDFUNDS
plot(0:24, irf_ci$irf$FEDFUNDS[, "FEDFUNDS"], type = "l",
     ylim = range(irf_ci$Lower$FEDFUNDS[, "FEDFUNDS"],
                  irf_ci$Upper$FEDFUNDS[, "FEDFUNDS"]),
     main = "FEDFUNDS → FEDFUNDS", xlab = "Months", ylab = "Response")
lines(0:24, irf_ci$Lower$FEDFUNDS[, "FEDFUNDS"], col = "red", lty = 2)
lines(0:24, irf_ci$Upper$FEDFUNDS[, "FEDFUNDS"], col = "red", lty = 2)
abline(h = 0)

#Returning to the original variables in levels

to_levels <- function(irf_matrix) {
  result <- irf_matrix
  result[, "INDPRO"] <- 100 * (exp(cumsum(irf_matrix[, "INDPRO"])) - 1)
  result[, "CPI"] <- 100 * (exp(cumsum(cumsum(irf_matrix[, "CPI"]))) - 1)
  result[, "FEDFUNDS"] <- cumsum(irf_matrix[, "FEDFUNDS"])
  result
}

irf_levels <- lapply(irf_points$irf, to_levels)

h <- 0:24

par(mfrow = c(3, 3), mar = c(4, 4, 2.5, 1))

plot(h, irf_levels$FEDFUNDS[, "INDPRO"],
     type = "l", xlab = "Months", ylab = "%", main = "INDPRO level")
abline(h = 0, col = "red")

plot(h, irf_levels$FEDFUNDS[, "CPI"],
     type = "l", xlab = "Months", ylab = "%", main = "CPI level")
abline(h = 0, col = "red")

plot(h, irf_levels$FEDFUNDS[, "FEDFUNDS"],
     type = "l", xlab = "Months", ylab = "Percentage points",
     main = "FEDFUNDS level")
abline(h = 0, col = "red")

h <- 0:24

plot(h, irf_levels$INDPRO[, "INDPRO"],
     type = "l", xlab = "Months", ylab = "%", main = "INDPRO level")
abline(h = 0, col = "red")

plot(h, irf_levels$INDPRO[, "CPI"],
     type = "l", xlab = "Months", ylab = "%", main = "CPI level")
abline(h = 0, col = "red")

plot(h, irf_levels$INDPRO[, "FEDFUNDS"],
     type = "l", xlab = "Months", ylab = "Percentage points",
     main = "FEDFUNDS level")
abline(h = 0, col = "red")

h <- 0:24

plot(h, irf_levels$CPI[, "INDPRO"],
     type = "l", xlab = "Months", ylab = "%", main = "INDPRO level")
abline(h = 0, col = "red")

plot(h, irf_levels$CPI[, "CPI"],
     type = "l", xlab = "Months", ylab = "%", main = "CPI level")
abline(h = 0, col = "red")

plot(h, irf_levels$CPI[, "FEDFUNDS"],
     type = "l", xlab = "Months", ylab = "Percentage points",
     main = "FEDFUNDS level")
abline(h = 0, col = "red")

#Reduced-form VARs versus Structural VARs. 

colnames(VAR_3$y)

irf_recursive <- irf(
  VAR_3,
  n.ahead = 24,
  ortho = TRUE,
  boot = TRUE,
  ci = 0.95,
  runs = 1000,
  seed = 123
)

plot(0:24, irf_recursive$irf$FEDFUNDS[, "INDPRO"],
     type = "l", main = "FEDFUNDS → INDPRO", xlab = "Months", ylab = "Response")
lines(0:24, irf_recursive$Lower$FEDFUNDS[, "INDPRO"], col = "red", lty = 2)
lines(0:24, irf_recursive$Upper$FEDFUNDS[, "INDPRO"], col = "red", lty = 2)
abline(h = 0)
plot(0:24, irf_recursive$irf$FEDFUNDS[, "CPI"],
     type = "l", main = "FEDFUNDS → CPI", xlab = "Months", ylab = "Response")
lines(0:24, irf_recursive$Lower$FEDFUNDS[, "CPI"], col = "red", lty = 2)
lines(0:24, irf_recursive$Upper$FEDFUNDS[, "CPI"], col = "red", lty = 2)
abline(h = 0)
plot(0:24, irf_recursive$irf$FEDFUNDS[, "FEDFUNDS"],
     type = "l", main = "FEDFUNDS → FEDFUNDS", xlab = "Months", ylab = "Response")
lines(0:24, irf_recursive$Lower$FEDFUNDS[, "FEDFUNDS"], col = "red", lty = 2)
lines(0:24, irf_recursive$Upper$FEDFUNDS[, "FEDFUNDS"], col = "red", lty = 2)
abline(h = 0)

plot(0:24, irf_recursive$irf$INDPRO[, "INDPRO"],
     type = "l", main = "INDPRO → INDPRO", xlab = "Months", ylab = "Response")
lines(0:24, irf_recursive$Lower$INDPRO[, "INDPRO"], col = "red", lty = 2)
lines(0:24, irf_recursive$Upper$INDPRO[, "INDPRO"], col = "red", lty = 2)
abline(h = 0)
plot(0:24, irf_recursive$irf$INDPRO[, "CPI"],
     type = "l", main = "INDPRO → CPI", xlab = "Months", ylab = "Response")
lines(0:24, irf_recursive$Lower$INDPRO[, "CPI"], col = "red", lty = 2)
lines(0:24, irf_recursive$Upper$INDPRO[, "CPI"], col = "red", lty = 2)
abline(h = 0)
plot(0:24, irf_recursive$irf$INDPRO[, "FEDFUNDS"],
     type = "l", main = "INDPRO → FEDFUNDS", xlab = "Months", ylab = "Response")
lines(0:24, irf_recursive$Lower$INDPRO[, "FEDFUNDS"], col = "red", lty = 2)
lines(0:24, irf_recursive$Upper$INDPRO[, "FEDFUNDS"], col = "red", lty = 2)
abline(h = 0)

plot(0:24, irf_recursive$irf$CPI[, "INDPRO"],
     type = "l", main = "CPI → INDPRO", xlab = "Months", ylab = "Response")
lines(0:24, irf_recursive$Lower$CPI[, "INDPRO"], col = "red", lty = 2)
lines(0:24, irf_recursive$Upper$CPI[, "INDPRO"], col = "red", lty = 2)
abline(h = 0)
plot(0:24, irf_recursive$irf$CPI[, "CPI"],
     type = "l", main = "CPI → CPI", xlab = "Months", ylab = "Response")
lines(0:24, irf_recursive$Lower$CPI[, "CPI"], col = "red", lty = 2)
lines(0:24, irf_recursive$Upper$CPI[, "CPI"], col = "red", lty = 2)
abline(h = 0)
plot(0:24, irf_recursive$irf$CPI[, "FEDFUNDS"],
     type = "l", main = "CPI → FEDFUNDS", xlab = "Months", ylab = "Response")
lines(0:24, irf_recursive$Lower$CPI[, "FEDFUNDS"], col = "red", lty = 2)
lines(0:24, irf_recursive$Upper$CPI[, "FEDFUNDS"], col = "red", lty = 2)
abline(h = 0)

#Different order:

var_data_other <- VAR_data[, c("FEDFUNDS", "INDPRO", "CPI")]

var_model_other <- VAR(var_data_other, p = p_final, type = "const")

irf_other <- irf(
  var_model_other,
  n.ahead = 24,
  ortho = TRUE,
  boot = TRUE,
  ci = 0.95,
  runs = 1000,
  seed = 123
)

plot(
  0:24, irf_other$irf$FEDFUNDS[, "INDPRO"],
  type = "l",
  ylim = range(
    irf_other$Lower$FEDFUNDS[, "INDPRO"],
    irf_other$Upper$FEDFUNDS[, "INDPRO"]
  ),
  main = "FEDFUNDS → INDPRO (FEDFUNDS first)",
  xlab = "Months", ylab = "Response"
)
lines(0:24, irf_other$Lower$FEDFUNDS[, "INDPRO"], col = "red", lty = 2)
lines(0:24, irf_other$Upper$FEDFUNDS[, "INDPRO"], col = "red", lty = 2)
abline(h = 0)

plot(
  0:24, irf_other$irf$FEDFUNDS[, "CPI"],
  type = "l",
  ylim = range(
    irf_other$Lower$FEDFUNDS[, "CPI"],
    irf_other$Upper$FEDFUNDS[, "CPI"]
  ),
  main = "FEDFUNDS → CPI (FEDFUNDS first)",
  xlab = "Months", ylab = "Response"
)
lines(0:24, irf_other$Lower$FEDFUNDS[, "CPI"], col = "red", lty = 2)
lines(0:24, irf_other$Upper$FEDFUNDS[, "CPI"], col = "red", lty = 2)
abline(h = 0)

plot(
  0:24, irf_other$irf$FEDFUNDS[, "FEDFUNDS"],
  type = "l",
  ylim = range(
    irf_other$Lower$FEDFUNDS[, "FEDFUNDS"],
    irf_other$Upper$FEDFUNDS[, "FEDFUNDS"]
  ),
  main = "FEDFUNDS → FEDFUNDS (FEDFUNDS first)",
  xlab = "Months", ylab = "Response"
)
lines(0:24, irf_other$Lower$FEDFUNDS[, "FEDFUNDS"], col = "red", lty = 2)
lines(0:24, irf_other$Upper$FEDFUNDS[, "FEDFUNDS"], col = "red", lty = 2)
abline(h = 0)


plot(
  0:24, irf_other$irf$INDPRO[, "INDPRO"],
  type = "l",
  ylim = range(
    irf_other$Lower$INDPRO[, "INDPRO"],
    irf_other$Upper$INDPRO[, "INDPRO"]
  ),
  main = "INDPRO → INDPRO (FEDFUNDS first)",
  xlab = "Months", ylab = "Response"
)
lines(0:24, irf_other$Lower$INDPRO[, "INDPRO"], col = "red", lty = 2)
lines(0:24, irf_other$Upper$INDPRO[, "INDPRO"], col = "red", lty = 2)
abline(h = 0)

plot(
  0:24, irf_other$irf$INDPRO[, "CPI"],
  type = "l",
  ylim = range(
    irf_other$Lower$INDPRO[, "CPI"],
    irf_other$Upper$INDPRO[, "CPI"]
  ),
  main = "INDPRO → CPI (FEDFUNDS first)",
  xlab = "Months", ylab = "Response"
)
lines(0:24, irf_other$Lower$INDPRO[, "CPI"], col = "red", lty = 2)
lines(0:24, irf_other$Upper$INDPRO[, "CPI"], col = "red", lty = 2)
abline(h = 0)

plot(
  0:24, irf_other$irf$INDPRO[, "FEDFUNDS"],
  type = "l",
  ylim = range(
    irf_other$Lower$INDPRO[, "FEDFUNDS"],
    irf_other$Upper$INDPRO[, "FEDFUNDS"]
  ),
  main = "INDPRO → FEDFUNDS (FEDFUNDS first)",
  xlab = "Months", ylab = "Response"
)
lines(0:24, irf_other$Lower$INDPRO[, "FEDFUNDS"], col = "red", lty = 2)
lines(0:24, irf_other$Upper$INDPRO[, "FEDFUNDS"], col = "red", lty = 2)
abline(h = 0)


plot(
  0:24, irf_other$irf$CPI[, "INDPRO"],
  type = "l",
  ylim = range(
    irf_other$Lower$CPI[, "INDPRO"],
    irf_other$Upper$CPI[, "INDPRO"]
  ),
  main = "CPI → INDPRO (FEDFUNDS first)",
  xlab = "Months", ylab = "Response"
)
lines(0:24, irf_other$Lower$CPI[, "INDPRO"], col = "red", lty = 2)
lines(0:24, irf_other$Upper$CPI[, "INDPRO"], col = "red", lty = 2)
abline(h = 0)

plot(
  0:24, irf_other$irf$CPI[, "CPI"],
  type = "l",
  ylim = range(
    irf_other$Lower$CPI[, "CPI"],
    irf_other$Upper$CPI[, "CPI"]
  ),
  main = "CPI → CPI (FEDFUNDS first)",
  xlab = "Months", ylab = "Response"
)
lines(0:24, irf_other$Lower$CPI[, "CPI"], col = "red", lty = 2)
lines(0:24, irf_other$Upper$CPI[, "CPI"], col = "red", lty = 2)
abline(h = 0)

plot(
  0:24, irf_other$irf$CPI[, "FEDFUNDS"],
  type = "l",
  ylim = range(
    irf_other$Lower$CPI[, "FEDFUNDS"],
    irf_other$Upper$CPI[, "FEDFUNDS"]
  ),
  main = "CPI → FEDFUNDS (FEDFUNDS first)",
  xlab = "Months", ylab = "Response"
)
lines(0:24, irf_other$Lower$CPI[, "FEDFUNDS"], col = "red", lty = 2)
lines(0:24, irf_other$Upper$CPI[, "FEDFUNDS"], col = "red", lty = 2)
abline(h = 0)

#IRF VIA VARS VERSUS IRF VIA LP 
# LOCAL PROJECTIONS 
H <- 24 
p <- p_final 
variables <- c("INDPRO", "CPI", "FEDFUNDS")
LP_IRF <- function(response, shock, data, H = 24, p = 3){
  estimates <- lower <- upper <- numeric(H + 1) 
  for (h in 0:H) {
    y <- data[(p + 1 + h):nrow(data), response] 
    shock_t <- data[(p + 1):(nrow(data) - h), shock] 
    X <- data.frame(shock = shock_t) 
    for (v in variables) 
      for (lag in 1:p) X[[paste0(v, "_lag", lag)]] <- data[(p + 1 - lag):(nrow(data) - h - lag), v] 
    model <- lm(y ~ ., data = X) 
    estimates[h + 1] <- coef(model)["shock"] 
    se <- summary(model)$coefficients["shock", "Std. Error"] 
    lower[h + 1] <- estimates[h + 1] - 1.96 * se 
    upper[h + 1] <- estimates[h + 1] + 1.96 * se 
  }
  list(estimate = estimates, lower = lower, upper = upper) } 

# Estimate all 9 LPs 
LP_IRFs <- list() 
for (shock in variables) for (response in variables) 
  LP_IRFs[[shock]][[response]] <- LP_IRF(response, shock, VAR_data, H, p) 
# Plot all 9 
h <- 0:H 
par(mfrow = c(3, 3)) 
for (shock in variables) for (response in variables) { 
  r <- LP_IRFs[[shock]][[response]] 
  plot(h, r$estimate, type = "l", xlab = "Horizon (months)", ylab = "Response", main = paste(shock, "→", response), ylim = range(c(r$lower, r$upper))) 
  lines(h, r$lower, lty = 2) 
  lines(h, r$upper, lty = 2)
  abline(h = 0, lty = 2) 
} 

#STRUCTURAL IRFs USING LOCAL PROJECTIONS 
LP_struct <- lp_lin( 
  endog_data = data.frame(VAR_data[, c("INDPRO", "CPI", "FEDFUNDS")]), 
  lags_endog_lin = p_final, 
  # 3 lags of each variable as controls 
  trend = 0, # no time trend 
  shock_type = 1, # unit shock: rate rises by 1 percentage point 
  confint = 1.96, # 95% confidence bands 
  hor = 24 # trace the response for 24 months 
  ) 
plot(LP_struct)

#SHOCK PROXY: LP-IV
mps <- read_csv("monetary-policy-surprises/mps.csv") 
#make the FMOC event monthly 
mps_monthly <- mps %>% group_by(month = floor_date(Date, "month")) %>% summarise(ME = sum(ME, na.rm = TRUE)) 
#match the data time period, so we put 0 for all the month we do not have 
VARproxy <- data.frame( date = case_data$sasdate[3:nrow(case_data)], VAR_data ) 
VARproxy <- VARproxy %>% left_join(mps_monthly, by = c("date" = "month")) %>% mutate(ME = coalesce(ME, 0)) 
#now let's reduced the data to 1994 so it matches 
LPIVdata <- VARproxy %>% filter(date >= "1994-01-01") 
LPIV <- lp_lin_iv( endog_data = LPIVdata[, c("INDPRO", "CPI", "FEDFUNDS")], 
                   shock = LPIVdata["FEDFUNDS"], 
                   instrum = LPIVdata["ME"], 
                   use_twosls = TRUE, 
                   lags_endog_lin = p_final, 
                   trend = 0, confint = 1.96, hor = 24 ) 
plot(LPIV)