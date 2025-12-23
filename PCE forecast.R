# Time Series Forecasting: U.S. Personal Consumption Expenditures (PCE)

# loading libraries 

suppressPackageStartupMessages({
  library(VIM)
  library(tseries)
  library(imputeTS)
  library(tsoutliers)
  library(forecast)
  library(ggplot2)
  library(zoo)
  library(knitr)
})

# Loading data
data <- read.csv("PCE.csv", encoding = "UTF-8", stringsAsFactors = FALSE)

# #Understanding data
head(data) 
str(data)

# Parsing date and sorting
data$observation_date <- as.Date(data$observation_date)
data <- data[order(data$observation_date), ]

# Visualisation of missing values
aggr(data, numbers = TRUE, prop = FALSE, cex.axis = 0.6)

# #Creating time-series 
pce_ts <- ts(data$PCE, start = c(1959, 1), frequency = 12)
print(pce_ts)

# Imputation 
train_end <- c(2023, 12)
pce_train_raw <- window(pce_ts, end = train_end)

# RMSE function
rmse <- function(true, predicted) {
  sqrt(mean((true - predicted)^2, na.rm = TRUE))
}

# Create artificial missing values
set.seed(689)
TRAIN_RATIO <- 0.8

observed_idx <- which(!is.na(pce_train_raw))
train_size <- floor(TRAIN_RATIO * length(observed_idx))

keep_idx <- sample(observed_idx, size = train_size, replace = FALSE)

train_masked <- pce_train_raw
# Mask the points not kept
mask_idx <- setdiff(observed_idx, keep_idx)
train_masked[mask_idx] <- NA

# Impute the training series 
imp_linear <- na_interpolation(train_masked)
imp_spline <- na_interpolation(train_masked, option = "spline")
imp_ma     <- na_ma(train_masked, k = 12, weighting = "exponential")
imp_kalman <- na_kalman(train_masked, model = "auto.arima")

# Evaluate RMSE
error_results <- data.frame(
  Method = c("Interpolation", "Spline", "MovingAvg", "Kalman_ARIMA"),
  RMSE = c(
    rmse(pce_train_raw[mask_idx], imp_linear[mask_idx]),
    rmse(pce_train_raw[mask_idx], imp_spline[mask_idx]),
    rmse(pce_train_raw[mask_idx], imp_ma[mask_idx]),
    rmse(pce_train_raw[mask_idx], imp_kalman[mask_idx])
  )
)

print(error_results)

best_method <- error_results$Method[which.min(error_results$RMSE)]
cat("\nBest imputation method (train-only RMSE):", best_method, "\n\n")

# Visualize imputations on FULL series
pce_imp_linear_full <- na_interpolation(pce_ts)
pce_imp_spline_full <- na_interpolation(pce_ts, option = "spline")
pce_imp_ma_full     <- na_ma(pce_ts, k = 12, weighting = "exponential")
pce_imp_kalman_full <- na_kalman(pce_ts, model = "auto.arima")

par(mfrow = c(2, 2))
plot(pce_imp_linear_full, main = "Interpolation", xlab = "Year", ylab = "PCE")
plot(pce_imp_ma_full,     main = "Moving average", xlab = "Year", ylab = "PCE")
plot(pce_imp_kalman_full, main = "Kalman ARIMA", xlab = "Year", ylab = "PCE")
plot(pce_imp_spline_full, main = "Spline", xlab = "Year", ylab = "PCE")
par(mfrow = c(1, 1))

# Choose imputed full series based on best_method
pce_imputed <- switch(
  best_method,
  "Interpolation" = pce_imp_linear_full,
  "Spline"        = pce_imp_spline_full,
  "MovingAvg"     = pce_imp_ma_full,
  "Kalman_ARIMA"  = pce_imp_kalman_full,
  pce_imp_spline_full
)

# Log transformation + outlier correction
log_pce <- log(pce_imputed)

par(mfrow = c(1, 2))
plot(pce_imputed, main = "Imputed PCE", ylab = "PCE", xlab = "Year")
plot(log_pce, main = "Log(PCE)", ylab = "log(PCE)", xlab = "Year")
par(mfrow = c(1, 1))

# Detect and replace outliers on log scale
out_log <- tsoutliers(log_pce)
if (!is.null(out_log$index) && length(out_log$index) > 0) {
  log_pce[out_log$index] <- out_log$replacements
  cat("Outliers replaced at dates:\n")
  print(as.yearmon(time(log_pce)[out_log$index]))
}

plot(log_pce, main = "Cleaned log(PCE)", ylab = "log(PCE)", xlab = "Year")

# Decompose 
plot(stl(log_pce, s.window = "periodic"))

# Suggested differencing checks 
cat("ndiffs(log_pce):", ndiffs(log_pce), "\n")
cat("nsdiffs(log_pce):", nsdiffs(log_pce), "\n")

# Train-test split

log_train <- window(log_pce, end = c(2023, 12))
log_test  <- window(log_pce, start = c(2024, 1))

# Fit models on train

# Drift model (random walk with drift)
fit_drift <- rwf(log_train, h = length(log_test), drift = TRUE)

# ETS
fit_ets <- ets(log_train)
print(fit_ets)
fc_ets <- forecast(fit_ets, h = length(log_test))

# ARIMA
fit_arima <- auto.arima(log_train, stepwise = FALSE, approximation = FALSE)
fc_arima <- forecast(fit_arima, h = length(log_test))

# Evaluation

checkresiduals(fit_drift)
checkresiduals(fit_ets)
checkresiduals(fit_arima)

# Ljung-Box tests 
ljung_results <- list(
  Drift = Box.test(residuals(fit_drift), lag = 24, type = "Ljung-Box", fitdf = length(fit_drift$coef)),
  ETS   = Box.test(residuals(fit_ets),   lag = 24, type = "Ljung-Box", fitdf = length(fit_ets$par)),
  ARIMA = Box.test(residuals(fit_arima), lag = 24, type = "Ljung-Box", fitdf = length(fit_arima$coef))
)
print(ljung_results)

# Accuracy 
acc_drift_log <- accuracy(fit_drift$mean, log_test)
acc_ets_log   <- accuracy(fc_ets$mean,    log_test)
acc_arima_log <- accuracy(fc_arima$mean,  log_test)

accuracy_results <- list(
  Drift = round(acc_drift_log, 2),
  ETS   = round(acc_ets_log,   2),
  ARIMA = round(acc_arima_log, 2)
)
print(accuracy_results)

rmse_results <- data.frame(
  Model = c("Drift", "ETS", "ARIMA"),
  RMSE  = c(
    acc_drift_log["Test set", "RMSE"],
    acc_ets_log["Test set", "RMSE"],
    acc_arima_log["Test set", "RMSE"]
  )
)

kable(rmse_results, caption = "Comparison of RMSE (log scale)")

best_model <- rmse_results[which.min(rmse_results$RMSE), , drop = FALSE]
cat("\nBest forecasting model (lowest log-RMSE):\n")
print(best_model)

# Back-transform forecasts 

orig_test <- exp(log_test)
drift_levels <- exp(fit_drift$mean)  
ets_levels   <- exp(fc_ets$mean)    
# bias correction 
arima_levels <- exp(fc_arima$mean + 0.5 * fit_arima$sigma2)

# Table of actual vs predicted 
table_PCE <- data.frame(
  Date = as.yearmon(time(orig_test)),
  Actual = as.numeric(orig_test),
  Predicted_ARIMA = as.numeric(arima_levels)
)
print(table_PCE)

# Plot actual vs model forecasts 
autoplot(orig_test, series = "Actual") +
  autolayer(drift_levels, series = "Drift") +
  autolayer(ets_levels, series = "ETS") +
  autolayer(arima_levels, series = "ARIMA") +
  xlab("Year") + ylab("PCE")

#  12-month forecast 
fit_arima_full <- auto.arima(log_pce, stepwise = FALSE, approximation = FALSE)
fc_final <- forecast(fit_arima_full, h = 12)

# Back-transform 
final_levels <- exp(fc_final$mean + 0.5 * fit_arima_full$sigma2)
pce_levels_full <- exp(log_pce)

# Residual diagnostics on final model
checkresiduals(fit_arima_full)

forecast_table <- data.frame(
  Date = as.yearmon(time(final_levels)),
  PCE  = as.numeric(final_levels)
)
print(forecast_table)

# Plot final 12-month forecast 
autoplot(pce_levels_full) +
  autolayer(final_levels, series = "ARIMA Forecast") +
  xlab("Year") + ylab("PCE") +
  ggtitle("PCE Forecast for the Next 12 Months (ARIMA)")
