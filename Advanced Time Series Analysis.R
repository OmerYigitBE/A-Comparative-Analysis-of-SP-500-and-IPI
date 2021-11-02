#Preparation
setwd("~/Desktop/Classes/Advanced Time Series Analysis/PROJECT/")
library(readr)
library(CADFtest)
library(forecast)
library(vars)
library(fGarch)

#Importing the raw data
#Datasets are modified manually to include only the values of interest.
trsp500 <- read_table("trsp500.txt", col_names = F)
indpro <- read_table("indpro.txt", col_names = F)

#Transforming the data as time series
ts_trsp500 <- ts(trsp500, frequency = 12, start = c(1970,1))
ts_indpro <- ts(indpro, frequency = 12, start = c(1970,1))

#Plots of the time series data
ts.plot(ts_trsp500)
ts.plot(ts_indpro)

#Correlograms of the time series data
acf(ts_trsp500)
pacf(ts_trsp500)
acf(ts_indpro)
pacf(ts_indpro)

#Tests on untransformed time series data
max.lag <- 24 #Sqrt of 588 (length of data), rounded to the nearest integer.
CADFtest(ts_trsp500, type = "trend", criterion = "BIC", max.lag.y = max.lag)
Box.test(ts_trsp500, lag = max.lag, type = "Ljung-Box")
CADFtest(ts_indpro, type = "trend", criterion = "BIC", max.lag.y = max.lag)
Box.test(ts_indpro, lag = max.lag, type = "Ljung-Box")
#Series are not stationary. Transformation is needed.

#Data transformation
ts_trsp500_dlog <- diff(log(ts_trsp500))
ts_indpro_dlog <- diff(log(ts_indpro))
ts.plot(ts_trsp500_dlog)
ts.plot(ts_indpro_dlog)

#Tests on transformed time series data
max.lag <- 24 #Sqrt of 588 (length of data), rounded to the nearest integer.
CADFtest(ts_trsp500_dlog, type = "trend", criterion = "BIC", max.lag.y = max.lag)
CADFtest(ts_trsp500_dlog, type = "drift", criterion = "BIC", max.lag.y = max.lag)
CADFtest(ts_indpro_dlog, type = "trend", criterion = "BIC", max.lag.y = max.lag)
CADFtest(ts_indpro_dlog, type = "drift", criterion = "BIC", max.lag.y = max.lag)
#Now, the series are stationary.

#Correlograms of the transformed data
acf(ts_trsp500_dlog); pacf(ts_trsp500_dlog) #No autocorrelation, no partial autocorrelation.
acf(ts_indpro_dlog); pacf(ts_indpro_dlog) #Seasonal autocorrelation with lag 2.

###########################
### UNIVARIATE ANALYSES ###
###########################

#ARIMA model attempts
model_trsp500_1 <- arima(ts_trsp500_dlog, order = c(1,0,1))
model_trsp500_2 <- auto.arima(ts_trsp500_dlog)
summary(model_trsp500_1) #All terms are significant.
summary(model_trsp500_2) #SAR1 term isn't significant.
acf(model_trsp500_1$residuals)
acf(model_trsp500_2$residuals)
Box.test(model_trsp500_1$residuals, lag = max.lag, type = "Ljung-Box")
Box.test(model_trsp500_2$residuals, lag = max.lag, type = "Ljung-Box")
#Models are valid.

model_indpro_1 <- arima(ts_indpro_dlog, order = c(1,0,1), seasonal = list(order = c(2,0,1)))
model_indpro_2 <- auto.arima(ts_indpro_dlog)
summary(model_indpro_1) #SAR1 and SMA1 aren't significant.
summary(model_indpro_2) #SAR1 isn't significant.
acf(model_indpro_1$residuals)
acf(model_indpro_2$residuals)
Box.test(model_indpro_1$residuals, lag = max.lag, type = "Ljung-Box")
Box.test(model_indpro_2$residuals, lag = max.lag, type = "Ljung-Box")
#Models are valid.

#Model comparison with squared forecast error
#For trsp500
y <- ts_trsp500_dlog
S <- round(0.75*length(y))
h <- 1
#Errors for model 1
error1_h <- c()
for (i in S:(length(y) - h))
{
        model_trsp500_1_sub <- arima(y[1:i], order = c(1,0,1))
        predict_h <- predict(model_trsp500_1_sub, n.ahead = h)$pred[h]
        error1_h <- c(error1_h, y[i+h] - predict_h)
}
#Errors for model 2
error2_h <- c()
for (i in S:(length(y) - h))
{
        model_trsp500_2_sub <- arima(y[1:i], order = c(0,0,0),seasonal=c(1,0,0))
        predict_h <- predict(model_trsp500_2_sub, n.ahead = h)$pred[h]
        error2_h <- c(error2_h, y[i+h] - predict_h)
}

#Diebold-Mariano test
dm.test(error1_h, error2_h, h = h, power = 2)

#For indpro
y <- ts_indpro_dlog
S <- round(0.75*length(y))
h <- 1
#Errors for model 1
error1_h <- c()
for (i in S:(length(y) - h))
{
        model_indpro_1_sub <- arima(y[1:i], order = c(1,0,1),seasonal=c(2,0,1))
        predict_h <- predict(model_indpro_1_sub, n.ahead = h)$pred[h]
        error1_h <- c(error1_h, y[i+h] - predict_h)
}
#Errors for model 2
error2_h <- c()
for (i in S:(length(y) - h))
{
        model_indpro_2_sub <- arima(y[1:i], order = c(1,0,1),seasonal=c(2,0,0))
        predict_h <- predict(model_indpro_2_sub, n.ahead = h)$pred[h]
        error2_h <- c(error2_h, y[i+h] - predict_h)
}

#Diebold-Mariano test
dm.test(error1_h, error2_h, h = h, power = 2)

#Forecasts
#For trsp500
forecast_trsp500_dlog <- forecast(model_trsp500_1, h = 12)
plot(forecast_trsp500_dlog, xlim = c(2015,2020))
forecast_trsp500_dlog_expected <- c(forecast_trsp500_dlog$mean)
forecast_trsp500_dlog_lower <- c(forecast_trsp500_dlog$lower)
forecast_trsp500_dlog_upper <- c(forecast_trsp500_dlog$upper)
#Back-transformation
forecast_trsp500_expected <- c(trsp500[588,1]$X1)
forecast_trsp500_lower <- c()
forecast_trsp500_upper <- c()
for (i in 1:12)
{
        forecast_trsp500_expected[i+1] <- exp(forecast_trsp500_dlog_expected[i])*forecast_trsp500_expected[i]
        forecast_trsp500_lower[i] <- exp(forecast_trsp500_dlog_lower[i])*forecast_trsp500_expected[i]
        forecast_trsp500_upper[i] <- exp(forecast_trsp500_dlog_upper[i])*forecast_trsp500_expected[i]
}
forecast_trsp500_expected <- forecast_trsp500_expected[-1]

#For indpro
forecast_indpro_dlog <- forecast(model_indpro_1, h = 12)
plot(forecast_indpro_dlog, xlim = c(2015,2020))
forecast_indpro_dlog_expected <- c(forecast_indpro_dlog$mean)
forecast_indpro_dlog_lower <- c(forecast_indpro_dlog$lower)
forecast_indpro_dlog_upper <- c(forecast_indpro_dlog$upper)
#Back transformation
forecast_indpro_expected <- c(indpro[588,1]$X1)
forecast_indpro_lower <- c()
forecast_indpro_upper <- c()
for (i in 1:12)
{
        forecast_indpro_expected[i+1] <- exp(forecast_indpro_dlog_expected[i])*forecast_indpro_expected[i]
        forecast_indpro_lower[i] <- exp(forecast_indpro_dlog_lower[i])*forecast_indpro_expected[i]
        forecast_indpro_upper[i] <- exp(forecast_indpro_dlog_upper[i])*forecast_indpro_expected[i]
}
forecast_indpro_expected <- forecast_indpro_expected[-1]

#Plots of forecasted values
#For trsp500
plot.ts(ts_trsp500,xlim=c(2015,2020), ylim=c(3000,6000))
lines(ts(forecast_trsp500_expected, frequency = 12, start = c(2019,1)), col="red")
lines(ts(forecast_trsp500_lower, frequency = 12, start = c(2019,1)),col="blue")
lines(ts(forecast_trsp500_upper, frequency = 12, start = c(2019,1)),col="blue")

#For indpro
plot.ts(ts_indpro,xlim=c(2015,2020), ylim=c(100,120))
lines(ts(forecast_indpro_expected, frequency = 12, start = c(2019,1)), col="red")
lines(ts(forecast_indpro_lower, frequency = 12, start = c(2019,1)),col="blue")
lines(ts(forecast_indpro_upper, frequency = 12, start = c(2019,1)),col="blue")

#############################
### MULTIVARIATE ANALYSES ###
#############################

#Cointegration (Engle - Granger approach) - Raw Data
lm_1 <- lm(ts_indpro ~ ts_trsp500); lm_2 <- lm(ts_trsp500 ~ ts_indpro)
summary(lm_1); summary(lm_2)
plot.ts(lm_1$residuals); plot.ts(lm_2$residuals)
acf(lm_1$residuals); acf(lm_2$residuals)
Box.test(lm_1$residuals, lag = max.lag, type = "Ljung-Box")
Box.test(lm_2$residuals, lag = max.lag, type = "Ljung-Box")
#Linear models aren't valid.

#Cointegration (Engle - Granger approach) - Differenced data
lm_3 <- lm(diff(ts_indpro) ~ diff(ts_trsp500)); lm_4 <- lm(diff(ts_trsp500) ~ diff(ts_indpro))
summary(lm_3); summary(lm_4)
plot.ts(lm_3$residuals); plot.ts(lm_4$residuals) #Heteroscedasticity on the residuals of lm_4
acf(lm_3$residuals); acf(lm_4$residuals) #Autocorrelation on the residuals of lm_3
Box.test(lm_3$residuals, lag = max.lag, type = "Ljung-Box")
Box.test(lm_4$residuals, lag = max.lag, type = "Ljung-Box")
#Linear models aren't valid.

#Cointegration (Engle - Granger approach) - Log-differenced data
lm_5 <- lm(ts_indpro_dlog ~ ts_trsp500_dlog); lm_6 <- lm(ts_trsp500_dlog ~ ts_indpro_dlog)
summary(lm_5); summary(lm_6) #High p-values. Models are not adequate.
plot.ts(lm_5$residuals); plot.ts(lm_6$residuals) #Homoscedastic residuals.
acf(lm_5$residuals); acf(lm_6$residuals) #Autocorrelation on the residuals of lm_5
Box.test(lm_5$residuals, lag = max.lag, type = "Ljung-Box")
Box.test(lm_6$residuals, lag = max.lag, type = "Ljung-Box")
#lm_6 is valid, but lm_5 isn't.
CADFtest(lm_6$residuals,type="drift",criterion="BIC",max.lag.y=max.lag)
#Cointegration occurs.

#Error correction model (ECM)
res_fit_lm_6 <- lm_6$residuals
ECT <- res_fit_lm_6[-length(res_fit_lm_6)]
fit_ecm <- lm(ts_trsp500_dlog[2:length(ts_trsp500_dlog)] ~ ts_indpro_dlog[2:length(ts_indpro_dlog)] + ECT) #Correction
summary(fit_ecm)

#DLM & ADLM
#DL(1)
lag <- 1
n <- length(ts_trsp500_dlog)
ts_trsp500_dlog.0 <- ts_trsp500_dlog[(lag+1):n]
ts_indpro_dlog.0 <- ts_indpro_dlog[(lag+1):n]
ts_indpro_dlog.1 <- ts_indpro_dlog[lag:(n-1)]
fit_dlm_1 <- lm(ts_trsp500_dlog.0 ~ ts_indpro_dlog.0 + ts_indpro_dlog.1)
summary(fit_dlm_1)
ts.plot(fit_dlm_1$residuals)
acf(fit_dlm_1$residuals)
Box.test(fit_dlm_1$residuals, lag = max.lag, type = "Ljung-Box")
#Model is valid.

#DL(2)
lag <- 2
n <- length(ts_trsp500_dlog)
ts_trsp500_dlog.0 <- ts_trsp500_dlog[(lag+1):n]
ts_indpro_dlog.0 <- ts_indpro_dlog[(lag+1):n]
ts_indpro_dlog.1 <- ts_indpro_dlog[lag:(n-1)]
ts_indpro_dlog.2 <- ts_indpro_dlog[(lag-1):(n-2)]
fit_dlm_2 <- lm(ts_trsp500_dlog.0 ~ ts_indpro_dlog.0 + ts_indpro_dlog.1 + ts_indpro_dlog.2)
summary(fit_dlm_2)
ts.plot(fit_dlm_2$residuals)
acf(fit_dlm_2$residuals)
Box.test(fit_dlm_2$residuals, lag = max.lag, type = "Ljung-Box")
#Model is valid.

#ADLM(1)
lag <- 1
n <- length(ts_trsp500_dlog)
ts_trsp500_dlog.0 <- ts_trsp500_dlog[(lag+1):n]
ts_trsp500_dlog.1 <- ts_trsp500_dlog[lag:(n-1)]
ts_indpro_dlog.0 <- ts_indpro_dlog[(lag+1):n]
fit_adlm_1 <- lm(ts_trsp500_dlog.0 ~ ts_trsp500_dlog.1 + ts_indpro_dlog.0)
summary(fit_adlm_1)
ts.plot(fit_adlm_1$residuals)
acf(fit_adlm_1$residuals)
Box.test(fit_adlm_1$residuals, lag = max.lag, type = "Ljung-Box")
#Model is valid.

#ADLM(2)
lag <- 2
n <- length(ts_trsp500_dlog)
ts_trsp500_dlog.0 <- ts_trsp500_dlog[(lag+1):n]
ts_trsp500_dlog.1 <- ts_trsp500_dlog[lag:(n-1)]
ts_trsp500_dlog.2 <- ts_trsp500_dlog[(lag-1):(n-2)]
ts_indpro_dlog.0 <- ts_indpro_dlog[(lag+1):n]
fit_adlm_2 <- lm(ts_trsp500_dlog.0 ~ ts_trsp500_dlog.1 + ts_trsp500_dlog.2 + ts_indpro_dlog.0)
summary(fit_adlm_2)
ts.plot(fit_adlm_2$residuals)
acf(fit_adlm_2$residuals)
Box.test(fit_adlm_2$residuals, lag = max.lag, type = "Ljung-Box")
#Model is valid.

#Granger-Causality
fit_dlm <- lm(ts_trsp500_dlog.0 ~ ts_trsp500_dlog.1 + ts_indpro_dlog.1)
fit_dlm_nox <- lm(ts_trsp500_dlog.0 ~ ts_trsp500_dlog.1)
anova(fit_dlm, fit_dlm_nox)
#No granger-causality occurs.

#VAR
investment <- ts_trsp500_dlog
production <- ts_indpro_dlog
mydata <- cbind(investment, production)
VARselect(mydata) #VAR(3) is selected
fit_var <- VAR(mydata, type = "const", p = 3)
summary(fit_var)
var_predictions <- predict(fit_var, n.ahead = 12)

#Forecasts
#For trsp500
forecast_investment_expected <- c(var_predictions$fcst$investment)
forecast_investment_lower <- c(var_predictions[[1]][[1]][,2])
forecast_investment_upper <- c(var_predictions[[1]][[1]][,3])
#Back-transformation
forecast_trsp500_expected_var <- c(trsp500[588,1]$X1)
forecast_trsp500_lower_var <- c()
forecast_trsp500_upper_var <- c()
for (i in 1:12)
{
        forecast_trsp500_expected_var[i+1] <- exp(forecast_investment_expected[i])*forecast_trsp500_expected_var[i]
        forecast_trsp500_lower_var[i] <- exp(forecast_investment_lower[i])*forecast_trsp500_expected_var[i]
        forecast_trsp500_upper_var[i] <- exp(forecast_investment_upper[i])*forecast_trsp500_expected_var[i]
}
forecast_trsp500_expected_var <- forecast_trsp500_expected_var[-1]

#For indpro
forecast_production_expected <- c(var_predictions$fcst$production)
forecast_production_lower <- c(var_predictions[[1]][[2]][,2])
forecast_production_upper <- c(var_predictions[[1]][[2]][,3])
#Back-transformation
forecast_indpro_expected_var <- c(indpro[588,1]$X1)
forecast_indpro_lower_var <- c()
forecast_indpro_upper_var <- c()
for (i in 1:12)
{
        forecast_indpro_expected_var[i+1] <- exp(forecast_production_expected[i])*forecast_indpro_expected_var[i]
        forecast_indpro_lower_var[i] <- exp(forecast_production_lower[i])*forecast_indpro_expected_var[i]
        forecast_indpro_upper_var[i] <- exp(forecast_production_upper[i])*forecast_indpro_expected_var[i]
}
forecast_indpro_expected_var <- forecast_indpro_expected_var[-1]

#Plots of forecasted values
#For trsp500
plot.ts(ts_trsp500,xlim=c(2015,2020), ylim=c(3000,6000))
lines(ts(forecast_trsp500_expected_var, frequency = 12, start = c(2019,1)), col="red")
lines(ts(forecast_trsp500_lower_var, frequency = 12, start = c(2019,1)),col="blue")
lines(ts(forecast_trsp500_upper_var, frequency = 12, start = c(2019,1)),col="blue")

#For indpro
plot.ts(ts_indpro,xlim=c(2015,2020), ylim=c(100,120))
lines(ts(forecast_indpro_expected_var, frequency = 12, start = c(2019,1)), col="red")
lines(ts(forecast_indpro_lower_var, frequency = 12, start = c(2019,1)),col="blue")
lines(ts(forecast_indpro_upper_var, frequency = 12, start = c(2019,1)),col="blue")

