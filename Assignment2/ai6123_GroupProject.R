library(forecast)
library(tseries)
library(urca)
library(stats4)
library(sarima)
library(ggplot2)
library(dplyr)
library(lubridate)
library(astsa)
library(gridExtra)

data = read.table("D:/NTU/Time Series/group_project/drug.txt", header=TRUE, sep = ",")
data <- data %>% mutate(date = ymd(date))
data.start_time <<- ymd(as.Date(data$date[[1]])) %>% print()
data.end_time <<- ymd(as.Date(data$date[[nrow(data)]])) %>% print()
data.frequence <<- 12
# convert data to time series
data_ts <- ts(data$value, start=c(year(data.start_time), month(data.start_time)),
               end=c(year(data.end_time), month(data.end_time)), frequency=data.frequence)
plot(data_ts, main="Monthly Anti-diabetic Drug Sales in Australia", xlab="Year", ylab="Drug Sales")
# analysis of data_ts
min(data_ts)  # 2.81452
max(data_ts)  # 29.66536
mean(data_ts) # 10.69443
data_ts %>% ur.kpss() %>% summary()
data_ts %>% ur.df() %>% summary()
# seasonal decomposition
autoplot(stl(data_ts, s.window="period"))
# ACF/PACF
ggtsdisplay(data_ts)

# split data_ts
train_len = round(length(data_ts)) * 0.8
data_train <- head(data_ts, train_len)
data_test <- tail(data_ts, round(length(data_ts) - train_len))
length(data_train) # 163
length(data_test)  # 41

# Find suitable SARIMA model
lambda_zero <- 0
lambda_auto <- BoxCox.lambda(data_ts, 'guerrero') %>% print() # 0.1313326
data_train_t_zero <- BoxCox(data_train,lambda_zero)
data_train_t_auto <- BoxCox(data_train,lambda_auto)
autoplot(stl(train_data_t_zero, s.window="period")) + ggtitle("BoxCox Transformed Data (lambda = 0)")
autoplot(stl(train_data_t_auto, s.window="period")) + ggtitle("BoxCox Transformed Data (auto lambda)")
# Choose transformation
data_train_t <- data_train_t_zero
#data_train_t <- data_train_t_auto
# Apply seasonal differencing
data_train_lag = diff(data_train_t)
ggtsdisplay(data_train_lag, main  = "Diff-1 Seasonal Differencing", ylab  = "Drug Sales", xlab  = "Year",lag.max = 60)
data_train_diff_lag = diff(data_train_lag, lag = data.frequence)
ggtsdisplay(data_train_diff_lag, main  = "Lag-12 + Diff-1 Seasonal Differencing", ylab  = "Drug Sales", xlab  = "Year",lag.max = 60)

# 017110 model
model.017110 = stats::arima(data_train_t, 
                            order=c(0,1,7), 
                            seasonal=list(order=c(1,1,0), 
                                          period=data.frequence))
checkresiduals(model.017110)
tsdiag(model.017110)
qqnorm(model.017110$residuals)
qqline(model.017110$residuals)
AIC(model.017110) # -416.971
BIC(model.017110) # -389.8753
summary(model.017110)
#    RMSE        MAE
#   0.05409898 0.04039567

pred.017110 <- InvBoxCox(forecast(model.017110, h = length(data_test)+36)$mean, lambda = lambda_zero)
plot(pred.017110,
     main="Prediction of Drug Sales of SARIMA(0,1,7)(1,1,0)", 
     xlab="Year", 
     ylab="Drug Sales",lwd=2,col="red",lty=3)
lines(data_test,lwd=2);
accuracy(pred.017110, data_test)
#    RMSE        MAE
#   2.229626 1.766395


# 510110 model
model.510110 = stats::arima(data_train_t, 
                            order=c(5,1,0), 
                            seasonal=list(order=c(1,1,0), 
                                          period=data.frequence))
checkresiduals(model.510110)
tsdiag(model.510110)
qqnorm(model.510110$residuals)
qqline(model.510110$residuals)
AIC(model.510110) # -422.1871
BIC(model.510110) # -401.1126
summary(model.510110)
#    RMSE        MAE
# 0.053885 0.04006019

pred.510110 <- InvBoxCox(forecast(model.510110, h = length(data_test)+36)$mean, lambda = lambda_zero)
plot(pred.510110,
     main="Prediction of Drug Sales of SARIMA(5,1,0)(1,1,0)", 
     xlab="Year", 
     ylab="Drug Sales",lwd=2,col="red",lty=3)
lines(data_test,lwd=2);
accuracy(pred.510110, data_test)
#    RMSE        MAE
#  2.289985 1.812177

# 510011 model
model.510011 = stats::arima(data_train_t, 
                            order=c(5,1,0), 
                            seasonal=list(order=c(0,1,1), 
                                          period=data.frequence))
checkresiduals(model.510011)
tsdiag(model.510011)
qqnorm(model.510011$residuals)
qqline(model.510011$residuals)
AIC(model.510011) # -434.6062
BIC(model.510011) # -413.5318
summary(model.510011)
#    RMSE        MAE
# 0.05098949 0.03778493

pred.510011 <- InvBoxCox(forecast(model.510011, h = length(data_test)+36)$mean, lambda = lambda_zero)
plot(pred.510011,
     main="Prediction of Drug Sales of SARIMA(5,1,0)(0,1,1)", 
     xlab="Year", 
     ylab="Drug Sales",lwd=2,col="red",lty=3)
lines(data_test,lwd=2);
accuracy(pred.510011, data_test)
#    RMSE        MAE
#  1.507028 1.198604


# 017011 model
model.017011 = stats::arima(data_train_t, 
                            order=c(0,1,7), 
                            seasonal=list(order=c(0,1,1), 
                                          period=data.frequence))
checkresiduals(model.017011)
tsdiag(model.017011)
qqnorm(model.017011$residuals)
qqline(model.017011$residuals)
AIC(model.017011) # -428.6995
BIC(model.017011) # -401.6038
summary(model.017011)
#    RMSE        MAE
# 0.05118094 0.03782738

pred.017011 <- InvBoxCox(forecast(model.017011, h = length(data_test)+36)$mean, lambda = lambda_zero)
plot(pred.017011,
     main="Prediction of Drug Sales of SARIMA(0,1,7)(0,1,1)", 
     xlab="Year", 
     ylab="Drug Sales",lwd=2,col="red",lty=3)
lines(data_test,lwd=2);
accuracy(pred.017011, data_test)
#    RMSE        MAE
# 1.483472 1.179295

# 517011 model
model.517011 = stats::arima(data_train_t, 
                            order=c(5,1,7), 
                            seasonal=list(order=c(0,1,1), 
                                          period=data.frequence))
checkresiduals(model.517011)
tsdiag(model.517011)
qqnorm(model.517011$residuals)
qqline(model.517011$residuals)
AIC(model.517011) # -437.5164
BIC(model.517011) # -395.3675
summary(model.517011)
#    RMSE        MAE
# 0.0471864 0.03341949

pred.517011 <- InvBoxCox(forecast(model.517011, h = length(data_test)+36)$mean, lambda = lambda_zero)
plot(pred.517011,
     main="Prediction of Drug Sales of SARIMA(5,1,7)(0,1,1)", 
     xlab="Year", 
     ylab="Drug Sales",lwd=2,col="red",lty=3)
lines(data_test,lwd=2);
accuracy(pred.517011, data_test)
#    RMSE        MAE
# 1.525499 1.238929

# 517111 model
model.517111 = stats::arima(data_train_t, 
                            order=c(5,1,7), 
                            seasonal=list(order=c(1,1,1), 
                                          period=data.frequence))
checkresiduals(model.517111)
tsdiag(model.517111)
qqnorm(model.517111$residuals)
qqline(model.517111$residuals)
AIC(model.517111) # -439.346
BIC(model.517111) # -394.1865
summary(model.517111)
#    RMSE        MAE
# 0.04605911 0.03299559

pred.517111 <- InvBoxCox(forecast(model.517111, h = length(data_test)+36)$mean, lambda = lambda_zero)
plot(pred.517111,
     main="Prediction of Drug Sales of SARIMA(5,1,7)(1,1,1)", 
     xlab="Year", 
     ylab="Drug Sales",lwd=2,col="red",lty=3)
lines(data_test,lwd=2);
accuracy(pred.517111, data_test)
#    RMSE        MAE
#  1.674455 1.376838

model.510010 = stats::arima(data_train_t, 
                            order=c(5,1,0), 
                            seasonal=list(order=c(0,1,0), 
                                          period=data.frequence))
checkresiduals(model.510010)
tsdiag(model.510010)
qqnorm(model.510010$residuals)
qqline(model.510010$residuals)
AIC(model.510010) # -416.2912
BIC(model.510010) # -398.2274
summary(model.510010)
#    RMSE        MAE
# 0.05547478 0.04048461

pred.510010 <- InvBoxCox(forecast(model.510010, h = length(data_test)+36)$mean, lambda = lambda_zero)
plot(pred.510010,
     main="Prediction of Drug Sales of SARIMA(5,1,0)(0,1,0)", 
     xlab="Year", 
     ylab="Drug Sales",lwd=2,col="red",lty=3)
lines(data_test,lwd=2);
accuracy(pred.510010, data_test)
#    RMSE        MAE
#  3.400101    2.956103

############ ARMA model with auto lambda ############
model.517011_auto = stats::arima(data_train_t_auto, 
                                 order=c(5,1,7), 
                                 seasonal=list(order=c(0,1,1), 
                                               period=data.frequence))
checkresiduals(model.517011_auto)
tsdiag(model.517011_auto)
qqnorm(model.517011_auto$residuals)
qqline(model.517011_auto$residuals)
AIC(model.517011_auto) 
BIC(model.517011_auto) 
summary(model.517011_auto)


model.517011_auto <- InvBoxCox(forecast(model.517011_auto, h = length(data_test)+36)$mean, lambda = lambda_auto)
plot(model.517011_auto,
     main="Prediction of Drug Sales of SARIMA(5,1,7)(0,1,1) with auto-lambda", 
     xlab="Year", 
     ylab="Drug Sales",lwd=2,col="red",lty=3)
lines(data_test,lwd=2);
accuracy(model.517011_auto, data_test)
#    RMSE        MAE
# 


# Brute-Force

# 13×13 matrix to store the AIC BIC RMSE MAE

len = 12

aic_mat_011 <- matrix(NA, len+1, len+1, dimnames = list(p = 0:len, q = 0:len)) 
bic_mat_011 <- matrix(NA, len+1, len+1, dimnames = list(p = 0:len, q = 0:len))
rmse_mat_011 <- matrix(NA, len+1, len+1, dimnames = list(p = 0:len, q = 0:len))
mae_mat_011 <- matrix(NA, len+1, len+1, dimnames = list(p = 0:len, q = 0:len))

aic_mat_110 <- matrix(NA, len+1, len+1, dimnames = list(p = 0:len, q = 0:len)) 
bic_mat_110 <- matrix(NA, len+1, len+1, dimnames = list(p = 0:len, q = 0:len))
rmse_mat_110 <- matrix(NA, len+1, len+1, dimnames = list(p = 0:len, q = 0:len))
mae_mat_110 <- matrix(NA, len+1, len+1, dimnames = list(p = 0:len, q = 0:len))

aic_mat_111 <- matrix(NA, len+1, len+1, dimnames = list(p = 0:len, q = 0:len)) 
bic_mat_111 <- matrix(NA, len+1, len+1, dimnames = list(p = 0:len, q = 0:len))
rmse_mat_111 <- matrix(NA, len+1, len+1, dimnames = list(p = 0:len, q = 0:len))
mae_mat_111 <- matrix(NA, len+1, len+1, dimnames = list(p = 0:len, q = 0:len))

for(p in 0:len)
  for(q in 0:len)
    if(!(p == 0 & q == 0)){
      tryCatch({
        model.cur = stats::arima(data_train_t,
                                 order=c(p,1,q),
                                 seasonal=list(order=c(0,1,1), 
                                               period=data.frequence), optim.method="BFGS")
        pred.cur <- InvBoxCox(forecast(model.cur, h = length(data_test)+36)$mean, lambda = lambda_zero)
        aic_mat_011[p+1, q+1] <- AIC(model.cur)
        bic_mat_011[p+1, q+1] <- BIC(model.cur)
        rmse_mat_011[p+1, q+1] <- accuracy(pred.cur, data_test)[2]
        mae_mat_011[p+1, q+1] <- accuracy(pred.cur, data_test)[3]},error=function(e){
          aic_mat_011[p+1, q+1] <- NA
          bic_mat_011[p+1, q+1] <- NA
          rmse_mat_011[p+1, q+1] <- NA
          mae_mat_011[p+1, q+1] <- NA
        })
      
      tryCatch({
        model.cur = stats::arima(data_train_t,
                                 order=c(p,1,q),
                                 seasonal=list(order=c(1,1,0), 
                                               period=data.frequence), optim.method="BFGS")
        pred.cur <- InvBoxCox(forecast(model.cur, h = length(data_test)+36)$mean, lambda = lambda_zero)
        aic_mat_110[p+1, q+1] <- AIC(model.cur)
        bic_mat_110[p+1, q+1] <- BIC(model.cur)
        rmse_mat_110[p+1, q+1] <- accuracy(pred.cur, data_test)[2]
        mae_mat_110[p+1, q+1] <- accuracy(pred.cur, data_test)[3]},error=function(e){
          aic_mat_110[p+1, q+1] <- NA
          bic_mat_110[p+1, q+1] <- NA
          rmse_mat_110[p+1, q+1] <- NA
          mae_mat_110[p+1, q+1] <- NA
        })
      
      tryCatch({
        model.cur = stats::arima(data_train_t,
                                 order=c(p,1,q),
                                 seasonal=list(order=c(1,1,1), 
                                               period=data.frequence), optim.method="BFGS")
        pred.cur <- InvBoxCox(forecast(model.cur, h = length(data_test)+36)$mean, lambda = lambda_zero)
        aic_mat_111[p+1, q+1] <- AIC(model.cur)
        bic_mat_111[p+1, q+1] <- BIC(model.cur)
        rmse_mat_111[p+1, q+1] <- accuracy(pred.cur, data_test)[2]
        mae_mat_111[p+1, q+1] <- accuracy(pred.cur, data_test)[3]},error=function(e){
          aic_mat_111[p+1, q+1] <- NA
          bic_mat_111[p+1, q+1] <- NA
          rmse_mat_111[p+1, q+1] <- NA
          mae_mat_111[p+1, q+1] <- NA
        })
      cat("p=", p, ", q=", q, "\n")
    }
# Find the location of minimum value
aic_mat - min(aic_mat_011, na.rm = TRUE)
bic_mat - min(bic_mat_011, na.rm = TRUE)
rmse_mat - min(rmse_mat_011, na.rm = TRUE)
mae_mat - min(mae_mat_011, na.rm = TRUE)

aic_mat - min(aic_mat_110, na.rm = TRUE)
bic_mat - min(bic_mat_110, na.rm = TRUE)
rmse_mat - min(rmse_mat_110, na.rm = TRUE)
mae_mat - min(mae_mat_110, na.rm = TRUE)

aic_mat - min(aic_mat_111, na.rm = TRUE)
bic_mat - min(bic_mat_111, na.rm = TRUE)
rmse_mat - min(rmse_mat_111, na.rm = TRUE)
mae_mat - min(mae_mat_111, na.rm = TRUE)


# best model according to AIC BIC
model.213011 = stats::arima(data_train_t, 
                            order=c(2,1,3), 
                            seasonal=list(order=c(0,1,1), 
                                          period=data.frequence))

pred.213011 <- InvBoxCox(forecast(model.213011, h = length(data_test)+36)$mean, lambda = lambda_zero)
plot(pred.213011,
     main="Prediction of Drug Sales of SARIMA(2,1,3)(0,1,1)", 
     xlab="Year", 
     ylab="Drug Sales",lwd=2,col="red",lty=3)
lines(data_test, lwd=2,col="black");
tsdiag(model.213011)


# Holt-Winter 

model.HW_add <- hw(data_train_t, seasonal = "additive")
model.HW_mul <- hw(data_train_t, seasonal = "multiplicative")
pred.HW_add <- InvBoxCox(forecast(model.HW_add$mean, h=length(data_test)+36)$mean, lambda = lambda_zero)
pred.HW_mul <- InvBoxCox(forecast(model.HW_mul$mean, h=length(data_test)+36)$mean, lambda = lambda_zero)
p_holtwinters <- ggplot() +
  geom_line(aes(x = 1:length(pred.HW_add), y = pred.HW_add),lwd=2,col="red",lty=3) +
  geom_line(aes(x = 1:length(pred.HW_mul), y = pred.HW_mul), col='purple',lwd=2,lty=3) +
  geom_line(aes(x = 1:length(data_test), y = data_test),lwd=2) +
  labs(title = "Prediction of Drug Sales of Holt-Winter Model", x = "Year", y = "Drug Sales")

accuracy(pred.HW_add, data_test) %>% print()
#    RMSE        MAE
#  1.857903   1.593702

accuracy(pred.HW_mul, data_test) %>% print()
#    RMSE        MAE
#  3.258563 2.777955

# store matrix to csv file
library(MASS)
write.matrix(aic_mat_011,file="aic_mat_011.csv",sep=",")
write.matrix(aic_mat_110,file="aic_mat_110.csv",sep=",")
write.matrix(aic_mat_111,file="aic_mat_111.csv",sep=",")
write.matrix(bic_mat_011,file="bic_mat_011.csv",sep=",")
write.matrix(bic_mat_110,file="bic_mat_110.csv",sep=",")
write.matrix(bic_mat_111,file="bic_mat_111.csv",sep=",")
write.matrix(mae_mat_011,file="mae_mat_011.csv",sep=",")
write.matrix(mae_mat_110,file="mae_mat_110.csv",sep=",")
write.matrix(mae_mat_111,file="mae_mat_111.csv",sep=",")
write.matrix(rmse_mat_011,file="rmse_mat_011.csv",sep=",")
write.matrix(rmse_mat_110,file="rmse_mat_110.csv",sep=",")
write.matrix(rmse_mat_111,file="rmse_mat_111.csv",sep=",")
# compaare different models
plot(pred.510011,
     main="Prediction of Drug Sales of different SARIMA", 
     xlab="Year", 
     ylab="Drug Sales",lwd=1,col="red",lty=1)
lines(pred.510110,lwd=1,col="green",lty=1)
lines(pred.017110,lwd=1,col="blue",lty=1)
lines(pred.017011,lwd=1,col="purple",lty=1)
lines(pred.517111,lwd=1,col="yellow",lty=1)
lines(data_test,lwd=2)
