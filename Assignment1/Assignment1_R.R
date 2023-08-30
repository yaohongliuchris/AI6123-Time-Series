library(forecast)
library(tseries)

location="C:/Users/liu92/Documents/MSAI/Time Series/Assignment1/wwwusage.txt"
y <- scan(location, skip = 1)
max(y)
min(y)
mean(y)

acf(y, lag.max = 100)
pacf(y, lag.max = 100)
ts.plot(y, gpars=list(xlab="minutes"
                      ,ylab="Number of user"
                      ,main='Numbers of users connected to the Internet every minute')
        
##D1        
d1<-diff(y)
max(d1)
min(d1)     

ts.plot(d1,xlim = c(0,99), ylim=c(-15, 15)
        ,gpars=list(xlab="1st differencing time point", ylab="Number of user"
                  ,main='Time Plot after apply 1st differencing'))
acf(d1, lag.max = 99,main='1st differencing ACF')
pacf(d1, lag.max = 99,main='1st differencing PACF')
ar.yw(d1, order.max = 10)
fit_d1 <- arima(y, order= c(3,1,0))
tsdiag(fit_d1)

trainSet <- y[1:90]
testSet <- y[91:100]
fit_d1_train = arima(trainSet, order= c(3,1,0))
forecast1 = predict(fit_d1_train, n.ahead=10)


plot(c(y),
      main = "Validation of ARIMA(3, 1, 0)",
      xlab = "minutes",
      ylab = "Number of users",
      type = "o",
      xlim = c(0,100),
      ylim = c(50,250))
lines(1:90, trainSet-fit_d1_train$residuals, type="o", col="red")
lines(91:100, forecast1$pred, type="o", col="red")
lines(91:100, forecast1$pred-1.96*forecast$se, col="blue")
lines(91:100, forecast1$pred+1.96*forecast$se, col="blue")

##D2
d2<-diff(d1)
max(d2)
min(d2) 

ts.plot(d2,xlim = c(0,98), ylim=c(-10, 10) 
        ,gpars=list(xlab="2nd differencing time point", ylab="Number of user"
                   ,main='Time Plot after apply 2nd differencing'))
acf(d2, lag.max = 98,main='2nd differencing ACF')
pacf(d2, lag.max = 98,main='2nd differencing PACF')
fit_d2 <- arima(y, order= c(2,2,0))

fit_d2_train = arima(trainSet, order= c(2,2,0))
forecast2 = predict(fit_d2_train, n.ahead=10)

plot(c(y),
     main = "Validation of ARIMA(2, 2, 0)",
     xlab = "minutes",
     ylab = "Number of users",
     type = "o",
     xlim = c(0,100),
     ylim = c(50,250))
lines(1:90, trainSet-fit_d2_train$residuals, type="o", col="red")
lines(91:100, forecast2$pred, type="o", col="red")
lines(91:100, forecast2$pred-1.96*forecast$se, col="blue")
lines(91:100, forecast2$pred+1.96*forecast$se, col="blue")


#sequential test
aicd1 <- matrix(NA, 7, 7, dimnames = list(p = 0:6, q = 0:6)) 
aicd2 <- matrix(NA, 7, 7, dimnames = list(p = 0:6, q = 0:6)) 
for(p in 0:6){
  for(q in 0:6){
    tryCatch({
      if(!(p == 0 & q == 0)){
        aicd1[p+1, q+1] <-AIC(arima(y, c(p, 1, q)))}
    }, error = function(e) {
      print(paste("An error occurred:", e$message))
      NaN
    }, warning = function(w) {
      print(paste("A warning occurred:", w$message))
      NaN
    }, finally = {
      print("Done.")
    })
  }}
for(p in 0:6){
  for(q in 0:6){
    tryCatch({
      if(!(p == 0 & q == 0)){
        aicd2[p+1, q+1] <-AIC(arima(y, c(p, 2, q)))}
    }, error = function(e) {
      print(paste("An error occurred:", e$message))
      NaN
    }, warning = function(w) {
      print(paste("A warning occurred:", w$message))
      NaN
    }, finally = {
      print("Done.")
    })
  }}



min(aics_d1, na.rm = TRUE)#511.1393
min(aics_d2, na.rm = TRUE)#509.8135
arima525= arima(y, order= c(5,2,5))
tsdiag(arima525)

arima525_train = arima(trainSet, order= c(5,2,5))
forecast3 = predict(arima525_train, n.ahead=10)

plot(c(y),
     main = "Validation of ARIMA(5, 2, 5)",
     xlab = "minutes",
     ylab = "Number of users",
     type = "o",
     xlim = c(0,100),
     ylim = c(50,250))
lines(1:90, trainSet-arima525_train$residuals, type="o", col="red")
lines(91:100, forecast3$pred, type="o", col="red")
lines(91:100, forecast3$pred-1.96*forecast$se, col="blue")
lines(91:100, forecast3$pred+1.96*forecast$se, col="blue")


##auto fit arima
autoarima <- auto.arima(y)
autoarima
autoarima_train = auto.arima(trainSet)
forecast4 = predict(autoarima_train, n.ahead=10)

plot(c(y),
     main = "Validation of ARIMA(1, 1, 1)",
     xlab = "minutes",
     ylab = "Number of users",
     type = "o",
     xlim = c(0,100),
     ylim = c(50,250))
lines(1:90, trainSet-autoarima_train$residuals, type="o", col="red")
lines(91:100, forecast4$pred, type="o", col="red")
lines(91:100, forecast4$pred-1.96*forecast$se, col="blue")
lines(91:100, forecast4$pred+1.96*forecast$se, col="blue")

##AIC Checking, Accuracy checking
AIC(fit_d1)
AIC(fit_d2)
AIC(arima525)
AIC(autoarima)
accuracy(testSet,forecast1$pred)
accuracy(testSet,forecast2$pred)
accuracy(testSet,forecast3$pred)
accuracy(testSet,forecast4$pred)
