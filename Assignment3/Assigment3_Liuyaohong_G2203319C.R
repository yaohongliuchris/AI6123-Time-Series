library(TSA)
library(astsa)
library(fBasics)
library(forecast)
library(ggplot2)
library(quantmod)
library(rugarch)
library(tseries)
library(zoo)

Sys.setlocale(category = "LC_TIME", locale = "C")
#dataset = getSymbols("AAPL", from='2002-02-01', to='2017-02-01'
#                     , src='yahoo', auto.assign = F) 
dataset =getSymbols('AAPL', from='2002-02-01', to='2017-02-01', src='yahoo', auto.assign = F) 
#write.zoo(dataset, file = "AAPL.csv", sep = ",", quote = FALSE)
dataset=na.omit(dataset)
ggtsdisplay(dataset$AAPL.Adjusted,main  = "AAPL Adjusted Close from 2002-02-01 to 2017-01-31"
            ,xlab="Date",ylab="Adj Close")

min(dataset$AAPL.Adjusted)#0.1994051
max(dataset$AAPL.Adjusted)#30.02323

adf.test(dataset$AAPL.Adjusted)
kpss.test(dataset$AAPL.Adjusted)
kpss.test(dataset$AAPL.Adjusted,null="Trend")

dataset.Weekly = ts(Ad(to.weekly(dataset)), frequency = 52)
autoplot(stl(dataset.Weekly[,1], s.window="period")
         ,main  = "Decomposition of Weekly Data",xlab  = "time")
test = stl(dataset.Weekly[,1], s.window="period")
datasetWithoutSeasonal = test$time.series[,2] + test$time.series[,3]
datasetWithoutSeasonal =na.omit(datasetWithoutSeasonal)


dataset.log =diff(log(dataset$AAPL.Adjusted))
dataset.log = dataset.log[!is.na(dataset.log)]



#dataset.log = diff(diff(log(dataset$AAPL.Adjusted))*100,lag=365)
#dataset.log = na.omit(dataset.log)
#qqnorm(dataset.log)
p1=ggtsdisplay(dataset.log,
         main  = sprintf("Log Returns"),
         xlab="Date",
         ylab="log Adj Close")
#p1=p1+xlab("Date")+ylab("log Adj Close")
#p1+geom_abline(slope = 0, intercept = 0, col = 'white')

adf.test(dataset.log)
kpss.test(dataset.log)

ggtsdisplay(abs(dataset.log),
            main  = sprintf("Absolute of Log Returns "),
            xlab="Date",
            ylab="log Adj Close")

ggtsdisplay(dataset.log^2,
            main  = sprintf("square of Log Returns "),
            xlab="Date",
            ylab="log Adj Close")
qqnorm(dataset.log)
qqline(dataset.log, col = "red") 
basicStats(dataset.log)


eacf(dataset.log)
eacf(abs(dataset.log))
eacf(dataset.log^2)


garch11=garch(dataset.log, order=c(1,1))
AIC(garch11)#-18556.44


garch02=garch(dataset.log, order=c(0,2))
AIC(garch02)#-18205.25

ggtsdisplay(residuals(garch11))
summary(garch11)
qqnorm(residuals(garch11))
qqline(residuals(garch11), col = "red")

ggtsdisplay(residuals(garch11)^2)
gBox(garch11,method='squared')

# sGARCH(1,1) model
# Normal Distribution
ugarchfit(spec = ugarchspec(
      variance.model=list(garchOrder=c(1,1)),
      mean.model=list(armaOrder = c(0,0))),
      data = dataset.log)


# Skew Normal Distribution
ugarchfit(spec = ugarchspec(
  variance.model=list(garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "snorm"),
  data = dataset.log)

# T-Distribution
ugarchfit(spec = ugarchspec(
  variance.model=list(garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "std"),
  data = dataset.log)

# Skew T-Distribution
ugarchfit(spec = ugarchspec(
  variance.model=list(garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "sstd"),
  data = dataset.log)

# Generalized Error Distribution
ugarchfit(spec = ugarchspec(
  variance.model=list(garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "ged"),
  data = dataset.log)

# Skew Generalized Error Distribution
ugarchfit(spec = ugarchspec(
  variance.model=list(garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "sged"),
  data = dataset.log)

# Normal Inverse Gaussian Distribution
ugarchfit(spec = ugarchspec(
  variance.model=list(garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "nig"),
  data = dataset.log)

# Generalized Hyperbolic Distribution
ugarchfit(spec = ugarchspec(
  variance.model=list(garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "ghyp"),
  data = dataset.log)

#model fitting
# fGARCH(1,1), GARCH
ugarchfit(spec = ugarchspec(
  variance.model=list(model = "fGARCH",
                      submodel = "GARCH",
                      garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "std"),
  data = dataset.log)

# fGARCH(1,1), TGARCH
ugarchfit(spec = ugarchspec(
  variance.model=list(model = "fGARCH",
                      submodel = "TGARCH",
                      garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "std"),
  data = dataset.log)

# fGARCH(1,1), AVGARCH
ugarchfit(spec = ugarchspec(
  variance.model=list(model = "fGARCH",
                      submodel = "AVGARCH",
                      garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "std"),
  data = dataset.log)

# fGARCH(1,1), NAGARCH
ugarchfit(spec = ugarchspec(
  variance.model=list(model = "fGARCH",
                      submodel = "NAGARCH",
                      garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "std"),
  data = dataset.log)

# fGARCH(1,1), APARCH
ugarchfit(spec = ugarchspec(
  variance.model=list(model = "fGARCH",
                      submodel = "APARCH",
                      garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "std"),
  data = dataset.log)

# fGARCH(1,1), GJRGARCH
ugarchfit(spec = ugarchspec(
  variance.model=list(model = "fGARCH",
                      submodel = "GJRGARCH",
                      garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "std"),
  data = dataset.log)

# fGARCH(1,1), ALLGARCH
ugarchfit(spec = ugarchspec(
  variance.model=list(model = "fGARCH",
                      submodel = "ALLGARCH",
                      garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "std"),
  data = dataset.log)

# eGARCH(1,1)
ugarchfit(spec = ugarchspec(
  variance.model=list(model = "eGARCH",
                      garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "std"),
  data = dataset.log)

# gjrGARCH(1,1)
ugarchfit(spec = ugarchspec(
  variance.model=list(model = "gjrGARCH",
                      garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "std"),
  data = dataset.log)

# apARCH(1,1)
ugarchfit(spec = ugarchspec(
  variance.model=list(model = "apARCH",
                      garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "std"),
  data = dataset.log)

# iGARCH(1,1)
ugarchfit(spec = ugarchspec(
  variance.model=list(model = "iGARCH",
                      garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "std"),
  data = dataset.log)

#csGARCH(1,1)
ugarchfit(spec = ugarchspec(
  variance.model=list(model = "csGARCH",
                      garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "std"),
  data = dataset.log)

#Forecasting
egarch11 = ugarchfit(spec = ugarchspec(
  variance.model=list(model = "eGARCH",
                      garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "std"),
  data = dataset.log,
  out.sample = 30)

egarch11.forecast = ugarchforecast(egarch11, n.ahead = 60, n.roll = 30)

plot(egarch11.forecast, which = 1)
plot(egarch11.forecast, which = 2)
plot(egarch11.forecast, which = 3)
plot(egarch11.forecast, which = 4)

