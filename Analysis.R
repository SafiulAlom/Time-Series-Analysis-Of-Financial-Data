
#-----------------------------------------------------------------------------------#
#*******************Time series analysis of financial data**************************#
#-----------------------------------------------------------------------------------#

#-----------------------------------------------------------------------------------#
#                       read data + install and load packages
#-----------------------------------------------------------------------------------#

#remove all objects
remove(list = ls())
dev.off()

#determine directory path and library path
dir.path = 'C:\\Users\\Himel\\OneDrive\\Studium\\M.Sc. Statistics\\1_Statistical Programming Language\\Lecture_Exercise\\Analysis\\Project_alomsafi'
libraryPath = 'C:\\Users\\Himel\\OneDrive\\Studium\\R\\Packages'

#read the Function.R, which includes all manually implemented functions
setwd(dir.path)
source(paste(getwd(),'\\R\\Function.R', sep = ""))

#list of the required packages
packages =  c( "ggplot2", "dplyr", 'gridExtra', 'forecast',
               'tseries', 'aTSA', 'ggExtra', 'aTSA', 'knitr', 'fGarch', 'LSTS')

#install and load packages
install.load.packages(packages, libraryPath)
#load dataset
SP500Series = read.csv(paste(dir.path, "\\Data\\SP500.csv",
                             sep = ""), sep = ",",  na.strings = "null")

#-----------------------------------------------------------------------------------#
#                       Data cleaning + discriptive analysis
#-----------------------------------------------------------------------------------#

#characteristics of dataset
str(SP500Series)
#Factor to Date
SP500Series$Date = as.Date(SP500Series$Date)
#summary of the data
summary(SP500Series[-c(1)])
#Timeseries plot of stock price(s&p 500)
ggplot2::ggplot(data = SP500Series, aes(x = Date, y = Close))+ geom_line() + 
  theme_bw() + ggExtra::removeGridX() +ggExtra::removeGridY()+ 
  ggtitle("Figure 1: Timeseries plot of stock price(s&p 500)")+
  theme(plot.title = element_text(hjust = 0.5))
#log-return plot of s&p 500
sP500.price = SP500Series$Close
Log.Return = diff(log(sP500.price))
Log.Return.df = data.frame(Date = SP500Series$Date[2:length(SP500Series$Date)],
                           Log.Return = Log.Return )
ggplot2::ggplot(data = Log.Return.df, aes(x = Date, y = Log.Return)) + geom_line() + 
  ggtitle("Figure 2: log-return plot of s&p 500")+ ylab('log-return')+
  theme_bw() + ggExtra::removeGridX() +ggExtra::removeGridY()+ 
  theme(plot.title = element_text(hjust = 0.5))

#Histogram, density and normal distribution\n of log-return of s&p 500
x.Value = seq(min(Log.Return), max(Log.Return), length.out = length(Log.Return))
y.Value = dnorm(x.Value, mean(Log.Return), sd(Log.Return))
df = data.frame(Log.Return, x.Value, y.Value )
colnames(df) = c("Log.Return", "x.Value", "y.Value")
ggplot2::ggplot(df, aes(x = Log.Return)) + 
  ggtitle("Figure 3.1: Histogram, density and 
          normal distribution of log-return of s&p 500")+
  geom_histogram(aes(x = Log.Return, y = ..density..),
                 binwidth = 0.005, fill = "white", color = "black") + 
  geom_density(size = 0.8)+ xlab('log-return')+
  geom_line(data = df, aes(x = x.Value, y = y.Value), col = "red", size = 0.8)+
  xlim(-0.08, 0.08) + theme_bw() + ggExtra::removeGridX() +ggExtra::removeGridY()+ 
  theme(plot.title = element_text(hjust = 0.5))

#Test of stationarity: dicky-fuller test
#(see manually implemented function test.ADF in Function.R)
test.ADF(Log.Return, lag.order = 0)

#ACF PACF plot of stock price (s&p 500)
#(see manually implemented function Auto.cf in Function.R)
acf.sp_price = Auto.cf(sP500.price, lag = 20, type = 'acf', title = "ACF: SP500 Price")
acf.sp_price$acf.plot
pacf = Auto.cf(sP500.price, lag = 20, type = "pacf", title = "PACF: SP500 Price" )
pacf$pacf.plot

#ACF PACF plot of stock return (s&p 500)
#(see manually implemented function Auto.cf in Function.R)
acf.sp_return = Auto.cf(Log.Return, lag = 20, type = 'acf', title = "ACF: SP500 return")
acf.sp_return$acf.plot
pacf.sp_return = Auto.cf(Log.Return, lag = 20, type = "pacf", title = "PACF: SP500 return" )
pacf.sp_return$pacf.plot

#-----------------------------------------------------------------------------------#
#                                ARIMA-Model
#-----------------------------------------------------------------------------------#

#finding best model by parameter tuning
#(see manually implemented function arima,model in Function.R)
fit = list()
p = c(0,1,2,3,4,5) #max. AR order= 5
d = 0              #lag difference = 0
q = p              #max. MA order = 5
order = expand.grid(p,d,q)
colnames(order) = c("p", "d", "q")
order = order %>% dplyr::filter(p != 0 | q != 0)
order$AIC = NA
order$BIC = NA
for(i in 1:nrow(order)){
  fit = arima.model(series = Log.Return,
                    order = as.integer(order[i,c(1,2,3)]),
                    add.mean = TRUE, method = 'ML')
  order$AIC[i] = fit$aic
  order$BIC[i] = fit$bic
}

#MA and AR order, which minimizes AIC and BIC
min.AIC_BIC = order[sapply(order[, c(4,5)], which.min),] # min.AIC_BIC = c(0,0,2)
#Best model
best.Fit = arima.model(series = Log.Return,
                       order = c(0,0,2), add.mean = TRUE, method = 'ML')
#estimated coeficients and standars error
best.Fit$coeficients

#Forecasting ARIMA(0,0,2) model with estimated parameters
best.Fit = stats::arima(Log.Return, order = c(0,0,2))
arma.forecast <- forecast::forecast(best.Fit, h = 1000, level = 0.95)
autoplot(arma.forecast, main = "3.2: ARIMA forecasts for S&P 500 returns",
         shadecols  = c('grey'), xlab = "days", ylab = "S&P 500 Return" ) + 
  xlab('Time') + theme_bw() + ggExtra::removeGridX() +ggExtra::removeGridY()+
  theme(plot.title = element_text(hjust = 0.5))

#Diagonestic checking
#ACF plot of residuals of identified model
plot.acf = Auto.cf(series = best.Fit$residuals, lag = 20, type = 'acf',
                   title = 'Figure 4.1: ACF plot of residuals')
plot.acf$acf.plot

#Ljung box test for independency of error term
Test.box.ljung(series = best.Fit$residuals, lag = 10,
               alpha = 0.05, output = TRUE,  plot = TRUE,
               title = "Figure 4.2: Ljung Box test of error term")

#Histogram, empirical and normal density of error term
par(mfrow = c(1,2))
plot(best.Fit$residuals, ylab = "Residuals", main = "Figure 5: Residual plot")
hist(best.Fit$residuals, prob=TRUE, breaks = 200, 
     ylim = c(0,60), xlim = c(-0.05, 0.05),
     main = "Figure 6: Histogram", xlab = "residuals")     
lines(density(best.Fit$residuals), col = 'red')
lines(seq(-0.04, 0.06, by= 0.0005), dnorm(seq(-0.04, 0.06, by= 0.0005),
                                          mean(best.Fit$residuals), sd(best.Fit$residuals)), col="blue")

#------------------------------------------------------------------------#
#                           t-GARCH-MODEL
#------------------------------------------------------------------------#

#plot of squared residuals
par(mfrow = c(1,1))
plot(best.Fit$residuals^2, ylab = "residuals^2", ylim = c(0,0.01))
#there exists some heterocedasticity effect
#finding best t-GARCH-model model by parameter tuning
#(see manually implemented function garch.model in Function.R)
lag = 4
Order_Garch = data.frame(p = integer(lag), q = integer(lag), Log.Likelihood = double(lag),
                         AIC = double(lag), BIC = double(lag))
Order_Garch$p = c(1,2,1,2) 
Order_Garch$q = c(1,1,2,2)
GARCH = list()
for(i in 1:lag){
  GARCH[[i]] =  garch.model(arma.order = c(0,2),
                            garch.order = c(Order_Garch$p[i],Order_Garch$q[i]),
                            series = Log.Return, cond.dist = 'std')
  Order_Garch$Log.Likelihood[i] = GARCH[[i]]@LogLikelihood
  Order_Garch$AIC[i] = GARCH[[i]]@measure$aic
  Order_Garch$BIC[i] = GARCH[[i]]@measure$bic
}

minlgk.GARCH= Order_Garch %>% dplyr::filter(Log.Likelihood == max(Log.Likelihood))
min.AIC_BIC_GaRCH = Order_Garch[unique(as.numeric(sapply(Order_Garch[, c(4,5)], which.min))),]
best.order = dplyr::bind_rows(minlgk.GARCH, min.AIC_BIC_GaRCH)#p = 1, q =1 are selected by 
                                                              #parameter tuning
#best t-GARCH model
best.Fit.tGArch = garch.model(arma.order = c(0,2), 
                              garch.order = c(1,1),
                              series = Log.Return, cond.dist = 'std')
#Parameter + significance
best.Fit.tGArch@model.coef$matrix.coef #all coeficients are significant
#residual check 
#ACF plot
plot.acf = Auto.cf(series = best.Fit.tGArch@residuals,
                   lag = 20, type = 'acf', title = 'Figure 7: ACF plot of residuals')
plot.acf$acf.plot
#Residuals plot of t-GARCH model
plot(best.Fit.tGArch@residuals, ylab = "Residuals", 
     main = "Figure 8: Residual plot", type = 'l')
#QQ-plot of Residuals
Resid.tGarch = data.frame(Residuals =  best.Fit.tGArch@residuals)
df = as.integer(best.Fit.tGArch@df)
#Q-Q Plot
qq.tGarch <- ggplot(Resid.tGarch, aes(sample = Residuals))
qq.tGarch + stat_qq(distribution = stats::qt, dparams = df) +
  stat_qq_line(distribution = stats::qt, dparams = df, size = 0.8, col = "red")+
  ylim(-0.1, 0.1) + theme_bw()+ ggExtra::removeGridX() +ggExtra::removeGridY()+
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle('Figure 9: QQ plot of Residuals of t-Garch-model')

#Prediction of log return with identified t-GARCH model
best.Fit.tGarch = fGarch::garchFit(as.formula(paste("~arma(0,2) +",
                   "garch(",1,",",1,")",sep="")),data = Log.Return,
                    cond.dist = "std",include.shape = TRUE,  trace = FALSE)
fGarch:: predict(best.Fit.tGarch,n.ahead=150,plot=TRUE,conf=.95,nx=500)
