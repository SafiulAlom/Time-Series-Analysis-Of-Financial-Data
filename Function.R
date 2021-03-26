

#*********************************Function.R******************************************#

#-------------------------------------------------------------------------------------#
#                     1. install.load.packages(pack, packagePath)
#-------------------------------------------------------------------------------------#

#function for installing and loading required packages
install.load.packages = function(pack, packagePath){
  all.packages = pack[!(pack %in% installed.packages(lib.loc = packagePath)[, 1])]
  if (length(all.packages)) 
    install.packages(all.packages, lib = packagePath, dependencies = TRUE)
  sapply(pack, require, lib.loc = packagePath, character.only = TRUE)
}


#-------------------------------------------------------------------------------------#
#                          2. NAtoNearestValue(a)
#-------------------------------------------------------------------------------------#

#function for replacing NA's by nearest value
NAtoNearestValue= function(a){
  
  if(is.na(a[1])) a[1] = a[min(which(!is.na(a)))] #if the 1st obs. = NA, take the nearest value != 0 
  for(i in 1:length(a)){
    if(is.na(a[i])) {a[i] = a[i-1]} #Last Observation Carried Forward
    else a[i] = a[i]
  }
  return(a) 
}


#-------------------------------------------------------------------------------------#
#   3.  Auto.cf(series, lag, type = c('acf', 'pacf'), title = c('Acf/pacf plot'))
#-------------------------------------------------------------------------------------#

#Autocorrelation function
Auto.cf = function(series, lag, type = c('acf', 'pacf'), title = c('Acf/pacf plot')){
  type = match.arg(type)
  if(lag >= length(series & lag < 0)){
    stop("Lag can not be larger than or equal to the length of the Series")
  }
  if(NCOL(series) > 1 || is.data.frame(series)){stop('series is not a vector')}
  if(any(is.na(series))){stop('series has NA values')}
  
  CI.upper = qnorm((1 + 0.95)/2)/sqrt(length(series))
  CI.lower =-qnorm((1 + 0.95)/2)/sqrt(length(series))
  #ACF
  if(type =='acf'){
    roh = vector(mode = 'numeric', length = lag)
    for(i in 1:lag)
      roh[i] = (1/((length(series)-1)*sd(series)^2))*
        sum((series[(1+i):length(series)] - mean(series))*(series[1:(length(series)-i)] - mean(series)))
    #ACF plot
    acf.df = data.frame(index = 1:length(roh), roh)
    acf.plot = ggplot(data = acf.df, aes(x = index, y = roh)) + xlab('lag') + ylab('acf')+
      geom_bar(stat="identity", position = 'dodge', width = 0.5)+
      theme_bw() + ggExtra::removeGridX() +ggExtra::removeGridY()+
      geom_hline(aes(yintercept = 0)) +
      theme(plot.title = element_text(hjust = 0.5))+
      ggtitle(title)+ geom_hline(aes(yintercept = CI.upper), color = 'blue',linetype="dashed")+
      geom_hline(aes(yintercept = CI.lower), color = 'blue',linetype="dashed")
    
    return(list(lag = lag, acf = roh, acf.plot = acf.plot,
                n.obs = paste("number of obs.:", length(series))))
  }
  
  else if(type == 'pacf'){
    #PACF
    pacf = stats::pacf(series,lag = lag, plot = FALSE)
    roh = pacf$acf
    acf.df = data.frame(index = 1:length(roh), roh)
    acf.plot = ggplot(data = acf.df, aes(x = index, y = roh)) + xlab('lag') + ylab('pacf')+
      geom_bar(stat="identity", position = 'dodge', width = 0.5)+
      theme_bw() + ggExtra::removeGridX() +ggExtra::removeGridY()+
      geom_hline(aes(yintercept = 0))+ theme(plot.title = element_text(hjust = 0.5))+
      ggtitle(title) + geom_hline(aes(yintercept = CI.upper), color = 'blue',linetype="dashed")+
      geom_hline(aes(yintercept = CI.lower), color = 'blue',linetype="dashed")
    
    return(list(lag = lag, pacf = roh, pacf.plot = acf.plot,
                n.obs = paste("number of obs.:", length(series))))
    }
}


#-------------------------------------------------------------------------------------#
#                      4. test.ADF(series, lag.order)
#-------------------------------------------------------------------------------------#

#Augmented dickey fuller test
test.ADF = function(series, lag.order){
  
  if(lag.order<0){stop('lag order must be non-negative')}
  if(NCOL(series) > 1 || is.data.frame(series)){stop('series is not a vector')}
  if(any(is.na(series))){stop('series has NA values')}
  
  lag.order.adj = lag.order + 1
  y = diff(series)
  n.obs = length(y)
  x1 = series[lag.order.adj:n.obs]
  x2 = rep(1, length(lag.order.adj:n.obs))
  x3 = lag.order.adj:n.obs
  y.series = embed(y, lag.order.adj)[,1]
  val.0.01 = -c(4.38, 4.15, 4.04, 3.99, 3.98, 3.96)
  val.0.025 = -c(3.95, 3.8, 3.73, 3.69, 3.68, 3.66)
  val.0.05 = -c(3.6, 3.5, 3.45, 3.43, 3.42, 3.41)
  val.0.1 = -c(3.24, 3.18, 3.15, 3.13, 3.13, 3.12)
  val.0.9 = -c(1.14, 1.19, 1.22, 1.23, 1.24, 1.25)
  val.0.95= -c(0.8, 0.87, 0.9, 0.92, 0.93, 0.94)
  val.0.975 = -c(0.5, 0.58, 0.62, 0.64, 0.65, 0.66)
  val.0.99 = -c(0.15, 0.24, 0.28, 0.31, 0.32, 0.33)
  crit.val = cbind(val.0.01, val.0.025, val.0.05,
                   val.0.1, val.0.9, val.0.95, val.0.975, val.0.99)
  col = ncol(crit.val)
  Obs.size = c(25, 50, 100, 250, 500, 10^5)
  p.val = c(0.01, 0.025, 0.05, 0.1, 0.9, 0.95, 0.975, 0.99)
  
  if(lag.order.adj == 1){
    mod = lm(y.series ~ x1 + x2 + x3)
  }
  else {
    x4 = embed(y, lag.order.adj)[,2:lag.order.adj]
    mod = lm(y.series ~ x1 + x2 + x3 + x4)}
  
  beta2 = coef(mod)[2]
  sd.beta2 = sqrt(diag(vcov(mod)))[2]
  test.stat = beta2/sd.beta2
  
  est.crit.val = lapply( as.data.frame(crit.val), stats::approx, 
                         x = Obs.size,  xout = n.obs, rule = 2)
  est.crit.val = as.numeric(unlist(lapply(est.crit.val, function(x) x[2])))
  est.p.val = approx(est.crit.val,p.val, test.stat, rule = 2)$y
  names(test.stat) <- "Dicky-Fuller Test-Statistic"
  if(est.p.val == 0.01){warning("p.value = 0.01 means p.value <= 0.01")}
  
  structure(list(data.name = deparse(substitute(series)),
                 statistic = test.stat, p.value = est.p.val,
                 method = "ADF-test:", alternative = "stationary"),
            class = "htest")
}

#-------------------------------------------------------------------------------#
#   5. arima.model(series, order = c(1,0,1),
#                 add.mean = TRUE, method = c("CSS-ML", "ML", "CSS"))
#-------------------------------------------------------------------------------#

#arima model
arima.model = function(series, order = c(1,0,1), add.mean = TRUE, method = c("CSS-ML", "ML", "CSS")){
  fit = stats::arima(x = series, order= order, include.mean = add.mean, method = method)
  n = length(series)
  k = length(coef(fit))
  rss = sum(fit$residuals^2)
  if(method %in% c("CSS-ML", "ML")){
    aic = n*(log(2*pi) + 1 + log((rss/n))) + ((k + 1)*2)
    bic = n + n*log(2*pi) + n*log(rss/n) + log(n)*(k + 1) #4 = k + 1
    loglik = -0.5*(aic - ((k + 1)*2))
  }
  sigma.sq = rss/(n-k)
  coef.fit =rbind(coef(fit), sqrt(diag(vcov(fit))))
  rownames(coef.fit) = c('coeficient', 's.e.')
  
  ME = mean(fit$residuals)
  MSE = rss/n
  RMSE = sqrt(MSE)
  MAE = sum(abs(fit$residuals))/n
  MAE.adj = (1/(n-1))*sum(abs(series[2:n] - series[1:(n-1)]))
  MASE = MAE/MAE.adj
  var.coef = fit$var.coef
  error.measure = c(ME, MSE, RMSE, MAE, MAE.adj, MASE)
  names(error.measure) = c('ME', 'MSE', 'RMSE', 'MAE', 'MAE.adj', 'MASE')
  
  output = structure(list(variance.coeficient = var.coef, series = series,
                coeficients = coef.fit, call = match.call(), 
                residuals = fit$residuals, loglik = loglik, error.measure = error.measure,
                nobs = n, sigma.squared = sigma.sq, aic = aic, bic = bic), class = "Arima")
  return(output)
}

#--------------------------------------------------------------------------------#
#        6. Test.box.ljung = function(series, lag, alpha = 0.05,
#                          output = TRUE, plot = FALSE, title = "Ljung Box test")
#--------------------------------------------------------------------------------#

#Ljung–Box test
Test.box.ljung = function(series, lag, alpha = 0.05, output = TRUE, plot = FALSE, title = "Ljung Box test"){
  n = length(series)
  k = lag
  p.value = vector()
  for(i in 1:k){
    Q = Auto.cf(series, i)
    Q = Q$acf^2
    d = (n-1):(n-i)
    t.stat = n*(n+2)*sum(Q/d)
    p.value[i] = 1 - pchisq(t.stat , df = i)
  }
  if(plot){
    plot(x = 1:k, y = p.value, xlab = "lag", ylim =c(0, 1),
         ylab = "p-value", las = 1, pch = 16,bty = "n", main = title)
    abline(h = alpha, col = 'red', lty = 2, lwd = 2)}
  df = data.frame(lag = 1:k, p.value = p.value)
  if(output){
  return(list(Box.Ljung.Test = df, data = series, alpha = alpha))}
}

#--------------------------------------------------------------------------------#
#     7. garch.model(series, arma.order = c(0,2), 
#               garch.order = c(1,0), cond.dist = 'norm',include.shape = NULL)
#--------------------------------------------------------------------------------#

#ARMA-GARCH model
garch.model = function(series, arma.order = c(0,2), garch.order = c(1,0), 
                       cond.dist = 'norm',
                       include.shape = NULL){
  
  fit = fGarch::garchFit(as.formula(paste("~arma(",arma.order[1],
                                          ",",arma.order[2],") +","garch(",garch.order[1],",",
                                          garch.order[2],")",sep="")),
                         data = series, trace = FALSE, include.shape = include.shape, cond.dist = cond.dist)
  coef = fit@fit$coef
  matrix.coef = fit@fit$matcoef
  covariance = fit@fit$cvar
  model.coef = list(coef, matrix.coef, covariance)
  names(model.coef) = c("coeficients", 'matrix.coef', 'covariance')
  residuals = fit@residuals
  method = fit@method
  call = match.call()
  title = fit@title
  description = fit@description
  data = series
  LogLikelihood = fit@fit$value
  n = length(series)
  rss = sum(fit@residuals^2)
  ME = mean(fit@residuals)
  MSE = rss/n
  RMSE = sqrt(MSE)
  MAE = sum(abs(fit@residuals))/n
  MAE.adj = (1/(n-1))*sum(abs(series[2:n] - series[1:(n-1)]))
  MASE = MAE/MAE.adj
  error.measure = c(ME, MSE, RMSE, MAE, MAE.adj, MASE)
  names(error.measure) = c('ME', 'MSE', 'RMSE', 'MAE', 'MAE.adj', 'MASE')
  aic = fit@fit$ics[1]
  bic = fit@fit$ics[2]
  df = fit@fit$coef["shape"]
  measure = list(error.measure, aic, bic)
  names(measure) = c("error.measures", "aic", 'bic')
  setClass("track", slots = c(model.coef="list", measure = "list", residuals ="numeric",
                              method = 'character', call = 'call',
                              title = 'character', description = 'character',
                              data = 'numeric', LogLikelihood = 'numeric', df = 'numeric'
                              ))
  model = new('track', model.coef = model.coef, measure = measure, residuals = residuals,
              method = method, call = call, title = title, description = description,
              data = data, LogLikelihood = LogLikelihood, df = df)
  return(model)
}
