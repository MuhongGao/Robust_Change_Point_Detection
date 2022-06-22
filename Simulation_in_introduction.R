###########################define some functions 
#Estimate the standard deviation of the intensities
estimateSigma<-function (Y, h = 10) {  
  n = length(Y)
  YBar = rep(0, n)
  for (i in 1:n) {
    a = min(n, i + h)
    b = max(1, i - h)
    YBar[i] = mean(Y[b:a])
  }
  return(sqrt(var(Y - YBar) * (2 * h + 1)/(2 * h)))
}

#Calculate the value for local diagnostic function
localDiagnostic<-function (y, h) { 
  yy = c(rep(0, h - 1), y, rep(0, h))
  n = length(y)
  z = rep(0, n)
  for(i in 1:n){
    z[i]=sum(yy[i:(h+i-1)])/h-sum(yy[(h+i):(2*h-1+i)])/h
  }
  return(z)
}

#Get the local maximizers of local diagnostic function
localMax<-function (y, span = 5) {  
  if (length(y) < span * 2 + 1) 
    return(NULL)
  n = length(y)
  index = NULL
  for (i in (span + 1):(n - span)) {
    if (y[i] == max(y[(i - span):(i + span)])) 
      index = c(index, i)
  }
  return(index)
}

clean<-function (LocalM, h) 
{
  len <- length(LocalM)
  rm.list <- NULL
  for (i in 1:(len - 1)) {
    if (LocalM[i] >= LocalM[i + 1] - h) {
      rm.list <- c(rm.list, i)
    }
  }
  if (length(rm.list) > 0) {
    LocalM <- LocalM[-as.vector(rm.list)]
  }
  return(LocalM = LocalM)
}

SARAp<-function (Y, h, hh = 2 * h, sigma = NULL) { 
  n = length(Y)
  LDF = localDiagnostic(Y, h)
  LDF.pos = LDF
  LDF.neg = -LDF
  if (is.null(sigma)) 
    sigma = estimateSigma(Y, h = max(3, 2 * floor(log(n))))
  pV.pos = 1 - 2 * pnorm(LDF.pos/(sqrt(2/h) * sigma))
  LocalMax = localMax(LDF.pos, span = hh)
  LocalMax = clean(LocalMax, h)
  LocalMaxValue = pV.pos[LocalMax]
  pV.neg = 1 - 2 * pnorm(LDF.neg/(sqrt(2/h) * sigma))
  LocalMin = localMax(LDF.neg, span = hh)
  LocalMin = clean(LocalMin, h)
  LocalMinValue = pV.neg[LocalMin]
  LocalExt <- c(LocalMax, LocalMin)
  LocalExtValue <- c(LocalMaxValue, LocalMinValue)
  LocalExtValue <- LocalExtValue[order(LocalExt)]
  LocalExt <- sort(LocalExt)
  return(list(index = LocalExt, pV = LocalExtValue))
}

#Get the inverse cumulative distribution function of local min p-values
fInverse<-function (n = 10000, h = 10, hh = 2 * h, precise = 10000, simT = 100) { 
  empirical = NULL
  for (i in 1:simT) {
    Y = rnorm(n)
    LDF = localDiagnostic(Y, h)
    LDF.pos = LDF
    LDF.neg = -LDF
    sigma = 1
    index.pos = localMax(y = LDF.pos, span = hh)
    pV.pos = 1 - 2 * pnorm(LDF.pos[index.pos]/(sqrt(2/h) * 
                                                 sigma))
    index.neg = localMax(y = LDF.neg, span = hh)
    pV.neg = 1 - 2 * pnorm(LDF.neg[index.neg]/(sqrt(2/h) * 
                                                 sigma))
    index <- c(index.pos, index.neg)
    pv <- c(pV.pos, pV.neg)
    pv <- pv[order(index)]
    index <- sort(index)
    len <- length(index)
    rm.list <- NULL
    for (j in 1:(len - 1)) {
      if (index[j] >= index[j + 1] - h) {
        rm.list <- c(rm.list, j)
      }
    }
    if (length(rm.list) > 0) {
      pv <- pv[-rm.list]
    }
    empirical <- c(empirical, pv)
    if (length(empirical) > 10 * precise) 
      break
  }
  return(quantile(empirical, probs = c(0:precise)/precise))
}

SARA<-function (Y, h = 10, hh = 2 * h, FINV = NULL, sigma = NULL, precise = 10000) {
  object = SARAp(Y = Y, h = h, hh = hh, sigma = sigma)
  index = object$index
  pV = object$pV
  if (is.null(FINV)) 
    FINV = fInverse(n = length(Y), h = h, hh = hh, precise = precise, 
                    simT = 100)
  pVcorr = pV
  for (i in 1:length(pV)) {
    pVcorr[i] = (length(which(FINV < pV[i])))/(precise + 
                                                 1)
  }
  return(list(index = index, pV = pVcorr))
}


#######################Simulation for Figure 1
set.seed(23)
n=200
x=1:200
signal=0.6*(x>=40)-0.4*(x>=160)
wave=0.05*(sin(2*pi*x/96)+2*sin(2*pi*x/240))
estimate1=NULL
estimate2=NULL
for(i in 1:500){
  e1=rnorm(n=200,mean=0,sd=0.2)
  y1=signal+e1
  e2=arima.sim(n=300, list(ar = 0.6), sd = 0.16)[101:300]
  y2=signal+wave+e2             
  model1=SARA(y1, h = 10)           #SaRa algorithm
  estimate1=c(estimate1, model1$index[which(p.adjust(model1$pV,"BH")<0.05)])      
  model2=SARA(y2, h = 10)
  estimate2=c(estimate2, model2$index[which(p.adjust(model2$pV,"BH")<0.05)])
}
D=abs(localDiagnostic(y1, h=10))     #scan statistic

par(mfrow=c(2,2))   #show plots
plot(x,signal,type="l", lwd=2, ylim=c(0,0.8), xlab="locations", ylab="signal", main="true signal")
abline(v=c(40,160), lty=2)
plot(x,D,type="l", xlab="locations", ylab="D", main="scan statistic")
abline(v=c(40,160), lty=2)
hist(estimate1,breaks=50,xlim=c(0,200),ylim=c(0,500),xlab="locations",main=paste("histogram of change points estimators","with i.i.d. errors",sep="\n"))
abline(v=c(40,160), lty=2)
hist(estimate2,breaks=50,xlim=c(0,200),ylim=c(0,500),xlab="locations",main=paste("histogram of change points estimators","with waves and correlated errors",sep="\n"))
abline(v=c(40,160), lty=2)