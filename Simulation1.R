library(KernSmooth)
library(Matrix)
library(MASS)

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

screening<-function(x, y, h){       #local linear smoothing (rightlimit-leftlimit)
  n=length(x)
  xx=1:(n+2*h)
  yy=c(rep(y[1],h),y,rep(y[n],h))       #create data outside the boundaries
  right=rep(0,n)     #rightlimit for xx[1:n]
  left=rep(0,n)       #leftlimit for xx[(1+2*h):(n+2*h)]
  for (i in 1:n){
    model=locpoly(xx[i:(i+2*h)],yy[i:(i+2*h)],kernel="epanech",bandwidth=h,gridsize=1+2*h)
    right[i]=model$y[1]
    left[i]=model$y[1+2*h]
  }
  L=c(rep(0,h),right[(2*h+1):n]-left[1:(n-2*h)],rep(0,h))
  return(L)
}

DWB<-function(x,l,B=500){    #dependent wild bootstrap 
  n=length(x)
  V=as(as(diag(0,n), "diagonalMatrix"), "CsparseMatrix")
  for (i in 1:n){
    for (j in 1:n){
      if (abs(i-j)<l)
        V[i,j]=(1-abs(i-j)/l)            #Bartlett kernel
    }
  }
  model=Cholesky(V, perm = TRUE, super = FALSE, Imult = 0)
  L=as(model,"Matrix")
  P=t(as(model, "pMatrix"))
  X=matrix(0,nrow=B,ncol=n)                 #show B bootstrap samples
  for(i in 1:B){
    X[i,]=x*(P%*%(L%*%rnorm(n,mean=0,sd=1)))
  }
  return(X)
}

Pvalue<-function(x, y, h, candidate, B=500){
  n=length(x)
  jump=0 
  for(k in 1:length(candidate)){
    jump<-jump+L[(candidate[k])]*(x>candidate[k])
  }
  hh=dpill(x,(y-jump))                #bandwidth selection
  yhat=locpoly(x,(y-jump), bandwidth=hh, gridsize=n)$y         #standard local linear regression
  residual=y-jump-yhat
  r=residual-mean(residual)       #centered residual
  rr=DWB(r,h,B)                 #bootstrap residual
  T=matrix(0,nrow=B,ncol=length(candidate))              #bootstrap statistics 
  segment=c(1,candidate[-length(candidate)]+floor(diff(candidate)/2),n)
  for(j in 1:B){
    yy=yhat+rr[j,]            #bootstrap y           
    LL=abs(screening(x,yy,h))
    for(k in 1:length(candidate)){
      T[j,k]=max(LL[(segment[k]:segment[k+1])])
    }
  }
  p=apply((T-matrix(rep(abs(L)[candidate],B),nrow=B,byrow=TRUE)>0),2,mean)
  return(p)
}

#######################################Simulation1
set.seed(18)
n=200
x=1:n
h=10
signal=1*(x>=40)-1*(x>=160)
wave=0.1*(sin(2*pi*x/100)+2*sin(2*pi*x/60))
############Scenario I: iid errors, no wave
estimate1=NULL       #SaRa algorithm
Estimate1=NULL       #proposed algorithm
n1=NULL        #number of change points by SaRa algorithm
N1=NULL       #number of change points by proposed algorithm
for(i in 1:200){
  e=rnorm(n=200,mean=0,sd=0.2)
  y=signal+e
  model=SARA(y, h)           #SaRa algorithm
  estimate=model$index[which(p.adjust(model$pV, "BH")<0.05)]
  estimate1=c(estimate1, estimate)  
  n1=c(n1,length(estimate))
  L=screening(x,y,h)
  candidate=localMax(abs(L),span=2*h)
  p=Pvalue(x,y,h,candidate)
  estimate=candidate[which(p.adjust(p,"BH")<0.05)]
  Estimate1=c(Estimate1, estimate) 
  N1=c(N1,length(estimate))
}
############Scenario II: correlated errors, no wave
estimate2=NULL
Estimate2=NULL
n2=NULL
N2=NULL
for(i in 1:200){
  e=arima.sim(n=300, list(ar = 0.6), sd = 0.16)[101:300]
  y=signal+e  
  model=SARA(y, h)           
  estimate=model$index[which(p.adjust(model$pV, "BH")<0.05)]
  estimate2=c(estimate2, estimate) 
  n2=c(n2,length(estimate)) 
  L=screening(x,y,h)
  candidate=localMax(abs(L),span=2*h)
  p=Pvalue(x,y,h,candidate)
  estimate=candidate[which(p.adjust(p,"BH")<0.05)]
  Estimate2=c(Estimate2, estimate) 
  N2=c(N2,length(estimate)) 
}
############Scenario III: iid errors, with waves
estimate3=NULL
Estimate3=NULL
n3=NULL
N3=NULL
for(i in 1:200){
  e=rnorm(n=200,mean=0,sd=0.2)
  y=signal+wave+e  
  model=SARA(y, h)           
  estimate=model$index[which(p.adjust(model$pV, "BH")<0.05)]
  estimate3=c(estimate3, estimate)
  n3=c(n3,length(estimate))  
  L=screening(x,y,h)
  candidate=localMax(abs(L),span=2*h)
  p=Pvalue(x,y,h,candidate)
  estimate=candidate[which(p.adjust(p,"BH")<0.05)]
  Estimate3=c(Estimate3, estimate) 
  N3=c(N3,length(estimate)) 
}
############Scenario IV: correlated errors, with waves
estimate4=NULL
Estimate4=NULL
n4=NULL
N4=NULL
for(i in 1:200){
  e=arima.sim(n=300, list(ar = 0.6), sd = 0.16)[101:300]
  y=signal+wave+e  
  model=SARA(y, h)           
  estimate=model$index[which(p.adjust(model$pV, "BH")<0.05)]
  estimate4=c(estimate4, estimate)  
  n4=c(n4,length(estimate)) 
  L=screening(x,y,h)
  candidate=localMax(abs(L),span=2*h)
  p=Pvalue(x,y,h,candidate)
  estimate=candidate[which(p.adjust(p,"BH")<0.05)]
  Estimate4=c(Estimate4, estimate) 
  N4=c(N4,length(estimate)) 
}

#######################################show table
rbind(c(sum(n1==0),sum(n1==1),sum(n1==2),sum(n1==3),sum(n1>3),mean(n1)),
      c(sum(N1==0),sum(N1==1),sum(N1==2),sum(N1==3),sum(N1>3),mean(N1)),
      c(sum(n2==0),sum(n2==1),sum(n2==2),sum(n2==3),sum(n2>3),mean(n2)),
      c(sum(N2==0),sum(N2==1),sum(N2==2),sum(N2==3),sum(N2>3),mean(N2)),
      c(sum(n3==0),sum(n3==1),sum(n3==2),sum(n3==3),sum(n3>3),mean(n3)),
      c(sum(N3==0),sum(N3==1),sum(N3==2),sum(N3==3),sum(N3>3),mean(N3)),
      c(sum(n4==0),sum(n4==1),sum(n4==2),sum(n4==3),sum(n4>3),mean(n4)),
      c(sum(N4==0),sum(N4==1),sum(N4==2),sum(N4==3),sum(N4>3),mean(N4)))
#######################################show figure
par(mfrow=c(4,2))
hist(estimate1,breaks=50,xlim=c(0,200),ylim=c(0,200),xlab="locations",main="SaRa (Scenario I)")
hist(Estimate1,breaks=50,xlim=c(0,200),ylim=c(0,200),xlab="locations",main="Proposed (Scenario I)")
hist(estimate2,breaks=50,xlim=c(0,200),ylim=c(0,200),xlab="locations",main="SaRa (Scenario II)")
hist(Estimate2,breaks=50,xlim=c(0,200),ylim=c(0,200),xlab="locations",main="Proposed (Scenario II)")
hist(estimate3,breaks=50,xlim=c(0,200),ylim=c(0,200),xlab="locations",main="SaRa (Scenario III)")
hist(Estimate3,breaks=50,xlim=c(0,200),ylim=c(0,200),xlab="locations",main="Proposed (Scenario III)")
hist(estimate4,breaks=50,xlim=c(0,200),ylim=c(0,200),xlab="locations",main="SaRa (Scenario IV)")
hist(Estimate4,breaks=50,xlim=c(0,200),ylim=c(0,200),xlab="locations",main="Proposed (Scenario IV)")