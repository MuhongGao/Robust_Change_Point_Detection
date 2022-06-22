library(KernSmooth)
library(MASS)
library(DNAcopy)
library(tilingArray)
library(cumSeg)

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
  right=rep(0,n+h)     #rightlimit for xx[1:(n+h)]
  for (i in 1:(n+h)){
    right[i]=locpoly(xx[i:(n+2*h)],yy[i:(n+2*h)],kernel="epanech",bandwidth=h,gridsize=length(xx[i:(n+2*h)]))$y[1]
  }
  left=rep(0,n+h)       #leftlimit for xx[(h+1):(n+2*h)]
  for (i in 1:(n+h)){
    left[i]=locpoly(xx[1:(h+i)],yy[1:(h+i)],kernel="epanech",bandwidth=h,gridsize=length(xx[1:(h+i)]))$y[length(xx[1:(h+i)])]
  }
  L=right[(h+1):(n+h)]-left[1:n]
  return(L)
}


DWB<-function(x,l,B){    #dependent wild bootstrap 
  n=length(x)
  V=matrix(0, nrow=n, ncol=n)
  for(i in 1:n){
    for(j in 1:n){
      V[i,j]=(1-abs(i-j)/l)*(abs(i-j)<l)                   #Bartlett kernel
    }
  }
  W=mvrnorm(n=B,mu=rep(0,n),Sigma=V)
  X=matrix(0,nrow=B,ncol=n)                 #show B bootstrap samples
  for(i in 1:B){
    X[i,]=x*W[i,]
  }
  return(X)
}

#######################################Simulation2
set.seed(18)
n=500
x=1:n
h=10
J=6
wave=0.1*(sin(2*pi*x/100)+2*sin(2*pi*x/60))
tau=(1:J)*80-30
########Case 1: iid errors without wave
N11=NULL            #estimated number of change points for method1 
N12=NULL
N13=NULL
N14=NULL
N15=NULL
FP11=NULL           #number of false positives for method1
FP12=NULL
FP13=NULL
FP14=NULL
FP15=NULL
for(j in 1:100){
  beta=sample(c(1,-1),size=J,replace=TRUE)
  signal=0
  for(i in 1:J) {signal<-signal+beta[i]*(x>tau[i])}
  y=signal+rnorm(n,mean=0,sd=0.2)
  ##############################CBS
  CBS=DNAcopy::segment(CNA(y, rep(1,n), 1:n))
  estimate=CBS$output[,4]
  estimate=estimate[-length(CBS$output[,4])]
  N11=c(N11,length(estimate))
  fp=0
  for(i in 1:length(estimate)){
    if(min(abs(estimate[i]-tau))>3)
      fp=fp+1
  }
  FP11=c(FP11,fp)
  ##############################DP
  model=tilingArray::segment(y, maxseg=40, maxk=n/2)
  estimate=model@breakpoints[[which.max(logLik(model, penalty="BIC"))]][,"estimate"]
  N12=c(N12,length(estimate))
  fp=0
  for(i in 1:length(estimate)){
    if(min(abs(estimate[i]-tau))>3)
      fp=fp+1
  }
  FP12=c(FP12,fp)
  ##############################cumSeg
  estimate=jumpoints(y,k=60,output="2")$psi
  N13=c(N13,length(estimate))
  fp=0
  for(i in 1:length(estimate)){
    if(min(abs(estimate[i]-tau))>3)
      fp=fp+1
  }
  FP13=c(FP13,fp)
  ##############################SaRa
  model=SARA(y,h=10)
  estimate=model$index[which(p.adjust(model$pV, "BH")<0.05)]
  N14=c(N14,length(estimate))
  fp=0
  for(i in 1:length(estimate)){
    if(min(abs(estimate[i]-tau))>3)
      fp=fp+1
  }
  FP14=c(FP14,fp)
  ##############################Proposed
  jump=0 
  L=screening(x,y,h)
  candidate=localMax(abs(L),span=2*h)
  for(k in 1:length(candidate)){
    jump<-jump+L[(candidate[k])]*(x>candidate[k])
  }
  hh=dpill(x,(y-jump))                #bandwidth selection
  yhat=locpoly(x,(y-jump), bandwidth=hh, gridsize=n)$y         #standard local linear regression
  residual=y-jump-yhat
  r=residual-mean(residual)       #centered residual
  rr=DWB(r,h,B=500)                 #bootstrap residual
  T=matrix(0,nrow=500,ncol=length(candidate))              #bootstrap statistics B=500
  segment=c(1,candidate[-length(candidate)]+floor(diff(candidate)/2),n)
  for(i in 1:500){
    yy=yhat+rr[i,]            #bootstrap y           
    LL=abs(screening(x,yy,h))
    for(k in 1:length(candidate)){
      T[i,k]=max(LL[(segment[k]:segment[k+1])])
    }
  }
  p=apply((T-matrix(rep(abs(L)[candidate],500),nrow=500,byrow=TRUE)>0),2,mean)
  estimate=candidate[which(p.adjust(p,"BH")<0.05)]
  N15=c(N15,length(estimate))
  fp=0
  for(i in 1:length(estimate)){
    if(min(abs(estimate[i]-tau))>3)
      fp=fp+1
  }
  FP15=c(FP15,fp)
}

########Case 2: correlated errors without wave
N21=NULL            #estimated number of change points for method1 
N22=NULL
N23=NULL
N24=NULL
N25=NULL
FP21=NULL           #number of false positives for method1
FP22=NULL
FP23=NULL
FP24=NULL
FP25=NULL
for(j in 1:100){
  beta=sample(c(1,-1),size=J,replace=TRUE)
  signal=0
  for(i in 1:J) {signal<-signal+beta[i]*(x>tau[i])}
  y=signal+arima.sim(n=600, list(ar = 0.6), sd = 0.16)[101:600]
  ##############################CBS
  CBS=DNAcopy::segment(CNA(y, rep(1,n), 1:n))
  estimate=CBS$output[,4]
  estimate=estimate[-length(CBS$output[,4])]
  N21=c(N21,length(estimate))
  fp=0
  for(i in 1:length(estimate)){
    if(min(abs(estimate[i]-tau))>3)
      fp=fp+1
  }
  FP21=c(FP21,fp)
  ##############################DP
  model=tilingArray::segment(y, maxseg=40, maxk=n/2)
  estimate=model@breakpoints[[which.max(logLik(model, penalty="BIC"))]][,"estimate"]
  N22=c(N22,length(estimate))
  fp=0
  for(i in 1:length(estimate)){
    if(min(abs(estimate[i]-tau))>3)
      fp=fp+1
  }
  FP22=c(FP22,fp)
  ##############################cumSeg
  estimate=jumpoints(y,k=60,output="2")$psi
  N23=c(N23,length(estimate))
  fp=0
  for(i in 1:length(estimate)){
    if(min(abs(estimate[i]-tau))>3)
      fp=fp+1
  }
  FP23=c(FP23,fp)
  ##############################SaRa
  model=SARA(y,h=10)
  estimate=model$index[which(p.adjust(model$pV, "BH")<0.05)]
  N24=c(N24,length(estimate))
  fp=0
  for(i in 1:length(estimate)){
    if(min(abs(estimate[i]-tau))>3)
      fp=fp+1
  }
  FP24=c(FP24,fp)
  ##############################Proposed
  jump=0 
  L=screening(x,y,h)
  candidate=localMax(abs(L),span=2*h)
  for(k in 1:length(candidate)){
    jump<-jump+L[(candidate[k])]*(x>candidate[k])
  }
  hh=dpill(x,(y-jump))                #bandwidth selection
  yhat=locpoly(x,(y-jump), bandwidth=hh, gridsize=n)$y         #standard local linear regression
  residual=y-jump-yhat
  r=residual-mean(residual)       #centered residual
  rr=DWB(r,h,B=500)                 #bootstrap residual
  T=matrix(0,nrow=500,ncol=length(candidate))              #bootstrap statistics B=500
  segment=c(1,candidate[-length(candidate)]+floor(diff(candidate)/2),n)
  for(i in 1:500){
    yy=yhat+rr[i,]            #bootstrap y           
    LL=abs(screening(x,yy,h))
    for(k in 1:length(candidate)){
      T[i,k]=max(LL[(segment[k]:segment[k+1])])
    }
  }
  p=apply((T-matrix(rep(abs(L)[candidate],500),nrow=500,byrow=TRUE)>0),2,mean)
  estimate=candidate[which(p.adjust(p,"BH")<0.05)]
  N25=c(N25,length(estimate))
  fp=0
  for(i in 1:length(estimate)){
    if(min(abs(estimate[i]-tau))>3)
      fp=fp+1
  }
  FP25=c(FP25,fp)
}

########Case 3: heteroscedastic errors+wave patterns
N31=NULL            #estimated number of change points for method1 
N32=NULL
N33=NULL
N34=NULL
N35=NULL
FP31=NULL           #number of false positives for method1
FP32=NULL
FP33=NULL
FP34=NULL
FP35=NULL
for(j in 1:100){
  beta=sample(c(1,-1),size=J,replace=TRUE)
  signal=0
  for(i in 1:J) {signal<-signal+beta[i]*(x>tau[i])}
  y=signal+wave+rnorm(n,mean=0,sd=sample(c(0.1,0.2,0.3),size=n,replace=TRUE))
  ##############################CBS
  CBS=DNAcopy::segment(CNA(y, rep(1,n), 1:n))
  estimate=CBS$output[,4]
  estimate=estimate[-length(CBS$output[,4])]
  N31=c(N31,length(estimate))
  fp=0
  for(i in 1:length(estimate)){
    if(min(abs(estimate[i]-tau))>3)
      fp=fp+1
  }
  FP31=c(FP31,fp)
  ##############################DP
  model=tilingArray::segment(y, maxseg=40, maxk=n/2)
  estimate=model@breakpoints[[which.max(logLik(model, penalty="BIC"))]][,"estimate"]
  N32=c(N32,length(estimate))
  fp=0
  for(i in 1:length(estimate)){
    if(min(abs(estimate[i]-tau))>3)
      fp=fp+1
  }
  FP32=c(FP32,fp)
  ##############################cumSeg
  estimate=jumpoints(y,k=60,output="2")$psi
  N33=c(N33,length(estimate))
  fp=0
  for(i in 1:length(estimate)){
    if(min(abs(estimate[i]-tau))>3)
      fp=fp+1
  }
  FP33=c(FP33,fp)
  ##############################SaRa
  model=SARA(y,h=10)
  estimate=model$index[which(p.adjust(model$pV, "BH")<0.05)]
  N34=c(N34,length(estimate))
  fp=0
  for(i in 1:length(estimate)){
    if(min(abs(estimate[i]-tau))>3)
      fp=fp+1
  }
  FP34=c(FP34,fp)
  ##############################Proposed
  jump=0 
  L=screening(x,y,h)
  candidate=localMax(abs(L),span=2*h)
  for(k in 1:length(candidate)){
    jump<-jump+L[(candidate[k])]*(x>candidate[k])
  }
  hh=dpill(x,(y-jump))                #bandwidth selection
  yhat=locpoly(x,(y-jump), bandwidth=hh, gridsize=n)$y         #standard local linear regression
  residual=y-jump-yhat
  r=residual-mean(residual)       #centered residual
  rr=DWB(r,h,B=500)                 #bootstrap residual
  T=matrix(0,nrow=500,ncol=length(candidate))              #bootstrap statistics B=500
  segment=c(1,candidate[-length(candidate)]+floor(diff(candidate)/2),n)
  for(i in 1:500){
    yy=yhat+rr[i,]            #bootstrap y           
    LL=abs(screening(x,yy,h))
    for(k in 1:length(candidate)){
      T[i,k]=max(LL[(segment[k]:segment[k+1])])
    }
  }
  p=apply((T-matrix(rep(abs(L)[candidate],500),nrow=500,byrow=TRUE)>0),2,mean)
  estimate=candidate[which(p.adjust(p,"BH")<0.05)]
  N35=c(N35,length(estimate))
  fp=0
  for(i in 1:length(estimate)){
    if(min(abs(estimate[i]-tau))>3)
      fp=fp+1
  }
  FP35=c(FP35,fp)
}

########Case 4: correlate errors+wave patterns
N41=NULL            #estimated number of change points for method1 
N42=NULL
N43=NULL
N44=NULL
N45=NULL
FP41=NULL           #number of false positives for method1
FP42=NULL
FP43=NULL
FP44=NULL
FP45=NULL
for(j in 1:100){
  beta=sample(c(1,-1),size=J,replace=TRUE)
  signal=0
  for(i in 1:J) {signal<-signal+beta[i]*(x>tau[i])}
  y=signal+wave+arima.sim(n=600, list(ar = 0.6), sd = 0.16)[101:600]
  ##############################CBS
  CBS=DNAcopy::segment(CNA(y, rep(1,n), 1:n))
  estimate=CBS$output[,4]
  estimate=estimate[-length(CBS$output[,4])]
  N41=c(N41,length(estimate))
  fp=0
  for(i in 1:length(estimate)){
    if(min(abs(estimate[i]-tau))>3)
      fp=fp+1
  }
  FP41=c(FP41,fp)
  ##############################DP
  model=tilingArray::segment(y, maxseg=40, maxk=n/2)
  estimate=model@breakpoints[[which.max(logLik(model, penalty="BIC"))]][,"estimate"]
  N42=c(N42,length(estimate))
  fp=0
  for(i in 1:length(estimate)){
    if(min(abs(estimate[i]-tau))>3)
      fp=fp+1
  }
  FP42=c(FP42,fp)
  ##############################cumSeg
  estimate=jumpoints(y,k=60,output="2")$psi
  N43=c(N43,length(estimate))
  fp=0
  for(i in 1:length(estimate)){
    if(min(abs(estimate[i]-tau))>3)
      fp=fp+1
  }
  FP43=c(FP43,fp)
  ##############################SaRa
  model=SARA(y,h=10)
  estimate=model$index[which(p.adjust(model$pV, "BH")<0.05)]
  N44=c(N44,length(estimate))
  fp=0
  for(i in 1:length(estimate)){
    if(min(abs(estimate[i]-tau))>3)
      fp=fp+1
  }
  FP44=c(FP44,fp)
  ##############################Proposed
  jump=0 
  L=screening(x,y,h)
  candidate=localMax(abs(L),span=2*h)
  for(k in 1:length(candidate)){
    jump<-jump+L[(candidate[k])]*(x>candidate[k])
  }
  hh=dpill(x,(y-jump))                #bandwidth selection
  yhat=locpoly(x,(y-jump), bandwidth=hh, gridsize=n)$y         #standard local linear regression
  residual=y-jump-yhat
  r=residual-mean(residual)       #centered residual
  rr=DWB(r,h,B=500)                 #bootstrap residual
  T=matrix(0,nrow=500,ncol=length(candidate))              #bootstrap statistics B=500
  segment=c(1,candidate[-length(candidate)]+floor(diff(candidate)/2),n)
  for(i in 1:500){
    yy=yhat+rr[i,]            #bootstrap y           
    LL=abs(screening(x,yy,h))
    for(k in 1:length(candidate)){
      T[i,k]=max(LL[(segment[k]:segment[k+1])])
    }
  }
  p=apply((T-matrix(rep(abs(L)[candidate],500),nrow=500,byrow=TRUE)>0),2,mean)
  estimate=candidate[which(p.adjust(p,"BH")<0.05)]
  N45=c(N45,length(estimate))
  fp=0
  for(i in 1:length(estimate)){
    if(min(abs(estimate[i]-tau))>3)
      fp=fp+1
  }
  FP45=c(FP45,fp)
}

####################################show table
cbind(
  c(mean(N11),mean(N12),mean(N13),mean(N14),mean(N15)),
  c(mean(FP11),mean(FP12),mean(FP13),mean(FP14),mean(FP15)),
  c(mean(N21),mean(N22),mean(N23),mean(N24),mean(N25)),
  c(mean(FP21),mean(FP22),mean(FP23),mean(FP24),mean(FP25)),
  c(mean(N31),mean(N32),mean(N33),mean(N34),mean(N35)),
  c(mean(FP31),mean(FP32),mean(FP33),mean(FP34),mean(FP35)),
  c(mean(N41),mean(N42),mean(N43),mean(N44),mean(N45)),
  c(mean(FP41),mean(FP42),mean(FP43),mean(FP44),mean(FP45)))