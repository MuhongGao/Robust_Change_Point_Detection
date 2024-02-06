library(KernSmooth)
library(MASS)
library(Matrix)
library(DNAcopy)
library(tilingArray)
library(cumSeg)
library(mgcv)
library(imputeTS)
library(mvtnorm)
library(np)

##############################Real Data
CGHdata <- read.csv("C:/Users/Acer/Desktop/My Document/Documents/Research/Projects/change points(bootstrap)/CGHdataset.csv", header=TRUE)
index=c(6,19,20,21)
data=CGHdata[2:2301,(1+3*index)]               
n=nrow(data)
d=ncol(data)
for (i in 1:d){
  data[,i]=na_ma(data[,i],k=5,weighting="linear")        #imputation for methods 1-3
}

###########################define some functions 
estimateSigma<-function (Y, h = 10) {       #Estimate the standard deviation of the intensities
  n = length(Y)
  YBar = rep(0, n)
  for (i in 1:n) {
    a = min(n, i + h)
    b = max(1, i - h)
    YBar[i] = mean(Y[b:a])
  }
  return(sqrt(var(Y - YBar) * (2 * h + 1)/(2 * h)))
}

localMax<-function (y, span = 5) {           #Get the local maximizers of local diagnostic function  
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

screening<-function(x, y, h){       #local linear smoothing (rightlimit-leftlimit)
  n=length(x)
  xx=c(2*x[1]-x[sort(which(x<(x[1]+h))[-1],decreasing=T)], x, 2*x[n]-x[sort(which(x>(x[n]-h)),decreasing=T)[-1]])
  yy=c(2*y[1]-y[sort(which(x<(x[1]+h))[-1],decreasing=T)], y, 2*y[n]-y[sort(which(x>(x[n]-h)),decreasing=T)[-1]])       #create data outside the boundaries
  rkernal<-function(t){
    0.75*(1-((xx-t)/h)^2)*(((xx-t)/h)>0)*(((xx-t)/h)<=1)
  }
  lkernal<-function(t){
    0.75*(1-((xx-t)/h)^2)*(((xx-t)/h)>=-1)*(((xx-t)/h)<=0)
  }
  rweight<-function(t){
    rkernal(t)*(sum(rkernal(t)*(xx-t)^2)-(xx-t)*sum(rkernal(t)*(xx-t)))
  }
  lweight<-function(t){
    lkernal(t)*(sum(lkernal(t)*(xx-t)^2)-(xx-t)*sum(lkernal(t)*(xx-t)))
  }
  rlimit<-function(t){
    sum(rweight(t)*yy)/sum(rweight(t))
  }
  llimit<-function(t){
    sum(lweight(t)*yy)/sum(lweight(t))
  }
  right=rep(0,n)     #rightlimit for x[1:n]
  left=rep(0,n)       #leftlimit for x[1:n]
  for (i in 1:n){
    right[i]=rlimit(x[i])
    left[i]=llimit(x[i])
  }
  D=abs(right-left)
  return(D)
}

refine1<-function(x,y,h,candidate){       #step 2 (proposed I)
  for(k in 1:5){
    Candidate=c(0,candidate,1)
    a=sample(1:length(candidate))
    for(j in 1:length(a)){
      index=which((x<Candidate[(a[j]+2)])*(x>Candidate[a[j]])==1)        #subsequence
      xsub=x[index]
      ysub=y[index]
      b=which((x<(candidate[a[j]]+h))*(x>(candidate[a[j]]-h))==1)          #neighborhood of candidate
      RSS=rep(0, length(b))
      for(i in 1:length(b)){
        z=(xsub>x[b[i]])+1-1
        model=gam(ysub~z+s(xsub), family=gaussian, method="GCV.Cp")
        RSS[i]=sum((model$residuals)^2)
      }
      candidate[a[j]]=x[b[which.min(RSS)]]
    }
  }
  return(candidate)
}

refine2<-function(x,y,h,candidate){       #step 2 (proposed II)
  for(k in 1:5){
    Candidate=c(0,candidate,1)
    a=sample(1:length(candidate))
    for(j in 1:length(a)){
      index=which((x<Candidate[(a[j]+2)])*(x>Candidate[a[j]])==1)        #subsequence
      xsub=x[index]
      ysub=y[index]
      b=which((x<(candidate[a[j]]+h))*(x>(candidate[a[j]]-h))==1)          #neighborhood of candidate
      estimator=rep(0, length(b))
      for(i in 1:length(b)){
        z=(xsub>x[b[i]])+1-1
        model=gam(ysub~z+s(xsub), family=gaussian, method="GCV.Cp")
        estimator[i]=(model$coefficients[2])^2
      }
      candidate[a[j]]=x[b[which.max(estimator)]]
    }
  }
  return(candidate)
}

Pvalue<-function(x,y,candidate,B=500){
  x2=c(0,(candidate[1:(length(candidate)-1)]+candidate[2:length(candidate)])/2,1.01)
  pvalue=1:length(candidate)
  for(k in 1:length(candidate)){
    subx=x[which((x>x2[k])&(x<x2[k+1]))]
    suby=y[which((x>x2[k])&(x<x2[k+1]))]
    z=(subx>candidate[k])+1-1
    model=gam(suby~z+s(subx),family=gaussian, method="GCV.Cp")           
    fitted=model$fitted.values-model$coefficients[2]*z
    residual=model$residuals
    #dependent wild bootstrap
    l=3/length(x)
    estimator=rep(0,B)
    Sigma=diag(length(subx))
    for(i in 1:(length(subx)-1)){ for(j in (i+1):(length(subx))){ Sigma[i,j]=(1-abs((subx[i]-subx[j])/l))*(abs((subx[i]-subx[j])/l)<1)}}
    for(i in 2:(length(subx))){ for(j in 1:(i-1)){ Sigma[i,j]=Sigma[j,i]}}
    E=rmvnorm(n=B, mean=rep(0, length(residual)), sigma=Sigma)
    for(b in 1:B){
      yy=fitted+residual*E[b,]
      estimator[b]=gam(yy~z+s(subx),family=gaussian, method="GCV.Cp")$coefficients[2]
    }
    pvalue[k]=sum(abs(estimator)>abs(model$coefficients[2]))/B
  }
  return(pvalue)
}

##############################sample X1333-4
#######CBS
set.seed(2)
num1=rep(0,5)
y1=data[,1]
CBS=DNAcopy::segment(CNA(y1, rep(1,n), 1:n))                
num1[1]=length(CBS$output[,4])-1
estimate1=CBS$output[1:num1[1],4]
######cumSeg
estimate2=jumpoints(y1,k=60,output="2")$psi
num1[2]=length(estimate2)
######DP
DP=tilingArray::segment(y1, maxseg=100, maxk=n/3)
num1[3]=which.max(logLik(DP, penalty="BIC"))-1
estimate3=DP@breakpoints[[which.max(logLik(DP, penalty="BIC"))]][,"estimate"]
######Proposed I
index=which(is.na(CGHdata[2:2301,19])==1)
x=((1:n)/n)[-index]
y=na.omit(CGHdata[2:2301,19])
h=20/length(x)
D=screening(x,y,h)
lambda=4*mad(D)
initial=x[localMax(D,span=(h*n))[which(D[localMax(D,span=(h*n))]>lambda)]]
candidate=refine1(x,y,h,initial)
pvalue=Pvalue(x,y,candidate)
estimate4=n*candidate[which(p.adjust(pvalue, "BH")<0.05)]
num1[4]=length(estimate4)
######Proposed II
candidate=refine2(x,y,h,initial)
pvalue=Pvalue(x,y,candidate)
estimate5=n*candidate[which(p.adjust(pvalue, "BH")<0.05)]
num1[5]=length(estimate5)
##############show figure
par(mfrow=c(2,3))
estimate11=c(0,estimate1,n)
fitted1=NULL
for(i in 1:(length(estimate1)+1)){
  m=mean(y1[(estimate11[i]+1):estimate11[i+1]])
  fitted1=c(fitted1,rep(m, estimate11[i+1]-estimate11[i]))
}
plot((1:n)/n,y1,xlab="locations",main="CBS",pch=20,col=8)
lines((1:n)/n,fitted1,col=2,lwd=2)
estimate22=c(0,estimate2,n)
fitted2=NULL
for(i in 1:(length(estimate2)+1)){
  m=mean(y1[(estimate22[i]+1):estimate22[i+1]])
  fitted2=c(fitted2,rep(m, estimate22[i+1]-estimate22[i]))
}
plot((1:n)/n,y1,xlab="locations",main="cumSeg",pch=20,col=8)
lines((1:n)/n,fitted2,col=2,lwd=2)
estimate33=c(0,estimate3,n)
fitted3=NULL
for(i in 1:(length(estimate3)+1)){
  m=mean(y1[(estimate33[i]+1):estimate33[i+1]])
  fitted3=c(fitted3,rep(m, estimate33[i+1]-estimate33[i]))
}
plot((1:n)/n,y1,xlab="locations",main="DP",pch=20,col=8)
lines((1:n)/n,fitted3,col=2,lwd=2)
estimate44=c(0,estimate4,n)
fitted4=NULL
for(i in 1:(length(estimate4)+1)){
  yy=y1[(estimate44[i]+1):estimate44[i+1]]
  xx=((estimate44[i]+1):estimate44[i+1])/n
  fitted4=c(fitted4, npreg(yy~xx, bws=0.2*length(xx)/n)$mean)
}
plot((1:n)/n,y1,xlab="locations",main="Proposed I",pch=20,col=8)
lines((1:n)/n,fitted4,col=2,lwd=2)
estimate55=c(0,estimate5,n)
fitted5=NULL
for(i in 1:(length(estimate5)+1)){
  yy=y1[(estimate55[i]+1):estimate55[i+1]]
  xx=((estimate55[i]+1):estimate55[i+1])/n
  fitted5=c(fitted5, npreg(yy~xx, bws=0.2*length(xx)/n)$mean)
}
plot((1:n)/n,y1,xlab="locations",main="Proposed II",pch=20,col=8)
lines((1:n)/n,fitted5,col=2,lwd=2)
residual=y1-fitted5
pacf(residual,lag.max=20,main="PACF")

##############################sample X1533-1
num2=rep(0,5)
y1=data[,2]
CBS=DNAcopy::segment(CNA(y1, rep(1,n), 1:n))                
num2[1]=length(CBS$output[,4])-1
estimate1=CBS$output[1:num1[1],4]
######cumSeg
estimate2=jumpoints(y1,k=60,output="2")$psi
num2[2]=length(estimate2)
######DP
DP=tilingArray::segment(y1, maxseg=100, maxk=n/3)
num2[3]=which.max(logLik(DP, penalty="BIC"))-1
estimate3=DP@breakpoints[[which.max(logLik(DP, penalty="BIC"))]][,"estimate"]
######Proposed I
index=which(is.na(CGHdata[2:2301,58])==1)
x=((1:n)/n)[-index]
y=na.omit(CGHdata[2:2301,58])
h=20/length(x)
D=screening(x,y,h)
lambda=4*mad(D)
initial=x[localMax(D,span=(h*n))[which(D[localMax(D,span=(h*n))]>lambda)]]
candidate=refine1(x,y,h,initial)
pvalue=Pvalue(x,y,candidate)
estimate4=n*candidate[which(p.adjust(pvalue, "BH")<0.05)]
num2[4]=length(estimate4)
######Proposed II
candidate=refine2(x,y,h,initial)
pvalue=Pvalue(x,y,candidate)
estimate5=n*candidate[which(p.adjust(pvalue, "BH")<0.05)]
num2[5]=length(estimate5)

##############################sample X1533-10
num3=rep(0,5)
y1=data[,3]
CBS=DNAcopy::segment(CNA(y1, rep(1,n), 1:n))                
num3[1]=length(CBS$output[,4])-1
estimate1=CBS$output[1:num1[1],4]
######cumSeg
estimate2=jumpoints(y1,k=60,output="2")$psi
num3[2]=length(estimate2)
######DP
DP=tilingArray::segment(y1, maxseg=100, maxk=n/3)
num3[3]=which.max(logLik(DP, penalty="BIC"))-1
estimate3=DP@breakpoints[[which.max(logLik(DP, penalty="BIC"))]][,"estimate"]
######Proposed I
index=which(is.na(CGHdata[2:2301,61])==1)
x=((1:n)/n)[-index]
y=na.omit(CGHdata[2:2301,61])
h=20/length(x)
D=screening(x,y,h)
lambda=4*mad(D)
initial=x[localMax(D,span=(h*n))[which(D[localMax(D,span=(h*n))]>lambda)]]
candidate=refine1(x,y,h,initial)
pvalue=Pvalue(x,y,candidate)
estimate4=n*candidate[which(p.adjust(pvalue, "BH")<0.05)]
num3[4]=length(estimate4)
######Proposed II
candidate=refine2(x,y,h,initial)
pvalue=Pvalue(x,y,candidate)
estimate5=n*candidate[which(p.adjust(pvalue, "BH")<0.05)]
num3[5]=length(estimate5)

##############################sample X1533-13
num4=rep(0,5)
y1=data[,4]
CBS=DNAcopy::segment(CNA(y1, rep(1,n), 1:n))                
num4[1]=length(CBS$output[,4])-1
estimate1=CBS$output[1:num1[1],4]
######cumSeg
estimate2=jumpoints(y1,k=60,output="2")$psi
num4[2]=length(estimate2)
######DP
DP=tilingArray::segment(y1, maxseg=100, maxk=n/3)
num4[3]=which.max(logLik(DP, penalty="BIC"))-1
estimate3=DP@breakpoints[[which.max(logLik(DP, penalty="BIC"))]][,"estimate"]
######Proposed I
index=which(is.na(CGHdata[2:2301,64])==1)
x=((1:n)/n)[-index]
y=na.omit(CGHdata[2:2301,64])
h=20/length(x)
D=screening(x,y,h)
lambda=4*mad(D)
initial=x[localMax(D,span=(h*n))[which(D[localMax(D,span=(h*n))]>lambda)]]
candidate=refine1(x,y,h,initial)
pvalue=Pvalue(x,y,candidate)
estimate4=n*candidate[which(p.adjust(pvalue, "BH")<0.05)]
num4[4]=length(estimate4)
######Proposed II
candidate=refine2(x,y,h,initial)
pvalue=Pvalue(x,y,candidate)
estimate5=n*candidate[which(p.adjust(pvalue, "BH")<0.05)]
num4[5]=length(estimate5)

############show table
rbind(num1,num2,num3,num4)
