library(KernSmooth)
library(MASS)
library(Matrix)
library(mgcv)
library(mvtnorm)
library(np)

##############################Get real data
load("C:/Users/Acer/Desktop/My Document/Documents/Research/Projects/change points(bootstrap)/data4.RData")
y=as.vector(t(data4[4,1:566,1+(1:4-1)*1200]))      #data20[4st stock, days 1-566, open prices]
n=length(y)
x=(1:n)/n

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

############################Real Data 2
set.seed(2)
h=20/n
D=screening(x,y,h)
lambda=4*mad(D)
initial=x[localMax(D,span=(h*n))[which(D[localMax(D,span=(h*n))]>lambda)]]
candidate=refine2(x,y,h,initial)
pvalue=Pvalue(x,y,candidate)
estimate=candidate[which(p.adjust(pvalue, "BH")<0.05)]

##############show figure
plot(y~x, type="l", main="hourly prices")
abline(v=estimate, lty=2)
