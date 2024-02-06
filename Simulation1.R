library(KernSmooth)
library(MASS)
library(Matrix)
library(gam)
library(mvtnorm)

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
        model=gam(ysub~z+s(xsub,5))
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
        model=gam(ysub~z+s(xsub,5))
        estimator[i]=(model$coefficients[2])^2
      }
      candidate[a[j]]=x[b[which.max(estimator)]]
    }
  }
  return(candidate)
}

#######################################Simulation1
###############Scenario I: regular design
set.seed(23)
SD1=NULL                #standard errors
n=400
x=(1:n)/n
signal=exp(x)+(x>0.2)-(x>0.5)+(x>0.8)
h=10/n
lambda=0.7
###############Case 1 iid error
initial1=NULL           #initial estimator of tau1
initial2=NULL
initial3=NULL
estimator1=NULL           #refined (proposed I) estimator of tau1 
estimator2=NULL
estimator3=NULL
Estimator1=NULL           #refined (proposed II) estimator of tau1 
Estimator2=NULL
Estimator3=NULL
for(j in 1:100){
  e=rnorm(n,mean=0,sd=0.2)
  y=signal+e
  D=screening(x,y,h)
  initial=x[localMax(D,span=(h*n))[which(D[localMax(D,span=(h*n))]>lambda)]]
  if(min(abs(initial-0.2))<=h) initial1=c(initial1, initial[which.min(abs(initial-0.2))])
  if(min(abs(initial-0.5))<=h) initial2=c(initial2, initial[which.min(abs(initial-0.5))])
  if(min(abs(initial-0.8))<=h) initial3=c(initial3, initial[which.min(abs(initial-0.8))])
  estimator=refine1(x,y,h,initial)
  if(min(abs(estimator-0.2))<=h) estimator1=c(estimator1, estimator[which.min(abs(estimator-0.2))])
  if(min(abs(estimator-0.5))<=h) estimator2=c(estimator2, estimator[which.min(abs(estimator-0.5))])
  if(min(abs(estimator-0.8))<=h) estimator3=c(estimator3, estimator[which.min(abs(estimator-0.8))])
  Estimator=refine2(x,y,h,initial)
  if(min(abs(Estimator-0.2))<=h) Estimator1=c(Estimator1, Estimator[which.min(abs(Estimator-0.2))])
  if(min(abs(Estimator-0.5))<=h) Estimator2=c(Estimator2, Estimator[which.min(abs(Estimator-0.5))])
  if(min(abs(Estimator-0.8))<=h) Estimator3=c(Estimator3, Estimator[which.min(abs(Estimator-0.8))])
}
SD1=rbind(SD1, c(sd(initial1),sd(initial2),sd(initial3)), c(sd(estimator1),sd(estimator2),sd(estimator3)), c(sd(Estimator1),sd(Estimator2),sd(Estimator3)))
###############Case 2  
initial1=NULL           #initial estimator of tau1
initial2=NULL
initial3=NULL
estimator1=NULL           #refined (proposed I) estimator of tau1 
estimator2=NULL
estimator3=NULL
Estimator1=NULL           #refined (proposed II) estimator of tau1 
Estimator2=NULL
Estimator3=NULL
for(j in 1:100){
  e=rnorm(n,mean=0,sd=0.1*(1+(1:400)/200))
  y=signal+e
  D=screening(x,y,h)
  initial=x[localMax(D,span=(h*n))[which(D[localMax(D,span=(h*n))]>lambda)]]
  if(min(abs(initial-0.2))<=h) initial1=c(initial1, initial[which.min(abs(initial-0.2))])
  if(min(abs(initial-0.5))<=h) initial2=c(initial2, initial[which.min(abs(initial-0.5))])
  if(min(abs(initial-0.8))<=h) initial3=c(initial3, initial[which.min(abs(initial-0.8))])
  estimator=refine1(x,y,h,initial)
  if(min(abs(estimator-0.2))<=h) estimator1=c(estimator1, estimator[which.min(abs(estimator-0.2))])
  if(min(abs(estimator-0.5))<=h) estimator2=c(estimator2, estimator[which.min(abs(estimator-0.5))])
  if(min(abs(estimator-0.8))<=h) estimator3=c(estimator3, estimator[which.min(abs(estimator-0.8))])
  Estimator=refine2(x,y,h,initial)
  if(min(abs(Estimator-0.2))<=h) Estimator1=c(Estimator1, Estimator[which.min(abs(Estimator-0.2))])
  if(min(abs(Estimator-0.5))<=h) Estimator2=c(Estimator2, Estimator[which.min(abs(Estimator-0.5))])
  if(min(abs(Estimator-0.8))<=h) Estimator3=c(Estimator3, Estimator[which.min(abs(Estimator-0.8))])
}
SD1=rbind(SD1, c(sd(initial1),sd(initial2),sd(initial3)), c(sd(estimator1),sd(estimator2),sd(estimator3)), c(sd(Estimator1),sd(Estimator2),sd(Estimator3)))
###############Case 3  
initial1=NULL           #initial estimator of tau1
initial2=NULL
initial3=NULL
estimator1=NULL           #refined (proposed I) estimator of tau1 
estimator2=NULL
estimator3=NULL
Estimator1=NULL           #refined (proposed II) estimator of tau1 
Estimator2=NULL
Estimator3=NULL
for(j in 1:100){
  e=arima.sim(n=500, list(ar = 0.3), sd = 0.2)[101:500]
  y=signal+e
  D=screening(x,y,h)
  initial=x[localMax(D,span=(h*n))[which(D[localMax(D,span=(h*n))]>lambda)]]
  if(min(abs(initial-0.2))<=h) initial1=c(initial1, initial[which.min(abs(initial-0.2))])
  if(min(abs(initial-0.5))<=h) initial2=c(initial2, initial[which.min(abs(initial-0.5))])
  if(min(abs(initial-0.8))<=h) initial3=c(initial3, initial[which.min(abs(initial-0.8))])
  estimator=refine1(x,y,h,initial)
  if(min(abs(estimator-0.2))<=h) estimator1=c(estimator1, estimator[which.min(abs(estimator-0.2))])
  if(min(abs(estimator-0.5))<=h) estimator2=c(estimator2, estimator[which.min(abs(estimator-0.5))])
  if(min(abs(estimator-0.8))<=h) estimator3=c(estimator3, estimator[which.min(abs(estimator-0.8))])
  Estimator=refine2(x,y,h,initial)
  if(min(abs(Estimator-0.2))<=h) Estimator1=c(Estimator1, Estimator[which.min(abs(Estimator-0.2))])
  if(min(abs(Estimator-0.5))<=h) Estimator2=c(Estimator2, Estimator[which.min(abs(Estimator-0.5))])
  if(min(abs(Estimator-0.8))<=h) Estimator3=c(Estimator3, Estimator[which.min(abs(Estimator-0.8))])
}
SD1=rbind(SD1, c(sd(initial1),sd(initial2),sd(initial3)), c(sd(estimator1),sd(estimator2),sd(estimator3)), c(sd(Estimator1),sd(Estimator2),sd(Estimator3)))
###############Case 4  
initial1=NULL           #initial estimator of tau1
initial2=NULL
initial3=NULL
estimator1=NULL           #refined (proposed I) estimator of tau1 
estimator2=NULL
estimator3=NULL
Estimator1=NULL           #refined (proposed II) estimator of tau1 
Estimator2=NULL
Estimator3=NULL
for(j in 1:100){
  e=0.7*(1+(1:400)/400)*arima.sim(n=500, list(ar = 0.3), sd = 0.2)[101:500]
  y=signal+e
  D=screening(x,y,h)
  initial=x[localMax(D,span=(h*n))[which(D[localMax(D,span=(h*n))]>lambda)]]
  if(min(abs(initial-0.2))<=h) initial1=c(initial1, initial[which.min(abs(initial-0.2))])
  if(min(abs(initial-0.5))<=h) initial2=c(initial2, initial[which.min(abs(initial-0.5))])
  if(min(abs(initial-0.8))<=h) initial3=c(initial3, initial[which.min(abs(initial-0.8))])
  estimator=refine1(x,y,h,initial)
  if(min(abs(estimator-0.2))<=h) estimator1=c(estimator1, estimator[which.min(abs(estimator-0.2))])
  if(min(abs(estimator-0.5))<=h) estimator2=c(estimator2, estimator[which.min(abs(estimator-0.5))])
  if(min(abs(estimator-0.8))<=h) estimator3=c(estimator3, estimator[which.min(abs(estimator-0.8))])
  Estimator=refine2(x,y,h,initial)
  if(min(abs(Estimator-0.2))<=h) Estimator1=c(Estimator1, Estimator[which.min(abs(Estimator-0.2))])
  if(min(abs(Estimator-0.5))<=h) Estimator2=c(Estimator2, Estimator[which.min(abs(Estimator-0.5))])
  if(min(abs(Estimator-0.8))<=h) Estimator3=c(Estimator3, Estimator[which.min(abs(Estimator-0.8))])
}
SD1=rbind(SD1, c(sd(initial1),sd(initial2),sd(initial3)), c(sd(estimator1),sd(estimator2),sd(estimator3)), c(sd(Estimator1),sd(Estimator2),sd(Estimator3)))

################Scenario II: irregular design
set.seed(33)
SD2=NULL                #standard errors
index=c(12,52,62,80,120,122,146,152,202,228,232,240,250,268,272,318,332,348,360,382)       #missing index
x=x[-index]
###############Case 1 iid error
initial1=NULL           #initial estimator of tau1
initial2=NULL
initial3=NULL
estimator1=NULL           #refined (proposed I) estimator of tau1 
estimator2=NULL
estimator3=NULL
Estimator1=NULL           #refined (proposed II) estimator of tau1 
Estimator2=NULL
Estimator3=NULL
total1=NULL
total2=NULL
total3=NULL
for(j in 1:100){
  e=rnorm(n,mean=0,sd=0.2)
  y=(signal+e)[-index]
  D=screening(x,y,h)
  initial=x[localMax(D,span=(h*n))[which(D[localMax(D,span=(h*n))]>lambda)]]
  if(min(abs(initial-0.2))<=h) initial1=c(initial1, initial[which.min(abs(initial-0.2))])
  if(min(abs(initial-0.5))<=h) initial2=c(initial2, initial[which.min(abs(initial-0.5))])
  if(min(abs(initial-0.8))<=h) initial3=c(initial3, initial[which.min(abs(initial-0.8))])
  estimator=refine1(x,y,h,initial)
  if(min(abs(estimator-0.2))<=h) estimator1=c(estimator1, estimator[which.min(abs(estimator-0.2))])
  if(min(abs(estimator-0.5))<=h) estimator2=c(estimator2, estimator[which.min(abs(estimator-0.5))])
  if(min(abs(estimator-0.8))<=h) estimator3=c(estimator3, estimator[which.min(abs(estimator-0.8))])
  Estimator=refine2(x,y,h,initial)
  if(min(abs(Estimator-0.2))<=h) Estimator1=c(Estimator1, Estimator[which.min(abs(Estimator-0.2))])
  if(min(abs(Estimator-0.5))<=h) Estimator2=c(Estimator2, Estimator[which.min(abs(Estimator-0.5))])
  if(min(abs(Estimator-0.8))<=h) Estimator3=c(Estimator3, Estimator[which.min(abs(Estimator-0.8))])
  total1=c(total1,initial)
  total2=c(total2,estimator)
  total3=c(total3,Estimator)
}
SD2=rbind(SD2, c(sd(initial1),sd(initial2),sd(initial3)), c(sd(estimator1),sd(estimator2),sd(estimator3)), c(sd(Estimator1),sd(Estimator2),sd(Estimator3)))
###############show plot
par(mfrow=c(4,3))
hist(total1, breaks=60, xlab="x", ylab="frequency", xlim=c(0,1), ylim=c(0,100), main="initial, Case 1")
hist(total2, breaks=60, xlab="x", ylab="frequency", xlim=c(0,1), ylim=c(0,100), main="proposed I, Case 1")
hist(total3, breaks=60, xlab="x", ylab="frequency", xlim=c(0,1), ylim=c(0,100), main="proposed II, Case 1")

###############Case 2  
initial1=NULL           #initial estimator of tau1
initial2=NULL
initial3=NULL
estimator1=NULL           #refined (proposed I) estimator of tau1 
estimator2=NULL
estimator3=NULL
Estimator1=NULL           #refined (proposed II) estimator of tau1 
Estimator2=NULL
Estimator3=NULL
total1=NULL
total2=NULL
total3=NULL
for(j in 1:100){
  e=rnorm(n,mean=0,sd=0.1*(1+(1:400)/200))
  y=(signal+e)[-index]
  D=screening(x,y,h)
  initial=x[localMax(D,span=(h*n))[which(D[localMax(D,span=(h*n))]>lambda)]]
  if(min(abs(initial-0.2))<=h) initial1=c(initial1, initial[which.min(abs(initial-0.2))])
  if(min(abs(initial-0.5))<=h) initial2=c(initial2, initial[which.min(abs(initial-0.5))])
  if(min(abs(initial-0.8))<=h) initial3=c(initial3, initial[which.min(abs(initial-0.8))])
  estimator=refine1(x,y,h,initial)
  if(min(abs(estimator-0.2))<=h) estimator1=c(estimator1, estimator[which.min(abs(estimator-0.2))])
  if(min(abs(estimator-0.5))<=h) estimator2=c(estimator2, estimator[which.min(abs(estimator-0.5))])
  if(min(abs(estimator-0.8))<=h) estimator3=c(estimator3, estimator[which.min(abs(estimator-0.8))])
  Estimator=refine2(x,y,h,initial)
  if(min(abs(Estimator-0.2))<=h) Estimator1=c(Estimator1, Estimator[which.min(abs(Estimator-0.2))])
  if(min(abs(Estimator-0.5))<=h) Estimator2=c(Estimator2, Estimator[which.min(abs(Estimator-0.5))])
  if(min(abs(Estimator-0.8))<=h) Estimator3=c(Estimator3, Estimator[which.min(abs(Estimator-0.8))])
  total1=c(total1,initial)
  total2=c(total2,estimator)
  total3=c(total3,Estimator)
}
SD2=rbind(SD2, c(sd(initial1),sd(initial2),sd(initial3)), c(sd(estimator1),sd(estimator2),sd(estimator3)), c(sd(Estimator1),sd(Estimator2),sd(Estimator3)))
###############show plot
hist(total1, breaks=60, xlab="x", ylab="frequency", xlim=c(0,1), ylim=c(0,100), main="initial, Case 2")
hist(total2, breaks=60, xlab="x", ylab="frequency", xlim=c(0,1), ylim=c(0,100), main="proposed I, Case 2")
hist(total3, breaks=60, xlab="x", ylab="frequency", xlim=c(0,1), ylim=c(0,100), main="proposed II, Case 2")

###############Case 3  
initial1=NULL           #initial estimator of tau1
initial2=NULL
initial3=NULL
estimator1=NULL           #refined (proposed I) estimator of tau1 
estimator2=NULL
estimator3=NULL
Estimator1=NULL           #refined (proposed II) estimator of tau1 
Estimator2=NULL
Estimator3=NULL
total1=NULL
total2=NULL
total3=NULL
for(j in 1:100){
  e=arima.sim(n=500, list(ar = 0.3), sd = 0.2)[101:500]
  y=(signal+e)[-index]
  D=screening(x,y,h)
  initial=x[localMax(D,span=(h*n))[which(D[localMax(D,span=(h*n))]>lambda)]]
  if(min(abs(initial-0.2))<=h) initial1=c(initial1, initial[which.min(abs(initial-0.2))])
  if(min(abs(initial-0.5))<=h) initial2=c(initial2, initial[which.min(abs(initial-0.5))])
  if(min(abs(initial-0.8))<=h) initial3=c(initial3, initial[which.min(abs(initial-0.8))])
  estimator=refine1(x,y,h,initial)
  if(min(abs(estimator-0.2))<=h) estimator1=c(estimator1, estimator[which.min(abs(estimator-0.2))])
  if(min(abs(estimator-0.5))<=h) estimator2=c(estimator2, estimator[which.min(abs(estimator-0.5))])
  if(min(abs(estimator-0.8))<=h) estimator3=c(estimator3, estimator[which.min(abs(estimator-0.8))])
  Estimator=refine2(x,y,h,initial)
  if(min(abs(Estimator-0.2))<=h) Estimator1=c(Estimator1, Estimator[which.min(abs(Estimator-0.2))])
  if(min(abs(Estimator-0.5))<=h) Estimator2=c(Estimator2, Estimator[which.min(abs(Estimator-0.5))])
  if(min(abs(Estimator-0.8))<=h) Estimator3=c(Estimator3, Estimator[which.min(abs(Estimator-0.8))])
  total1=c(total1,initial)
  total2=c(total2,estimator)
  total3=c(total3,Estimator)
}
SD2=rbind(SD2, c(sd(initial1),sd(initial2),sd(initial3)), c(sd(estimator1),sd(estimator2),sd(estimator3)), c(sd(Estimator1),sd(Estimator2),sd(Estimator3)))
###############show plot
hist(total1, breaks=60, xlab="x", ylab="frequency", xlim=c(0,1), ylim=c(0,100), main="initial, Case 3")
hist(total2, breaks=60, xlab="x", ylab="frequency", xlim=c(0,1), ylim=c(0,100), main="proposed I, Case 3")
hist(total3, breaks=60, xlab="x", ylab="frequency", xlim=c(0,1), ylim=c(0,100), main="proposed II, Case 3")

###############Case 4  
initial1=NULL           #initial estimator of tau1
initial2=NULL
initial3=NULL
estimator1=NULL           #refined (proposed I) estimator of tau1 
estimator2=NULL
estimator3=NULL
Estimator1=NULL           #refined (proposed II) estimator of tau1 
Estimator2=NULL
Estimator3=NULL
total1=NULL
total2=NULL
total3=NULL
for(j in 1:100){
  e=0.7*(1+(1:400)/400)*arima.sim(n=500, list(ar = 0.3), sd = 0.2)[101:500]
  y=(signal+e)[-index]
  D=screening(x,y,h)
  initial=x[localMax(D,span=(h*n))[which(D[localMax(D,span=(h*n))]>lambda)]]
  if(min(abs(initial-0.2))<=h) initial1=c(initial1, initial[which.min(abs(initial-0.2))])
  if(min(abs(initial-0.5))<=h) initial2=c(initial2, initial[which.min(abs(initial-0.5))])
  if(min(abs(initial-0.8))<=h) initial3=c(initial3, initial[which.min(abs(initial-0.8))])
  estimator=refine1(x,y,h,initial)
  if(min(abs(estimator-0.2))<=h) estimator1=c(estimator1, estimator[which.min(abs(estimator-0.2))])
  if(min(abs(estimator-0.5))<=h) estimator2=c(estimator2, estimator[which.min(abs(estimator-0.5))])
  if(min(abs(estimator-0.8))<=h) estimator3=c(estimator3, estimator[which.min(abs(estimator-0.8))])
  Estimator=refine2(x,y,h,initial)
  if(min(abs(Estimator-0.2))<=h) Estimator1=c(Estimator1, Estimator[which.min(abs(Estimator-0.2))])
  if(min(abs(Estimator-0.5))<=h) Estimator2=c(Estimator2, Estimator[which.min(abs(Estimator-0.5))])
  if(min(abs(Estimator-0.8))<=h) Estimator3=c(Estimator3, Estimator[which.min(abs(Estimator-0.8))])
  total1=c(total1,initial)
  total2=c(total2,estimator)
  total3=c(total3,Estimator)
}
SD2=rbind(SD2, c(sd(initial1),sd(initial2),sd(initial3)), c(sd(estimator1),sd(estimator2),sd(estimator3)), c(sd(Estimator1),sd(Estimator2),sd(Estimator3)))
###############show plot
hist(total1, breaks=60, xlab="x", ylab="frequency", xlim=c(0,1), ylim=c(0,100), main="initial, Case 4")
hist(total2, breaks=60, xlab="x", ylab="frequency", xlim=c(0,1), ylim=c(0,100), main="proposed I, Case 4")
hist(total3, breaks=60, xlab="x", ylab="frequency", xlim=c(0,1), ylim=c(0,100), main="proposed II, Case 4")

##############show table
round(cbind(SD1,SD2)*1000,2)