library(KernSmooth)
library(MASS)
library(gam)
library(mvtnorm)


#######################################Simulation2
set.seed(23)
index=c(3,26,31,40,60,61,73,76,101,114,116,118,125,134,136,159,166,174,180,191)       #missing index
B=500
h=3/200
alpha=0.05
beta=0.2*(1:5)-0.2

############Case 1, regular design
x=(1:200)/200
n=length(x)
Sigma=diag(n)
for(i in 1:(n-1)){ 
  for(j in (i+1):n){ 
    Sigma[i,j]=(1-abs((x[i]-x[j])/h))*(abs((x[i]-x[j])/h)<1)
  }
}
for(i in 2:n){ 
  for(j in 1:(i-1)){ 
    Sigma[i,j]=Sigma[j,i]
  }
}
Power1=matrix(0,nrow=5,ncol=100)                 #power function
for(k in 1:5){
  for(l in 1:100){
    signal=exp(x)+beta[k]*(x>0.5)
    e=rnorm(200,mean=0,sd=0.2)
    y=signal+e
    z=(x>0.5)+1-1
    model=gam(y~z+s(x,5))           
    fitted=model$fitted.values-model$coefficients[2]*z
    residual=model$residuals
    estimator=rep(0,B)           #dependent wild bootstrap
    E=rmvnorm(n=B, mean=rep(0, length(residual)), sigma=Sigma)
    for(b in 1:B){
      y1=fitted+residual*E[b,]
      estimator[b]=gam(y1~z+s(x,5))$coefficients[2]
    }
    Power1[k,l]=((sum(abs(estimator)>abs(model$coefficients[2]))/B)<alpha)
  }
}
############Case 1, irregular design
x=x[-index]
n=length(x)
sigma=diag(n)
for(i in 1:(n-1)){ 
  for(j in (i+1):n){ 
    sigma[i,j]=(1-abs((x[i]-x[j])/h))*(abs((x[i]-x[j])/h)<1)
  }
}
for(i in 2:n){ 
  for(j in 1:(i-1)){ 
    sigma[i,j]=sigma[j,i]
  }
}
power1=matrix(0,nrow=5,ncol=100)                        #dependent wild bootstrap
for(k in 1:5){
  for(l in 1:100){
    signal=exp(x)+beta[k]*(x>0.5)
    e=rnorm(200,mean=0,sd=0.2)[-index]
    y=signal+e
    z=(x>0.5)+1-1
    model=gam(y~z+s(x,5))           
    fitted=model$fitted.values-model$coefficients[2]*z
    residual=model$residuals
    estimator=rep(0,B)
    E=rmvnorm(n=B, mean=rep(0, length(residual)), sigma=sigma)
    for(b in 1:B){
      y1=fitted+residual*E[b,]
      estimator[b]=gam(y1~z+s(x,5))$coefficients[2]
    }
    power1[k,l]=((sum(abs(estimator)>abs(model$coefficients[2]))/B)<alpha)
  }
}

############Case 2, regular design
x=(1:200)/200
Power2=matrix(0,nrow=5,ncol=100)                 #power function
for(k in 1:5){
  for(l in 1:100){
    signal=exp(x)+beta[k]*(x>0.5)
    e=rnorm(200,mean=0,sd=0.1*(1+(1:200)/100))
    y=signal+e
    z=(x>0.5)+1-1
    model=gam(y~z+s(x,5))           
    fitted=model$fitted.values-model$coefficients[2]*z
    residual=model$residuals
    estimator=rep(0,B)           #dependent wild bootstrap
    E=rmvnorm(n=B, mean=rep(0, length(residual)), sigma=Sigma)
    for(b in 1:B){
      y1=fitted+residual*E[b,]
      estimator[b]=gam(y1~z+s(x,5))$coefficients[2]
    }
    Power2[k,l]=((sum(abs(estimator)>abs(model$coefficients[2]))/B)<alpha)
  }
}
############Case 2, irregular design
x=x[-index]
power2=matrix(0,nrow=5,ncol=100)                        #dependent wild bootstrap
for(k in 1:5){
  for(l in 1:100){
    signal=exp(x)+beta[k]*(x>0.5)
    e=rnorm(200,mean=0,sd=0.1*(1+(1:200)/100))[-index]
    y=signal+e
    z=(x>0.5)+1-1
    model=gam(y~z+s(x,5))           
    fitted=model$fitted.values-model$coefficients[2]*z
    residual=model$residuals
    estimator=rep(0,B)
    E=rmvnorm(n=B, mean=rep(0, length(residual)), sigma=sigma)
    for(b in 1:B){
      y1=fitted+residual*E[b,]
      estimator[b]=gam(y1~z+s(x,5))$coefficients[2]
    }
    power2[k,l]=((sum(abs(estimator)>abs(model$coefficients[2]))/B)<alpha)
  }
}

############Case 3, regular design
x=(1:200)/200
Power3=matrix(0,nrow=5,ncol=100)                 #power function
for(k in 1:5){
  for(l in 1:100){
    signal=exp(x)+beta[k]*(x>0.5)
    e=arima.sim(n=300, list(ar = 0.3), sd = 0.2)[101:300]
    y=signal+e
    z=(x>0.5)+1-1
    model=gam(y~z+s(x,5))           
    fitted=model$fitted.values-model$coefficients[2]*z
    residual=model$residuals
    estimator=rep(0,B)           #dependent wild bootstrap
    E=rmvnorm(n=B, mean=rep(0, length(residual)), sigma=Sigma)
    for(b in 1:B){
      y1=fitted+residual*E[b,]
      estimator[b]=gam(y1~z+s(x,5))$coefficients[2]
    }
    Power3[k,l]=((sum(abs(estimator)>abs(model$coefficients[2]))/B)<alpha)
  }
}
############Case 3, irregular design
x=x[-index]
power3=matrix(0,nrow=5,ncol=100)                        #dependent wild bootstrap
for(k in 1:5){
  for(l in 1:100){
    signal=exp(x)+beta[k]*(x>0.5)
    e=arima.sim(n=300, list(ar = 0.3), sd = 0.2)[101:300][-index]
    y=signal+e
    z=(x>0.5)+1-1
    model=gam(y~z+s(x,5))           
    fitted=model$fitted.values-model$coefficients[2]*z
    residual=model$residuals
    estimator=rep(0,B)
    E=rmvnorm(n=B, mean=rep(0, length(residual)), sigma=sigma)
    for(b in 1:B){
      y1=fitted+residual*E[b,]
      estimator[b]=gam(y1~z+s(x,5))$coefficients[2]
    }
    power3[k,l]=((sum(abs(estimator)>abs(model$coefficients[2]))/B)<alpha)
  }
}

############Case 4, regular design
x=(1:200)/200
Power4=matrix(0,nrow=5,ncol=100)                 #power function
for(k in 1:5){
  for(l in 1:100){
    signal=exp(x)+beta[k]*(x>0.5)
    e=0.7*(1+(1:200)/200)*arima.sim(n=300, list(ar = 0.3), sd = 0.2)[101:300]
    y=signal+e
    z=(x>0.5)+1-1
    model=gam(y~z+s(x,5))           
    fitted=model$fitted.values-model$coefficients[2]*z
    residual=model$residuals
    estimator=rep(0,B)           #dependent wild bootstrap
    E=rmvnorm(n=B, mean=rep(0, length(residual)), sigma=Sigma)
    for(b in 1:B){
      y1=fitted+residual*E[b,]
      estimator[b]=gam(y1~z+s(x,5))$coefficients[2]
    }
    Power4[k,l]=((sum(abs(estimator)>abs(model$coefficients[2]))/B)<alpha)
  }
}
############Case 4, irregular design
x=x[-index]
power4=matrix(0,nrow=5,ncol=100)                        #dependent wild bootstrap
for(k in 1:5){
  for(l in 1:100){
    signal=exp(x)+beta[k]*(x>0.5)
    e=(0.7*(1+(1:200)/200)*arima.sim(n=300, list(ar = 0.3), sd = 0.2)[101:300])[-index]
    y=signal+e
    z=(x>0.5)+1-1
    model=gam(y~z+s(x,5))           
    fitted=model$fitted.values-model$coefficients[2]*z
    residual=model$residuals
    estimator=rep(0,B)
    E=rmvnorm(n=B, mean=rep(0, length(residual)), sigma=sigma)
    for(b in 1:B){
      y1=fitted+residual*E[b,]
      estimator[b]=gam(y1~z+s(x,5))$coefficients[2]
    }
    power4[k,l]=((sum(abs(estimator)>abs(model$coefficients[2]))/B)<alpha)
  }
}

##############show plots
par(mfrow=c(2,2))
plot(x=beta, y=apply(Power1,1,mean), type="b", pch=1, ylim=c(0,1), xlab="beta", ylab ="power", main="Case 1")
lines(x=beta, y=apply(power1,1,mean), type="b", pch=4)
legend("bottomright", legend=c("regular","irregular"),pch=c(1,4),bty="n")
plot(x=beta, y=apply(Power2,1,mean), type="b", pch=1, ylim=c(0,1), xlab="beta", ylab ="power", main="Case 2")
lines(x=beta, y=apply(power2,1,mean), type="b", pch=4)
legend("bottomright", legend=c("regular","irregular"),pch=c(1,4),bty="n")
plot(x=beta, y=apply(Power3,1,mean), type="b", pch=1, ylim=c(0,1), xlab="beta", ylab ="power", main="Case 3")
lines(x=beta, y=apply(power3,1,mean), type="b", pch=4)
legend("bottomright", legend=c("regular","irregular"),pch=c(1,4),bty="n")
plot(x=beta, y=apply(Power4,1,mean), type="b", pch=1, ylim=c(0,1), xlab="beta", ylab ="power", main="Case 4")
lines(x=beta, y=apply(power4,1,mean), type="b", pch=4)
legend("bottomright", legend=c("regular","irregular"),pch=c(1,4),bty="n")
