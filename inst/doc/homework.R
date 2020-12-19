## -----------------------------------------------------------------------------
library(StatComp20009)

## -----------------------------------------------------------------------------
x <- 0:64/64
y1 <- sin(3*pi*x)
y2 <- cos(3*pi*x)

plot(x, y1, type = "l", col = "red", main = "sin(ax)")
legend("topright", "sin(ax)", lty = 1, col = "red")
abline(h=0 ,lty=3)
plot(x, y2, type = "l", col = "blue", main = "cos(ax)")
legend("topright", "cos(ax)", lty = 1, col = "blue")
abline(h=0 ,lty=3)

## ----results='asis'-----------------------------------------------------------
library(xtable)
data(tli)
tli.table <- xtable(tli[1:15, ])
print(tli.table, type = "html")

## ----results='asis'-----------------------------------------------------------
## Demonstrate aov
fm1 <- aov(tlimth ~ sex + ethnicty + grade + disadvg, data = tli)
fm1t <- xtable(fm1)
print(fm1t, type = "html")

## -----------------------------------------------------------------------------
set.seed(123)
n <- 1000
u <- runif(n) 
x <- 2/(u^(1/2)) 
hist(x,breaks = 40, prob = TRUE, main = expression(f(x)==8/x^3),col = "lightblue") #density histogram of sample 
y <- seq(min(x), max(x), .01) 
lines(y, 8/y^3, col="red") #density curve f(x)


## -----------------------------------------------------------------------------
set.seed(123)
rEpanechnikov<-function(n){
  x0<-c()
  for (i in 1:n) {
    U<-runif(3,min = -1,max = 1)
    if((abs(U[1])<=abs(U[3]))&(abs(U[2])<=abs(U[3])))
      x0<-c(x0,U[2])
    else
      x0<-c(x0,U[3])
  }
  x0
}
x<-rEpanechnikov(1000)
hist(x, prob = TRUE, main = expression(f(x)==3/4*(1-x^2)),col = "lightblue") #density histogram of x
y <- seq(min(x), max(x), .01) 
lines(y, 3/4*(1-y^2), col="red") #density curve f(x)


## -----------------------------------------------------------------------------
set.seed(123)
n <- 1000
rPareto<-function(n,beta,r){
  u <- runif(n) 
  x <- beta/(u^(1/r))-beta
  return(x)
}
x<-rPareto(n,2,4)
hist(x, breaks=40,prob = TRUE, main = expression(f(x)==64/(2+x)^5),col = "lightblue") #density histogram of sample 
y <- seq(min(x), max(x), .01) 
lines(y, 64*(2+y)^(-5), col="red") #density curve f(x)


## -----------------------------------------------------------------------------
n<-10000
set.seed(123)
t<-runif(n,min=0,max=pi/3)
a_hat<-mean(sin(t))
a<-cos(0)-cos(pi/3)
cbind(a_hat,a)

## -----------------------------------------------------------------------------
set.seed(123)
MC.theta<-function(n,anti){
  x1<-runif(n/2,min=0,max=1)
  if(!anti) {
    x2<-runif(n/2,min=0,max=1)#simple Monte Carlo method
    x<-exp(c(x1,x2))
  }
  else {
    x2<-1-x1#antithetic variate approach
    x<-(exp(x1)+exp(x2))/2
  }
  x
}
m<-5000
mc1<-MC.theta(m,anti = FALSE)
theta_hat<-mean(mc1)#Simple MC
mc2<-MC.theta(m,anti = TRUE)
theta1_hat<-mean(mc2)#Antithetic approach


## -----------------------------------------------------------------------------
rbind("Simple MC"=theta_hat,"Antithetic approach"=theta1_hat,"True value"=exp(1)-exp(0))

## -----------------------------------------------------------------------------
rbind("var of simple MC"=var(mc1)/m,"var of antithetic approach"=var(mc2)/(m/2),"the percent reduction in var"=(var(mc1)/m-var(mc2)/(m/2))/(var(mc1)/m))

## -----------------------------------------------------------------------------
rbind("theoretical variance of simple MC"=(2*exp(1)-exp(2)/2-3/2)/m,"theoretical variance of antithetic approach"=(5*exp(1)-3*exp(2)/2-5/2)/m)

## -----------------------------------------------------------------------------
set.seed(123)
m<-5000
g<-function(x) x^2/(sqrt(2*pi)*exp(x^2/2))
#inverse transform method#inverse transform method
y<-rexp(m,rate = 1)
x1<-sqrt(2*y+1)
f1g<-g(x1)/(x1*exp((1-x1^2)/2))
theta_hat1<-mean(f1g)

## -----------------------------------------------------------------------------
#inverse transform method
set.seed(123)
x2<-rnorm(10*m)
x2<-x2[x2>1]
f2g<-g(x2)/(exp(-x2^2/2)/sqrt(2*pi)/(1-pnorm(1)))
theta_hat2<-mean(f2g)
true.value=integrate(g,1,Inf)
cat("The estimated value is:\n")
cbind(theta_hat1,theta_hat2,TRUE.VALUE=true.value$value)
cat("The variance using f1 and f2: \n",c(var(f1g),var(f2g)))

## -----------------------------------------------------------------------------
x <- seq(1, 10, .1)
w <- 2
f1 <- x*exp((1-x^2)/2)
f2 <- exp(-x^2/2)/sqrt(2*pi)/(1-pnorm(1))
g0 <- x^2/(sqrt(2*pi)*exp(x^2/2))
#figure (b)
plot(x, g0, type = "n", main = "",ylim = c(0,10), ylab = "",lwd = w, lty = 2)
lines(x, g0/f1, col='red', lwd = w)
lines(x, g0/f2, col='blue', lwd = w)
legend("topright", legend = c('g/f1','g/f2'),col=c('red','blue'), lwd = w, inset = 0.02)

## -----------------------------------------------------------------------------
set.seed(123)
M <- 5000 #number of replicates
k <- 5 #number of strata
r <- M / k #replicates per stratum
N <- 50 #number of times to repeat the estimation
T2 <- numeric(k)
estimates <- matrix(0,N,2)
g <- function(x) {
  exp(-x)/(1+x^2)
}
for (i in 1:N) {
for (j in 1:k){
  u<-runif(M/k)
  x<--log(exp(-(j-1)/k)-u*(exp(-(j-1)/k)-exp(-j/k)))
  fg<-g(x)/(exp(-x)/(exp(-(j-1)/k)-exp(-j/k)))
  T2[j] <- mean(fg)
}
estimates[i,1] <- sum(T2)
}

## -----------------------------------------------------------------------------
g <- function(x) {
  exp(-x)/(1+x^2)*(x>0)*(x<1)
}
for (i in 1:N) {
for (j in 1:k){
  u<-runif(M/k,(j-1)/k,j/k)
  x<- -log(1-u*(1-exp(-1)))
  fg<-g(x)/(exp(-x)/(1-exp(-1))*k)
  T2[j] <- mean(fg)
}
estimates[i,2] <- sum(T2)
}
cat("The estimated value are:\n")
apply(estimates,2,mean)
cat("The ses are:\n")
apply(estimates,2,sd)

## -----------------------------------------------------------------------------
set.seed(123)
n<-2000
alpha <- .05
Up <- replicate(1000, expr = {
x <- rlnorm(n,1,2)
w<-sd(log(x))
mean(log(x))+w*qt(1-alpha/2, df = n-1)/sqrt(n)
} )
Low <- replicate(1000, expr = {
x <- rlnorm(n,1,2)
w<-sd(log(x))
mean(log(x))-w*qt(1-alpha/2, df = n-1)/sqrt(n)
} )
cat("one emprical estimate of 95% CI for mean is :","\n",'[',mean(Low),mean(Up),']')

## -----------------------------------------------------------------------------
cat("an empirical estimate of the confidence level:\n",mean((Up>1)&(Low<1)))

## -----------------------------------------------------------------------------
set.seed(123)
n<-20
alpha <- .05
Up <- replicate(1000, expr = {
x <- rchisq(n,df=2)
mean(x)+sd(x)*qt(1-alpha/2, df = n-1)/sqrt(n)
} )
Low <- replicate(1000, expr = {
x <- rchisq(n,df=2)
mean(x)-sd(x)*qt(1-alpha/2, df = n-1)/sqrt(n)
} )

## -----------------------------------------------------------------------------
cat('The probability that the confidence interval covers the mean is:\n',mean((Low<2)&(Up>2)))

## -----------------------------------------------------------------------------
Up <- replicate(1000, expr = {
x <- rchisq(n,df=2)
var(x)*(n-1)/qchisq(alpha, df = n-1)
} )
cat("The simulation results in Example 6.4\n",mean(Up>4))

## -----------------------------------------------------------------------------
sks <- function(x) {
  m <- mean(x)
  a <- mean((x - m)^3)
  b <- mean((x - m)^2)
  return(a/b^1.5)
}
#For this experiment, the significance level is 0.1 and the sample size is n = 30.
sl <- .1
n <- 30
m <- 2000
alpha <- c(seq(0.5, 50, 0.5))
N <- length(alpha)#不同的参数值个数
pw1 <- pw2 <- numeric(N)
#critical value for the skewness test
cv <- qnorm(1-sl/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) { #for each epsilon
e <- alpha[j]
sksts1 <- sksts2 <- numeric(m)
for (i in 1:m) { #for each replicate
x <- rbeta(n, e, e)#Beta distribution
sksts1[i] <- as.integer(abs(sks(x)) >= cv)
y <- rt(n, e)#t distribution
sksts2[i] <- as.integer(abs(sks(y)) >= cv)
}
pw1[j] <- mean(sksts1)
pw2[j] <- mean(sksts2)
}

## -----------------------------------------------------------------------------
plot(alpha, pw1, type = "l", xlab = bquote(alpha), ylim = c(0,1))
abline(h = .1, lty = 3)
se1 <- sqrt(pw1 * (1-pw1) / m) #add standard errors
lines(alpha, pw1+se1, lty = 3,col="red")
lines(alpha, pw1-se1, lty = 3,col="blue")

## -----------------------------------------------------------------------------
plot(alpha, pw2, type = "l", xlab = bquote(alpha), ylim = c(0,1))
abline(h = .1, lty = 3)
se2 <- sqrt(pw2 * (1-pw2) / m) #add standard errors
lines(alpha, pw2+se2, lty = 3,col="red")
lines(alpha, pw2-se2, lty = 3,col="blue")

## -----------------------------------------------------------------------------
c5t <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  return(as.integer(max(c(outx, outy)) > 5))
}
Ft <- function(x, y, p) {
  a<-var.test(x, y, conf.level = 1-p)
  return(as.integer(a$p.value <= p))
}

## -----------------------------------------------------------------------------
# generate samples under H1 to estimate power
s1 <- 1
s2 <- 1.5
N<-c(15,100,1000)
m<-1000
p5<-pf<-numeric(3)
for (i in 1:3) {
  n=N[i]
  #Count Five test
  power5 <- mean(replicate(m, expr={
    x <- rnorm(n, 0, s1)
    y <- rnorm(n, 0, s2)
    c5t(x, y)
    }))
  p5[i]<-power5
  #F test
  powerF <- mean(replicate(m, expr={
    x <- rnorm(n, 0, s1)
    y <- rnorm(n, 0, s2)
    Ft(x, y,0.055)
    }))
  pf[i]<-powerF
}
P<-rbind(p5,pf)

## -----------------------------------------------------------------------------
colnames(P)<-c('small size','medium size','large size')
rownames(P)<-c('Count Five test','F test')
knitr::kable(P,align = "c",caption = "the power of the Count Five test vs F test")

## -----------------------------------------------------------------------------
library(MASS)
Mardia<-function(mydata){
  n=nrow(mydata)
  c=ncol(mydata)
  central<-mydata
  for(i in 1:c){
    central[,i]<-mydata[,i]-mean(mydata[,i])
  }
  sigmah<-t(central)%*%central/n
  a<-central%*%solve(sigmah)%*%t(central)
  b<-sum(colSums(a^{3}))/(n*n)
  test<-n*b/6
  chi<-qchisq(0.95,c*(c+1)*(c+2)/6)
  as.integer(test>chi)
}

set.seed(1234)
mu <- c(0,0,0)
sigma <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
m=1000
n<-c(10, 20, 30, 50, 100, 500)
#m: number of replicates; n: sample size
a=numeric(length(n))
for(i in 1:length(n)){
  a[i]=mean(replicate(m, expr={
    mydata <- mvrnorm(n[i],mu,sigma) 
    Mardia(mydata)
  }))
}

## -----------------------------------------------------------------------------
print(a)

## -----------------------------------------------------------------------------
library(MASS)
set.seed(1234)
set.seed(1234)
mu1 <- mu2 <- c(0,0,0)
sigma1 <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
sigma2 <- matrix(c(100,0,0,0,100,0,0,0,100),nrow=3,ncol=3)
sigma=list(sigma1,sigma2)
m=1000
n=50
#m: number of replicates; n: sample size
epsilon <- c(seq(0, .06, .01), seq(.1, 1, .05))
N <- length(epsilon)
pwr <- numeric(N)
for (j in 1:N) { #for each epsilon
  e <- epsilon[j]
  sktests <- numeric(m)
  for (i in 1:m) { #for each replicate
    index=sample(c(1, 2), replace = TRUE, size = n, prob = c(1-e, e))
    mydata<-matrix(0,nrow=n,ncol=3)
    for(t in 1:n){
      if(index[t]==1) mydata[t,]=mvrnorm(1,mu1,sigma1) 
      else mydata[t,]=mvrnorm(1,mu2,sigma2)
    }
    sktests[i] <- Mardia(mydata)
  }
  pwr[j] <- mean(sktests)
}
plot(epsilon, pwr, type = "b",
     xlab = bquote(epsilon), ylim = c(0,1))
abline(h = .05, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(epsilon, pwr+se, lty = 3)
lines(epsilon, pwr-se, lty = 3)

## -----------------------------------------------------------------------------
library(bootstrap) #for the law data
attach(law)
set.seed(123)
R.hat<-cor(LSAT, GPA)
n<-dim(law)[1]
R.jack <- numeric(n)
for (i in 1:n){
  R.jack[i]<-cor(LSAT[-i],GPA[-i])
}
R.bias <- (n - 1) * (mean(R.jack) - R.hat)
#jackknife estimate of bias
R.se <- sqrt((n-1) * mean((R.jack - mean(R.jack))^2))
cbind('jackknife estimate of bias'=R.bias,'jackknife estimate of se'=R.se)
detach(law)

## -----------------------------------------------------------------------------
library(boot)
attach(aircondit)
set.seed(123)
lambda.boot <- function(x, i) {
  #function to compute the statistic
  mean(x[i])
}
boot.result <- boot(hours, statistic = lambda.boot, R = 2000)
print(boot.ci(boot.result,
              type = c("norm","basic", "perc", "bca")))

detach(aircondit)

## -----------------------------------------------------------------------------
set.seed(123)
#compute theta_hat
lambda.hat<-eigen(cov(scor))$values
theta.hat<-lambda.hat[1]/sum(lambda.hat)
#number of rows
n<-nrow(scor)
#Jackknife method
theta.j<-numeric(n)
for (i in 1:n) {
  x<-scor[-i,]
  lambda<-eigen(cov(x))$values
  theta.j[i]<-lambda[1]/sum(lambda)
}
bias.j<-(n-1)*(mean(theta.j)-theta.hat)
se.j<-sqrt((n-1)*mean((theta.j-mean(theta.j))^2))
cbind(bias.j,se.j,theta.hat)

## -----------------------------------------------------------------------------
library(DAAG)
attach(ironslag)
set.seed(1234)
n <- length(magnetic) #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- numeric(n*(n-2)/2)
i<-0
# for n-fold cross validation
# fit models on leave-one-out samples
for (k in 1:(n-1)) {
  for (j in (k+1):n) {
    i<-i+1
    y<-magnetic[c(-k,-j)]
    x <- chemical[c(-k,-j)]
    J1 <- lm(y ~ x)
    yhat1 <- J1$coef[1] + J1$coef[2] * chemical[c(k,j)]
    e1[i] <- sum((magnetic[c(k,j)] - yhat1)^2)
    J2 <- lm(y ~ x + I(x^2))
    yhat2 <- J2$coef[1] + J2$coef[2] * chemical[c(k,j)] +
      J2$coef[3] * chemical[c(k,j)]^2
    e2[i] <- sum((magnetic[c(k,j)] - yhat2)^2)
    J3 <- lm(log(y) ~ x)
    logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[c(k,j)]
    yhat3 <- exp(logyhat3)
    e3[i] <- sum((magnetic[c(k,j)] - yhat3)^2)
    J4 <- lm(log(y) ~ log(x))
    logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[c(k,j)])
    yhat4 <- exp(logyhat4)
    e4[i] <- sum((magnetic[c(k,j)] - yhat4)^2)
  }
}
c(mean(e1),mean(e2),mean(e3),mean(e4))
detach(ironslag)

## ----echo=FALSE---------------------------------------------------------------
library(RANN)
library(boot)
library(energy)
library(Ball)
library(ggplot2)

## -----------------------------------------------------------------------------
#the function computes the maximum number of extreme points
maxout <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  return(max(c(outx, outy)))
}

#the function computes the statistics 
st<-function(z,i,n){
  x<-z[i][1:n]
  y<-z[i][-(1:n)]
  maxout(x,y)
}
# the function calculates p value
c5.p<-function(n,mu,sd){
  x<-rnorm(n[1],mu[1],sd[1])
  y<-rnorm(n[2],mu[2],sd[2])
  z<-c(x,y)
  R=999
  boot.obj<-boot(z,statistic = st,R=R,sim='permutation',n=n[1])
  count<-c(boot.obj$t0, boot.obj$t)
  p.value<-mean(count>=count[1])
  return(p.value)
}
set.seed(123)
n<-c(20,30)#different sample sizes
mu<-c(0,0)
N<-1000
p<-numeric(N)

# calculate the empirical type I error rate
for(i in 1:N) p[i]<-c5.p(n,mu,sd=c(1,1))
type.1.error<-mean(p<0.05)

# calculate the power
for(i in 1:N) p[i]<-c5.p(n,mu,sd=c(1,2))
e.power<-mean(p<0.05)

cbind(type.1.error,e.power)

## -----------------------------------------------------------------------------
# NN method
Tn <- function(z, ix, sizes,k) {
  n1 <- sizes[1]
  n2 <- sizes[2]
  n <- n1 + n2
  if(is.vector(z)) z <- data.frame(z,0);
  z <- z[ix, ];
  NN <- nn2(data=z, k=k+1) # what's the first column?
  block1 <- NN$nn.idx[1:n1,-1]
  block2 <- NN$nn.idx[(n1+1):n,-1]
  i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
  (i1 + i2) / (k * n)
}

eqdist.nn <- function(z,sizes,k){
  boot.obj <- boot(data=z,statistic=Tn,R=R,
                   sim = "permutation", sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}


m <- 500 #permutation
p <- 2 # dimension of data
n1 <- n2 <- 50 #the sample sizes
R<-999 #boot parameter
k<-3
N <- c(n1,n2)

# the function compares 3 methods 
meth3<-function(n,mu,sd,m){
  p.values <- matrix(NA,m,3)
  for(i in 1:m){
    x <- matrix(rnorm(n[1]*p,mu[1],sd[1]),ncol=p);
    y <- matrix(rnorm(n[2]*p,mu[2],sd[2]),ncol=p);
    z <- rbind(x,y)
    #NN method
    p.values[i,1] <- eqdist.nn(z,N,k)$p.value
    #energy method
    p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
    #ball method
    p.values[i,3] <- bd.test(x=x,y=y,R=999,seed=i*12345)$p.value
    }
  alpha <- 0.05
  colMeans(p.values<alpha)
}


## -----------------------------------------------------------------------------
set.seed(12345)
pow1<-meth3(N,mu=c(0,0),sd=c(1,1.5),m)
power1 <- data.frame(methods = c('NN','Energy','Ball'),pow1)
power1

## -----------------------------------------------------------------------------
set.seed(12345)
pow2<-meth3(N,mu=c(0,0.5),sd=c(1,1.5),m)
power2 <- data.frame(methods = c('NN','Energy','Ball'),pow2)
power2

## -----------------------------------------------------------------------------
set.seed(12345)
mu <- 0.5
sd <- 2
N<-c(n1,n2)
p.values <- matrix(NA,m,3)
#for mixture of two normal distributions VS normal distribution
for(i in 1:m){
  x <- matrix(rnorm(n1*p),ncol=p)#t distribution with 1 df
  y1 <- rnorm(n2*p)
  y2 <- rnorm(n2*p,mean=mu,sd=sd)
  a <- rbinom(n, 1, .5) 
  y <- matrix(a*y1 + (1-a)*y2,ncol=p)# normal distributions mixture
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,R=999,seed=i*12345)$p.value
}
alpha <- 0.05
pow31 <- colMeans(p.values<alpha)
power31 <- data.frame(methods = c('NN','Energy','Ball'),pow31)
power31

## -----------------------------------------------------------------------------
set.seed(12345)
N<-c(n1,n2)
p.values <- matrix(NA,m,3)
#t distribution with 1 df VS normal distribution
for(i in 1:m){
  x <- matrix(rt(n1*p,df=1),ncol=p)
  y <- matrix(rnorm(n2*p,sd=1.5),ncol=p)
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,R=999,seed=i*12345)$p.value
}
alpha <- 0.05
pow32 <- colMeans(p.values<alpha)
power32 <- data.frame(methods = c('NN','Energy','Ball'),pow32)
power32

## -----------------------------------------------------------------------------
set.seed(12345)
N<-c(n1,n2)
p.values <- matrix(NA,m,3)
mu<-0.5
sd<-2
#t distribution with 1 df VS mixture of two normal distributions
for(i in 1:m){
  x <- matrix(rt(n1*p,df=1),ncol=p)
  y1 <- rnorm(n2*p)
  y2 <- rnorm(n2*p,mean=mu,sd=sd)
  a <- rbinom(n, 1, .5) 
  y <- matrix(a*y1 + (1-a)*y2,ncol=p)# normal distributions mixture
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,R=999,seed=i*12345)$p.value
}
alpha <- 0.05
pow32 <- colMeans(p.values<alpha)
power32 <- data.frame(methods = c('NN','Energy','Ball'),pow32)
power32

## -----------------------------------------------------------------------------
set.seed(123)
mu <- 0.5
sd <- 1
N<-c(n1,n2*10)
p.values <- matrix(NA,m,3) 
#sample size of Y is 10 times of X
for(i in 1:m){
  x <- matrix(rnorm(n1*p),ncol=p);
  y <- matrix(rnorm(n2*10*p,mean=mu,sd=sd),ncol=p);
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value#NN method
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value#energy methods
  p.values[i,3] <- bd.test(x=x,y=y,R=999,seed=i*12345)$p.value# ball method
}
alpha <- 0.05;
pow4 <- colMeans(p.values<alpha)
power4 <- data.frame(methods = c('NN','Energy','Ball'),pow4)
power4

## -----------------------------------------------------------------------------
rw.M <- function(sigma, x0, N) {
x <- numeric(N)
x[1] <- x0
u <- runif(N)
k <- 0
for (i in 2:N) {
y <- rnorm(1, x[i-1], sigma)
if (u[i] <= (exp(abs(x[i-1])-abs(y)))) x[i] <- y 
else {
x[i] <- x[i-1]
k <- k + 1
}
}
return(list(x=x, k=k))
}

## -----------------------------------------------------------------------------
N <- 4500
sigma <- c(.05, .5, 2, 16)#different variances
x0 <- 25
k <- length(sigma)
rw<-matrix(NA,ncol = N,nrow = k)# random walk chains
accept.rate<-numeric(k)# acceptance rate
for (i in 1:length(sigma)) {
  rw[i,] <- rw.M(sigma[i], x0, N)$x
  accept.rate[i] <- (N-rw.M(sigma[i], x0, N)$k)/N
}
rbind('sigma'=sigma,'acceptance rate'=accept.rate)
for (i in 1:length(sigma)) {
  plot(1:length(rw[i,]),rw[i,],"l",ylab = 'x',xlab = 'sigma')
}

## -----------------------------------------------------------------------------
l.chain <- function(sigma, N, X1) {
#generates a Metropolis chain for the standard Laplace distribution with Normal(X[t], sigma) proposal distribution and starting value X1
x <- rep(0, N)
x[1] <- X1
u <- runif(N)
for (i in 2:N) {
xt <- x[i-1]
y <- rnorm(1, xt, sigma) #candidate point
if (u[i] <= (exp(abs(x[i-1])-abs(y))))
x[i] <- y 
else {
x[i] <- x[i-1]
k <- k + 1
}
}
return(x)
}

Gelman.Rubin <- function(psi) {
# psi[i,j] is the statistic psi(X[i,1:j])
# for chain in i-th row of X
psi <- as.matrix(psi)
n <- ncol(psi)
k <- nrow(psi)
psi.means <- rowMeans(psi) #row means
B <- n * var(psi.means) #between variance est.
psi.w <- apply(psi, 1, "var") #within variances
W <- mean(psi.w) #within est.
v.hat <- W*(n-1)/n + (B/n) #upper variance est.
r.hat <- v.hat / W #G-R statistic
return(r.hat)
}



## -----------------------------------------------------------------------------
sigma <- 1 #parameter of proposal distribution
k <- 4 #number of chains to generate
n <- 3500 #length of chains
b <- 1000 #burn-in length
#choose over dispersed initial values
x0 <- c(-10, -5, 5, 10)
#generate the chains
X <- matrix(0, nrow=k, ncol=n)
for (i in 1:k) X[i, ] <- l.chain(sigma, n, x0[i])
#compute diagnostic statistics
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi)) psi[i,] <- psi[i,] / (1:ncol(psi))
print(Gelman.Rubin(psi))
#plot psi for the four chains
for (i in 1:k)
plot(psi[i, (b+1):n], type="l",xlab=i, ylab=bquote(psi))
par(mfrow=c(1,1)) #restore default
#plot the sequence of R-hat statistics
rhat <- rep(0, n)
for (j in (b+1):n)
rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):n], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
k<-c(4:25,100,500,1000)
sk<-function(a,ki){
  pt(sqrt((a^2*(ki-1))/(ki-a^2)),df=ki-1)-pt(sqrt((a^2*ki)/(ki+1-a^2)),df=ki)
  }
out<-numeric(length(k))
for (i in 1:length(k)) {
   out[i]<-uniroot(sk, lower = 0.5, upper = 2,ki = k[i])$root
}
cbind(k,out)

## -----------------------------------------------------------------------------
set.seed(123)
N <- 5000 #max. number of iterations
pqr <- c(.5, .4, .1) #initial est. for lambdas
tol <- .Machine$double.eps^0.5
na.<-444
nb.<-132
noo<-361
nab<-63
n<-na.+nb.+noo+nab
pqr<-c(0.5,0.4,0.1)
pqrn<-pqr+1
l.pqr<-function(x){
  p<-x[1];q<-x[2];r<-x[3]
  na.*log(p^2+2*p*r)+nb.*log(q^2+2*q*r)+2*noo*log(r)+nab*log(2*p*q)
}
PQR<-c('p'=pqr[1],'q'=pqr[2],llhd=l.pqr(pqr))
for(i in 1:N){
  nao.t<-na.*(2*pqr[3]/(2*pqr[3]+pqr[1]))
  nbo.t<-nb.*(2*pqr[3]/(2*pqr[3]+pqr[2]))
  naa.t<-na.*(pqr[1]/(2*pqr[3]+pqr[1]))
  nbb.t<-nb.*(pqr[2]/(2*pqr[3]+pqr[2]))
  na.<-nao.t+naa.t
  nb.<-nbo.t+nbb.t
  pqrn[1]<-(2*na.-nao.t+nab)/(2*n)
  pqrn[2]<-(2*nb.-nbo.t+nab)/(2*n)
  pqrn[3]<-1-pqr[1]-pqr[2]
  if (sum(abs(pqr-pqrn)/pqr) < tol) break
  pqr <- pqrn
  PQR<-rbind(PQR,c(pqr[1:2],l.pqr(pqr)))
}
PQR<-cbind(1-PQR[,1]-PQR[,2],PQR)
rownames(PQR)<-c(1:i)
colnames(PQR)<-c('r','p','q','log likelihood value')
knitr::kable(PQR)

## -----------------------------------------------------------------------------
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)

for(i in 1:length(formulas)){
  out<-lm(formulas[[i]],mtcars)
  print(out)
}

lapply(formulas,function(l) lm(l,mtcars))

## -----------------------------------------------------------------------------
#using sapply
trials <- replicate(100,
                    t.test(rpois(10, 10), rpois(7, 10)),
                    simplify = FALSE)
sapply(trials, function(x) x$p.value)
#using [[ directly
sapply(trials, '[[' ,'p.value')

## -----------------------------------------------------------------------------
testset <- list(mtcars, cars)
lapply(testset, function(x) vapply(x, mean, numeric(1)))

## -----------------------------------------------------------------------------
mvapply <- function(X, f, f.value, simplify = FALSE){
  out <- Map(function(x) vapply(x, f, f.value), X)
  if(simplify == TRUE){return(simplify2array(out))}
  out<-unlist(out)
  return(out)
}
mvapply(testset, mean, numeric(1))

## -----------------------------------------------------------------------------
library(Rcpp)

## -----------------------------------------------------------------------------
rw.R <- function(sigma, x0, N) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= (exp(abs(x[i-1])-abs(y)))) x[i] <- y 
    else {
      x[i] <- x[i-1]
      k <- k + 1
    }
  }
  return(list(x=x, k=k))
}

## -----------------------------------------------------------------------------
library(microbenchmark)
N <- 8000
x0 <- 25
set.seed(123)
sigma <- .05
ts <- microbenchmark(rwR=rw.R(sigma, x0, N),rwC=rwRcpp(sigma, x0, N))
summary(ts)[,c(1,3,5,6)]
rwR=rw.R(sigma, x0, N)
rwC=rwRcpp(sigma, x0, N)
qqplot(rwR$x,rwC$x,xlab = 'rwR',ylab = 'rwRcpp',pch = 20)
abline(a=0,b=1)

## -----------------------------------------------------------------------------
set.seed(123)
sigma <- .5
ts <- microbenchmark(rwR=rw.R(sigma, x0, N),rwC=rwRcpp(sigma, x0, N))
summary(ts)[,c(1,3,5,6)]
rwR=rw.R(sigma, x0, N);rwC=rwRcpp(sigma, x0, N)
qqplot(rwR$x,rwC$x,xlab = 'rwR',ylab = 'rwRcpp',pch = 20)
abline(a=0,b=1)

## -----------------------------------------------------------------------------
set.seed(123)
sigma <- 2
ts <- microbenchmark(rwR=rw.R(sigma, x0, N),rwC=rwRcpp(sigma, x0, N))
summary(ts)[,c(1,3,5,6)]
rwR=rw.R(sigma, x0, N);rwC=rwRcpp(sigma, x0, N)
qqplot(rwR$x,rwC$x,xlab = 'rwR',ylab = 'rwRcpp',pch = 20)
abline(a=0,b=1)

## -----------------------------------------------------------------------------
set.seed(123)
sigma <- 16
ts <- microbenchmark(rwR=rw.R(sigma, x0, N),rwC=rwRcpp(sigma, x0, N))
summary(ts)[,c(1,3,5,6)]
rwR=rw.R(sigma, x0, N);rwC=rwRcpp(sigma, x0, N)
qqplot(rwR$x,rwC$x,xlab = 'rwR',ylab = 'rwRcpp',pch = 20)
abline(a=0,b=1)

