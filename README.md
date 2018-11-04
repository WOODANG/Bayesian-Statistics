# Bayesian-Statistics
This repository is made for studying Bayesian Statistics.

I want to upload my Bayesian R code.
a graduate student, majoring in Satistics at Chung-Ang University.
Interested in studying Machine Learning, Beginner of Github. have a lot of things to learn about github, statistics, and machine learning.

#1. Normal distribution model with Gibbs sampling, assuming mu and sigma prior are independent

library(MCMCpack)

ybar=200
s2=5
n=200
B=1000
mu<-sigma2<-rep(NA,10000)
sigma2[1]<-7
for (i in 2:10000){
  mu[i] = rnorm(1,mean=ybar,sd=sqrt(sigma2[i-1]/n))
  sigma2[i]=rinvgamma(1,shape=n/2,scale=((n-1)*s2+n*(mu[i]-ybar)^2)/2)
}

mu <- mu[-(1:B)]
sigma2<-sigma2[-(1:B)]
par(mfrow=c(1,1))
plot(mu,type="l")
plot(sigma2,type="l")
mean(mu)
mean(sigma2)
