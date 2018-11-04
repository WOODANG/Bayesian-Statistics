#gibbs.R 파일의 코드를 이용한 깁스 샘플링 연습
#정규분포 말고 다른 분포, 이항분포나 포아송분포 깁스 샘플링도 연습해보기

library(MCMCpack)

#joint prior with mu proportional to 1, sigma^2 proportional to 1/sigma^2
#mu prior and sigma prior are independent
#data follows normal distribution with ybar and sigma^2/n
#in this case, we can't conduct Bayesian inference only using joint posterior because it is too complex to use
#so we're gonna use conditional posteior distributions of each parameters.
#originally we have to use marginal posterior distributions but marginal posterior is too hard to calculate integrals

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

######another gibbs sampling
#bwt in baby data, joint prior, parameters are independent
#using Normal priors for mean and Inverse Gamma priors for sd

baby=read.table("http://www.stat.berkeley.edu/~statlabs/data/babies.data",head=TRUE)
attach(baby)

hist(bwt)
mean(bwt)
sd(bwt)

ybar <- mean(bwt)
s2<-(sd(bwt))^2
n<-length(bwt)
B<-1000 #burnin

mu <- sigma2 <- rep(NA,5000)
sigma2[1]<-0.00001
for ( i in 2:5000){
  mu[i] <- rnorm(1,mean=ybar,sd=sigma2[i-1]/sqrt(n))
  sigma2[i]<-rinvgamma(1,shape=n/2,scale=((n-1)*s2+n*(ybar-mu[i])^2)/2)
}

mu2<-mu[-(1:B)]
sigma22<-sigma2[-(1:B)]
plot(mu2,type="l")

mu[1]<-1000
for ( i in 2:5000){
  sigma2[i]<-rinvgamma(1,shape=n/2,scale=((n-1)*s2+n*(ybar-mu[i-1])^2)/2)
  mu[i] <- rnorm(1,mean=ybar,sd=sigma2[i]/sqrt(n))
}

###tried to converge mu and sigma2 but it didn't work
#each iteration value is going bigger, reaching Infinite at some iterations.
#Thus we can't use MCMC because mu and sigma2 is NaN due to Inf.

#Data in Binomial with Prior in Beta -> Posterior in Beta. let's do sampling from posterior distribution.
install.packages("statip")
library(statip)

a <- 10
N <- 10
theta <- 0.3

x <- rbinom(a,N,theta)
alpha <- 2
beta <- 2
post.alpha <- alpha+sum(x)
post.beta <- beta-sum(x)+a*N
grid2<-seq(0,1,by=0.001)
theta.norm <- rbeta(grid2,post.alpha,post.beta)
hist(theta.norm,freq=F)
lines(density(theta.norm),col="red",lwd=2) ##histrogram and density plot
mean(theta.norm) #similar to theta=0.3

#Data in Poisoon with Prior in Gamma -> Posterior in Gamma. let's do sampling from posterior distribution.

lambda <- 2
x2 <- rpois(100,lambda)
alpha<-3
beta<-3
post.alpha2 <- alpha+sum(x2)
post.beta2 <- 1/(a+(1/beta))
lambda.gamma<-rgamma(1000,post.alpha2,post.beta2)
hist(lambda.gamma,freq=F)
lines(density(lambda.gamma)) #if you can't draw density distribution directly on a histogram with lines function, run plot(density()) first and then run lines(density())
mean(lambda.gamma) #why mean seems so large? because posterior mean includes n, wich is so large
gamma23<-rgamma(1000,2,2)
plot(density(gamma23))
plot(density(x2))
plot(density(lambda.gamma)) #parameter distribution moves a little bit to the center.

#I need to explain exactly about posterior normal, beta, gamma.
#For Binomial and Poisson case, posterior distribution is about data's parameter, and parameter is single value.
#so it is not complicated as we think, because we can just use one parameter to interpret posterior distribution.
#but when data follows multiparameter distribution such as Normal, we need joint posterior distribution, and interpretation becomes more difficult.
#since data follows multiparameter distribution we need each parameter's marginal distribution. however, some marginal distribution are mathmatically difficult to calculate.
#in this case we use parameter's conditional distribution instead of marginal distribution.
#when conditional(or marginal) distribution follows certain known distribution, we can use Gibbs Sampling to display the characteristics of parameter.
#but if even conditional parameter doesn't follow certain known distribution, which means there is no closed form in joint posterior distribution, we can use Metropolis-Hastings Algorithm to get parameter's characteristics.

