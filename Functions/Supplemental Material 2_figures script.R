##### Function to load packages

for (package in c("bda","moments","onewaytests","fGarch", "dplyr")) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}

# moments = package to compute kurtosis and skewness
# smoothmest = package to generate data from double exponential distribution
# fGarch = package to generate data from normal skewed distribution

##### Function for data generation (for different distributions)

get_sample     <- function(distName,n,              # n = sample size; distName = distribution underlying the data
                           lambda,                  # arguments for double exponential distribution
                           m, sds,                  # arguments for normal, normal skewed or double exponential distribution
                           p1,p2,sd1,sd2,           # arguments for mixed normal
                           min, max,                # arguments for unif distribution
                           df)                      # argument for chi square distribution
{
  
  if (distName=="normal"){ out <- rnorm(n, mean=m, sd=sds)}
  else if (distName=="doublex"){out <- rlap(n, mu=m, rate=1/lambda)} # transformation in order that lambda = sd
  else if (distName=="skewpos"){out <- rsnorm(n, mean=m, sd=sds,xi=10)}
  else if (distName=="skewneg"){out <- rsnorm(n, mean=m, sd=sds,xi=-10)}
  else if (distName=="unif"){out <- runif(n, min=min, max=max)}
  else if (distName=="chi2"){out <- rchisq(n, df=df)}
  else if (distName=="mixed"){out <- rmixnorm(n,p=c(p1,p2),mean=rep(m,2),sd=c(sd1,sd2))}
  return (out)
  
}

setwd("C:/Users/Marie/Dropbox/ANOVA's Welch/Appendix/Figures used in Appendixes/Appendix 3_ type 1 error rate")


## Normal distributions
A=get_sample(n=1000000,distName="normal",m=0,sds=1)
B=get_sample(n=1000000,distName="normal",m=0,sds=2)
C=get_sample(n=1000000,distName="normal",m=0,sds=4)
D=get_sample(n=1000000,distName="normal",m=0,sds=8)

png("Normal.png",width=3000,height=1200, res = 300)
  par(mai=c(.5,1,.5,1))
  plot(density(A),lty=1,xlim=c(-15,15),ylim=c(0,.7),lwd=2,main="",xlab="",cex.lab=1.2)
  lines(density(B),lty=2,lwd=2)
  lines(density(C),lty=3,lwd=2)
  lines(density(D),lty=6,lwd=2)
  legend(0,.83,legend=c("sd=1","sd=2","sd=4","sd=8"),lty=c(1,2,3,6),lwd=c(2,2,2,2),bty="n",xjust=0.5,yjust=1,horiz=TRUE,xpd=TRUE,cex=1.2)
dev.off() 

## double exponential distrbutions
A=get_sample(n=1000000,distName="doublex",m=0,lambda=1/sqrt(2))
B=get_sample(n=1000000,distName="doublex",m=0,lambda=2/sqrt(2))
C=get_sample(n=1000000,distName="doublex",m=0,lambda=4/sqrt(2))
D=get_sample(n=1000000,distName="doublex",m=0,lambda=8/sqrt(2))

png("doublex.png",width=3000,height=1200, res = 300)
par(mai=c(.5,1,.5,1))
plot(density(A),lty=1,xlim=c(-15,15),ylim=c(0,.7),lwd=2,main="",xlab="",cex.lab=1.2)
lines(density(B),lty=2,lwd=2)
lines(density(C),lty=3,lwd=2)
lines(density(D),lty=6,lwd=2)
legend(0,.83,legend=c("sd=1","sd=2","sd=4","sd=8"),lty=c(1,2,3,6),lwd=c(2,2,2,2),bty="n",xjust=0.5,yjust=1,horiz=TRUE,xpd=TRUE,cex=1.2)
dev.off() 

## mixed distrbutions
A=get_sample(n=1000000,distName="mixed",m=0,p1=.1,p2=.9,sd1=2.53,sd2=.6325)
B=get_sample(n=1000000,distName="mixed",m=0,p1=.1,p2=.9,sd1=5.06,sd2=1.265)
C=get_sample(n=1000000,distName="mixed",m=0,p1=.1,p2=.9,sd1=10.119,sd2=2.53)
D=get_sample(n=1000000,distName="mixed",m=0,p1=.1,p2=.9,sd1=20.239,sd2=5.06)

png("mixed.png",width=3000,height=1200, res = 300)
par(mai=c(.5,1,.5,1))
plot(density(A),lty=1,xlim=c(-15,15),ylim=c(0,.7),lwd=2,main="",xlab="",cex.lab=1.2)
lines(density(B),lty=2,lwd=2)
lines(density(C),lty=3,lwd=2)
lines(density(D),lty=6,lwd=2)
legend(0,.83,legend=c("sd=1","sd=2","sd=4","sd=8"),lty=c(1,2,3,6),lwd=c(2,2,2,2),bty="n",xjust=0.5,yjust=1,horiz=TRUE,xpd=TRUE,cex=1.2)
dev.off() 

## positively normal skewed
A=get_sample(n=1000000,distName="skewpos",m=0,sds=1)
B=get_sample(n=1000000,distName="skewpos",m=0,sds=2)
C=get_sample(n=1000000,distName="skewpos",m=0,sds=4)
D=get_sample(n=1000000,distName="skewpos",m=0,sds=8)

png("skewpos.png",width=3000,height=1200, res = 300)
par(mai=c(.5,1,.5,1))
plot(density(A),lty=1,xlim=c(-15,15),ylim=c(0,.7),lwd=2,main="",xlab="",cex.lab=1.2)
lines(density(B),lty=2,lwd=2)
lines(density(C),lty=3,lwd=2)
lines(density(D),lty=6,lwd=2)
legend(0,.83,legend=c("sd=1","sd=2","sd=4","sd=8"),lty=c(1,2,3,6),lwd=c(2,2,2,2),bty="n",xjust=0.5,yjust=1,horiz=TRUE,xpd=TRUE,cex=1.2)
dev.off() 

## negatively normal skewed
A=get_sample(n=1000000,distName="skewneg",m=0,sds=1)
B=get_sample(n=1000000,distName="skewneg",m=0,sds=2)
C=get_sample(n=1000000,distName="skewneg",m=0,sds=4)
D=get_sample(n=1000000,distName="skewneg",m=0,sds=8)

png("skewneg.png",width=3000,height=1200, res = 300)
par(mai=c(.5,1,.5,1))
plot(density(A),lty=1,xlim=c(-15,15),ylim=c(0,.7),lwd=2,main="",xlab="",cex.lab=1.2)
lines(density(B),lty=2,lwd=2)
lines(density(C),lty=3,lwd=2)
lines(density(D),lty=6,lwd=2)
legend(0,.83,legend=c("sd=1","sd=2","sd=4","sd=8"),lty=c(1,2,3,6),lwd=c(2,2,2,2),bty="n",xjust=0.5,yjust=1,horiz=TRUE,xpd=TRUE,cex=1.2)
dev.off() 

## chi square distribution
A=get_sample(n=1000000,distName="chi2",df=2)

png("chi2.png",width=3000,height=1200, res = 300)
par(mai=c(.5,1,.5,1))
plot(density(A),lty=1,xlim=c(-15,15),ylim=c(0,.7),lwd=2,main="",xlab="",cex.lab=1.2)
dev.off() 
