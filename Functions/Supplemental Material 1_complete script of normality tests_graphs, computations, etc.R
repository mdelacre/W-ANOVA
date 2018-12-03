for (package in c("bda","moments","fGarch","onewaytests", "dplyr","tseries","nortest","smoothmest")) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}

get_sample     <- function(distName,n,              # n = sample size; distName = distribution underlying the data
                           lambda,                  # arguments for double exponential distribution
                           m, sds,                  # arguments for normal, double exponential distribution
                           p1,p2,sd1,sd2,           # arguments for mixed normal
                           min, max,                # arguments for unif distribution
                           df)                      # argument for chi square distribution
{

  if (distName=="normal"){ out <- rnorm(n, mean=m, sd=sds)}
  else if (distName=="mixed"){out <- rmixnorm(n,p=c(p1,p2),mean=rep(m,2),sd=c(sd1,sd2))}
  else if (distName=="doublex"){out <- rdoublex(n, mu=m, lambda=lambda/sqrt(2))} # transformation in order that lambda = sd
  else if (distName=="skewpos"){out <- rsnorm(n, mean=m, sd=sds,xi=10)}
  else if (distName=="unif"){out <- runif(n, min=min, max=max)}
  else if (distName=="chi2"){out <- rchisq(n, df=df)}

  return (out)

}

get_write <- function(object,distName,n,     # object = what we want to export, k=number of groups, n = sample size; distName = distribution underlying the data
                      lambda,
                      m,sds,                  # arguments for normal distribution
                      p1,p2,sd1,sd2,          # arguments for mixed normal distribution
                      min,max,
                      df)
       {

# compute mean and standard deviation

        distr_y <- paste0("Distr=",distName,", ")
        nobs_y <- paste0("nobs=",n,", ")
        mu_y <- "mu="
        std_y <- "std="

	  if (distName=="normal"){ 
     		 mu=m
        	 std=sds
             }
	  else if (distName=="mixed"){
             mu=m
             std=round(sqrt(p1*sd1^2+p2*sd2^2),2)
     		 }
	  else if (distName=="skewpos"){
             mu=m
             std=sds
     		 }
	  else if (distName=="unif"){
      	mu=round((min+max)/2,2)
      	std=round(sqrt((min-max)^2/12),2)
     		}
  	  else if (distName=="chi2"){
      	mu=df
      	std=round(sqrt(2*df),2)
     		}

        mu_y <- paste0(mu_y,mu,", ")
        std_y <- paste0(std_y,std)

        fname <-  paste0(distr_y, nobs_y, mu_y,std_y,".rds")
        saveRDS(object, file = fname)


  }

# get_write(5,distName="normal",n=100,m=0,sds=2)
# get_write(5,distName="mixed",n=100,m=0,p1=.1,p2=.9,sd1=5.06,sd2=1.265)
# get_write(5,distName="doublex",n=100,m=0,lambda=2)
# get_write(5,distName="unif",n=100,min=-3,max=3)
# get_write(5,distName="chi2",n=100,df=2)
 

get_simu     <- function(nSims=1000000,k,distName,n,# k=number of groups, n = sample size; distName = distribution underlying the data
                           lambda,                  # arguments for double exponential distribution
                           m, sds,                  # arguments for normal, double exponential distribution
                           p1,p2,sd1,sd2,           # arguments for mixed normal distribution
                           min, max,                # arguments for uniform distribution
                           df)                      # degrees of freedom (for chi square distribution
{

    # set up empty container for all estimated parameters
    p_val<-matrix(0,nSims,5+4) # 3 colums to store the p-value.
                               # 4 columns to store mean, sd, kurtosis and skewness

    colnames(p_val) <- c("ks-testSPSS","Lillefors KS test","shapiro-test", "moments tests","Jarque-Bera test","mean","sd","skewness","kurtosis")

    # Variables generation (dependant variable and factor)
    for (i in 1:nSims){

      # Sample as a function of the sample size and distribution
      y <- get_sample(distName,n,lambda,m,sds,p1,p2,sd1,sd2,min,max,df)

      # normality tests
      
      p_val[i,1] <- ks.test(y,"pnorm",mean(y),sd(y))$p.value
      p_val[i,2] <- lillie.test(y)$p.value # lillefors correction (as in SPSS)
      p_val[i,3] <- shapiro.test(y)$p.value

      # Computing moments and SE of moments, as in SPSS
      y_skewness <- (n^2*moment(y, order=3,central=T))/((n-1)*(n-2)*sd(y)^3)
      y_kurtosis <- (n^2*(n+1)*(moment(y, order=4,central=T)))/((n-1)*(n-2)*(n-3)*sd(y)^4)-(3*(n-1)^2)/((n-2)*(n-3))

      SE_skewness <- sqrt((6*(n)*(n-1))/((n-2)*(n+1)*(n+3)))
      SE_kurtosis <- sqrt((4*(n^2-1)*SE_skewness^2)/((n-3)*(n+5)))

      z_skewness <- y_skewness/SE_skewness
      z_kurtosis <- y_kurtosis/SE_kurtosis

      if(abs(z_skewness)>1.96 | abs(z_kurtosis)>1.96) {p_val[i,4] <-.01}  else {p_val[i,4] <-.99}        

      p_val[i,5] <- jarque.bera.test(y)$p.value
      
          ### descriptives
      p_val[i,6:9]=c(mean(y),sd(y),y_skewness,y_kurtosis)
                          }

  # p-values extraction in a specified file (thx to ? setwd ?)
    setwd("C:/Users/mdela/Dropbox/ANOVA's Welch/Appendix/Appendix 1_KS test/KS test")
    get_write(p_val,distName,n,lambda,m,sds,p1,p2,sd1,sd2,min,max,df)

}

#----------------------------------------------------------------------------------------

read_sample=function(ssdossier="",distName,n,m,sds){

  # compute mean and standard deviation
  distr <-  paste0("Distr=",distName,", ")
  nobs <- paste0("nobs=",n,", ")
  means <-  paste0("mu=",m,", ")
  stdevs <- paste0("std=", sds)
   
  fname <-  paste0(distr, nobs, means, stdevs)

  setwd(dir="E:/Welch's W ANOVA/Outputs Appendix/Appendix 1_normality tests/") # destination file

  data=readRDS(file = paste0(fname,".rds"))
  return(data)
}



#----------------------------------------------------------------------------------------

graphique=function(variable1,variable2,variable3,variable4,variable5,name){
setwd("C:/Users/mdela/Dropbox/ANOVA's Welch/Appendix/KS test")
png(paste0(name,".png"),width=800,height=1500, res = 300)
par(mfrow=c(5,1))
hist(variable1[,1],main="N=10",xlab="p-value",cex.lab=1.5,cex.main=1.5)
hist(variable2[,1],main="N=20",xlab="p-value",cex.lab=1.5,cex.main=1.5)
hist(variable3[,1],main="N=30",xlab="p-value",cex.lab=1.5,cex.main=1.5)
hist(variable4[,1],main="N=50",xlab="p-value",cex.lab=1.5,cex.main=1.5)
hist(variable5[,1],main="N=100",xlab="p-value",cex.lab=1.5,cex.main=1.5)
dev.off() 
}

#----------------------------------------------------------------------------------------
alpha_fct=function(variable1,variable2,variable3,variable4,variable5){
Alpha=matrix(0,5,4)
Alpha[1,1]<-round(sum(variable1[,1]<.05)/length(variable1[,1]),2)
Alpha[2,1]<-round(sum(variable2[,1]<.05)/length(variable2[,1]),2)
Alpha[3,1]<-round(sum(variable3[,1]<.05)/length(variable3[,1]),2)
Alpha[4,1]<-round(sum(variable4[,1]<.05)/length(variable4[,1]),2)
Alpha[5,1]<-round(sum(variable5[,1]<.05)/length(variable5[,1]),2)
Alpha[1,2]<-round(sum(variable1[,2]<.05)/length(variable1[,2]),2)
Alpha[2,2]<-round(sum(variable2[,2]<.05)/length(variable2[,2]),2)
Alpha[3,2]<-round(sum(variable3[,2]<.05)/length(variable3[,2]),2)
Alpha[4,2]<-round(sum(variable4[,2]<.05)/length(variable4[,2]),2)
Alpha[5,2]<-round(sum(variable5[,2]<.05)/length(variable5[,2]),2)
Alpha[1,3]<-round(sum(variable1[,3]<.05)/length(variable1[,3]),2)
Alpha[2,3]<-round(sum(variable2[,3]<.05)/length(variable2[,3]),2)
Alpha[3,3]<-round(sum(variable3[,3]<.05)/length(variable3[,3]),2)
Alpha[4,3]<-round(sum(variable4[,3]<.05)/length(variable4[,3]),2)
Alpha[5,3]<-round(sum(variable5[,3]<.05)/length(variable5[,3]),2)
Alpha[1,4]<-round(sum(variable1[,5]<.05)/length(variable1[,5]),2)
Alpha[2,4]<-round(sum(variable2[,5]<.05)/length(variable2[,5]),2)
Alpha[3,4]<-round(sum(variable3[,5]<.05)/length(variable3[,5]),2)
Alpha[4,4]<-round(sum(variable4[,5]<.05)/length(variable4[,5]),2)
Alpha[5,4]<-round(sum(variable5[,5]<.05)/length(variable5[,5]),2)

return(Alpha)
}

#----------------------------------------------------------------------------------------
#                                        APPLICATIONS

get_simu(distName="normal",n=10,m=0,sds=1)
get_simu(distName="normal",n=20,m=0,sds=1)
get_simu(distName="normal",n=30,m=0,sds=1)
get_simu(distName="normal",n=50,m=0,sds=1)
get_simu(distName="normal",n=100,m=0,sds=1)
get_simu(distName="normal",n=1000,m=0,sds=1)

get_simu(distName="mixed",n=10,m=0,p1=.1,p2=.9,sd1=2.53,sd2=.6325)
get_simu(distName="mixed",n=20,m=0,p1=.1,p2=.9,sd1=2.53,sd2=.6325)
get_simu(distName="mixed",n=30,m=0,p1=.1,p2=.9,sd1=2.53,sd2=.6325)
get_simu(distName="mixed",n=50,m=0,p1=.1,p2=.9,sd1=2.53,sd2=.6325)
get_simu(distName="mixed",n=100,m=0,p1=.1,p2=.9,sd1=2.53,sd2=.6325)

get_simu(distName="unif",n=10,min=-1.7325,max=1.7325)
get_simu(distName="unif",n=20,min=-1.7325,max=1.7325)
get_simu(distName="unif",n=30,min=-1.7325,max=1.7325)
get_simu(distName="unif",n=50,min=-1.7325,max=1.7325)
get_simu(distName="unif",n=100,min=-1.7325,max=1.7325)

get_simu(distName="skewpos",n=10,m=0,sds=1)
get_simu(distName="skewpos",n=20,m=0,sds=1)
get_simu(distName="skewpos",n=30,m=0,sds=1)
get_simu(distName="skewpos",n=50,m=0,sds=1)
get_simu(distName="skewpos",n=100,m=0,sds=1)

get_simu(distName="chi2",n=10,df=0.5)
get_simu(distName="chi2",n=20,df=0.5)
get_simu(distName="chi2",n=30,df=0.5)
get_simu(distName="chi2",n=50,df=0.5)
get_simu(distName="chi2",n=100,df=0.5)

N_n10_sd1<-read_sample(distName="normal",n=10,m=0,sds=1)
N_n20_sd1<-read_sample(distName="normal",n=20,m=0,sds=1)
N_n30_sd1<-read_sample(distName="normal",n=30,m=0,sds=1)
N_n50_sd1<-read_sample(distName="normal",n=50,m=0,sds=1)
N_n100_sd1<-read_sample(distName="normal",n=100,m=0,sds=1)
N_n1000_sd1<-read_sample(distName="normal",n=1000,m=0,sds=1)


MX_n10_sd1<-read_sample(distName="mixed",n=10,m=0,sds=1)
MX_n20_sd1<-read_sample(distName="mixed",n=20,m=0,sds=1)
MX_n30_sd1<-read_sample(distName="mixed",n=30,m=0,sds=1)
MX_n50_sd1<-read_sample(distName="mixed",n=50,m=0,sds=1)
MX_n100_sd1<-read_sample(distName="mixed",n=100,m=0,sds=1)

U_n10_sd1<-read_sample(distName="unif",n=10,m=0,sds=1)
U_n20_sd1<-read_sample(distName="unif",n=20,m=0,sds=1)
U_n30_sd1<-read_sample(distName="unif",n=30,m=0,sds=1)
U_n50_sd1<-read_sample(distName="unif",n=50,m=0,sds=1)
U_n100_sd1<-read_sample(distName="unif",n=100,m=0,sds=1)

NS_n10_sd1<-read_sample(distName="skewpos",n=10,m=0,sds=1)
NS_n20_sd1<-read_sample(distName="skewpos",n=20,m=0,sds=1)
NS_n30_sd1<-read_sample(distName="skewpos",n=30,m=0,sds=1)
NS_n50_sd1<-read_sample(distName="skewpos",n=50,m=0,sds=1)
NS_n100_sd1<-read_sample(distName="skewpos",n=100,m=0,sds=1)

K_n10_sd1<-read_sample(distName="chi2",n=10,m=0.5,sds=1)
K_n20_sd1<-read_sample(distName="chi2",n=20,m=0.5,sds=1)
K_n30_sd1<-read_sample(distName="chi2",n=30,m=0.5,sds=1)
K_n50_sd1<-read_sample(distName="chi2",n=50,m=0.5,sds=1)
K_n100_sd1<-read_sample(distName="chi2",n=100,m=0.5,sds=1)

skew=function(variable){
plot(density(variable[,5]),main=paste0("kurtosis=",round(kurtosis(variable[,6]),2),", skewness=",round(skewness(variable[,6]),2)),xlab=round(sd(variable[,6]),3))
}

kurt=function(variable){
par(mfrow=c(2,1))
plot(density(variable[,7]),main=paste0("kurtosis=",round(kurtosis(variable[,7]),2),", skewness=",round(skewness(variable[,7]),2)))
new=variable[,7]/SE_kurtosis
plot(density(new),main=paste0("kurtosis=",round(kurtosis(new),2),", skewness=",round(skewness(new),2)))
  }


moy=function(variable){
  par(mfrow=c(2,1))
  plot(density(variable[,4]),main=paste0("kurtosis=",round(kurtosis(variable[,4]),2),", skewness=",round(skewness(variable[,4]),2)))
}


moy(N_n10_sd1)
kurt(N_n1000_sd1)

View(N_n100_sd1)
n=10
SE_skewness <- round(sqrt((6*(n)*(n-1))/((n-2)*(n+1)*(n+3))),3)
SE_kurtosis <- round(sqrt((4*(n^2-1)*SE_skewness^2)/((n-3)*(n+5))),3)

variable=N_n10_sd1[,7]
new=variable/SE_kurtosis
plot(density(new))
kurtosis(new)
skewness(new)

kurt(N_n10_sd1)

round(mean(variable[,3]),3)
sd(variable[,3])

round(mean(variable[,4]),3)

variable=N_n1000_sd1
View(variable)

mean(variable[,5])
mean(variable[,6])

sd(variable[,5])
sd(variable[,6])


skew(N_n30_sd1)
skew(N_n40_sd1)
skew(N_n50_sd1)
skew(N_n100_sd1)

kurt(N_n20_sd1) 
kurt(N_n30_sd1)
kurt(N_n40_sd1)
kurt(N_n50_sd1)
kurt(N_n100_sd1)


graphique(N_n10_sd1,N_n20_sd1,N_n30_sd1,N_n50_sd1,N_n100_sd1,"normal")
graphique(MX_n10_sd1,MX_n20_sd1,MX_n30_sd1,MX_n50_sd1,MX_n100_sd1,"mixed")
graphique(U_n10_sd1,U_n20_sd1,U_n30_sd1,U_n50_sd1,U_n100_sd1,"uniform")
graphique(NS_n10_sd1,NS_n20_sd1,NS_n30_sd1,NS_n50_sd1,NS_n100_sd1,"normal skewed")
graphique((K_n10_sd1,K_n20_sd1,K_n30_sd1,K_n50_sd1,K_n100_sd1,"chi²")
          
  
alpha_fct(N_n10_sd1,N_n20_sd1,N_n30_sd1,N_n50_sd1,N_n100_sd1)
alpha_fct(MX_n10_sd1,MX_n20_sd1,MX_n30_sd1,MX_n50_sd1,MX_n100_sd1)
alpha_fct(U_n10_sd1,U_n20_sd1,U_n30_sd1,U_n50_sd1,U_n100_sd1)
alpha_fct(NS_n10_sd1,NS_n20_sd1,NS_n30_sd1,NS_n50_sd1,NS_n100_sd1)
alpha_fct(K_n10_sd1,K_n20_sd1,K_n30_sd1,K_n50_sd1,K_n100_sd1)

setwd("C:/Users/mdela/Dropbox/ANOVA's Welch/Appendix/Appendix 1_KS test/KS test")
A=rnorm(1000000,0,1)
png("stdnormal.png",width=2000,height=1200, res = 300)
plot(density(A),xlab="",main="")
dev.off()

B=get_sample(distName="mixed",n=1000000,m=0,p1=.1,p2=.9,sd1=2.53,sd2=.6325)
png("mixed.png",width=2000,height=1200, res = 300)
par(xpd=T)
plot(density(B),xlab="",main="")
lines(density(A), lty=2)
legend(-6,.7,lty=c(2,1),legend=c("normal","mixed normal"),horiz=T,bty="n")
dev.off()

C=get_sample(distName="unif",n=10000000,min=-1.7325,max=1.7325)
png("unif.png",width=2000,height=1200, res = 300)
par(xpd=T)
plot(density(A),xlab="",main="", lty=2)
lines(density(C))
legend(-3,.5,lty=c(2,1),legend=c("normal","uniform"),horiz=T,bty="n")
dev.off()


D=get_sample(distName="skewpos",n=10000000,m=0,sds=1)
png("skewpos.png",width=2000,height=1200, res = 300)
par(xpd=T)
plot(density(D),xlab="",main="",xlim=c(-6,6))
lines(density(A), lty=2)
legend(-4.4,.6,lty=c(2,1),legend=c("normal","right skewed normal"),horiz=T,bty="n")
dev.off()

E=get_sample(distName="chi2",n=1000000,df=0.5)
png("chi2.png",width=2000,height=1200, res = 300)
par(xpd=F)
plot(density(E),xlab="",main="",ylim=c(0,.962),xlim=c(-5,5))
lines(density(A), lty=2)
par(xpd=T)
legend(-4,1.2,lty=c(2,1),legend=c("normal","right skewed normal"),horiz=T,bty="n")
dev.off()
