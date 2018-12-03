##### Function to load packages

for (package in c("bda","moments","onewaytests","smoothmest","fGarch", "dplyr")) {
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
  else if (distName=="doublex"){out <- rdoublex(n, mu=m, lambda=lambda/sqrt(2))} # transformation in order that lambda = sd
  else if (distName=="doublex_SP"){out <- rdoublex(n, mu=m, lambda=lambda)} # transformation in order that lambda = sd
  else if (distName=="skewpos"){out <- rsnorm(n, mean=m, sd=sds,xi=10)}
  else if (distName=="skewneg"){out <- rsnorm(n, mean=m, sd=sds,xi=-10)}
  else if (distName=="unif"){out <- runif(n, min=min, max=max)}
  else if (distName=="chi2"){out <- rchisq(n, df=df)}
  else if (distName=="mixed"){out <- rmixnorm(n,p=c(p1,p2),mean=rep(m,2),sd=c(sd1,sd2))}
  return (out)

}

# example: get_sample(n=1000, distName="normal", m=0, sds=2) will generate 1000 observations extracted from a N(0,2)

##### Function to export results in a .rds file with customized name
#### name containing informations about n, sd and mean of each group, and distribution underlying data for each group

get_write <- function(object,k,distName,n,    # object = what we want to export, k=number of groups, n = sample size; distName = distribution underlying the data
                      lambda,                 # arguments for double exponential distribution
                      m, sds,                 # arguments for normal, normal skewed or double exponential distribution
                      p1,p2,sd1,sd2,          # arguments for mixed normal distribution
                      min, max,               # arguments for unif distribution
                      df)                     # degrees of freedom (argument for chi square distribution)

       {

# compute mean and standard deviation
  mu <- rep(0, k)
  std <- rep(0,k)
  distr <-  "Distr=["
  nobs <- "n=["
  means <-  "means=["
  stdevs <- "sds=["
  
  for (j in 1:k){
	  if (distName[j]=="normal"){ 
     		 mu[j]=m[j]
        	 std[j]=sds[j]
             }
	  else if (distName[j]=="doublex"){
             mu[j]=m[j]
             std[j]=lambda[j]
     		 }
	  else if (distName[j]=="doublex_SP"){
             mu[j]=m[j]
             std[j]=lambda[j]
     		 }
	  else if (distName[j]=="skewpos"){
      	 mu[j]=m[j]
      	 std[j]=sds[j]
     		}
	  else if (distName[j]=="skewneg"){
      	 mu[j]=m[j]
      	 std[j]=sds[j]
     		}
	  else if (distName[j]=="unif"){
      	mu[j]=round((min[j]+max[j])/2,2)
      	std[j]=round(sqrt((min[j]-max[j])^2/12),2)
     		}
  	  else if (distName[j]=="chi2"){
      	mu[j]=df[j]
      	std[j]=round(sqrt(2*df[j]),2)
     		}
	  else if (distName[j]=="mixed"){
             mu[j]=m[j]
             std[j]=round(sqrt(p1[j]*sd1[j]^2+p2[j]*sd2[j]^2),2)
     		 }
	
      if (j == k){
          distr <- paste0(distr, distName[j], "]")
          nobs <- paste0(nobs, n[j], "]")
          means <- paste0(means, mu[j], "]")
          stdevs <- paste0(stdevs, std[j], "]")
      } else {
          distr <- paste0(distr, distName[j], ",")
          nobs <- paste0(nobs, n[j], ",")
          means <- paste0(means, mu[j], ",")
          stdevs <- paste0(stdevs, std[j], ",")
      }
  }
  
  fname <-  paste(distr, nobs, means, stdevs, sep=",")
  fname <-  paste0(fname, ".rds")
  saveRDS(object, file = fname)
}
  

##### Function to generate dataset,realize statistical test and compute (and extract) p-value

get_simu     <- function(nSims=1000000,k,distName,n,        # k=number of groups, n = sample size; distName = distribution underlying the data
                           lambda,                  # arguments for double exponential distribution
                           m, sds,                # arguments for normal, normal skewed or double exponential distribution
                           p1,p2,sd1,sd2,           # arguments for mixed normal distribution
                           min, max,                # arguments for uniform distribution
                           df)                      # degrees of freedom (for chi square distribution

{

     # set up empty container for all estimated parameters
     p_val<-matrix(0,nSims,3+2*k) # Three colums to store the p-value of ANOVA F-test, F*-test and W-test
                                  # 2*k columns to store mean and sd of each sample
                                  # example: if there are 3 groups, the function will store mean_sample1 in c4, mean_sample2 in c5, and so on, until sd_sample3 in c9.

    descr_name=expand.grid(paste("ech",1:k),c("mean","sd"))
    desc_name=paste(descr_name[,2],descr_name[,1])
    colnames(p_val) <- c("Oneway", "Welch", "Anova B-F",desc_name)

    # Variables generation (dependant variable and factor)
    for (i in 1:nSims){

      # Dataset as a function of the number of groups (k)
      y <- c(); factor <- c()
      for (s in 1:k){
          y <- c(y, 
                 get_sample(distName[s],n[s],lambda[s],m[s],sds[s],p1[s],p2[s],sd1[s],sd2[s],min[s], max[s],df[s])
          )
          factor <- c(factor, rep(s, n[s]))
      }
      
      # Data construction
      data<-cbind(factor = factor(factor), y)


       # between group tests
       
       p_val[i,1] <- oneway.test(y ~ factor, data=data, var.equal=TRUE)$p.value # extract the p-value of the classical ANOVA F-test
       p_val[i,2] <- oneway.test(y ~ factor, data=data, var.equal=FALSE)$p.value # extract the p-value of the Welch's W ANOVA test
       p_val[i,3] <- bf.test(y ~ factor,data=data, verbose = F)$p.value # extract the p-value of the Brown-Forsythe F*-test
   
       p_val[i,4:(4+2*k)] <- c(tapply(mean,factor),tapply(sd,factor))
 
                    }

   # p-values extraction in a specified file 
    setwd(dir="C:/Users/Marie/Documents/SIMU/Manquant annexes") # destination file  
    get_write(p_val,k,distName,n,lambda,m,sds,p1,p2,sd1,sd2,min,max,df) # generate the .rds file (see function above)

}

#####################################################################################
#############                       APPLICATIONS                        #############
#####################################################################################

#############                    Type 1 error rate                      #############

### normal distributions

get_simu(k=2,distName=rep("normal",2),n=c(20,10),m=c(0,0),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(20,10),m=c(0,0),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(20,10),m=c(0,0),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(20,10),m=c(0,0),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(20,20),m=c(0,0),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(20,20),m=c(0,0),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(20,20),m=c(0,0),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(20,20),m=c(0,0),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(20,30),m=c(0,0),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(20,30),m=c(0,0),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(20,30),m=c(0,0),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(20,30),m=c(0,0),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(20,40),m=c(0,0),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(20,40),m=c(0,0),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(20,40),m=c(0,0),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(20,40),m=c(0,0),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(30,15),m=c(0,0),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(30,15),m=c(0,0),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(30,15),m=c(0,0),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(30,15),m=c(0,0),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(30,30),m=c(0,0),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(30,30),m=c(0,0),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(30,30),m=c(0,0),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(30,30),m=c(0,0),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(30,45),m=c(0,0),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(30,45),m=c(0,0),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(30,45),m=c(0,0),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(30,45),m=c(0,0),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(30,60),m=c(0,0),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(30,60),m=c(0,0),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(30,60),m=c(0,0),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(30,60),m=c(0,0),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(40,20),m=c(0,0),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(40,20),m=c(0,0),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(40,20),m=c(0,0),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(40,20),m=c(0,0),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(40,40),m=c(0,0),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(40,40),m=c(0,0),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(40,40),m=c(0,0),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(40,40),m=c(0,0),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(40,60),m=c(0,0),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(40,60),m=c(0,0),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(40,60),m=c(0,0),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(40,60),m=c(0,0),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(40,80),m=c(0,0),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(40,80),m=c(0,0),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(40,80),m=c(0,0),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(40,80),m=c(0,0),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(50,25),m=c(0,0),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(50,25),m=c(0,0),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(50,25),m=c(0,0),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(50,25),m=c(0,0),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(50,50),m=c(0,0),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(50,50),m=c(0,0),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(50,50),m=c(0,0),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(50,50),m=c(0,0),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(50,75),m=c(0,0),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(50,75),m=c(0,0),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(50,75),m=c(0,0),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(50,75),m=c(0,0),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(50,100),m=c(0,0),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(50,100),m=c(0,0),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(50,100),m=c(0,0),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(50,100),m=c(0,0),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(100,50),m=c(0,0),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(100,50),m=c(0,0),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(100,50),m=c(0,0),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(100,50),m=c(0,0),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(100,100),m=c(0,0),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(100,100),m=c(0,0),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(100,100),m=c(0,0),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(100,100),m=c(0,0),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(100,150),m=c(0,0),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(100,150),m=c(0,0),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(100,150),m=c(0,0),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(100,150),m=c(0,0),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(100,200),m=c(0,0),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(100,200),m=c(0,0),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(100,200),m=c(0,0),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(100,200),m=c(0,0),sds=c(2,8))

get_simu(k=3,distName=rep("normal",3),n=c(20,20,10),m=c(0,0,0),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(20,20,10),m=c(0,0,0),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(20,20,10),m=c(0,0,0),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(20,20,10),m=c(0,0,0),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(20,20,20),m=c(0,0,0),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(20,20,20),m=c(0,0,0),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(20,20,20),m=c(0,0,0),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(20,20,20),m=c(0,0,0),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(20,20,30),m=c(0,0,0),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(20,20,30),m=c(0,0,0),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(20,20,30),m=c(0,0,0),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(20,20,30),m=c(0,0,0),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(20,20,40),m=c(0,0,0),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(20,20,40),m=c(0,0,0),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(20,20,40),m=c(0,0,0),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(20,20,40),m=c(0,0,0),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(30,30,15),m=c(0,0,0),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(30,30,15),m=c(0,0,0),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(30,30,15),m=c(0,0,0),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(30,30,15),m=c(0,0,0),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(30,30,30),m=c(0,0,0),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(30,30,30),m=c(0,0,0),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(30,30,30),m=c(0,0,0),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(30,30,30),m=c(0,0,0),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(30,30,45),m=c(0,0,0),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(30,30,45),m=c(0,0,0),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(30,30,45),m=c(0,0,0),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(30,30,45),m=c(0,0,0),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(30,30,60),m=c(0,0,0),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(30,30,60),m=c(0,0,0),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(30,30,60),m=c(0,0,0),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(30,30,60),m=c(0,0,0),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(40,40,20),m=c(0,0,0),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(40,40,20),m=c(0,0,0),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(40,40,20),m=c(0,0,0),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(40,40,20),m=c(0,0,0),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(40,40,40),m=c(0,0,0),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(40,40,40),m=c(0,0,0),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(40,40,40),m=c(0,0,0),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(40,40,40),m=c(0,0,0),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(40,40,60),m=c(0,0,0),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(40,40,60),m=c(0,0,0),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(40,40,60),m=c(0,0,0),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(40,40,60),m=c(0,0,0),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(40,40,80),m=c(0,0,0),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(40,40,80),m=c(0,0,0),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(40,40,80),m=c(0,0,0),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(40,40,80),m=c(0,0,0),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(50,50,25),m=c(0,0,0),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(50,50,25),m=c(0,0,0),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(50,50,25),m=c(0,0,0),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(50,50,25),m=c(0,0,0),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(50,50,50),m=c(0,0,0),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(50,50,50),m=c(0,0,0),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(50,50,50),m=c(0,0,0),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(50,50,50),m=c(0,0,0),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(50,50,75),m=c(0,0,0),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(50,50,75),m=c(0,0,0),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(50,50,75),m=c(0,0,0),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(50,50,75),m=c(0,0,0),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(50,50,100),m=c(0,0,0),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(50,50,100),m=c(0,0,0),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(50,50,100),m=c(0,0,0),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(50,50,100),m=c(0,0,0),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(100,100,50),m=c(0,0,0),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(100,100,50),m=c(0,0,0),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(100,100,50),m=c(0,0,0),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(100,100,50),m=c(0,0,0),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(100,100,100),m=c(0,0,0),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(100,100,100),m=c(0,0,0),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(100,100,100),m=c(0,0,0),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(100,100,100),m=c(0,0,0),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(100,100,150),m=c(0,0,0),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(100,100,150),m=c(0,0,0),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(100,100,150),m=c(0,0,0),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(100,100,150),m=c(0,0,0),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(100,100,200),m=c(0,0,0),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(100,100,200),m=c(0,0,0),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(100,100,200),m=c(0,0,0),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(100,100,200),m=c(0,0,0),sds=c(2,2,8))

get_simu(k=4,distName=rep("normal",4),n=c(20,20,20,10),m=c(0,0,0,0),sds=c(2,2,2,1))
get_simu(k=4,distName=rep("normal",4),n=c(20,20,20,10),m=c(0,0,0,0),sds=c(2,2,2,2))
get_simu(k=4,distName=rep("normal",4),n=c(20,20,20,10),m=c(0,0,0,0),sds=c(2,2,2,4))
get_simu(k=4,distName=rep("normal",4),n=c(20,20,20,10),m=c(0,0,0,0),sds=c(2,2,2,8))
get_simu(k=4,distName=rep("normal",4),n=c(20,20,20,20),m=c(0,0,0,0),sds=c(2,2,2,1))
get_simu(k=4,distName=rep("normal",4),n=c(20,20,20,20),m=c(0,0,0,0),sds=c(2,2,2,2))
get_simu(k=4,distName=rep("normal",4),n=c(20,20,20,20),m=c(0,0,0,0),sds=c(2,2,2,4))
get_simu(k=4,distName=rep("normal",4),n=c(20,20,20,20),m=c(0,0,0,0),sds=c(2,2,2,8))
get_simu(k=4,distName=rep("normal",4),n=c(20,20,20,30),m=c(0,0,0,0),sds=c(2,2,2,1))
get_simu(k=4,distName=rep("normal",4),n=c(20,20,20,30),m=c(0,0,0,0),sds=c(2,2,2,2))
get_simu(k=4,distName=rep("normal",4),n=c(20,20,20,30),m=c(0,0,0,0),sds=c(2,2,2,4))
get_simu(k=4,distName=rep("normal",4),n=c(20,20,20,30),m=c(0,0,0,0),sds=c(2,2,2,8))
get_simu(k=4,distName=rep("normal",4),n=c(20,20,20,40),m=c(0,0,0,0),sds=c(2,2,2,1))
get_simu(k=4,distName=rep("normal",4),n=c(20,20,20,40),m=c(0,0,0,0),sds=c(2,2,2,2))
get_simu(k=4,distName=rep("normal",4),n=c(20,20,20,40),m=c(0,0,0,0),sds=c(2,2,2,4))
get_simu(k=4,distName=rep("normal",4),n=c(20,20,20,40),m=c(0,0,0,0),sds=c(2,2,2,8))
get_simu(k=4,distName=rep("normal",4),n=c(30,30,30,15),m=c(0,0,0,0),sds=c(2,2,2,1))
get_simu(k=4,distName=rep("normal",4),n=c(30,30,30,15),m=c(0,0,0,0),sds=c(2,2,2,2))
get_simu(k=4,distName=rep("normal",4),n=c(30,30,30,15),m=c(0,0,0,0),sds=c(2,2,2,4))
get_simu(k=4,distName=rep("normal",4),n=c(30,30,30,15),m=c(0,0,0,0),sds=c(2,2,2,8))
get_simu(k=4,distName=rep("normal",4),n=c(30,30,30,30),m=c(0,0,0,0),sds=c(2,2,2,1))
get_simu(k=4,distName=rep("normal",4),n=c(30,30,30,30),m=c(0,0,0,0),sds=c(2,2,2,2))
get_simu(k=4,distName=rep("normal",4),n=c(30,30,30,30),m=c(0,0,0,0),sds=c(2,2,2,4))
get_simu(k=4,distName=rep("normal",4),n=c(30,30,30,30),m=c(0,0,0,0),sds=c(2,2,2,8))
get_simu(k=4,distName=rep("normal",4),n=c(30,30,30,45),m=c(0,0,0,0),sds=c(2,2,2,1))
get_simu(k=4,distName=rep("normal",4),n=c(30,30,30,45),m=c(0,0,0,0),sds=c(2,2,2,2))
get_simu(k=4,distName=rep("normal",4),n=c(30,30,30,45),m=c(0,0,0,0),sds=c(2,2,2,4))
get_simu(k=4,distName=rep("normal",4),n=c(30,30,30,45),m=c(0,0,0,0),sds=c(2,2,2,8))
get_simu(k=4,distName=rep("normal",4),n=c(30,30,30,60),m=c(0,0,0,0),sds=c(2,2,2,1))
get_simu(k=4,distName=rep("normal",4),n=c(30,30,30,60),m=c(0,0,0,0),sds=c(2,2,2,2))
get_simu(k=4,distName=rep("normal",4),n=c(30,30,30,60),m=c(0,0,0,0),sds=c(2,2,2,4))
get_simu(k=4,distName=rep("normal",4),n=c(30,30,30,60),m=c(0,0,0,0),sds=c(2,2,2,8))
get_simu(k=4,distName=rep("normal",4),n=c(40,40,40,20),m=c(0,0,0,0),sds=c(2,2,2,1))
get_simu(k=4,distName=rep("normal",4),n=c(40,40,40,20),m=c(0,0,0,0),sds=c(2,2,2,2))
get_simu(k=4,distName=rep("normal",4),n=c(40,40,40,20),m=c(0,0,0,0),sds=c(2,2,2,4))
get_simu(k=4,distName=rep("normal",4),n=c(40,40,40,20),m=c(0,0,0,0),sds=c(2,2,2,8))
get_simu(k=4,distName=rep("normal",4),n=c(40,40,40,40),m=c(0,0,0,0),sds=c(2,2,2,1))
get_simu(k=4,distName=rep("normal",4),n=c(40,40,40,40),m=c(0,0,0,0),sds=c(2,2,2,2))
get_simu(k=4,distName=rep("normal",4),n=c(40,40,40,40),m=c(0,0,0,0),sds=c(2,2,2,4))
get_simu(k=4,distName=rep("normal",4),n=c(40,40,40,40),m=c(0,0,0,0),sds=c(2,2,2,8))
get_simu(k=4,distName=rep("normal",4),n=c(40,40,40,60),m=c(0,0,0,0),sds=c(2,2,2,1))
get_simu(k=4,distName=rep("normal",4),n=c(40,40,40,60),m=c(0,0,0,0),sds=c(2,2,2,2))
get_simu(k=4,distName=rep("normal",4),n=c(40,40,40,60),m=c(0,0,0,0),sds=c(2,2,2,4))
get_simu(k=4,distName=rep("normal",4),n=c(40,40,40,60),m=c(0,0,0,0),sds=c(2,2,2,8))
get_simu(k=4,distName=rep("normal",4),n=c(40,40,40,80),m=c(0,0,0,0),sds=c(2,2,2,1))
get_simu(k=4,distName=rep("normal",4),n=c(40,40,40,80),m=c(0,0,0,0),sds=c(2,2,2,2))
get_simu(k=4,distName=rep("normal",4),n=c(40,40,40,80),m=c(0,0,0,0),sds=c(2,2,2,4))
get_simu(k=4,distName=rep("normal",4),n=c(40,40,40,80),m=c(0,0,0,0),sds=c(2,2,2,8))
get_simu(k=4,distName=rep("normal",4),n=c(50,50,50,25),m=c(0,0,0,0),sds=c(2,2,2,1))
get_simu(k=4,distName=rep("normal",4),n=c(50,50,50,25),m=c(0,0,0,0),sds=c(2,2,2,2))
get_simu(k=4,distName=rep("normal",4),n=c(50,50,50,25),m=c(0,0,0,0),sds=c(2,2,2,4))
get_simu(k=4,distName=rep("normal",4),n=c(50,50,50,25),m=c(0,0,0,0),sds=c(2,2,2,8))
get_simu(k=4,distName=rep("normal",4),n=c(50,50,50,50),m=c(0,0,0,0),sds=c(2,2,2,1))
get_simu(k=4,distName=rep("normal",4),n=c(50,50,50,50),m=c(0,0,0,0),sds=c(2,2,2,2))
get_simu(k=4,distName=rep("normal",4),n=c(50,50,50,50),m=c(0,0,0,0),sds=c(2,2,2,4))
get_simu(k=4,distName=rep("normal",4),n=c(50,50,50,50),m=c(0,0,0,0),sds=c(2,2,2,8))
get_simu(k=4,distName=rep("normal",4),n=c(50,50,50,75),m=c(0,0,0,0),sds=c(2,2,2,1))
get_simu(k=4,distName=rep("normal",4),n=c(50,50,50,75),m=c(0,0,0,0),sds=c(2,2,2,2))
get_simu(k=4,distName=rep("normal",4),n=c(50,50,50,75),m=c(0,0,0,0),sds=c(2,2,2,4))
get_simu(k=4,distName=rep("normal",4),n=c(50,50,50,75),m=c(0,0,0,0),sds=c(2,2,2,8))
get_simu(k=4,distName=rep("normal",4),n=c(50,50,50,100),m=c(0,0,0,0),sds=c(2,2,2,1))
get_simu(k=4,distName=rep("normal",4),n=c(50,50,50,100),m=c(0,0,0,0),sds=c(2,2,2,2))
get_simu(k=4,distName=rep("normal",4),n=c(50,50,50,100),m=c(0,0,0,0),sds=c(2,2,2,4))
get_simu(k=4,distName=rep("normal",4),n=c(50,50,50,100),m=c(0,0,0,0),sds=c(2,2,2,8))
get_simu(k=4,distName=rep("normal",4),n=c(100,100,100,50),m=c(0,0,0,0),sds=c(2,2,2,1))
get_simu(k=4,distName=rep("normal",4),n=c(100,100,100,50),m=c(0,0,0,0),sds=c(2,2,2,2))
get_simu(k=4,distName=rep("normal",4),n=c(100,100,100,50),m=c(0,0,0,0),sds=c(2,2,2,4))
get_simu(k=4,distName=rep("normal",4),n=c(100,100,100,50),m=c(0,0,0,0),sds=c(2,2,2,8))
get_simu(k=4,distName=rep("normal",4),n=c(100,100,100,100),m=c(0,0,0,0),sds=c(2,2,2,1))
get_simu(k=4,distName=rep("normal",4),n=c(100,100,100,100),m=c(0,0,0,0),sds=c(2,2,2,2))
get_simu(k=4,distName=rep("normal",4),n=c(100,100,100,100),m=c(0,0,0,0),sds=c(2,2,2,4))
get_simu(k=4,distName=rep("normal",4),n=c(100,100,100,100),m=c(0,0,0,0),sds=c(2,2,2,8))
get_simu(k=4,distName=rep("normal",4),n=c(100,100,100,150),m=c(0,0,0,0),sds=c(2,2,2,1))
get_simu(k=4,distName=rep("normal",4),n=c(100,100,100,150),m=c(0,0,0,0),sds=c(2,2,2,2))
get_simu(k=4,distName=rep("normal",4),n=c(100,100,100,150),m=c(0,0,0,0),sds=c(2,2,2,4))
get_simu(k=4,distName=rep("normal",4),n=c(100,100,100,150),m=c(0,0,0,0),sds=c(2,2,2,8))
get_simu(k=4,distName=rep("normal",4),n=c(100,100,100,200),m=c(0,0,0,0),sds=c(2,2,2,1))
get_simu(k=4,distName=rep("normal",4),n=c(100,100,100,200),m=c(0,0,0,0),sds=c(2,2,2,2))
get_simu(k=4,distName=rep("normal",4),n=c(100,100,100,200),m=c(0,0,0,0),sds=c(2,2,2,4))
get_simu(k=4,distName=rep("normal",4),n=c(100,100,100,200),m=c(0,0,0,0),sds=c(2,2,2,8))


get_simu(k=5,distName=rep("normal",5),n=c(20,20,20,20,10),m=c(0,0,0,0,0),sds=c(2,2,2,2,1))
get_simu(k=5,distName=rep("normal",5),n=c(20,20,20,20,10),m=c(0,0,0,0,0),sds=c(2,2,2,2,2))
get_simu(k=5,distName=rep("normal",5),n=c(20,20,20,20,10),m=c(0,0,0,0,0),sds=c(2,2,2,2,4))
get_simu(k=5,distName=rep("normal",5),n=c(20,20,20,20,10),m=c(0,0,0,0,0),sds=c(2,2,2,2,8))
get_simu(k=5,distName=rep("normal",5),n=c(20,20,20,20,20),m=c(0,0,0,0,0),sds=c(2,2,2,2,1))
get_simu(k=5,distName=rep("normal",5),n=c(20,20,20,20,20),m=c(0,0,0,0,0),sds=c(2,2,2,2,2))
get_simu(k=5,distName=rep("normal",5),n=c(20,20,20,20,20),m=c(0,0,0,0,0),sds=c(2,2,2,2,4))
get_simu(k=5,distName=rep("normal",5),n=c(20,20,20,20,20),m=c(0,0,0,0,0),sds=c(2,2,2,2,8))
get_simu(k=5,distName=rep("normal",5),n=c(20,20,20,20,30),m=c(0,0,0,0,0),sds=c(2,2,2,2,1))
get_simu(k=5,distName=rep("normal",5),n=c(20,20,20,20,30),m=c(0,0,0,0,0),sds=c(2,2,2,2,2))
get_simu(k=5,distName=rep("normal",5),n=c(20,20,20,20,30),m=c(0,0,0,0,0),sds=c(2,2,2,2,4))
get_simu(k=5,distName=rep("normal",5),n=c(20,20,20,20,30),m=c(0,0,0,0,0),sds=c(2,2,2,2,8))
get_simu(k=5,distName=rep("normal",5),n=c(20,20,20,20,40),m=c(0,0,0,0,0),sds=c(2,2,2,2,1))
get_simu(k=5,distName=rep("normal",5),n=c(20,20,20,20,40),m=c(0,0,0,0,0),sds=c(2,2,2,2,2))
get_simu(k=5,distName=rep("normal",5),n=c(20,20,20,20,40),m=c(0,0,0,0,0),sds=c(2,2,2,2,4))
get_simu(k=5,distName=rep("normal",5),n=c(20,20,20,20,40),m=c(0,0,0,0,0),sds=c(2,2,2,2,8))
get_simu(k=5,distName=rep("normal",5),n=c(30,30,30,30,15),m=c(0,0,0,0,0),sds=c(2,2,2,2,1))
get_simu(k=5,distName=rep("normal",5),n=c(30,30,30,30,15),m=c(0,0,0,0,0),sds=c(2,2,2,2,2))
get_simu(k=5,distName=rep("normal",5),n=c(30,30,30,30,15),m=c(0,0,0,0,0),sds=c(2,2,2,2,4))
get_simu(k=5,distName=rep("normal",5),n=c(30,30,30,30,15),m=c(0,0,0,0,0),sds=c(2,2,2,2,8))
get_simu(k=5,distName=rep("normal",5),n=c(30,30,30,30,30),m=c(0,0,0,0,0),sds=c(2,2,2,2,1))
get_simu(k=5,distName=rep("normal",5),n=c(30,30,30,30,30),m=c(0,0,0,0,0),sds=c(2,2,2,2,2))
get_simu(k=5,distName=rep("normal",5),n=c(30,30,30,30,30),m=c(0,0,0,0,0),sds=c(2,2,2,2,4))
get_simu(k=5,distName=rep("normal",5),n=c(30,30,30,30,30),m=c(0,0,0,0,0),sds=c(2,2,2,2,8))
get_simu(k=5,distName=rep("normal",5),n=c(30,30,30,30,45),m=c(0,0,0,0,0),sds=c(2,2,2,2,1))
get_simu(k=5,distName=rep("normal",5),n=c(30,30,30,30,45),m=c(0,0,0,0,0),sds=c(2,2,2,2,2))
get_simu(k=5,distName=rep("normal",5),n=c(30,30,30,30,45),m=c(0,0,0,0,0),sds=c(2,2,2,2,4))
get_simu(k=5,distName=rep("normal",5),n=c(30,30,30,30,45),m=c(0,0,0,0,0),sds=c(2,2,2,2,8))
get_simu(k=5,distName=rep("normal",5),n=c(30,30,30,30,60),m=c(0,0,0,0,0),sds=c(2,2,2,2,1))
get_simu(k=5,distName=rep("normal",5),n=c(30,30,30,30,60),m=c(0,0,0,0,0),sds=c(2,2,2,2,2))
get_simu(k=5,distName=rep("normal",5),n=c(30,30,30,30,60),m=c(0,0,0,0,0),sds=c(2,2,2,2,4))
get_simu(k=5,distName=rep("normal",5),n=c(30,30,30,30,60),m=c(0,0,0,0,0),sds=c(2,2,2,2,8))
get_simu(k=5,distName=rep("normal",5),n=c(40,40,40,40,20),m=c(0,0,0,0,0),sds=c(2,2,2,2,1))
get_simu(k=5,distName=rep("normal",5),n=c(40,40,40,40,20),m=c(0,0,0,0,0),sds=c(2,2,2,2,2))
get_simu(k=5,distName=rep("normal",5),n=c(40,40,40,40,20),m=c(0,0,0,0,0),sds=c(2,2,2,2,4))
get_simu(k=5,distName=rep("normal",5),n=c(40,40,40,40,20),m=c(0,0,0,0,0),sds=c(2,2,2,2,8))
get_simu(k=5,distName=rep("normal",5),n=c(40,40,40,40,40),m=c(0,0,0,0,0),sds=c(2,2,2,2,1))
get_simu(k=5,distName=rep("normal",5),n=c(40,40,40,40,40),m=c(0,0,0,0,0),sds=c(2,2,2,2,2))
get_simu(k=5,distName=rep("normal",5),n=c(40,40,40,40,40),m=c(0,0,0,0,0),sds=c(2,2,2,2,4))
get_simu(k=5,distName=rep("normal",5),n=c(40,40,40,40,40),m=c(0,0,0,0,0),sds=c(2,2,2,2,8))
get_simu(k=5,distName=rep("normal",5),n=c(40,40,40,40,60),m=c(0,0,0,0,0),sds=c(2,2,2,2,1))
get_simu(k=5,distName=rep("normal",5),n=c(40,40,40,40,60),m=c(0,0,0,0,0),sds=c(2,2,2,2,2))
get_simu(k=5,distName=rep("normal",5),n=c(40,40,40,40,60),m=c(0,0,0,0,0),sds=c(2,2,2,2,4))
get_simu(k=5,distName=rep("normal",5),n=c(40,40,40,40,60),m=c(0,0,0,0,0),sds=c(2,2,2,2,8))
get_simu(k=5,distName=rep("normal",5),n=c(40,40,40,40,80),m=c(0,0,0,0,0),sds=c(2,2,2,2,1))
get_simu(k=5,distName=rep("normal",5),n=c(40,40,40,40,80),m=c(0,0,0,0,0),sds=c(2,2,2,2,2))
get_simu(k=5,distName=rep("normal",5),n=c(40,40,40,40,80),m=c(0,0,0,0,0),sds=c(2,2,2,2,4))
get_simu(k=5,distName=rep("normal",5),n=c(40,40,40,40,80),m=c(0,0,0,0,0),sds=c(2,2,2,2,8))
get_simu(k=5,distName=rep("normal",5),n=c(50,50,50,50,25),m=c(0,0,0,0,0),sds=c(2,2,2,2,1))
get_simu(k=5,distName=rep("normal",5),n=c(50,50,50,50,25),m=c(0,0,0,0,0),sds=c(2,2,2,2,2))
get_simu(k=5,distName=rep("normal",5),n=c(50,50,50,50,25),m=c(0,0,0,0,0),sds=c(2,2,2,2,4))
get_simu(k=5,distName=rep("normal",5),n=c(50,50,50,50,25),m=c(0,0,0,0,0),sds=c(2,2,2,2,8))
get_simu(k=5,distName=rep("normal",5),n=c(50,50,50,50,50),m=c(0,0,0,0,0),sds=c(2,2,2,2,1))
get_simu(k=5,distName=rep("normal",5),n=c(50,50,50,50,50),m=c(0,0,0,0,0),sds=c(2,2,2,2,2))
get_simu(k=5,distName=rep("normal",5),n=c(50,50,50,50,50),m=c(0,0,0,0,0),sds=c(2,2,2,2,4))
get_simu(k=5,distName=rep("normal",5),n=c(50,50,50,50,50),m=c(0,0,0,0,0),sds=c(2,2,2,2,8))
get_simu(k=5,distName=rep("normal",5),n=c(50,50,50,50,75),m=c(0,0,0,0,0),sds=c(2,2,2,2,1))
get_simu(k=5,distName=rep("normal",5),n=c(50,50,50,50,75),m=c(0,0,0,0,0),sds=c(2,2,2,2,2))
get_simu(k=5,distName=rep("normal",5),n=c(50,50,50,50,75),m=c(0,0,0,0,0),sds=c(2,2,2,2,4))
get_simu(k=5,distName=rep("normal",5),n=c(50,50,50,50,75),m=c(0,0,0,0,0),sds=c(2,2,2,2,8))
get_simu(k=5,distName=rep("normal",5),n=c(50,50,50,50,100),m=c(0,0,0,0,0),sds=c(2,2,2,2,1))
get_simu(k=5,distName=rep("normal",5),n=c(50,50,50,50,100),m=c(0,0,0,0,0),sds=c(2,2,2,2,2))
get_simu(k=5,distName=rep("normal",5),n=c(50,50,50,50,100),m=c(0,0,0,0,0),sds=c(2,2,2,2,4))
get_simu(k=5,distName=rep("normal",5),n=c(50,50,50,50,100),m=c(0,0,0,0,0),sds=c(2,2,2,2,8))
get_simu(k=5,distName=rep("normal",5),n=c(100,100,100,100,50),m=c(0,0,0,0,0),sds=c(2,2,2,2,1))
get_simu(k=5,distName=rep("normal",5),n=c(100,100,100,100,50),m=c(0,0,0,0,0),sds=c(2,2,2,2,2))
get_simu(k=5,distName=rep("normal",5),n=c(100,100,100,100,50),m=c(0,0,0,0,0),sds=c(2,2,2,2,4))
get_simu(k=5,distName=rep("normal",5),n=c(100,100,100,100,50),m=c(0,0,0,0,0),sds=c(2,2,2,2,8))
get_simu(k=5,distName=rep("normal",5),n=c(100,100,100,100,100),m=c(0,0,0,0,0),sds=c(2,2,2,2,1))
get_simu(k=5,distName=rep("normal",5),n=c(100,100,100,100,100),m=c(0,0,0,0,0),sds=c(2,2,2,2,2))
get_simu(k=5,distName=rep("normal",5),n=c(100,100,100,100,100),m=c(0,0,0,0,0),sds=c(2,2,2,2,4))
get_simu(k=5,distName=rep("normal",5),n=c(100,100,100,100,100),m=c(0,0,0,0,0),sds=c(2,2,2,2,8))
get_simu(k=5,distName=rep("normal",5),n=c(100,100,100,100,150),m=c(0,0,0,0,0),sds=c(2,2,2,2,1))
get_simu(k=5,distName=rep("normal",5),n=c(100,100,100,100,150),m=c(0,0,0,0,0),sds=c(2,2,2,2,2))
get_simu(k=5,distName=rep("normal",5),n=c(100,100,100,100,150),m=c(0,0,0,0,0),sds=c(2,2,2,2,4))
get_simu(k=5,distName=rep("normal",5),n=c(100,100,100,100,150),m=c(0,0,0,0,0),sds=c(2,2,2,2,8))
get_simu(k=5,distName=rep("normal",5),n=c(100,100,100,100,200),m=c(0,0,0,0,0),sds=c(2,2,2,2,1))
get_simu(k=5,distName=rep("normal",5),n=c(100,100,100,100,200),m=c(0,0,0,0,0),sds=c(2,2,2,2,2))
get_simu(k=5,distName=rep("normal",5),n=c(100,100,100,100,200),m=c(0,0,0,0,0),sds=c(2,2,2,2,4))
get_simu(k=5,distName=rep("normal",5),n=c(100,100,100,100,200),m=c(0,0,0,0,0),sds=c(2,2,2,2,8))

### double exponential distributions

get_simu(k=2,distName=rep("doublex",2),n=c(20,10),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(20,10),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(20,10),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(20,10),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(20,20),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(20,20),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(20,20),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(20,20),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(20,30),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(20,30),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(20,30),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(20,30),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(20,40),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(20,40),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(20,40),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(20,40),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(30,15),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(30,15),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(30,15),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(30,15),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(30,30),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(30,30),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(30,30),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(30,30),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(30,45),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(30,45),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(30,45),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(30,45),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(30,60),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(30,60),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(30,60),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(30,60),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(40,20),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(40,20),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(40,20),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(40,20),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(40,40),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(40,40),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(40,40),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(40,40),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(40,60),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(40,60),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(40,60),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(40,60),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(40,80),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(40,80),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(40,80),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(40,80),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(50,25),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(50,25),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(50,25),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(50,25),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(50,50),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(50,50),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(50,50),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(50,50),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(50,75),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(50,75),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(50,75),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(50,75),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(50,100),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(50,100),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(50,100),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(50,100),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(100,50),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(100,50),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(100,50),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(100,50),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(100,100),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(100,100),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(100,100),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(100,100),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(100,150),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(100,150),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(100,150),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(100,150),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(100,200),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(100,200),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(100,200),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(100,200),m=c(0,0),lambda=c(2,8))

get_simu(k=3,distName=rep("doublex",3),n=c(20,20,10),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(20,20,10),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(20,20,10),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(20,20,10),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(20,20,20),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(20,20,20),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(20,20,20),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(20,20,20),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(20,20,30),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(20,20,30),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(20,20,30),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(20,20,30),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(20,20,40),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(20,20,40),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(20,20,40),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(20,20,40),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(30,30,15),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(30,30,15),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(30,30,15),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(30,30,15),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(30,30,30),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(30,30,30),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(30,30,30),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(30,30,30),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(30,30,45),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(30,30,45),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(30,30,45),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(30,30,45),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(30,30,60),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(30,30,60),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(30,30,60),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(30,30,60),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(40,40,20),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(40,40,20),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(40,40,20),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(40,40,20),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(40,40,40),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(40,40,40),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(40,40,40),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(40,40,40),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(40,40,60),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(40,40,60),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(40,40,60),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(40,40,60),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(40,40,80),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(40,40,80),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(40,40,80),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(40,40,80),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(50,50,25),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(50,50,25),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(50,50,25),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(50,50,25),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(50,50,50),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(50,50,50),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(50,50,50),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(50,50,50),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(50,50,75),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(50,50,75),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(50,50,75),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(50,50,75),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(50,50,100),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(50,50,100),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(50,50,100),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(50,50,100),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(100,100,50),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(100,100,50),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(100,100,50),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(100,100,50),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(100,100,100),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(100,100,100),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(100,100,100),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(100,100,100),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(100,100,150),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(100,100,150),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(100,100,150),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(100,100,150),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(100,100,200),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(100,100,200),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(100,100,200),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(100,100,200),m=c(0,0,0),lambda=c(2,2,8))

get_simu(k=4,distName=rep("doublex",4),n=c(20,20,20,10),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex",4),n=c(20,20,20,10),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex",4),n=c(20,20,20,10),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex",4),n=c(20,20,20,10),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex",4),n=c(20,20,20,20),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex",4),n=c(20,20,20,20),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex",4),n=c(20,20,20,20),m=c(0,0,0,0),lambda=c(2,2,2,4))sd2
get_simu(k=4,distName=rep("doublex",4),n=c(20,20,20,20),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex",4),n=c(20,20,20,30),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex",4),n=c(20,20,20,30),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex",4),n=c(20,20,20,30),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex",4),n=c(20,20,20,30),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex",4),n=c(20,20,20,40),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex",4),n=c(20,20,20,40),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex",4),n=c(20,20,20,40),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex",4),n=c(20,20,20,40),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex",4),n=c(30,30,30,15),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex",4),n=c(30,30,30,15),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex",4),n=c(30,30,30,15),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex",4),n=c(30,30,30,15),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex",4),n=c(30,30,30,30),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex",4),n=c(30,30,30,30),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex",4),n=c(30,30,30,30),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex",4),n=c(30,30,30,30),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex",4),n=c(30,30,30,45),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex",4),n=c(30,30,30,45),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex",4),n=c(30,30,30,45),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex",4),n=c(30,30,30,45),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex",4),n=c(30,30,30,60),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex",4),n=c(30,30,30,60),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex",4),n=c(30,30,30,60),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex",4),n=c(30,30,30,60),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex",4),n=c(40,40,40,20),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex",4),n=c(40,40,40,20),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex",4),n=c(40,40,40,20),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex",4),n=c(40,40,40,20),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex",4),n=c(40,40,40,40),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex",4),n=c(40,40,40,40),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex",4),n=c(40,40,40,40),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex",4),n=c(40,40,40,40),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex",4),n=c(40,40,40,60),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex",4),n=c(40,40,40,60),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex",4),n=c(40,40,40,60),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex",4),n=c(40,40,40,60),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex",4),n=c(40,40,40,80),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex",4),n=c(40,40,40,80),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex",4),n=c(40,40,40,80),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex",4),n=c(40,40,40,80),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex",4),n=c(50,50,50,25),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex",4),n=c(50,50,50,25),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex",4),n=c(50,50,50,25),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex",4),n=c(50,50,50,25),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex",4),n=c(50,50,50,50),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex",4),n=c(50,50,50,50),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex",4),n=c(50,50,50,50),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex",4),n=c(50,50,50,50),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex",4),n=c(50,50,50,75),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex",4),n=c(50,50,50,75),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex",4),n=c(50,50,50,75),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex",4),n=c(50,50,50,75),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex",4),n=c(50,50,50,100),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex",4),n=c(50,50,50,100),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex",4),n=c(50,50,50,100),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex",4),n=c(50,50,50,100),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex",4),n=c(100,100,100,50),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex",4),n=c(100,100,100,50),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex",4),n=c(100,100,100,50),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex",4),n=c(100,100,100,50),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex",4),n=c(100,100,100,100),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex",4),n=c(100,100,100,100),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex",4),n=c(100,100,100,100),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex",4),n=c(100,100,100,100),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex",4),n=c(100,100,100,150),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex",4),n=c(100,100,100,150),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex",4),n=c(100,100,100,150),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex",4),n=c(100,100,100,150),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex",4),n=c(100,100,100,200),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex",4),n=c(100,100,100,200),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex",4),n=c(100,100,100,200),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex",4),n=c(100,100,100,200),m=c(0,0,0,0),lambda=c(2,2,2,8))

get_simu(k=5,distName=rep("doublex",5),n=c(20,20,20,20,10),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex",5),n=c(20,20,20,20,10),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex",5),n=c(20,20,20,20,10),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex",5),n=c(20,20,20,20,10),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex",5),n=c(20,20,20,20,20),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex",5),n=c(20,20,20,20,20),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex",5),n=c(20,20,20,20,20),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex",5),n=c(20,20,20,20,20),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex",5),n=c(20,20,20,20,30),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex",5),n=c(20,20,20,20,30),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex",5),n=c(20,20,20,20,30),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex",5),n=c(20,20,20,20,30),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex",5),n=c(20,20,20,20,40),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex",5),n=c(20,20,20,20,40),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex",5),n=c(20,20,20,20,40),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex",5),n=c(20,20,20,20,40),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex",5),n=c(30,30,30,30,15),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex",5),n=c(30,30,30,30,15),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex",5),n=c(30,30,30,30,15),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex",5),n=c(30,30,30,30,15),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex",5),n=c(30,30,30,30,30),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex",5),n=c(30,30,30,30,30),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex",5),n=c(30,30,30,30,30),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex",5),n=c(30,30,30,30,30),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex",5),n=c(30,30,30,30,45),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex",5),n=c(30,30,30,30,45),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex",5),n=c(30,30,30,30,45),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex",5),n=c(30,30,30,30,45),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex",5),n=c(30,30,30,30,60),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex",5),n=c(30,30,30,30,60),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex",5),n=c(30,30,30,30,60),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex",5),n=c(30,30,30,30,60),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex",5),n=c(40,40,40,40,20),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex",5),n=c(40,40,40,40,20),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex",5),n=c(40,40,40,40,20),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex",5),n=c(40,40,40,40,20),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex",5),n=c(40,40,40,40,40),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex",5),n=c(40,40,40,40,40),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex",5),n=c(40,40,40,40,40),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex",5),n=c(40,40,40,40,40),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex",5),n=c(40,40,40,40,60),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex",5),n=c(40,40,40,40,60),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex",5),n=c(40,40,40,40,60),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex",5),n=c(40,40,40,40,60),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex",5),n=c(40,40,40,40,80),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex",5),n=c(40,40,40,40,80),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex",5),n=c(40,40,40,40,80),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex",5),n=c(40,40,40,40,80),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex",5),n=c(50,50,50,50,25),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex",5),n=c(50,50,50,50,25),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex",5),n=c(50,50,50,50,25),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex",5),n=c(50,50,50,50,25),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex",5),n=c(50,50,50,50,50),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex",5),n=c(50,50,50,50,50),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex",5),n=c(50,50,50,50,50),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex",5),n=c(50,50,50,50,50),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex",5),n=c(50,50,50,50,75),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex",5),n=c(50,50,50,50,75),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex",5),n=c(50,50,50,50,75),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex",5),n=c(50,50,50,50,75),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex",5),n=c(50,50,50,50,100),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex",5),n=c(50,50,50,50,100),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex",5),n=c(50,50,50,50,100),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex",5),n=c(50,50,50,50,100),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex",5),n=c(100,100,100,100,50),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex",5),n=c(100,100,100,100,50),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex",5),n=c(100,100,100,100,50),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex",5),n=c(100,100,100,100,50),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex",5),n=c(100,100,100,100,100),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex",5),n=c(100,100,100,100,100),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex",5),n=c(100,100,100,100,100),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex",5),n=c(100,100,100,100,100),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex",5),n=c(100,100,100,100,150),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex",5),n=c(100,100,100,100,150),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex",5),n=c(100,100,100,100,150),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex",5),n=c(100,100,100,100,150),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex",5),n=c(100,100,100,100,200),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex",5),n=c(100,100,100,100,200),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex",5),n=c(100,100,100,100,200),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex",5),n=c(100,100,100,100,200),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))


### double exponential distribution (using default scale parameter of the distribution)

get_simu(k=2,distName=rep("doublex_SP",2),n=c(20,10),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(20,10),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(20,10),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(20,10),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(20,20),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(20,20),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(20,20),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(20,20),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(20,30),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(20,30),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(20,30),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(20,30),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(20,40),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(20,40),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(20,40),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(20,40),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(30,15),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(30,15),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(30,15),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(30,15),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(30,30),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(30,30),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(30,30),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(30,30),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(30,45),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(30,45),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(30,45),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(30,45),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(30,60),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(30,60),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(30,60),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(30,60),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(40,20),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(40,20),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(40,20),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(40,20),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(40,40),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(40,40),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(40,40),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(40,40),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(40,60),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(40,60),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(40,60),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(40,60),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(40,80),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(40,80),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(40,80),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(40,80),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(50,25),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(50,25),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(50,25),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(50,25),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(50,50),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(50,50),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(50,50),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(50,50),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(50,75),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(50,75),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(50,75),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(50,75),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(50,100),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(50,100),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(50,100),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(50,100),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(100,50),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(100,50),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(100,50),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(100,50),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(100,100),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(100,100),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(100,100),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(100,100),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(100,150),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(100,150),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(100,150),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(100,150),m=c(0,0),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(100,200),m=c(0,0),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(100,200),m=c(0,0),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(100,200),m=c(0,0),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(100,200),m=c(0,0),lambda=c(2,8))

get_simu(k=3,distName=rep("doublex_SP",3),n=c(20,20,10),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(20,20,10),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(20,20,10),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(20,20,10),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(20,20,20),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(20,20,20),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(20,20,20),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(20,20,20),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(20,20,30),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(20,20,30),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(20,20,30),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(20,20,30),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(20,20,40),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(20,20,40),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(20,20,40),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(20,20,40),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(30,30,15),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(30,30,15),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(30,30,15),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(30,30,15),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(30,30,30),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(30,30,30),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(30,30,30),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(30,30,30),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(30,30,45),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(30,30,45),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(30,30,45),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(30,30,45),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(30,30,60),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(30,30,60),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(30,30,60),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(30,30,60),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(40,40,20),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(40,40,20),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(40,40,20),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(40,40,20),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(40,40,40),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(40,40,40),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(40,40,40),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(40,40,40),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(40,40,60),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(40,40,60),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(40,40,60),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(40,40,60),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(40,40,80),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(40,40,80),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(40,40,80),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(40,40,80),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(50,50,25),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(50,50,25),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(50,50,25),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(50,50,25),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(50,50,50),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(50,50,50),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(50,50,50),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(50,50,50),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(50,50,75),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(50,50,75),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(50,50,75),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(50,50,75),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(50,50,100),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(50,50,100),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(50,50,100),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(50,50,100),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(100,100,50),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(100,100,50),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(100,100,50),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(100,100,50),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(100,100,100),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(100,100,100),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(100,100,100),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(100,100,100),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(100,100,150),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(100,100,150),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(100,100,150),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(100,100,150),m=c(0,0,0),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(100,100,200),m=c(0,0,0),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(100,100,200),m=c(0,0,0),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(100,100,200),m=c(0,0,0),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(100,100,200),m=c(0,0,0),lambda=c(2,2,8))

get_simu(k=4,distName=rep("doublex_SP",4),n=c(20,20,20,10),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(20,20,20,10),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(20,20,20,10),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(20,20,20,10),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(20,20,20,20),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(20,20,20,20),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(20,20,20,20),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(20,20,20,20),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(20,20,20,30),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(20,20,20,30),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(20,20,20,30),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(20,20,20,30),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(20,20,20,40),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(20,20,20,40),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(20,20,20,40),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(20,20,20,40),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(30,30,30,15),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(30,30,30,15),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(30,30,30,15),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(30,30,30,15),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(30,30,30,30),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(30,30,30,30),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(30,30,30,30),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(30,30,30,30),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(30,30,30,45),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(30,30,30,45),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(30,30,30,45),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(30,30,30,45),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(30,30,30,60),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(30,30,30,60),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(30,30,30,60),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(30,30,30,60),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(40,40,40,20),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(40,40,40,20),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(40,40,40,20),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(40,40,40,20),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(40,40,40,40),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(40,40,40,40),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(40,40,40,40),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(40,40,40,40),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(40,40,40,60),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(40,40,40,60),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(40,40,40,60),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(40,40,40,60),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(40,40,40,80),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(40,40,40,80),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(40,40,40,80),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(40,40,40,80),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(50,50,50,25),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(50,50,50,25),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(50,50,50,25),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(50,50,50,25),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(50,50,50,50),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(50,50,50,50),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(50,50,50,50),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(50,50,50,50),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(50,50,50,75),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(50,50,50,75),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(50,50,50,75),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(50,50,50,75),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(50,50,50,100),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(50,50,50,100),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(50,50,50,100),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(50,50,50,100),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(100,100,100,50),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(100,100,100,50),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(100,100,100,50),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(100,100,100,50),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(100,100,100,100),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(100,100,100,100),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(100,100,100,100),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(100,100,100,100),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(100,100,100,150),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(100,100,100,150),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(100,100,100,150),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(100,100,100,150),m=c(0,0,0,0),lambda=c(2,2,2,8))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(100,100,100,200),m=c(0,0,0,0),lambda=c(2,2,2,1))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(100,100,100,200),m=c(0,0,0,0),lambda=c(2,2,2,2))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(100,100,100,200),m=c(0,0,0,0),lambda=c(2,2,2,4))
get_simu(k=4,distName=rep("doublex_SP",4),n=c(100,100,100,200),m=c(0,0,0,0),lambda=c(2,2,2,8))

get_simu(k=5,distName=rep("doublex_SP",5),n=c(20,20,20,20,10),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(20,20,20,20,10),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(20,20,20,20,10),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(20,20,20,20,10),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(20,20,20,20,20),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(20,20,20,20,20),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(20,20,20,20,20),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(20,20,20,20,20),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(20,20,20,20,30),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(20,20,20,20,30),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(20,20,20,20,30),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(20,20,20,20,30),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(20,20,20,20,40),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(20,20,20,20,40),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(20,20,20,20,40),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(20,20,20,20,40),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(30,30,30,30,15),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(30,30,30,30,15),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(30,30,30,30,15),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(30,30,30,30,15),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(30,30,30,30,30),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(30,30,30,30,30),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(30,30,30,30,30),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(30,30,30,30,30),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(30,30,30,30,45),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(30,30,30,30,45),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(30,30,30,30,45),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(30,30,30,30,45),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(30,30,30,30,60),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(30,30,30,30,60),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(30,30,30,30,60),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(30,30,30,30,60),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(40,40,40,40,20),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(40,40,40,40,20),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(40,40,40,40,20),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(40,40,40,40,20),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(40,40,40,40,40),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(40,40,40,40,40),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(40,40,40,40,40),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(40,40,40,40,40),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(40,40,40,40,60),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(40,40,40,40,60),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(40,40,40,40,60),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(40,40,40,40,60),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(40,40,40,40,80),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(40,40,40,40,80),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(40,40,40,40,80),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(40,40,40,40,80),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(50,50,50,50,25),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(50,50,50,50,25),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(50,50,50,50,25),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(50,50,50,50,25),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(50,50,50,50,50),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(50,50,50,50,50),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(50,50,50,50,50),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(50,50,50,50,50),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(50,50,50,50,75),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(50,50,50,50,75),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(50,50,50,50,75),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(50,50,50,50,75),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(50,50,50,50,100),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(50,50,50,50,100),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(50,50,50,50,100),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(50,50,50,50,100),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(100,100,100,100,50),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(100,100,100,100,50),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(100,100,100,100,50),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(100,100,100,100,50),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(100,100,100,100,100),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(100,100,100,100,100),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(100,100,100,100,100),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(100,100,100,100,100),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(100,100,100,100,150),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(100,100,100,100,150),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(100,100,100,100,150),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(100,100,100,100,150),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(100,100,100,100,200),m=c(0,0,0,0,0),lambda=c(2,2,2,2,1))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(100,100,100,100,200),m=c(0,0,0,0,0),lambda=c(2,2,2,2,2))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(100,100,100,100,200),m=c(0,0,0,0,0),lambda=c(2,2,2,2,4))
get_simu(k=5,distName=rep("doublex_SP",5),n=c(100,100,100,100,200),m=c(0,0,0,0,0),lambda=c(2,2,2,2,8))

# mixed normal distributions

get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(20,10),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,2.53),sd2=c(1.265,.6325))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(20,10),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,5.06),sd2=c(1.265,1.265))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(20,10),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,10.119),sd2=c(1.265,2.53))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(20,10),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,20.239),sd2=c(1.265,5.06))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(20,20),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,2.53),sd2=c(1.265,.6325))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(20,20),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,5.06),sd2=c(1.265,1.265))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(20,20),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,10.119),sd2=c(1.265,2.53))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(20,20),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,20.239),sd2=c(1.265,5.06))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(20,30),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,2.53),sd2=c(1.265,.6325))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(20,30),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,5.06),sd2=c(1.265,1.265))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(20,30),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,10.119),sd2=c(1.265,2.53))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(20,30),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,20.239),sd2=c(1.265,5.06))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(20,40),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,2.53),sd2=c(1.265,.6325))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(20,40),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,5.06),sd2=c(1.265,1.265))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(20,40),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,10.119),sd2=c(1.265,2.53))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(20,40),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,20.239),sd2=c(1.265,5.06))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(30,15),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,2.53),sd2=c(1.265,.6325))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(30,15),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,5.06),sd2=c(1.265,1.265))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(30,15),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,10.119),=c(1.265,2.53))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(30,15),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,20.239),sd2=c(1.265,5.06))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(30,30),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,2.53),sd2=c(1.265,.6325))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(30,30),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,5.06),sd2=c(1.265,1.265))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(30,30),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,10.119),sd2=c(1.265,2.53))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(30,30),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,20.239),sd2=c(1.265,5.06))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(30,45),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,2.53),sd2=c(1.265,.6325))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(30,45),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,5.06),sd2=c(1.265,1.265))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(30,45),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,10.119),sd2=c(1.265,2.53))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(30,45),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,20.239),sd2=c(1.265,5.06))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(30,60),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,2.53),sd2=c(1.265,.6325))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(30,60),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,5.06),sd2=c(1.265,1.265))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(30,60),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,10.119),sd2=c(1.265,2.53))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(30,60),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,20.239),sd2=c(1.265,5.06))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(40,20),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,2.53),sd2=c(1.265,.6325))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(40,20),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,5.06),sd2=c(1.265,1.265))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(40,20),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,10.119),sd2=c(1.265,2.53))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(40,20),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,20.239),sd2=c(1.265,5.06))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(40,40),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,2.53),sd2=c(1.265,.6325))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(40,40),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,5.06),sd2=c(1.265,1.265))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(40,40),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,10.119),sd2=c(1.265,2.53))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(40,40),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,20.239),sd2=c(1.265,5.06))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(40,60),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,2.53),sd2=c(1.265,.6325))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(40,60),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,5.06),sd2=c(1.265,1.265))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(40,60),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,10.119),sd2=c(1.265,2.53))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(40,60),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,20.239),sd2=c(1.265,5.06))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(40,80),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,2.53),sd2=c(1.265,.6325))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(40,80),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,5.06),sd2=c(1.265,1.265))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(40,80),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,10.119),sd2=c(1.265,2.53))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(40,80),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,20.239),sd2=c(1.265,5.06))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(50,25),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,2.53),sd2=c(1.265,.6325))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(50,25),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,5.06),sd2=c(1.265,1.265))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(50,25),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,10.119),sd2=c(1.265,2.53))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(50,25),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,20.239),sd2=c(1.265,5.06))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(50,50),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,2.53),sd2=c(1.265,.6325))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(50,50),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,5.06),sd2=c(1.265,1.265))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(50,50),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,10.119),sd2=c(1.265,2.53))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(50,50),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,20.239),sd2=c(1.265,5.06))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(50,75),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,2.53),sd2=c(1.265,.6325))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(50,75),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,5.06),sd2=c(1.265,1.265))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(50,75),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,10.119),sd2=c(1.265,2.53))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(50,75),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,20.239),sd2=c(1.265,5.06))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(50,100),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,2.53),sd2=c(1.265,.6325))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(50,100),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,5.06),sd2=c(1.265,1.265))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(50,100),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,10.119),sd2=c(1.265,2.53))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(50,100),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,20.239),sd2=c(1.265,5.06))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(100,50),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,2.53),sd2=c(1.265,.6325))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(100,50),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,5.06),sd2=c(1.265,1.265))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(100,50),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,10.119),sd2=c(1.265,2.53))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(100,50),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,20.239),sd2=c(1.265,5.06))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(100,100),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,2.53),sd2=c(1.265,.6325))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(100,100),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,5.06),sd2=c(1.265,1.265))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(100,100),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,10.119),sd2=c(1.265,2.53))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(100,100),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,20.239),sd2=c(1.265,5.06))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(100,150),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,2.53),sd2=c(1.265,.6325))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(100,150),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,5.06),sd2=c(1.265,1.265))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(100,150),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,10.119),sd2=c(1.265,2.53))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(100,150),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,20.239),sd2=c(1.265,5.06))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(100,200),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,2.53),sd2=c(1.265,.6325))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(100,200),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,5.06),sd2=c(1.265,1.265))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(100,200),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,10.119),sd2=c(1.265,2.53))
get_simu(nSims=1,k=2,distName=rep("mixed",2),n=c(100,200),m=c(0,0),p1=c(.1,.1),p2=c(.9,.9),sd1=c(5.06,20.239),sd2=c(1.265,5.06))

get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(20,20,10),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,2.53),sd2=c(1.265,1.265,.6325))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(20,20,10),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,5.06),sd2=c(1.265,1.265,1.265))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(20,20,10),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,10.119),sd2=c(1.265,1.265,2.53))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(20,20,10),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,20.239),sd2=c(1.265,1.265,5.06))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(20,20,20),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,2.53),sd2=c(1.265,1.265,.6325))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(20,20,20),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,5.06),sd2=c(1.265,1.265,1.265))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(20,20,20),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,10.119),sd2=c(1.265,1.265,2.53))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(20,20,20),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,20.239),sd2=c(1.265,1.265,5.06))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(20,20,30),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,2.53),sd2=c(1.265,1.265,.6325))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(20,20,30),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,5.06),sd2=c(1.265,1.265,1.265))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(20,20,30),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,10.119),sd2=c(1.265,1.265,2.53))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(20,20,30),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,20.239),sd2=c(1.265,1.265,5.06))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(20,20,40),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,2.53),sd2=c(1.265,1.265,.6325))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(20,20,40),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,5.06),sd2=c(1.265,1.265,1.265))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(20,20,40),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,10.119),sd2=c(1.265,1.265,2.53))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(20,20,40),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,20.239),sd2=c(1.265,1.265,5.06))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(30,30,15),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,2.53),sd2=c(1.265,1.265,.6325))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(30,30,15),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,5.06),sd2=c(1.265,1.265,1.265))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(30,30,15),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,10.119),sd2=c(1.265,1.265,2.53))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(30,30,15),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,20.239),sd2=c(1.265,1.265,5.06))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(30,30,30),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,2.53),sd2=c(1.265,1.265,.6325))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(30,30,30),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,5.06),sd2=c(1.265,1.265,1.265))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(30,30,30),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,10.119),sd2=c(1.265,1.265,2.53))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(30,30,30),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,20.239),sd2=c(1.265,1.265,5.06))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(30,30,45),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,2.53),sd2=c(1.265,1.265,.6325))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(30,30,45),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,5.06),sd2=c(1.265,1.265,1.265))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(30,30,45),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,10.119),sd2=c(1.265,1.265,2.53))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(30,30,45),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,20.239),sd2=c(1.265,1.265,5.06))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(30,30,60),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,2.53),sd2=c(1.265,1.265,.6325))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(30,30,60),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,5.06),sd2=c(1.265,1.265,1.265))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(30,30,60),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,10.119),sd2=c(1.265,1.265,2.53))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(30,30,60),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,20.239),sd2=c(1.265,1.265,5.06))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(40,40,20),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,2.53),sd2=c(1.265,1.265,.6325))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(40,40,20),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,5.06),sd2=c(1.265,1.265,1.265))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(40,40,20),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,10.119),sd2=c(1.265,1.265,2.53))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(40,40,20),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,20.239),sd2=c(1.265,1.265,5.06))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(40,40,40),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,2.53),sd2=c(1.265,1.265,.6325))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(40,40,40),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,5.06),sd2=c(1.265,1.265,1.265))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(40,40,40),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,10.119),sd2=c(1.265,1.265,2.53))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(40,40,40),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,20.239),sd2=c(1.265,1.265,5.06))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(40,40,60),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,2.53),sd2=c(1.265,1.265,.6325))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(40,40,60),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,5.06),sd2=c(1.265,1.265,1.265))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(40,40,60),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,10.119),sd2=c(1.265,1.265,2.53))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(40,40,60),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,20.239),sd2=c(1.265,1.265,5.06))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(40,40,80),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,2.53),sd2=c(1.265,1.265,.6325))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(40,40,80),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,5.06),sd2=c(1.265,1.265,1.265))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(40,40,80),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,10.119),sd2=c(1.265,1.265,2.53))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(40,40,80),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,20.239),sd2=c(1.265,1.265,5.06))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(50,50,25),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,2.53),sd2=c(1.265,1.265,.6325))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(50,50,25),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,5.06),sd2=c(1.265,1.265,1.265))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(50,50,25),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,10.119),sd2=c(1.265,1.265,2.53))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(50,50,25),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,20.239),sd2=c(1.265,1.265,5.06))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(50,50,50),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,2.53),sd2=c(1.265,1.265,.6325))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(50,50,50),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,5.06),sd2=c(1.265,1.265,1.265))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(50,50,50),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,10.119),sd2=c(1.265,1.265,2.53))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(50,50,50),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,20.239),sd2=c(1.265,1.265,5.06))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(50,50,75),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,2.53),sd2=c(1.265,1.265,.6325))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(50,50,75),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,5.06),sd2=c(1.265,1.265,1.265))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(50,50,75),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,10.119),sd2=c(1.265,1.265,2.53))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(50,50,75),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,20.239),sd2=c(1.265,1.265,5.06))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(50,50,100),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,2.53),sd2=c(1.265,1.265,.6325))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(50,50,100),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,5.06),sd2=c(1.265,1.265,1.265))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(50,50,100),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,10.119),sd2=c(1.265,1.265,2.53))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(50,50,100),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,20.239),sd2=c(1.265,1.265,5.06))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(100,100,50),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,2.53),sd2=c(1.265,1.265,.6325))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(100,100,50),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,5.06),sd2=c(1.265,1.265,1.265))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(100,100,50),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,10.119),sd2=c(1.265,1.265,2.53))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(100,100,50),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,20.239),sd2=c(1.265,1.265,5.06))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(100,100,100),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,2.53),sd2=c(1.265,1.265,.6325))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(100,100,100),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,5.06),sd2=c(1.265,1.265,1.265))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(100,100,100),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,10.119),sd2=c(1.265,1.265,2.53))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(100,100,100),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,20.239),sd2=c(1.265,1.265,5.06))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(100,100,150),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,2.53),sd2=c(1.265,1.265,.6325))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(100,100,150),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,5.06),sd2=c(1.265,1.265,1.265))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(100,100,150),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,10.119),sd2=c(1.265,1.265,2.53))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(100,100,150),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,20.239),sd2=c(1.265,1.265,5.06))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(100,100,200),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,2.53),sd2=c(1.265,1.265,.6325))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(100,100,200),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,5.06),sd2=c(1.265,1.265,1.265))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(100,100,200),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,10.119),sd2=c(1.265,1.265,2.53))
get_simu(nSims=1,k=3,distName=rep("mixed",3),n=c(100,100,200),m=c(0,0,0),p1=c(.1,.1,.1),p2=c(.9,.9,.9),sd1=c(5.06,5.06,20.239),sd2=c(1.265,1.265,5.06))

get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(20,20,20,10),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(20,20,20,10),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(20,20,20,10),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(20,20,20,10),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(20,20,20,20),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(20,20,20,20),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(20,20,20,20),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(20,20,20,20),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(20,20,20,30),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(20,20,20,30),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(20,20,20,30),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(20,20,20,30),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(20,20,20,40),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(20,20,20,40),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(20,20,20,40),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(20,20,20,40),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(30,30,30,15),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(30,30,30,15),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(30,30,30,15),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(30,30,30,15),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(30,30,30,30),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(30,30,30,30),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(30,30,30,30),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(30,30,30,30),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(30,30,30,45),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(30,30,30,45),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(30,30,30,45),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(30,30,30,45),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(30,30,30,60),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(30,30,30,60),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(30,30,30,60),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(30,30,30,60),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(40,40,40,20),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(40,40,40,20),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(40,40,40,20),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(40,40,40,20),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(40,40,40,40),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(40,40,40,40),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(40,40,40,40),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(40,40,40,40),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(40,40,40,60),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(40,40,40,60),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(40,40,40,60),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(40,40,40,60),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(40,40,40,80),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(40,40,40,80),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(40,40,40,80),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(40,40,40,80),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(50,50,50,25),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(50,50,50,25),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(50,50,50,25),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(50,50,50,25),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(50,50,50,50),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(50,50,50,50),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(50,50,50,50),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(50,50,50,50),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(50,50,50,75),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(50,50,50,75),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(50,50,50,75),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(50,50,50,75),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(50,50,50,100),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(50,50,50,100),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(50,50,50,100),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(50,50,50,100),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(100,100,100,50),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(100,100,100,50),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(100,100,100,50),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(100,100,100,50),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(100,100,100,100),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(100,100,100,100),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(100,100,100,100),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(100,100,100,100),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(100,100,100,150),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(100,100,100,150),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(100,100,100,150),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(100,100,100,150),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(100,100,100,200),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(100,100,100,200),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(100,100,100,200),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=4,distName=rep("mixed",4),n=c(100,100,100,200),m=c(0,0,0,0),p1=c(.1,.1,.1,.1),p2=c(.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,5.06))

get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(20,20,20,20,10),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(20,20,20,20,10),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(20,20,20,20,10),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(20,20,20,20,10),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(20,20,20,20,20),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(20,20,20,20,20),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(20,20,20,20,20),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(20,20,20,20,20),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(20,20,20,20,30),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(20,20,20,20,30),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(20,20,20,20,30),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(20,20,20,20,30),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(20,20,20,20,40),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(20,20,20,20,40),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(20,20,20,20,40),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(20,20,20,20,40),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(30,30,30,30,15),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(30,30,30,30,15),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(30,30,30,30,15),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(30,30,30,30,15),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(30,30,30,30,30),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(30,30,30,30,30),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(30,30,30,30,30),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(30,30,30,30,30),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(30,30,30,30,45),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(30,30,30,30,45),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(30,30,30,30,45),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(30,30,30,30,45),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(30,30,30,30,60),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(30,30,30,30,60),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(30,30,30,30,60),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(30,30,30,30,60),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(40,40,40,40,20),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(40,40,40,40,20),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(40,40,40,40,20),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(40,40,40,40,20),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(40,40,40,40,40),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(40,40,40,40,40),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(40,40,40,40,40),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(40,40,40,40,40),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(40,40,40,40,60),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(40,40,40,40,60),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(40,40,40,40,60),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(40,40,40,40,60),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(40,40,40,40,80),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(40,40,40,40,80),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(40,40,40,40,80),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(40,40,40,40,80),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(50,50,50,50,25),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(50,50,50,50,25),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(50,50,50,50,25),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(50,50,50,50,25),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(50,50,50,50,50),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(50,50,50,50,50),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(50,50,50,50,50),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(50,50,50,50,50),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(50,50,50,50,75),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(50,50,50,50,75),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(50,50,50,50,75),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(50,50,50,50,75),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(50,50,50,50,100),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(50,50,50,50,100),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(50,50,50,50,100),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(50,50,50,50,100),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(100,100,100,100,50),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(100,100,100,100,50),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(100,100,100,100,50),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(100,100,100,100,50),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(100,100,100,100,100),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(100,100,100,100,100),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(100,100,100,100,100),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(100,100,100,100,100),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(100,100,100,100,150),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(100,100,100,100,150),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(100,100,100,100,150),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(100,100,100,100,150),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,1.265,5.06))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(100,100,100,100,200),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,2.53),sd2=c(1.265,1.265,1.265,1.265,.6325))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(100,100,100,100,200),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,5.06),sd2=c(1.265,1.265,1.265,1.265,1.265))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(100,100,100,100,200),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,10.119),sd2=c(1.265,1.265,1.265,1.265,2.53))
get_simu(nSims=1,k=5,distName=rep("mixed",5),n=c(100,100,100,100,200),m=c(0,0,0,0,0),p1=c(.1,.1,.1,.1,.1),p2=c(.9,.9,.9,.9,.9),sd1=c(5.06,5.06,5.06,5.06,20.239),sd2=c(1.265,1.265,1.265,1.265,5.06))

### chi square and normal right-skewed distributions

get_simu(k=2,distName=c("chi2","skewpos"),n=c(20,10),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(20,10),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(20,10),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(20,10),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(20,20),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(20,20),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(20,20),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(20,20),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(20,30),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(20,30),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(20,30),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(20,30),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(20,40),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(20,40),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(20,40),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(20,40),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(30,15),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(30,15),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(30,15),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(30,15),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(30,30),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(30,30),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(30,30),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(30,30),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(30,45),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(30,45),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(30,45),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(30,45),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(30,60),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(30,60),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(30,60),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(30,60),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(40,20),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(40,20),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(40,20),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(40,20),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(40,40),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(40,40),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(40,40),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(40,40),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(40,60),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(40,60),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(40,60),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(40,60),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(40,80),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(40,80),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(40,80),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(40,80),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(50,25),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(50,25),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(50,25),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(50,25),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(50,50),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(50,50),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(50,50),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(50,50),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(50,75),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(50,75),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(50,75),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(50,75),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(50,100),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(50,100),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(50,100),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(50,100),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(100,50),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(100,50),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(100,50),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(100,50),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(100,100),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(100,100),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(100,100),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(100,100),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(100,150),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(100,150),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(100,150),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(100,150),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(100,200),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(100,200),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(100,200),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(100,200),df=c(2,2),m=rep(2,2),sds=c(2,8))

get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(20,20,10),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(20,20,10),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(20,20,10),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(20,20,10),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(20,20,20),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(20,20,20),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(20,20,20),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(20,20,20),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(20,20,30),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(20,20,30),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(20,20,30),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(20,20,30),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(20,20,40),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(20,20,40),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(20,20,40),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(20,20,40),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(30,30,15),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(30,30,15),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(30,30,15),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(30,30,15),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(30,30,30),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(30,30,30),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(30,30,30),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(30,30,30),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(30,30,45),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(30,30,45),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(30,30,45),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(30,30,45),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(30,30,60),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(30,30,60),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(30,30,60),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(30,30,60),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(40,40,20),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(40,40,20),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(40,40,20),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(40,40,20),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(40,40,40),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(40,40,40),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(40,40,40),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(40,40,40),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(40,40,60),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(40,40,60),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(40,40,60),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(40,40,60),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(40,40,80),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(40,40,80),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(40,40,80),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(40,40,80),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(50,50,25),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(50,50,25),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(50,50,25),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(50,50,25),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(50,50,50),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(50,50,50),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(50,50,50),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(50,50,50),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(50,50,75),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(50,50,75),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(50,50,75),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(50,50,75),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(50,50,100),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(50,50,100),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(50,50,100),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(50,50,100),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(100,100,50),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(100,100,50),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(100,100,50),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(100,100,50),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(100,100,100),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(100,100,100),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(100,100,100),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(100,100,100),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(100,100,150),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(100,100,150),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(100,100,150),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(100,100,150),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(100,100,200),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(100,100,200),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(100,100,200),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(100,100,200),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))

get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(20,20,20,10),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(20,20,20,10),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(20,20,20,10),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(20,20,20,10),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(20,20,20,20),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(20,20,20,20),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(20,20,20,20),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(20,20,20,20),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(20,20,20,30),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(20,20,20,30),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(20,20,20,30),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(20,20,20,30),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(20,20,20,40),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(20,20,20,40),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(20,20,20,40),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(20,20,20,40),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(30,30,30,15),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(30,30,30,15),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(30,30,30,15),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(30,30,30,15),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(30,30,30,30),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(30,30,30,30),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(30,30,30,30),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(30,30,30,30),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(30,30,30,45),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(30,30,30,45),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(30,30,30,45),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(30,30,30,45),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(30,30,30,60),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(30,30,30,60),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(30,30,30,60),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(30,30,30,60),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(40,40,40,20),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(40,40,40,20),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(40,40,40,20),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(40,40,40,20),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(40,40,40,40),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(40,40,40,40),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(40,40,40,40),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(40,40,40,40),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(40,40,40,60),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(40,40,40,60),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(40,40,40,60),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(40,40,40,60),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(40,40,40,80),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(40,40,40,80),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(40,40,40,80),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(40,40,40,80),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(50,50,50,25),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(50,50,50,25),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(50,50,50,25),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(50,50,50,25),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(50,50,50,50),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(50,50,50,50),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(50,50,50,50),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(50,50,50,50),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(50,50,50,75),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(50,50,50,75),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(50,50,50,75),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(50,50,50,75),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(50,50,50,100),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(50,50,50,100),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(50,50,50,100),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(50,50,50,100),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(100,100,100,50),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(100,100,100,50),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(100,100,100,50),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(100,100,100,50),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(100,100,100,100),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(100,100,100,100),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(100,100,100,100),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(100,100,100,100),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(100,100,100,150),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(100,100,100,150),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(100,100,100,150),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(100,100,100,150),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(100,100,100,200),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(100,100,100,200),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(100,100,100,200),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewpos"),n=c(100,100,100,200),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))

get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(20,20,20,20,10),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(20,20,20,20,10),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(20,20,20,20,10),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(20,20,20,20,10),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(20,20,20,20,20),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(20,20,20,20,20),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(20,20,20,20,20),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(20,20,20,20,20),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(20,20,20,20,30),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(20,20,20,20,30),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(20,20,20,20,30),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(20,20,20,20,30),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(20,20,20,20,40),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(20,20,20,20,40),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(20,20,20,20,40),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(20,20,20,20,40),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(30,30,30,30,15),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(30,30,30,30,15),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(30,30,30,30,15),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(30,30,30,30,15),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(30,30,30,30,30),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(30,30,30,30,30),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(30,30,30,30,30),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(30,30,30,30,30),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(30,30,30,30,45),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(30,30,30,30,45),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(30,30,30,30,45),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(30,30,30,30,45),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(30,30,30,30,60),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(30,30,30,30,60),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(30,30,30,30,60),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(30,30,30,30,60),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(40,40,40,40,20),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(40,40,40,40,20),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(40,40,40,40,20),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(40,40,40,40,20),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(40,40,40,40,40),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(40,40,40,40,40),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(40,40,40,40,40),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(40,40,40,40,40),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(40,40,40,40,60),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(40,40,40,40,60),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(40,40,40,40,60),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(40,40,40,40,60),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(40,40,40,40,80),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(40,40,40,40,80),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(40,40,40,40,80),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(40,40,40,40,80),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(50,50,50,50,25),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(50,50,50,50,25),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(50,50,50,50,25),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(50,50,50,50,25),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(50,50,50,50,50),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(50,50,50,50,50),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(50,50,50,50,50),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(50,50,50,50,50),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(50,50,50,50,75),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(50,50,50,50,75),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(50,50,50,50,75),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(50,50,50,50,75),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(50,50,50,50,100),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(50,50,50,50,100),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(50,50,50,50,100),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(50,50,50,50,100),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(100,100,100,100,50),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(100,100,100,100,50),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(100,100,100,100,50),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(100,100,100,100,50),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(100,100,100,100,100),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(100,100,100,100,100),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(100,100,100,100,100),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(100,100,100,100,100),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(100,100,100,100,150),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(100,100,100,100,150),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(100,100,100,100,150),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(100,100,100,100,150),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(100,100,100,100,200),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(100,100,100,100,200),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(100,100,100,100,200),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewpos"),n=c(100,100,100,100,200),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))

### chi square and normal left-skewed distributions

get_simu(k=2,distName=c("chi2","skewneg"),n=c(20,10),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(20,10),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(20,10),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(20,10),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(20,20),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(20,20),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(20,20),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(20,20),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(20,30),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(20,30),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(20,30),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(20,30),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(20,40),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(20,40),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(20,40),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(20,40),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(30,15),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(30,15),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(30,15),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(30,15),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(30,30),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(30,30),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(30,30),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(30,30),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(30,45),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(30,45),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(30,45),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(30,45),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(30,60),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(30,60),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(30,60),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(30,60),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(40,20),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(40,20),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(40,20),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(40,20),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(40,40),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(40,40),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(40,40),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(40,40),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(40,60),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(40,60),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(40,60),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(40,60),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(40,80),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(40,80),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(40,80),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(40,80),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(50,25),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(50,25),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(50,25),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(50,25),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(50,50),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(50,50),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(50,50),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(50,50),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(50,75),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(50,75),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(50,75),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(50,75),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(50,100),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(50,100),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(50,100),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(50,100),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(100,50),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(100,50),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(100,50),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(100,50),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(100,100),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(100,100),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(100,100),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(100,100),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(100,150),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(100,150),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(100,150),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(100,150),df=c(2,2),m=rep(2,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(100,200),df=c(2,2),m=rep(2,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(100,200),df=c(2,2),m=rep(2,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(100,200),df=c(2,2),m=rep(2,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(100,200),df=c(2,2),m=rep(2,2),sds=c(2,8))

get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(20,20,10),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(20,20,10),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(20,20,10),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(20,20,10),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(20,20,20),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(20,20,20),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(20,20,20),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(20,20,20),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(20,20,30),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(20,20,30),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(20,20,30),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(20,20,30),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(20,20,40),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(20,20,40),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(20,20,40),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(20,20,40),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(30,30,15),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(30,30,15),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(30,30,15),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(30,30,15),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(30,30,30),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(30,30,30),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(30,30,30),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(30,30,30),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(30,30,45),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(30,30,45),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(30,30,45),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(30,30,45),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(30,30,60),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(30,30,60),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(30,30,60),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(30,30,60),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(40,40,20),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(40,40,20),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(40,40,20),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(40,40,20),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(40,40,40),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(40,40,40),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(40,40,40),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(40,40,40),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(40,40,60),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(40,40,60),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(40,40,60),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(40,40,60),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(40,40,80),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(40,40,80),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(40,40,80),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(40,40,80),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(50,50,25),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(50,50,25),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(50,50,25),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(50,50,25),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(50,50,50),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(50,50,50),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(50,50,50),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(50,50,50),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(50,50,75),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(50,50,75),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(50,50,75),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(50,50,75),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(50,50,100),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(50,50,100),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(50,50,100),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(50,50,100),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(100,100,50),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(100,100,50),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(100,100,50),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(100,100,50),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(100,100,100),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(100,100,100),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(100,100,100),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(100,100,100),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(100,100,150),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(100,100,150),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(100,100,150),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(100,100,150),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(100,100,200),df=c(2,2,2),m=rep(2,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(100,100,200),df=c(2,2,2),m=rep(2,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(100,100,200),df=c(2,2,2),m=rep(2,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(100,100,200),df=c(2,2,2),m=rep(2,3),sds=c(2,2,8))

get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(20,20,20,10),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(20,20,20,10),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(20,20,20,10),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(20,20,20,10),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(20,20,20,20),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(20,20,20,20),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(20,20,20,20),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(20,20,20,20),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(20,20,20,30),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(20,20,20,30),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(20,20,20,30),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(20,20,20,30),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(20,20,20,40),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(20,20,20,40),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(20,20,20,40),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(20,20,20,40),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(30,30,30,15),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(30,30,30,15),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(30,30,30,15),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(30,30,30,15),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(30,30,30,30),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(30,30,30,30),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(30,30,30,30),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(30,30,30,30),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(30,30,30,45),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(30,30,30,45),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(30,30,30,45),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(30,30,30,45),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(30,30,30,60),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(30,30,30,60),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(30,30,30,60),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(30,30,30,60),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(40,40,40,20),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(40,40,40,20),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(40,40,40,20),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(40,40,40,20),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(40,40,40,40),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(40,40,40,40),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(40,40,40,40),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(40,40,40,40),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(40,40,40,60),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(40,40,40,60),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(40,40,40,60),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(40,40,40,60),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(40,40,40,80),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(40,40,40,80),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(40,40,40,80),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(40,40,40,80),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(50,50,50,25),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(50,50,50,25),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(50,50,50,25),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(50,50,50,25),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(50,50,50,50),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(50,50,50,50),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(50,50,50,50),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(50,50,50,50),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(50,50,50,75),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(50,50,50,75),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(50,50,50,75),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(50,50,50,75),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(50,50,50,100),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(50,50,50,100),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(50,50,50,100),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(50,50,50,100),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(100,100,100,50),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(100,100,100,50),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(100,100,100,50),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(100,100,100,50),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(100,100,100,100),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(100,100,100,100),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(100,100,100,100),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(100,100,100,100),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(100,100,100,150),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(100,100,100,150),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(100,100,100,150),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(100,100,100,150),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(100,100,100,200),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,1))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(100,100,100,200),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,2))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(100,100,100,200),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,4))
get_simu(k=4,distName=c("chi2","chi2","chi2","skewneg"),n=c(100,100,100,200),df=c(2,2,2,2),m=rep(2,4),sds=c(2,2,2,8))

get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(20,20,20,20,10),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(20,20,20,20,10),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(20,20,20,20,10),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(20,20,20,20,10),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(20,20,20,20,20),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(20,20,20,20,20),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(20,20,20,20,20),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(20,20,20,20,20),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(20,20,20,20,30),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(20,20,20,20,30),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(20,20,20,20,30),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(20,20,20,20,30),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(20,20,20,20,40),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(20,20,20,20,40),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(20,20,20,20,40),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(20,20,20,20,40),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(30,30,30,30,15),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(30,30,30,30,15),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(30,30,30,30,15),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(30,30,30,30,15),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(30,30,30,30,30),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(30,30,30,30,30),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(30,30,30,30,30),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(30,30,30,30,30),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(30,30,30,30,45),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(30,30,30,30,45),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(30,30,30,30,45),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(30,30,30,30,45),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(30,30,30,30,60),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(30,30,30,30,60),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(30,30,30,30,60),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(30,30,30,30,60),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(40,40,40,40,20),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(40,40,40,40,20),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(40,40,40,40,20),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(40,40,40,40,20),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(40,40,40,40,40),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(40,40,40,40,40),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(40,40,40,40,40),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(40,40,40,40,40),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(40,40,40,40,60),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(40,40,40,40,60),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(40,40,40,40,60),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(40,40,40,40,60),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(40,40,40,40,80),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(40,40,40,40,80),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(40,40,40,40,80),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(40,40,40,40,80),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(50,50,50,50,25),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(50,50,50,50,25),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(50,50,50,50,25),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(50,50,50,50,25),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(50,50,50,50,50),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(50,50,50,50,50),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(50,50,50,50,50),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(50,50,50,50,50),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(50,50,50,50,75),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(50,50,50,50,75),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(50,50,50,50,75),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(50,50,50,50,75),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(50,50,50,50,100),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(50,50,50,50,100),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(50,50,50,50,100),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(50,50,50,50,100),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(100,100,100,100,50),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(100,100,100,100,50),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(100,100,100,100,50),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(100,100,100,100,50),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(100,100,100,100,100),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(100,100,100,100,100),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(100,100,100,100,100),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(100,100,100,100,100),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(100,100,100,100,150),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(100,100,100,100,150),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(100,100,100,100,150),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(100,100,100,100,150),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(100,100,100,100,200),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,1))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(100,100,100,100,200),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,2))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(100,100,100,100,200),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,4))
get_simu(k=5,distName=c("chi2","chi2","chi2","chi2","skewneg"),n=c(100,100,100,100,200),df=c(2,2,2,2,2),m=rep(2,5),sds=c(2,2,2,2,8))

### normal right-skewed and normal left-skewed distributions

get_simu(k=2,distName=c("skewpos","skewneg"),n=c(20,10),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(20,10),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(20,10),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(20,10),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(20,20),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(20,20),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(20,20),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(20,20),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(20,30),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(20,30),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(20,30),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(20,30),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(20,40),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(20,40),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(20,40),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(20,40),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(30,15),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(30,15),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(30,15),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(30,15),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(30,30),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(30,30),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(30,30),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(30,30),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(30,45),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(30,45),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(30,45),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(30,45),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(30,60),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(30,60),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(30,60),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(30,60),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(40,20),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(40,20),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(40,20),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(40,20),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(40,40),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(40,40),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(40,40),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(40,40),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(40,60),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(40,60),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(40,60),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(40,60),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(40,80),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(40,80),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(40,80),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(40,80),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(50,25),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(50,25),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(50,25),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(50,25),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(50,50),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(50,50),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(50,50),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(50,50),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(50,75),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(50,75),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(50,75),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(50,75),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(50,100),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(50,100),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(50,100),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(50,100),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(100,50),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(100,50),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(100,50),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(100,50),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(100,100),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(100,100),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(100,100),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(100,100),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(100,150),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(100,150),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(100,150),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(100,150),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(100,200),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(100,200),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(100,200),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(100,200),m=c(0,0),sds=c(2,8)) # inequal skewness

get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(20,20,10),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(20,20,10),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(20,20,10),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(20,20,10),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(20,20,20),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(20,20,20),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(20,20,20),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(20,20,20),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(20,20,30),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(20,20,30),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(20,20,30),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(20,20,30),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(20,20,40),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(20,20,40),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(20,20,40),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(20,20,40),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(30,30,15),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(30,30,15),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(30,30,15),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(30,30,15),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(30,30,30),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(30,30,30),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(30,30,30),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(30,30,30),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(30,30,45),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(30,30,45),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(30,30,45),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(30,30,45),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(30,30,60),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(30,30,60),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(30,30,60),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(30,30,60),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(40,40,20),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(40,40,20),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(40,40,20),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(40,40,20),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(40,40,40),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(40,40,40),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(40,40,40),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(40,40,40),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(40,40,60),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(40,40,60),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(40,40,60),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(40,40,60),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(40,40,80),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(40,40,80),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(40,40,80),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(40,40,80),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(50,50,25),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(50,50,25),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(50,50,25),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(50,50,25),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(50,50,50),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(50,50,50),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(50,50,50),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(50,50,50),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(50,50,75),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(50,50,75),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(50,50,75),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(50,50,75),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(50,50,100),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(50,50,100),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(50,50,100),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(50,50,100),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(100,100,50),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(100,100,50),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(100,100,50),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(100,100,50),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(100,100,100),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(100,100,100),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(100,100,100),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(100,100,100),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(100,100,150),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(100,100,150),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(100,100,150),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(100,100,150),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(100,100,200),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(100,100,200),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(100,100,200),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(100,100,200),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness

get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(20,20,20,10),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(20,20,20,10),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(20,20,20,10),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(20,20,20,10),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(20,20,20,20),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(20,20,20,20),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(20,20,20,20),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(20,20,20,20),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(20,20,20,30),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(20,20,20,30),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(20,20,20,30),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(20,20,20,30),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(20,20,20,40),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(20,20,20,40),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(20,20,20,40),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(20,20,20,40),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(30,30,30,15),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(30,30,30,15),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(30,30,30,15),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(30,30,30,15),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(30,30,30,30),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(30,30,30,30),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(30,30,30,30),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(30,30,30,30),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(30,30,30,45),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(30,30,30,45),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(30,30,30,45),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(30,30,30,45),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(30,30,30,60),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(30,30,30,60),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(30,30,30,60),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(30,30,30,60),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(40,40,40,20),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(40,40,40,20),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(40,40,40,20),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(40,40,40,20),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(40,40,40,40),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(40,40,40,40),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(40,40,40,40),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(40,40,40,40),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(40,40,40,60),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(40,40,40,60),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(40,40,40,60),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(40,40,40,60),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(40,40,40,80),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(40,40,40,80),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(40,40,40,80),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(40,40,40,80),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(50,50,50,25),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(50,50,50,25),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(50,50,50,25),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(50,50,50,25),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(50,50,50,50),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(50,50,50,50),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(50,50,50,50),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(50,50,50,50),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(50,50,50,75),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(50,50,50,75),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(50,50,50,75),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(50,50,50,75),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(50,50,50,100),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(50,50,50,100),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(50,50,50,100),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(50,50,50,100),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(100,100,100,50),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(100,100,100,50),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(100,100,100,50),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(100,100,100,50),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(100,100,100,100),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(100,100,100,100),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(100,100,100,100),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(100,100,100,100),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(100,100,100,150),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(100,100,100,150),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(100,100,100,150),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(100,100,100,150),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(100,100,100,200),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(100,100,100,200),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(100,100,100,200),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewneg"),n=c(100,100,100,200),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness

get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(20,20,20,20,10),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(20,20,20,20,10),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(20,20,20,20,10),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(20,20,20,20,10),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(20,20,20,20,20),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(20,20,20,20,20),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(20,20,20,20,20),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(20,20,20,20,20),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(20,20,20,20,30),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(20,20,20,20,30),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(20,20,20,20,30),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(20,20,20,20,30),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(20,20,20,20,40),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(20,20,20,20,40),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(20,20,20,20,40),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(20,20,20,20,40),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(30,30,30,30,15),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(30,30,30,30,15),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(30,30,30,30,15),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(30,30,30,30,15),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(30,30,30,30,30),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(30,30,30,30,30),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(30,30,30,30,30),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(30,30,30,30,30),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(30,30,30,30,45),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(30,30,30,30,45),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(30,30,30,30,45),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(30,30,30,30,45),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(30,30,30,30,60),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(30,30,30,30,60),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(30,30,30,30,60),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(30,30,30,30,60),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(40,40,40,40,20),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(40,40,40,40,20),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(40,40,40,40,20),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(40,40,40,40,20),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(40,40,40,40,40),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(40,40,40,40,40),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(40,40,40,40,40),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(40,40,40,40,40),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(40,40,40,40,60),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(40,40,40,40,60),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(40,40,40,40,60),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(40,40,40,40,60),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(40,40,40,40,80),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(40,40,40,40,80),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(40,40,40,40,80),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(40,40,40,40,80),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(50,50,50,50,25),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(50,50,50,50,25),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(50,50,50,50,25),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(50,50,50,50,25),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(50,50,50,50,50),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(50,50,50,50,50),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(50,50,50,50,50),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(50,50,50,50,50),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(50,50,50,50,75),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(50,50,50,50,75),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(50,50,50,50,75),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(50,50,50,50,75),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(50,50,50,50,100),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(50,50,50,50,100),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(50,50,50,50,100),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(50,50,50,50,100),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(100,100,100,100,50),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(100,100,100,100,50),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(100,100,100,100,50),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(100,100,100,100,50),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(100,100,100,100,100),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(100,100,100,100,100),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(100,100,100,100,100),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(100,100,100,100,100),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(100,100,100,100,150),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(100,100,100,100,150),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(100,100,100,100,150),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(100,100,100,100,150),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(100,100,100,100,200),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(100,100,100,100,200),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(100,100,100,100,200),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewneg"),n=c(100,100,100,100,200),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness

### normal right-skewed distributions

get_simu(k=2,distName=c("skewpos","skewpos"),n=c(20,10),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(20,10),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(20,10),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(20,10),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(20,20),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(20,20),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(20,20),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(20,20),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(20,30),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(20,30),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(20,30),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(20,30),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(20,40),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(20,40),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(20,40),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(20,40),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(30,15),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(30,15),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(30,15),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(30,15),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(30,30),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(30,30),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(30,30),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(30,30),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(30,45),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(30,45),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(30,45),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(30,45),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(30,60),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(30,60),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(30,60),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(30,60),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(40,20),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(40,20),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(40,20),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(40,20),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(40,40),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(40,40),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(40,40),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(40,40),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(40,60),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(40,60),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(40,60),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(40,60),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(40,80),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(40,80),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(40,80),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(40,80),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(50,25),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(50,25),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(50,25),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(50,25),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(50,50),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(50,50),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(50,50),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(50,50),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(50,75),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(50,75),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(50,75),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(50,75),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(50,100),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(50,100),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(50,100),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(50,100),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(100,50),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(100,50),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(100,50),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(100,50),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(100,100),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(100,100),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(100,100),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(100,100),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(100,150),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(100,150),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(100,150),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(100,150),m=c(0,0),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(100,200),m=c(0,0),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(100,200),m=c(0,0),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(100,200),m=c(0,0),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(100,200),m=c(0,0),sds=c(2,8)) # inequal skewness

get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(20,20,10),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(20,20,10),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(20,20,10),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(20,20,10),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(20,20,20),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(20,20,20),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(20,20,20),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(20,20,20),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(20,20,30),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(20,20,30),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(20,20,30),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(20,20,30),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(20,20,40),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(20,20,40),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(20,20,40),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(20,20,40),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(30,30,15),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(30,30,15),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(30,30,15),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(30,30,15),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(30,30,30),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(30,30,30),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(30,30,30),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(30,30,30),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(30,30,45),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(30,30,45),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(30,30,45),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(30,30,45),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(30,30,60),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(30,30,60),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(30,30,60),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(30,30,60),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(40,40,20),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(40,40,20),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(40,40,20),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(40,40,20),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(40,40,40),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(40,40,40),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(40,40,40),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(40,40,40),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(40,40,60),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(40,40,60),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(40,40,60),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(40,40,60),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(40,40,80),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(40,40,80),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(40,40,80),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(40,40,80),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(50,50,25),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(50,50,25),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(50,50,25),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(50,50,25),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(50,50,50),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(50,50,50),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(50,50,50),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(50,50,50),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(50,50,75),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(50,50,75),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(50,50,75),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(50,50,75),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(50,50,100),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(50,50,100),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(50,50,100),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(50,50,100),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(100,100,50),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(100,100,50),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(100,100,50),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(100,100,50),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(100,100,100),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(100,100,100),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(100,100,100),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(100,100,100),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(100,100,150),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(100,100,150),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(100,100,150),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(100,100,150),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(100,100,200),m=c(0,0,0),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(100,100,200),m=c(0,0,0),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(100,100,200),m=c(0,0,0),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(100,100,200),m=c(0,0,0),sds=c(2,2,8)) # inequal skewness

get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(20,20,20,10),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(20,20,20,10),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(20,20,20,10),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(20,20,20,10),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(20,20,20,20),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(20,20,20,20),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(20,20,20,20),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(20,20,20,20),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(20,20,20,30),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(20,20,20,30),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(20,20,20,30),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(20,20,20,30),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(20,20,20,40),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(20,20,20,40),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(20,20,20,40),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(20,20,20,40),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(30,30,30,15),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(30,30,30,15),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(30,30,30,15),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(30,30,30,15),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(30,30,30,30),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(30,30,30,30),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(30,30,30,30),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(30,30,30,30),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(30,30,30,45),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(30,30,30,45),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(30,30,30,45),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(30,30,30,45),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(30,30,30,60),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(30,30,30,60),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(30,30,30,60),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(30,30,30,60),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(40,40,40,20),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(40,40,40,20),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(40,40,40,20),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(40,40,40,20),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(40,40,40,40),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(40,40,40,40),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(40,40,40,40),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(40,40,40,40),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(40,40,40,60),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(40,40,40,60),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(40,40,40,60),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(40,40,40,60),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(40,40,40,80),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(40,40,40,80),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(40,40,40,80),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(40,40,40,80),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(50,50,50,25),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(50,50,50,25),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(50,50,50,25),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(50,50,50,25),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(50,50,50,50),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(50,50,50,50),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(50,50,50,50),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(50,50,50,50),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(50,50,50,75),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(50,50,50,75),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(50,50,50,75),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(50,50,50,75),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(50,50,50,100),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(50,50,50,100),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(50,50,50,100),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(50,50,50,100),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(100,100,100,50),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(100,100,100,50),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(100,100,100,50),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(100,100,100,50),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(100,100,100,100),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(100,100,100,100),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(100,100,100,100),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(100,100,100,100),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(100,100,100,150),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(100,100,100,150),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(100,100,100,150),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(100,100,100,150),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(100,100,100,200),m=c(0,0,0,0),sds=c(2,2,2,1)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(100,100,100,200),m=c(0,0,0,0),sds=c(2,2,2,2)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(100,100,100,200),m=c(0,0,0,0),sds=c(2,2,2,4)) # inequal skewness
get_simu(k=4,distName=c("skewpos","skewpos","skewpos","skewpos"),n=c(100,100,100,200),m=c(0,0,0,0),sds=c(2,2,2,8)) # inequal skewness

get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(20,20,20,20,10),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(20,20,20,20,10),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(20,20,20,20,10),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(20,20,20,20,10),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(20,20,20,20,20),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(20,20,20,20,20),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(20,20,20,20,20),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(20,20,20,20,20),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(20,20,20,20,30),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(20,20,20,20,30),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(20,20,20,20,30),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(20,20,20,20,30),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(20,20,20,20,40),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(20,20,20,20,40),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(20,20,20,20,40),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(20,20,20,20,40),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(30,30,30,30,15),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(30,30,30,30,15),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(30,30,30,30,15),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(30,30,30,30,15),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(30,30,30,30,30),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(30,30,30,30,30),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(30,30,30,30,30),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(30,30,30,30,30),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(30,30,30,30,45),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(30,30,30,30,45),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(30,30,30,30,45),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(30,30,30,30,45),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(30,30,30,30,60),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(30,30,30,30,60),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(30,30,30,30,60),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(30,30,30,30,60),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(40,40,40,40,20),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(40,40,40,40,20),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(40,40,40,40,20),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(40,40,40,40,20),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(40,40,40,40,40),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(40,40,40,40,40),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(40,40,40,40,40),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(40,40,40,40,40),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(40,40,40,40,60),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(40,40,40,40,60),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(40,40,40,40,60),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(40,40,40,40,60),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(40,40,40,40,80),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(40,40,40,40,80),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(40,40,40,40,80),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(40,40,40,40,80),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(50,50,50,50,25),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(50,50,50,50,25),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(50,50,50,50,25),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(50,50,50,50,25),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(50,50,50,50,50),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(50,50,50,50,50),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(50,50,50,50,50),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(50,50,50,50,50),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(50,50,50,50,75),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(50,50,50,50,75),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(50,50,50,50,75),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(50,50,50,50,75),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(50,50,50,50,100),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(50,50,50,50,100),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(50,50,50,50,100),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(50,50,50,50,100),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(100,100,100,100,50),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(100,100,100,100,50),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(100,100,100,100,50),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(100,100,100,100,50),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(100,100,100,100,100),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(100,100,100,100,100),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(100,100,100,100,100),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(100,100,100,100,100),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(100,100,100,100,150),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(100,100,100,100,150),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(100,100,100,100,150),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(100,100,100,100,150),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(100,100,100,100,200),m=c(0,0,0,0,0),sds=c(2,2,2,2,1)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(100,100,100,100,200),m=c(0,0,0,0,0),sds=c(2,2,2,2,2)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(100,100,100,100,200),m=c(0,0,0,0,0),sds=c(2,2,2,2,4)) # inequal skewness
get_simu(k=5,distName=c("skewpos","skewpos","skewpos","skewpos","skewpos"),n=c(100,100,100,100,200),m=c(0,0,0,0,0),sds=c(2,2,2,2,8)) # inequal skewness

#############                       Power                         #############

### normal distributions

get_simu(k=2,distName=rep("normal",2),n=c(20,10),m=c(0,1),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(20,10),m=c(0,1),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(20,10),m=c(0,1),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(20,10),m=c(0,1),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(20,20),m=c(0,1),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(20,20),m=c(0,1),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(20,20),m=c(0,1),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(20,20),m=c(0,1),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(20,30),m=c(0,1),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(20,30),m=c(0,1),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(20,30),m=c(0,1),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(20,30),m=c(0,1),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(20,40),m=c(0,1),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(20,40),m=c(0,1),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(20,40),m=c(0,1),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(20,40),m=c(0,1),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(30,15),m=c(0,1),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(30,15),m=c(0,1),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(30,15),m=c(0,1),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(30,15),m=c(0,1),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(30,30),m=c(0,1),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(30,30),m=c(0,1),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(30,30),m=c(0,1),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(30,30),m=c(0,1),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(30,45),m=c(0,1),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(30,45),m=c(0,1),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(30,45),m=c(0,1),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(30,45),m=c(0,1),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(30,60),m=c(0,1),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(30,60),m=c(0,1),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(30,60),m=c(0,1),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(30,60),m=c(0,1),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(40,20),m=c(0,1),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(40,20),m=c(0,1),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(40,20),m=c(0,1),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(40,20),m=c(0,1),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(40,40),m=c(0,1),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(40,40),m=c(0,1),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(40,40),m=c(0,1),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(40,40),m=c(0,1),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(40,60),m=c(0,1),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(40,60),m=c(0,1),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(40,60),m=c(0,1),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(40,60),m=c(0,1),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(40,80),m=c(0,1),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(40,80),m=c(0,1),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(40,80),m=c(0,1),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(40,80),m=c(0,1),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(50,25),m=c(0,1),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(50,25),m=c(0,1),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(50,25),m=c(0,1),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(50,25),m=c(0,1),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(50,50),m=c(0,1),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(50,50),m=c(0,1),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(50,50),m=c(0,1),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(50,50),m=c(0,1),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(50,75),m=c(0,1),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(50,75),m=c(0,1),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(50,75),m=c(0,1),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(50,75),m=c(0,1),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(50,100),m=c(0,1),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(50,100),m=c(0,1),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(50,100),m=c(0,1),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(50,100),m=c(0,1),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(100,50),m=c(0,1),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(100,50),m=c(0,1),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(100,50),m=c(0,1),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(100,50),m=c(0,1),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(100,100),m=c(0,1),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(100,100),m=c(0,1),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(100,100),m=c(0,1),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(100,100),m=c(0,1),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(100,150),m=c(0,1),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(100,150),m=c(0,1),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(100,150),m=c(0,1),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(100,150),m=c(0,1),sds=c(2,8))
get_simu(k=2,distName=rep("normal",2),n=c(100,200),m=c(0,1),sds=c(2,1))
get_simu(k=2,distName=rep("normal",2),n=c(100,200),m=c(0,1),sds=c(2,2))
get_simu(k=2,distName=rep("normal",2),n=c(100,200),m=c(0,1),sds=c(2,4))
get_simu(k=2,distName=rep("normal",2),n=c(100,200),m=c(0,1),sds=c(2,8))

get_simu(k=3,distName=rep("normal",3),n=c(20,20,10),m=c(0,0,1),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(20,20,10),m=c(0,0,1),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(20,20,10),m=c(0,0,1),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(20,20,10),m=c(0,0,1),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(20,20,20),m=c(0,0,1),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(20,20,20),m=c(0,0,1),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(20,20,20),m=c(0,0,1),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(20,20,20),m=c(0,0,1),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(20,20,30),m=c(0,0,1),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(20,20,30),m=c(0,0,1),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(20,20,30),m=c(0,0,1),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(20,20,30),m=c(0,0,1),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(20,20,40),m=c(0,0,1),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(20,20,40),m=c(0,0,1),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(20,20,40),m=c(0,0,1),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(20,20,40),m=c(0,0,1),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(30,30,15),m=c(0,0,1),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(30,30,15),m=c(0,0,1),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(30,30,15),m=c(0,0,1),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(30,30,15),m=c(0,0,1),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(30,30,30),m=c(0,0,1),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(30,30,30),m=c(0,0,1),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(30,30,30),m=c(0,0,1),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(30,30,30),m=c(0,0,1),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(30,30,45),m=c(0,0,1),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(30,30,45),m=c(0,0,1),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(30,30,45),m=c(0,0,1),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(30,30,45),m=c(0,0,1),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(30,30,60),m=c(0,0,1),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(30,30,60),m=c(0,0,1),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(30,30,60),m=c(0,0,1),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(30,30,60),m=c(0,0,1),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(40,40,20),m=c(0,0,1),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(40,40,20),m=c(0,0,1),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(40,40,20),m=c(0,0,1),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(40,40,20),m=c(0,0,1),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(40,40,40),m=c(0,0,1),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(40,40,40),m=c(0,0,1),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(40,40,40),m=c(0,0,1),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(40,40,40),m=c(0,0,1),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(40,40,60),m=c(0,0,1),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(40,40,60),m=c(0,0,1),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(40,40,60),m=c(0,0,1),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(40,40,60),m=c(0,0,1),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(40,40,80),m=c(0,0,1),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(40,40,80),m=c(0,0,1),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(40,40,80),m=c(0,0,1),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(40,40,80),m=c(0,0,1),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(50,50,25),m=c(0,0,1),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(50,50,25),m=c(0,0,1),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(50,50,25),m=c(0,0,1),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(50,50,25),m=c(0,0,1),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(50,50,50),m=c(0,0,1),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(50,50,50),m=c(0,0,1),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(50,50,50),m=c(0,0,1),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(50,50,50),m=c(0,0,1),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(50,50,75),m=c(0,0,1),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(50,50,75),m=c(0,0,1),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(50,50,75),m=c(0,0,1),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(50,50,75),m=c(0,0,1),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(50,50,100),m=c(0,0,1),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(50,50,100),m=c(0,0,1),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(50,50,100),m=c(0,0,1),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(50,50,100),m=c(0,0,1),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(100,100,50),m=c(0,0,1),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(100,100,50),m=c(0,0,1),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(100,100,50),m=c(0,0,1),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(100,100,50),m=c(0,0,1),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(100,100,100),m=c(0,0,1),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(100,100,100),m=c(0,0,1),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(100,100,100),m=c(0,0,1),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(100,100,100),m=c(0,0,1),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(100,100,150),m=c(0,0,1),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(100,100,150),m=c(0,0,1),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(100,100,150),m=c(0,0,1),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(100,100,150),m=c(0,0,1),sds=c(2,2,8))
get_simu(k=3,distName=rep("normal",3),n=c(100,100,200),m=c(0,0,1),sds=c(2,2,1))
get_simu(k=3,distName=rep("normal",3),n=c(100,100,200),m=c(0,0,1),sds=c(2,2,2))
get_simu(k=3,distName=rep("normal",3),n=c(100,100,200),m=c(0,0,1),sds=c(2,2,4))
get_simu(k=3,distName=rep("normal",3),n=c(100,100,200),m=c(0,0,1),sds=c(2,2,8))

### double exponential distributions

get_simu(k=2,distName=rep("doublex",2),n=c(20,10),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(20,10),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(20,10),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(20,10),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(20,20),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(20,20),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(20,20),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(20,20),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(20,30),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(20,30),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(20,30),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(20,30),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(20,40),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(20,40),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(20,40),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(20,40),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(30,15),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(30,15),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(30,15),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(30,15),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(30,30),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(30,30),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(30,30),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(30,30),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(30,45),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(30,45),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(30,45),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(30,45),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(30,60),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(30,60),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(30,60),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(30,60),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(40,20),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(40,20),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(40,20),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(40,20),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(40,40),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(40,40),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(40,40),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(40,40),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(40,60),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(40,60),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(40,60),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(40,60),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(40,80),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(40,80),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(40,80),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(40,80),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(50,25),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(50,25),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(50,25),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(50,25),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(50,50),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(50,50),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(50,50),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(50,50),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(50,75),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(50,75),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(50,75),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(50,75),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(50,100),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(50,100),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(50,100),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(50,100),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(100,50),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(100,50),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(100,50),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(100,50),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(100,100),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(100,100),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(100,100),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(100,100),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(100,150),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(100,150),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(100,150),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(100,150),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex",2),n=c(100,200),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex",2),n=c(100,200),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex",2),n=c(100,200),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex",2),n=c(100,200),m=c(0,1),lambda=c(2,8))

get_simu(k=3,distName=rep("doublex",3),n=c(20,20,10),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(20,20,10),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(20,20,10),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(20,20,10),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(20,20,20),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(20,20,20),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(20,20,20),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(20,20,20),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(20,20,30),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(20,20,30),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(20,20,30),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(20,20,30),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(20,20,40),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(20,20,40),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(20,20,40),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(20,20,40),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(30,30,15),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(30,30,15),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(30,30,15),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(30,30,15),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(30,30,30),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(30,30,30),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(30,30,30),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(30,30,30),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(30,30,45),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(30,30,45),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(30,30,45),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(30,30,45),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(30,30,60),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(30,30,60),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(30,30,60),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(30,30,60),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(40,40,20),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(40,40,20),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(40,40,20),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(40,40,20),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(40,40,40),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(40,40,40),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(40,40,40),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(40,40,40),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(40,40,60),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(40,40,60),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(40,40,60),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(40,40,60),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(40,40,80),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(40,40,80),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(40,40,80),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(40,40,80),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(50,50,25),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(50,50,25),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(50,50,25),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(50,50,25),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(50,50,50),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(50,50,50),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(50,50,50),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(50,50,50),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(50,50,75),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(50,50,75),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(50,50,75),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(50,50,75),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(50,50,100),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(50,50,100),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(50,50,100),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(50,50,100),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(100,100,50),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(100,100,50),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(100,100,50),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(100,100,50),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(100,100,100),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(100,100,100),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(100,100,100),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(100,100,100),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(100,100,150),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(100,100,150),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(100,100,150),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(100,100,150),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex",3),n=c(100,100,200),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex",3),n=c(100,100,200),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex",3),n=c(100,100,200),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex",3),n=c(100,100,200),m=c(0,0,1),lambda=c(2,2,8))

### double exponential distribution (using default scale parameter of the distribution)

get_simu(k=2,distName=rep("doublex_SP",2),n=c(20,10),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(20,10),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(20,10),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(20,10),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(20,20),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(20,20),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(20,20),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(20,20),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(20,30),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(20,30),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(20,30),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(20,30),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(20,40),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(20,40),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(20,40),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(20,40),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(30,15),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(30,15),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(30,15),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(30,15),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(30,30),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(30,30),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(30,30),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(30,30),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(30,45),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(30,45),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(30,45),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(30,45),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(30,60),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(30,60),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(30,60),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(30,60),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(40,20),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(40,20),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(40,20),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(40,20),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(40,40),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(40,40),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(40,40),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(40,40),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(40,60),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(40,60),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(40,60),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(40,60),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(40,80),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(40,80),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(40,80),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(40,80),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(50,25),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(50,25),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(50,25),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(50,25),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(50,50),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(50,50),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(50,50),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(50,50),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(50,75),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(50,75),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(50,75),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(50,75),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(50,100),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(50,100),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(50,100),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(50,100),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(100,50),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(100,50),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(100,50),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(100,50),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(100,100),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(100,100),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(100,100),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(100,100),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(100,150),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(100,150),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(100,150),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(100,150),m=c(0,1),lambda=c(2,8))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(100,200),m=c(0,1),lambda=c(2,1))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(100,200),m=c(0,1),lambda=c(2,2))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(100,200),m=c(0,1),lambda=c(2,4))
get_simu(k=2,distName=rep("doublex_SP",2),n=c(100,200),m=c(0,1),lambda=c(2,8))

get_simu(k=3,distName=rep("doublex_SP",3),n=c(20,20,10),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(20,20,10),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(20,20,10),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(20,20,10),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(20,20,20),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(20,20,20),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(20,20,20),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(20,20,20),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(20,20,30),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(20,20,30),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(20,20,30),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(20,20,30),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(20,20,40),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(20,20,40),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(20,20,40),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(20,20,40),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(30,30,15),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(30,30,15),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(30,30,15),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(30,30,15),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(30,30,30),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(30,30,30),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(30,30,30),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(30,30,30),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(30,30,45),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(30,30,45),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(30,30,45),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(30,30,45),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(30,30,60),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(30,30,60),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(30,30,60),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(30,30,60),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(40,40,20),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(40,40,20),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(40,40,20),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(40,40,20),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(40,40,40),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(40,40,40),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(40,40,40),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(40,40,40),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(40,40,60),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(40,40,60),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(40,40,60),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(40,40,60),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(40,40,80),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(40,40,80),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(40,40,80),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(40,40,80),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(50,50,25),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(50,50,25),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(50,50,25),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(50,50,25),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(50,50,50),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(50,50,50),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(50,50,50),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(50,50,50),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(50,50,75),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(50,50,75),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(50,50,75),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(50,50,75),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(50,50,100),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(50,50,100),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(50,50,100),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(50,50,100),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(100,100,50),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(100,100,50),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(100,100,50),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(100,100,50),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(100,100,100),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(100,100,100),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(100,100,100),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(100,100,100),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(100,100,150),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(100,100,150),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(100,100,150),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(100,100,150),m=c(0,0,1),lambda=c(2,2,8))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(100,100,200),m=c(0,0,1),lambda=c(2,2,1))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(100,100,200),m=c(0,0,1),lambda=c(2,2,2))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(100,100,200),m=c(0,0,1),lambda=c(2,2,4))
get_simu(k=3,distName=rep("doublex_SP",3),n=c(100,100,200),m=c(0,0,1),lambda=c(2,2,8))

### chi square and normal right-skewed distributions

get_simu(k=2,distName=c("chi2","skewpos"),n=c(20,10),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(20,10),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(20,10),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(20,10),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(20,20),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(20,20),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(20,20),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(20,20),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(20,30),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(20,30),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(20,30),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(20,30),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(20,40),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(20,40),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(20,40),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(20,40),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(30,15),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(30,15),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(30,15),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(30,15),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(30,30),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(30,30),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(30,30),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(30,30),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(30,45),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(30,45),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(30,45),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(30,45),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(30,60),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(30,60),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(30,60),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(30,60),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(40,20),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(40,20),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(40,20),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(40,20),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(40,40),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(40,40),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(40,40),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(40,40),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(40,60),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(40,60),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(40,60),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(40,60),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(40,80),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(40,80),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(40,80),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(40,80),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(50,25),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(50,25),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(50,25),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(50,25),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(50,50),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(50,50),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(50,50),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(50,50),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(50,75),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(50,75),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(50,75),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(50,75),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(50,100),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(50,100),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(50,100),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(50,100),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(100,50),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(100,50),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(100,50),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(100,50),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(100,100),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(100,100),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(100,100),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(100,100),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(100,150),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(100,150),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(100,150),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(100,150),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(100,200),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(100,200),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(100,200),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewpos"),n=c(100,200),df=c(2,2),m=rep(3,2),sds=c(2,8))

get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(20,20,10),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(20,20,10),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(20,20,10),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(20,20,10),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(20,20,20),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(20,20,20),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(20,20,20),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(20,20,20),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(20,20,30),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(20,20,30),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(20,20,30),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(20,20,30),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(20,20,40),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(20,20,40),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(20,20,40),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(20,20,40),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(30,30,15),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(30,30,15),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(30,30,15),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(30,30,15),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(30,30,30),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(30,30,30),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(30,30,30),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(30,30,30),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(30,30,45),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(30,30,45),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(30,30,45),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(30,30,45),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(30,30,60),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(30,30,60),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(30,30,60),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(30,30,60),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(40,40,20),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(40,40,20),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(40,40,20),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(40,40,20),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(40,40,40),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(40,40,40),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(40,40,40),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(40,40,40),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(40,40,60),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(40,40,60),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(40,40,60),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(40,40,60),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(40,40,80),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(40,40,80),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(40,40,80),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(40,40,80),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(50,50,25),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(50,50,25),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(50,50,25),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(50,50,25),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(50,50,50),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(50,50,50),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(50,50,50),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(50,50,50),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(50,50,75),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(50,50,75),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(50,50,75),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(50,50,75),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(50,50,100),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(50,50,100),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(50,50,100),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(50,50,100),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(100,100,50),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(100,100,50),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(100,100,50),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(100,100,50),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(100,100,100),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(100,100,100),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(100,100,100),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(100,100,100),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(100,100,150),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(100,100,150),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(100,100,150),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(100,100,150),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(100,100,200),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(100,100,200),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(100,100,200),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewpos"),n=c(100,100,200),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))

### chi square and normal left-skewed distributions

get_simu(k=2,distName=c("chi2","skewneg"),n=c(20,10),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(20,10),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(20,10),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(20,10),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(20,20),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(20,20),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(20,20),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(20,20),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(20,30),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(20,30),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(20,30),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(20,30),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(20,40),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(20,40),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(20,40),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(20,40),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(30,15),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(30,15),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(30,15),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(30,15),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(30,30),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(30,30),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(30,30),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(30,30),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(30,45),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(30,45),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(30,45),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(30,45),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(30,60),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(30,60),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(30,60),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(30,60),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(40,20),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(40,20),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(40,20),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(40,20),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(40,40),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(40,40),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(40,40),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(40,40),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(40,60),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(40,60),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(40,60),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(40,60),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(40,80),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(40,80),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(40,80),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(40,80),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(50,25),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(50,25),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(50,25),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(50,25),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(50,50),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(50,50),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(50,50),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(50,50),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(50,75),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(50,75),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(50,75),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(50,75),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(50,100),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(50,100),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(50,100),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(50,100),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(100,50),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(100,50),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(100,50),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(100,50),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(100,100),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(100,100),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(100,100),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(100,100),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(100,150),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(100,150),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(100,150),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(100,150),df=c(2,2),m=rep(3,2),sds=c(2,8))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(100,200),df=c(2,2),m=rep(3,2),sds=c(2,1))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(100,200),df=c(2,2),m=rep(3,2),sds=c(2,2))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(100,200),df=c(2,2),m=rep(3,2),sds=c(2,4))
get_simu(k=2,distName=c("chi2","skewneg"),n=c(100,200),df=c(2,2),m=rep(3,2),sds=c(2,8))

get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(20,20,10),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(20,20,10),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(20,20,10),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(20,20,10),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(20,20,20),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(20,20,20),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(20,20,20),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(20,20,20),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(20,20,30),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(20,20,30),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(20,20,30),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(20,20,30),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(20,20,40),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(20,20,40),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(20,20,40),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(20,20,40),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(30,30,15),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(30,30,15),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(30,30,15),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(30,30,15),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(30,30,30),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(30,30,30),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(30,30,30),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(30,30,30),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(30,30,45),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(30,30,45),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(30,30,45),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(30,30,45),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(30,30,60),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(30,30,60),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(30,30,60),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(30,30,60),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(40,40,20),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(40,40,20),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(40,40,20),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(40,40,20),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(40,40,40),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(40,40,40),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(40,40,40),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(40,40,40),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(40,40,60),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(40,40,60),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(40,40,60),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(40,40,60),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(40,40,80),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(40,40,80),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(40,40,80),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(40,40,80),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(50,50,25),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(50,50,25),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(50,50,25),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(50,50,25),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(50,50,50),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(50,50,50),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(50,50,50),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(50,50,50),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(50,50,75),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(50,50,75),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(50,50,75),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(50,50,75),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(50,50,100),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(50,50,100),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(50,50,100),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(50,50,100),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(100,100,50),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(100,100,50),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(100,100,50),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(100,100,50),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(100,100,100),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(100,100,100),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(100,100,100),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(100,100,100),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(100,100,150),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(100,100,150),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(100,100,150),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(100,100,150),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(100,100,200),df=c(2,2,2),m=rep(3,3),sds=c(2,2,1))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(100,100,200),df=c(2,2,2),m=rep(3,3),sds=c(2,2,2))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(100,100,200),df=c(2,2,2),m=rep(3,3),sds=c(2,2,4))
get_simu(k=3,distName=c("chi2","chi2","skewneg"),n=c(100,100,200),df=c(2,2,2),m=rep(3,3),sds=c(2,2,8))


### uniform distributions

get_simu(k=3,distName=rep("unif",3),n=c(20,20,10),min=c(-3.465,-3.465,-0.7325),max=c(3.465,3.465,2.7325))
get_simu(k=3,distName=rep("unif",3),n=c(20,20,10),min=c(-3.465,-3.465,-2.465),max=c(3.465,3.465,4.465))
get_simu(k=3,distName=rep("unif",3),n=c(20,20,10),min=c(-3.465,-3.465,-5.93),max=c(3.465,3.465,7.93))
get_simu(k=3,distName=rep("unif",3),n=c(20,20,10),min=c(-3.465,-3.465,-12.86),max=c(3.465,3.465,14.86))
get_simu(k=3,distName=rep("unif",3),n=c(20,20,20),min=c(-3.465,-3.465,-0.7325),max=c(3.465,3.465,2.7325))
get_simu(k=3,distName=rep("unif",3),n=c(20,20,20),min=c(-3.465,-3.465,-2.465),max=c(3.465,3.465,4.465))
get_simu(k=3,distName=rep("unif",3),n=c(20,20,20),min=c(-3.465,-3.465,-5.93),max=c(3.465,3.465,7.93))
get_simu(k=3,distName=rep("unif",3),n=c(20,20,20),min=c(-3.465,-3.465,-12.86),max=c(3.465,3.465,14.86))
get_simu(k=3,distName=rep("unif",3),n=c(20,20,30),min=c(-3.465,-3.465,-0.7325),max=c(3.465,3.465,2.7325))
get_simu(k=3,distName=rep("unif",3),n=c(20,20,30),min=c(-3.465,-3.465,-2.465),max=c(3.465,3.465,4.465))
get_simu(k=3,distName=rep("unif",3),n=c(20,20,30),min=c(-3.465,-3.465,-5.93),max=c(3.465,3.465,7.93))
get_simu(k=3,distName=rep("unif",3),n=c(20,20,30),min=c(-3.465,-3.465,-12.86),max=c(3.465,3.465,14.86))
get_simu(k=3,distName=rep("unif",3),n=c(20,20,40),min=c(-3.465,-3.465,-0.7325),max=c(3.465,3.465,2.7325))
get_simu(k=3,distName=rep("unif",3),n=c(20,20,40),min=c(-3.465,-3.465,-2.465),max=c(3.465,3.465,4.465))
get_simu(k=3,distName=rep("unif",3),n=c(20,20,40),min=c(-3.465,-3.465,-5.93),max=c(3.465,3.465,7.93))
get_simu(k=3,distName=rep("unif",3),n=c(20,20,40),min=c(-3.465,-3.465,-12.86),max=c(3.465,3.465,14.86))
get_simu(k=3,distName=rep("unif",3),n=c(30,30,15),min=c(-3.465,-3.465,-0.7325),max=c(3.465,3.465,2.7325))
get_simu(k=3,distName=rep("unif",3),n=c(30,30,15),min=c(-3.465,-3.465,-2.465),max=c(3.465,3.465,4.465))
get_simu(k=3,distName=rep("unif",3),n=c(30,30,15),min=c(-3.465,-3.465,-5.93),max=c(3.465,3.465,7.93))
get_simu(k=3,distName=rep("unif",3),n=c(30,30,15),min=c(-3.465,-3.465,-12.86),max=c(3.465,3.465,14.86))
get_simu(k=3,distName=rep("unif",3),n=c(30,30,30),min=c(-3.465,-3.465,-0.7325),max=c(3.465,3.465,2.7325))
get_simu(k=3,distName=rep("unif",3),n=c(30,30,30),min=c(-3.465,-3.465,-2.465),max=c(3.465,3.465,4.465))
get_simu(k=3,distName=rep("unif",3),n=c(30,30,30),min=c(-3.465,-3.465,-5.93),max=c(3.465,3.465,7.93))
get_simu(k=3,distName=rep("unif",3),n=c(30,30,30),min=c(-3.465,-3.465,-12.86),max=c(3.465,3.465,14.86))
get_simu(k=3,distName=rep("unif",3),n=c(30,30,45),min=c(-3.465,-3.465,-0.7325),max=c(3.465,3.465,2.7325))
get_simu(k=3,distName=rep("unif",3),n=c(30,30,45),min=c(-3.465,-3.465,-2.465),max=c(3.465,3.465,4.465))
get_simu(k=3,distName=rep("unif",3),n=c(30,30,45),min=c(-3.465,-3.465,-5.93),max=c(3.465,3.465,7.93))
get_simu(k=3,distName=rep("unif",3),n=c(30,30,45),min=c(-3.465,-3.465,-12.86),max=c(3.465,3.465,14.86))
get_simu(k=3,distName=rep("unif",3),n=c(30,30,60),min=c(-3.465,-3.465,-0.7325),max=c(3.465,3.465,2.7325))
get_simu(k=3,distName=rep("unif",3),n=c(30,30,60),min=c(-3.465,-3.465,-2.465),max=c(3.465,3.465,4.465))
get_simu(k=3,distName=rep("unif",3),n=c(30,30,60),min=c(-3.465,-3.465,-5.93),max=c(3.465,3.465,7.93))
get_simu(k=3,distName=rep("unif",3),n=c(30,30,60),min=c(-3.465,-3.465,-12.86),max=c(3.465,3.465,14.86))
get_simu(k=3,distName=rep("unif",3),n=c(40,40,20),min=c(-3.465,-3.465,-0.7325),max=c(3.465,3.465,2.7325))
get_simu(k=3,distName=rep("unif",3),n=c(40,40,20),min=c(-3.465,-3.465,-2.465),max=c(3.465,3.465,4.465))
get_simu(k=3,distName=rep("unif",3),n=c(40,40,20),min=c(-3.465,-3.465,-5.93),max=c(3.465,3.465,7.93))
get_simu(k=3,distName=rep("unif",3),n=c(40,40,20),min=c(-3.465,-3.465,-12.86),max=c(3.465,3.465,14.86))
get_simu(k=3,distName=rep("unif",3),n=c(40,40,40),min=c(-3.465,-3.465,-0.7325),max=c(3.465,3.465,2.7325))
get_simu(k=3,distName=rep("unif",3),n=c(40,40,40),min=c(-3.465,-3.465,-2.465),max=c(3.465,3.465,4.465))
get_simu(k=3,distName=rep("unif",3),n=c(40,40,40),min=c(-3.465,-3.465,-5.93),max=c(3.465,3.465,7.93))
get_simu(k=3,distName=rep("unif",3),n=c(40,40,40),min=c(-3.465,-3.465,-12.86),max=c(3.465,3.465,14.86))
get_simu(k=3,distName=rep("unif",3),n=c(40,40,60),min=c(-3.465,-3.465,-0.7325),max=c(3.465,3.465,2.7325))
get_simu(k=3,distName=rep("unif",3),n=c(40,40,60),min=c(-3.465,-3.465,-2.465),max=c(3.465,3.465,4.465))
get_simu(k=3,distName=rep("unif",3),n=c(40,40,60),min=c(-3.465,-3.465,-5.93),max=c(3.465,3.465,7.93))
get_simu(k=3,distName=rep("unif",3),n=c(40,40,60),min=c(-3.465,-3.465,-12.86),max=c(3.465,3.465,14.86))
get_simu(k=3,distName=rep("unif",3),n=c(40,40,80),min=c(-3.465,-3.465,-0.7325),max=c(3.465,3.465,2.7325))
get_simu(k=3,distName=rep("unif",3),n=c(40,40,80),min=c(-3.465,-3.465,-2.465),max=c(3.465,3.465,4.465))
get_simu(k=3,distName=rep("unif",3),n=c(40,40,80),min=c(-3.465,-3.465,-5.93),max=c(3.465,3.465,7.93))
get_simu(k=3,distName=rep("unif",3),n=c(40,40,80),min=c(-3.465,-3.465,-12.86),max=c(3.465,3.465,14.86))
get_simu(k=3,distName=rep("unif",3),n=c(50,50,25),min=c(-3.465,-3.465,-0.7325),max=c(3.465,3.465,2.7325))
get_simu(k=3,distName=rep("unif",3),n=c(50,50,25),min=c(-3.465,-3.465,-2.465),max=c(3.465,3.465,4.465))
get_simu(k=3,distName=rep("unif",3),n=c(50,50,25),min=c(-3.465,-3.465,-5.93),max=c(3.465,3.465,7.93))
get_simu(k=3,distName=rep("unif",3),n=c(50,50,25),min=c(-3.465,-3.465,-12.86),max=c(3.465,3.465,14.86))
get_simu(k=3,distName=rep("unif",3),n=c(50,50,50),min=c(-3.465,-3.465,-0.7325),max=c(3.465,3.465,2.7325))
get_simu(k=3,distName=rep("unif",3),n=c(50,50,50),min=c(-3.465,-3.465,-2.465),max=c(3.465,3.465,4.465))
get_simu(k=3,distName=rep("unif",3),n=c(50,50,50),min=c(-3.465,-3.465,-5.93),max=c(3.465,3.465,7.93))
get_simu(k=3,distName=rep("unif",3),n=c(50,50,50),min=c(-3.465,-3.465,-12.86),max=c(3.465,3.465,14.86))
get_simu(k=3,distName=rep("unif",3),n=c(50,50,75),min=c(-3.465,-3.465,-0.7325),max=c(3.465,3.465,2.7325))
get_simu(k=3,distName=rep("unif",3),n=c(50,50,75),min=c(-3.465,-3.465,-2.465),max=c(3.465,3.465,4.465))
get_simu(k=3,distName=rep("unif",3),n=c(50,50,75),min=c(-3.465,-3.465,-5.93),max=c(3.465,3.465,7.93))
get_simu(k=3,distName=rep("unif",3),n=c(50,50,75),min=c(-3.465,-3.465,-12.86),max=c(3.465,3.465,14.86))
get_simu(k=3,distName=rep("unif",3),n=c(50,50,100),min=c(-3.465,-3.465,-0.7325),max=c(3.465,3.465,2.7325))
get_simu(k=3,distName=rep("unif",3),n=c(50,50,100),min=c(-3.465,-3.465,-2.465),max=c(3.465,3.465,4.465))
get_simu(k=3,distName=rep("unif",3),n=c(50,50,100),min=c(-3.465,-3.465,-5.93),max=c(3.465,3.465,7.93))
get_simu(k=3,distName=rep("unif",3),n=c(50,50,100),min=c(-3.465,-3.465,-12.86),max=c(3.465,3.465,14.86))
get_simu(k=3,distName=rep("unif",3),n=c(100,100,50),min=c(-3.465,-3.465,-0.7325),max=c(3.465,3.465,2.7325))
get_simu(k=3,distName=rep("unif",3),n=c(100,100,50),min=c(-3.465,-3.465,-2.465),max=c(3.465,3.465,4.465))
get_simu(k=3,distName=rep("unif",3),n=c(100,100,50),min=c(-3.465,-3.465,-5.93),max=c(3.465,3.465,7.93))
get_simu(k=3,distName=rep("unif",3),n=c(100,100,50),min=c(-3.465,-3.465,-12.86),max=c(3.465,3.465,14.86))
get_simu(k=3,distName=rep("unif",3),n=c(100,100,100),min=c(-3.465,-3.465,-0.7325),max=c(3.465,3.465,2.7325))
get_simu(k=3,distName=rep("unif",3),n=c(100,100,100),min=c(-3.465,-3.465,-2.465),max=c(3.465,3.465,4.465))
get_simu(k=3,distName=rep("unif",3),n=c(100,100,100),min=c(-3.465,-3.465,-5.93),max=c(3.465,3.465,7.93))
get_simu(k=3,distName=rep("unif",3),n=c(100,100,100),min=c(-3.465,-3.465,-12.86),max=c(3.465,3.465,14.86))
get_simu(k=3,distName=rep("unif",3),n=c(100,100,150),min=c(-3.465,-3.465,-0.7325),max=c(3.465,3.465,2.7325))
get_simu(k=3,distName=rep("unif",3),n=c(100,100,150),min=c(-3.465,-3.465,-2.465),max=c(3.465,3.465,4.465))
get_simu(k=3,distName=rep("unif",3),n=c(100,100,150),min=c(-3.465,-3.465,-5.93),max=c(3.465,3.465,7.93))
get_simu(k=3,distName=rep("unif",3),n=c(100,100,150),min=c(-3.465,-3.465,-12.86),max=c(3.465,3.465,14.86))
get_simu(k=3,distName=rep("unif",3),n=c(100,100,200),min=c(-3.465,-3.465,-0.7325),max=c(3.465,3.465,2.7325))
get_simu(k=3,distName=rep("unif",3),n=c(100,100,200),min=c(-3.465,-3.465,-2.465),max=c(3.465,3.465,4.465))
get_simu(k=3,distName=rep("unif",3),n=c(100,100,200),min=c(-3.465,-3.465,-5.93),max=c(3.465,3.465,7.93))
get_simu(k=3,distName=rep("unif",3),n=c(100,100,200),min=c(-3.465,-3.465,-12.86),max=c(3.465,3.465,14.86))

### normal right-skewed and normal left-skewed distributions

get_simu(k=2,distName=c("skewpos","skewneg"),n=c(20,10),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(20,10),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(20,10),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(20,10),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(20,20),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(20,20),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(20,20),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(20,20),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(20,30),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(20,30),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(20,30),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(20,30),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(20,40),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(20,40),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(20,40),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(20,40),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(30,15),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(30,15),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(30,15),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(30,15),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(30,30),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(30,30),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(30,30),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(30,30),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(30,45),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(30,45),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(30,45),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(30,45),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(30,60),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(30,60),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(30,60),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(30,60),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(40,20),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(40,20),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(40,20),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(40,20),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(40,40),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(40,40),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(40,40),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(40,40),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(40,60),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(40,60),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(40,60),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(40,60),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(40,80),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(40,80),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(40,80),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(40,80),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(50,25),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(50,25),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(50,25),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(50,25),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(50,50),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(50,50),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(50,50),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(50,50),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(50,75),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(50,75),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(50,75),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(50,75),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(50,100),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(50,100),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(50,100),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(50,100),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(100,50),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(100,50),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(100,50),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(100,50),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(100,100),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(100,100),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(100,100),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(100,100),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(100,150),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(100,150),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(100,150),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(100,150),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(100,200),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(100,200),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(100,200),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewneg"),n=c(100,200),m=c(0,1),sds=c(2,8)) # inequal skewness

get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(20,20,10),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(20,20,10),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(20,20,10),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(20,20,10),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(20,20,20),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(20,20,20),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(20,20,20),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(20,20,20),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(20,20,30),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(20,20,30),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(20,20,30),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(20,20,30),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(20,20,40),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(20,20,40),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(20,20,40),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(20,20,40),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(30,30,15),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(30,30,15),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(30,30,15),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(30,30,15),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(30,30,30),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(30,30,30),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(30,30,30),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(30,30,30),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(30,30,45),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(30,30,45),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(30,30,45),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(30,30,45),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(30,30,60),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(30,30,60),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(30,30,60),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(30,30,60),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(40,40,20),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(40,40,20),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(40,40,20),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(40,40,20),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(40,40,40),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(40,40,40),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(40,40,40),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(40,40,40),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(40,40,60),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(40,40,60),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(40,40,60),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(40,40,60),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(40,40,80),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(40,40,80),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(40,40,80),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(40,40,80),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(50,50,25),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(50,50,25),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(50,50,25),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(50,50,25),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(50,50,50),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(50,50,50),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(50,50,50),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(50,50,50),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(50,50,75),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(50,50,75),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(50,50,75),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(50,50,75),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(50,50,100),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(50,50,100),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(50,50,100),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(50,50,100),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(100,100,50),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(100,100,50),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(100,100,50),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(100,100,50),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(100,100,100),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(100,100,100),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(100,100,100),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(100,100,100),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(100,100,150),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(100,100,150),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(100,100,150),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(100,100,150),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(100,100,200),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(100,100,200),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(100,100,200),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewneg"),n=c(100,100,200),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness

### normal right-skewed distributions

get_simu(k=2,distName=c("skewpos","skewpos"),n=c(20,10),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(20,10),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(20,10),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(20,10),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(20,20),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(20,20),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(20,20),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(20,20),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(20,30),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(20,30),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(20,30),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(20,30),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(20,40),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(20,40),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(20,40),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(20,40),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(30,15),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(30,15),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(30,15),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(30,15),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(30,30),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(30,30),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(30,30),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(30,30),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(30,45),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(30,45),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(30,45),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(30,45),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(30,60),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(30,60),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(30,60),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(30,60),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(40,20),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(40,20),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(40,20),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(40,20),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(40,40),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(40,40),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(40,40),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(40,40),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(40,60),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(40,60),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(40,60),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(40,60),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(40,80),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(40,80),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(40,80),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(40,80),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(50,25),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(50,25),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(50,25),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(50,25),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(50,50),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(50,50),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(50,50),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(50,50),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(50,75),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(50,75),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(50,75),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(50,75),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(50,100),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(50,100),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(50,100),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(50,100),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(100,50),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(100,50),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(100,50),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(100,50),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(100,100),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(100,100),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(100,100),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(100,100),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(100,150),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(100,150),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(100,150),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(100,150),m=c(0,1),sds=c(2,8)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(100,200),m=c(0,1),sds=c(2,1)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(100,200),m=c(0,1),sds=c(2,2)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(100,200),m=c(0,1),sds=c(2,4)) # inequal skewness
get_simu(k=2,distName=c("skewpos","skewpos"),n=c(100,200),m=c(0,1),sds=c(2,8)) # inequal skewness

get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(20,20,10),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(20,20,10),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(20,20,10),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(20,20,10),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(20,20,20),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(20,20,20),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(20,20,20),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(20,20,20),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(20,20,30),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(20,20,30),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(20,20,30),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(20,20,30),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(20,20,40),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(20,20,40),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(20,20,40),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(20,20,40),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(30,30,15),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(30,30,15),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(30,30,15),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(30,30,15),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(30,30,30),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(30,30,30),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(30,30,30),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(30,30,30),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(30,30,45),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(30,30,45),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(30,30,45),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(30,30,45),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(30,30,60),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(30,30,60),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(30,30,60),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(30,30,60),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(40,40,20),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(40,40,20),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(40,40,20),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(40,40,20),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(40,40,40),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(40,40,40),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(40,40,40),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(40,40,40),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(40,40,60),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(40,40,60),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(40,40,60),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(40,40,60),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(40,40,80),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(40,40,80),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(40,40,80),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(40,40,80),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(50,50,25),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(50,50,25),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(50,50,25),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(50,50,25),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(50,50,50),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(50,50,50),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(50,50,50),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(50,50,50),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(50,50,75),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(50,50,75),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(50,50,75),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(50,50,75),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(50,50,100),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(50,50,100),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(50,50,100),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(50,50,100),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(100,100,50),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(100,100,50),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(100,100,50),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(100,100,50),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(100,100,100),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(100,100,100),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(100,100,100),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(100,100,100),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(100,100,150),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(100,100,150),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(100,100,150),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(100,100,150),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(100,100,200),m=c(0,0,1),sds=c(2,2,1)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(100,100,200),m=c(0,0,1),sds=c(2,2,2)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(100,100,200),m=c(0,0,1),sds=c(2,2,4)) # inequal skewness
get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(100,100,200),m=c(0,0,1),sds=c(2,2,8)) # inequal skewness
