##### Function to load required packages

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

# example: get_simu(k=3,distName=c("skewpos","skewpos","skewpos"),n=c(100,100,200),m=c(0,0,1),sds=c(2,2,8)) 
# will generate 1000000 datasets where there are 3 samples extracted frop normal right-skewed distribution
# n1 and n2 = 100, n3 = 200
# SD-ratio = 4

