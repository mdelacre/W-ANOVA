library(stringr)

power=function(alpha){

  # Defining all subsections containing RDS files
  Mainfolder="F:/Welch's W ANOVA/ANOVA's Welch/Outputs of simulations/Data files/All simulations/With effect/"
  subfolder1=list.files(Mainfolder)
  subfolder2=list.files(paste0(Mainfolder,subfolder1[1]))
  subfolders=expand.grid(subfolder1,subfolder2)
  subsection=paste0(subfolders[,1],"/",subfolders[,2])
  
  for (j in 1:length(subsection)){ 
  target_section=paste0(Mainfolder,subsection[j])
  setwd(target_section)
  files=list.files()# listing all files in the subsection

  # organize files in the desired order
  k=length(as.numeric(str_extract_all(files[1], "[0-9]+")[[1]]))/3
  filenames_n=matrix(0,length(files),k)
  for (a in 1:length(files)){
    filenames_n[a,]=as.numeric(str_extract_all(files[a], "[0-9]+")[[1]])[1:k]
  }
  sorting=filenames_n[order(filenames_n[,1],filenames_n[,ncol(filenames_n)]),]
  ordre=NULL 

  for (b in 1:(length(files)/4)){
    liste=files[grep(paste0(sorting[(4*b),(k-1)],",",sorting[(4*b),k],"]"),files)]
    ordre=c(ordre,liste)
    }
    files=ordre

    # Computing elements of power for each file 
    results=data.frame(matrix(0,length(files),9))
    colnames(results)=c("Distributions","Nb_groups","Condition","F_power","W_power","BF_power","F_power_expected","W_power_expected","BF_power_expected")

    for (i in 1:length(files)){ 
    A=readRDS(files[i])

    # Extract informations about the condition in the name file
    Distributions=sub("Distr=", "\\",gsub("\\..*","",gsub(",n=", ".", files[i])))    

    filenames_numbers=as.numeric(str_extract_all(files[i], "[0-9]+")[[1]])
    k=length(filenames_numbers)/3 # because all file names contain 3k numbers
    
    Condition=i
    
    colnames(results)=c("Distributions","Nb_groups","Condition","F_power","W_power","BF_power","F_power_expected","W_power_expected","BF_power_expected")
    
    # Compute the observed power of each condition (= proportion of pvalues under alpha, for each test)
    power_F = sum(A[,1]<.05)/length(A[,1]) # proportion of pvalues under alpha, for the F-test
    power_W = sum(A[,2]<.05)/length(A[,2]) # proportion of pvalues under alpha, for the W-test
    power_BF = sum(A[,3]<.05)/length(A[,3])# proportion of pvalues under alpha, for the F*-test
    
    # Compute the expected power
    
    # All input parameters are specified in the file names
    # always in the same order: All n, All means, All sds
    # in other words, file names always contain 3*k numbers

    n=filenames_numbers[1:k]
    m=filenames_numbers[(k+1):(2*k)]
    sd=filenames_numbers[(2*k+1):(3*k)]
    
    sd_pooled<-sqrt(sum((n-1)*sd^2)/(sum(n)-k)) #pooled standard deviation
    sd_BF<-sqrt(sum((1-n/sum(n))*sd^2)) #sd of B-F
    mu <- sum(n*m)/sum(n)
    mean_SS <- sum(n*(m-mu)^2)
    mean_sd <- sqrt(sum(n*(m-mu)^2)/sum(n))
    f <- mean_sd/sd_pooled
    
    #Expected power for ANOVA
    df1 <- k - 1 
    df2_anova <- sum(n) - k
    crit_f <- qf(1 - alpha, df1, df2_anova)
    ncp_anova <- mean_SS/(sd_pooled^2) #ncp_anova <- f^2 * sum(n)
    power_F_theo <- pf(crit_f, df1, df2_anova, ncp = ncp_anova, lower.tail = FALSE)
    
    #Expected power for Welch's ANOVA
    df1 <- k - 1 
    w_j <- n/sd^2 
    mu <- sum((w_j*m)/sum(w_j))
    df2_welch <- (k^2-1)/(3*sum((1-w_j/sum(w_j))^2/(n-1)))
    crit_w <- qf(1 - alpha, df1, df2_welch)
    ncp_w <- sum(w_j*(m-mu)^2)
    power_W_theo <- pf(crit_w, df1, df2_welch, ncp = ncp_w, lower.tail = FALSE)
    
    #Expected power for Brown-Forsythe's test of comparison of means
    df1 <- k - 1 
    df2_BF <- 1/sum((((1-n/sum(n))*sd^2)/sum((1-n/sum(n))*sd^2))^2/(n-1))
    crit_BF <- qf(1 - alpha, df1, df2_BF)
    ncp_BF <- mean_SS/(sd_BF^2)*(k-1) 
    power_BF_theo <- pf(crit_BF, df1, df2_BF, ncp = ncp_BF, lower.tail = FALSE)
    results[i,]=c(Distributions, k, Condition,power_F,power_W,power_BF,power_F_theo,power_W_theo,power_BF_theo)
  
    #rm(list=setdiff(ls(), list("Mainfolder","subsection","results","files")))
    rm(A)
      }

setwd("F:/Welch's W ANOVA/ANOVA's Welch/Outputs of simulations/statistics_power and type 1 error rate")
write.table(results,paste0(results[1,1]," when k= ",results[1,2]," and alpha= ",alpha,".txt"),sep=";",dec=",")  
rm(results)
    }
}

power(alpha=.05)
