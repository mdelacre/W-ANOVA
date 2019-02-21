##########################################################
####                 COMPUTING ALPHA                  ####
##########################################################
library(stringr)

typeIerrorrate=function(alpha){
  
  # Defining all subsections containing RDS files
  Mainfolder="G:/Welch's W ANOVA/ANOVA's Welch/Outputs of simulations/Data files/All simulations/Without effect/"
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
    results=data.frame(matrix(0,length(files),6)) # How many conditions?
    colnames(results)=c("Distributions","Nb_groups","Condition","F_alpha","W_alpha","BF_alpha")
    
    for (i in 1:length(files)){ 
      A=readRDS(files[i])
      
      # Extract informations about the condition in the name file
      Distributions=sub("Distr=", "\\",gsub("\\..*","",gsub(",n=", ".", files[i]))) 
      
      filenames_numbers=as.numeric(str_extract_all(files[i], "[0-9]+")[[1]])
      k=length(filenames_numbers)/3 # because all file names contain 3k numbers
      
      Condition=i
      
      # Compute the observed power of each condition (= proportion of pvalues under alpha, for each test)
      alpha_F = sum(A[,1]<alpha)/length(A[,1]) # proportion of pvalues under alpha, for the F-test
      alpha_W = sum(A[,2]<alpha)/length(A[,2]) # proportion of pvalues under alpha, for the W-test
      alpha_BF = sum(A[,3]<alpha)/length(A[,3])# proportion of pvalues under alpha, for the F*-test
      
      results[i,]=c(Distributions, k, Condition,alpha_F,alpha_W,alpha_BF)
      
      rm(A)
    }
    
    setwd("G:/Welch's W ANOVA/ANOVA's Welch/Outputs of simulations/statistics_power and type 1 error rate/alpha/")
    write.table(results,paste0(results[1,1]," when k= ",results[1,2]," and alpha= ",alpha,".txt"),sep=";",dec=",")  
    rm(results)
  }
}

typeIerrorrate(alpha=.01)
typeIerrorrate(alpha=.05)
typeIerrorrate(alpha=.10)