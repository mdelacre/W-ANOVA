# Codes for subcategories
# Hom vs. sdSD vs. SDsd 
##### Hom = homoscedasticity 
##### sdSD = the last group has the biggest sd 
##### SDsd = the last group has the smallest sd

# bal vs. nN vs. Nn 
##### bal = equal n across groups 
##### nN = the last group has the biggest sample size 
##### Nn = the last group has the smallest sample size

RECAP=NULL

summary=function(K,alpha){
  
  files_path="G:/Welch's W ANOVA/ANOVA's Welch/Outputs of simulations/statistics_power and type 1 error rate/alpha/"
  setwd(files_path)
  allfiles=list.files(files_path)
  
  if (K==2){
    Files=allfiles[grepl("k= 2",allfiles)==TRUE]
    if(alpha==0.05){
      files=Files[grepl("alpha= 0.05",Files)==TRUE]
    } else if (alpha==0.01){
      files=Files[grepl("alpha= 0.01",Files)==TRUE]
    } else if (alpha==0.1){
      files=Files[grepl("alpha= 0.1",Files)==TRUE]}
  } else if (K==3){
    Files=allfiles[grepl("k= 3",allfiles)==TRUE]
    if(alpha==0.05){
      files=Files[grepl("alpha= 0.05",Files)==TRUE]
    } else if (alpha==0.01){
      files=Files[grepl("alpha= 0.01",Files)==TRUE]
    } else if (alpha==0.1){
      files=Files[grepl("alpha= 0.1",Files)==TRUE]}
  }  else if (K==4){
    Files=allfiles[grepl("k= 4",allfiles)==TRUE]
    if(alpha==0.05){
      files=Files[grepl("alpha= 0.05",Files)==TRUE]
    } else if (alpha==0.01){
      files=Files[grepl("alpha= 0.01",Files)==TRUE]
    } else if (alpha==0.1){
      files=Files[grepl("alpha= 0.1",Files)==TRUE]}
  }  else if (K==5){
    Files=allfiles[grepl("k= 5",allfiles)==TRUE]
    if(alpha==0.05){
      files=Files[grepl("alpha= 0.05",Files)==TRUE]
    } else if (alpha==0.01){
      files=Files[grepl("alpha= 0.01",Files)==TRUE]
    } else if (alpha==0.1){
      files=Files[grepl("alpha= 0.1",Files)==TRUE]}
  }
  
  for (i in 1:length(files)){
    
    File=read.table(files[i],sep=";",dec=".")
    
    # Subdivise datasets into 9 categories
    # Depending on sd-ratio, sample sizes ratio, 
    # On the fact that there is a positive, negative or no correlation between n and sd
    # And on the fact that the last group (where mean = 1) has biggest sd/sample size that other groups (where mean = 0) or not
    # see codes line 5:14 for subcategories)
    
    # Conditions id 
    id_Hom_bal=c(6,22,38,54,70) #
    id_Hom_unbal=c(10,14,26,30,42,46,58,62,74,78,2,18,34,50,66)
    id_Het_bal=c(7,8,23,24,39,40,55,56,71,72,5,21,37,53,69)
    id_Het_positivecor=c(11,12,15,16,27,28,31,32,43,44,47,48,59,60,63,64,75,76,79,80,1,17,33,49,65) 
    id_Het_negativecor=c(9,13,25,29,41,45,57,61,73,77,3,4,19,20,35,36,51,52,67,68)

    Hom_bal=apply(File[id_Hom_bal,c(4:6)],2,mean)
    Hom_unbal=apply(File[id_Hom_unbal,c(4:6)],2,mean)
    Het_bal=apply(File[id_Het_bal,c(4:6)],2,mean)
    Het_positivecor=apply(File[id_Het_positivecor,c(4:6)],2,mean)
    Het_negativecor=apply(File[id_Het_negativecor,c(4:6)],2,mean)

    subcateg=list(Het_bal=Het_bal,Het_negativecor=Het_negativecor,Het_positivecor=Het_positivecor,Hom_bal=Hom_bal,Hom_unbal=Hom_unbal)

    Hom_bal_sd=apply(File[id_Hom_bal,c(4:6)],2,sd)
    Hom_unbal_sd=apply(File[id_Hom_unbal,c(4:6)],2,sd)
    Het_bal_sd=apply(File[id_Het_bal,c(4:6)],2,sd)
    Het_positivecor_sd=apply(File[id_Het_positivecor,c(4:6)],2,sd)
    Het_negativecor_sd=apply(File[id_Het_negativecor,c(4:6)],2,sd)
    
    subcateg_sd=list(Het_bal_sd=Het_bal_sd,Het_negativecor_sd=Het_negativecor_sd,Het_positivecor_sd=Het_positivecor_sd,Hom_bal_sd=Hom_bal_sd,Hom_unbal_sd=Hom_unbal_sd)
    
    # Print results
    
    alpha_results=data.frame(matrix(0, ncol = 9, nrow = length(subcateg)))
    colnames(alpha_results)=c("Distribution","K","subcategory","alpha_F","alpha_W","alpha_BF","stderror_F","stderror_W","stderror_BF")
    
    alpha_results[,1]=File[1,1] # Distributions of simulations
    alpha_results[,2]=File[1,2] #  Number of groups
    for (j in 1:length(subcateg)){
      alpha_results[j,3]=ls(subcateg,sorted=TRUE)[j]
      alpha_results[j,4:6]=round(subcateg[[j]],3)  
      alpha_results[j,7:9]=round(subcateg_sd[[j]],3)
    }
    
    RECAP[[i]]=alpha_results
    
  }
  return(RECAP)
  cat(capture.output(print(RECAP), file=paste0("alpha RECAP, K=",K,".txt")))
}

RECAP_K2_alpha.01=summary(K=2,alpha=.01)
RECAP_K2_alpha.05=summary(K=2,alpha=.05)
RECAP_K2_alpha.10=summary(K=2,alpha=.1)

RECAP_K3_alpha.01=summary(K=3,alpha=.01)
RECAP_K3_alpha.05=summary(K=3,alpha=.05)
RECAP_K3_alpha.10=summary(K=3,alpha=.1)

RECAP_K4_alpha.01=summary(K=4,alpha=.01)
RECAP_K4_alpha.05=summary(K=4,alpha=.05)
RECAP_K4_alpha.10=summary(K=4,alpha=.1)

RECAP_K5_alpha.01=summary(K=5,alpha=.01)
RECAP_K5_alpha.05=summary(K=5,alpha=.05)
RECAP_K5_alpha.10=summary(K=5,alpha=.1)

############################## Graphics ############################## 

# Legend

plot(1:3,RECAP[[1]][1,4:6],bty="n",xaxt="n",yaxt="n",ylim=c(.62,.67),main="",xlab="",ylab="",pch=19,type="o")
legend("center", legend=c("Chi-square and normal Left-skewed","Chi-square and normal Righ-skewed","Double exponential","Mixed normal","Normal","Normal Right-skewed and Normal Left-skewed","Normal right-skewed"),
       lty=1:7,pch=1:7,cex=1.1)

graphs=function(K,yliminf,ylimsup,alpha){
  
  if (K==2){
    if(alpha==.01){RECAP=RECAP_K2_alpha.01
    } else if (alpha==.05){RECAP=RECAP_K2_alpha.05
    } else if (alpha==.1){RECAP=RECAP_K2_alpha.10}
  } else if (K==3){
    if(alpha==.01){RECAP=RECAP_K3_alpha.01
    } else if (alpha==.05){RECAP=RECAP_K3_alpha.05
    } else if (alpha==.1){RECAP=RECAP_K3_alpha.10}
  } else if (K==4){
    if(alpha==.01){RECAP=RECAP_K4_alpha.01
    } else if (alpha==.05){RECAP=RECAP_K4_alpha.05
    } else if (alpha==.1){RECAP=RECAP_K4_alpha.10}
  } else if (K==5){
    if(alpha==.01){RECAP=RECAP_K5_alpha.01
    } else if (alpha==.05){RECAP=RECAP_K5_alpha.05
    } else if (alpha==.1){RECAP=RECAP_K5_alpha.10}
  }
  
  subcategory=sapply(RECAP, '[[', 3)[,1]
  
  index=c("C","E","D","A","B")
  for (S in 1:length(subcategory)){
    setwd("C:/Users/Administrateur/Desktop/Plots W-test")
    png(file=paste0("Fig1",index[S],", K=",K,", alpha=",alpha,".png"),width=2000,height=1700, units = "px", res = 300)  
    par(xpd=FALSE,mar=c(3,4,4,1))  
    plot(1:3,NULL,bty="n",ylim=c(yliminf,ylimsup),xaxt="n",xlab="",ylab="averaged Type I error rate",pch=7,type="o")
    axis(side=1,1:3,c("F-test","W-test","F*-test"))
    for (j in 1:length(RECAP)){ 
      lines(1:3,RECAP[[j]][S,4:6],bty="n",xaxt="n",main="Averaged alpha of 3 tests when n and sd are equal across groups",pch=j,type="o",lty=j)}
    abline(h=alpha,lty=2,lwd=2,col="red")
    
    rect(.5,.5*alpha,3.5,1.5*alpha, col= rgb(0, 0, 0, alpha=.05),border=NA)
    rect(.5,.9*alpha,3.5,1.1*alpha, col= rgb(0, 0, 0, alpha=.25),border=NA)
    dev.off()
  }  
} 

getwd()
graphs(K=3,.02,.12,alpha=.05)

(K,yliminf,ylimsup,alpha)


graphs(K=4,.02,.15,alpha=.05)
graphs(K=5,.02,.15,alpha=.05)
