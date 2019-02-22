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

  files_path="G:/Welch's W ANOVA/ANOVA's Welch/Outputs of simulations/statistics_power and type 1 error rate/Power/"
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
  }
  
for (i in 1:length(files)){

  File=read.table(files[i],sep=";",dec=".")

# Computing [O-E]/E in order to check consistency between observed and expected power (for F, F* and W tests)

  OE_deviation_F=(as.numeric(File[,4])-as.numeric(File[,7]))/as.numeric(File[,7])
  OE_deviation_W=(as.numeric(File[,5])-as.numeric(File[,8]))/as.numeric(File[,8])
  OE_deviation_BF=(as.numeric(File[,6])-as.numeric(File[,9]))/as.numeric(File[,9])
  File[,10]=OE_deviation_F
  File[,11]=OE_deviation_W
  File[,12]=OE_deviation_BF

# Subdivise datasets into 9 categories
# Depending on sd-ratio, sample sizes ratio, 
# On the fact that there is a positive, negative or no correlation between n and sd
# And on the fact that the last group (where mean = 1) has biggest sd/sample size that other groups (where mean = 0) or not
# see codes line 5:14 for subcategories)
  
  # Conditions id 
  id_Hom_bal=c(6,22,38,54,70)
  id_Hom_nN=c(10,14,26,30,42,46,58,62,74,78)
  id_Hom_Nn=c(2,18,34,50,66)
  id_sdSD_bal=c(7,8,23,24,39,40,55,56,71,72)
  id_SDsd_bal=c(5,21,37,53,69)
  id_sdSD_nN=c(11,12,15,16,27,28,31,32,43,44,47,48,59,60,63,64,75,76,79,80)
  id_SDsd_Nn=c(1,17,33,49,65)
  id_SDsd_nN=c(9,13,25,29,41,45,57,61,73,77)
  id_sdSD_Nn=c(3,4,19,20,35,36,51,52,67,68)
  
  Hom_bal=apply(File[id_Hom_bal,c(4:6,10:12)],2,mean)
  Hom_nN=apply(File[id_Hom_nN,c(4:6,10:12)],2,mean)
  Hom_Nn=apply(File[id_Hom_Nn,c(4:6,10:12)],2,mean)
  sdSD_bal=apply(File[id_sdSD_bal,c(4:6,10:12)],2,mean)
  SDsd_bal=apply(File[id_SDsd_bal,c(4:6,10:12)],2,mean)
  sdSD_nN=apply(File[id_sdSD_nN,c(4:6,10:12)],2,mean)
  SDsd_Nn=apply(File[id_SDsd_Nn,c(4:6,10:12)],2,mean)
  SDsd_nN=apply(File[id_SDsd_nN,c(4:6,10:12)],2,mean)
  sdSD_Nn=apply(File[id_sdSD_Nn,c(4:6,10:12)],2,mean)
  
  subcateg=list(Hom_bal=Hom_bal,Hom_nN=Hom_nN,Hom_Nn=Hom_Nn,sdSD_bal=sdSD_bal,SDsd_bal=SDsd_bal,sdSD_nN=sdSD_nN,sdSD_Nn=sdSD_Nn,SDsd_nN=SDsd_nN,SDsd_Nn=SDsd_Nn)
  
  # Print results
  
  power_results=data.frame(matrix(0, ncol = 9, nrow = 9))
  colnames(power_results)=c("Distribution","K","subcategory","power_F","power_W","power_BF","consistency_F","consistency_W","consistency_BF")

  power_results[,1]=File[1,1] # Distributions of simulations
  power_results[,2]=File[1,2] #  Number of groups
  for (j in 1:length(subcateg)){
    power_results[j,3]=ls(subcateg)[j]
    power_results[j,4:9]=round(subcateg[[j]],3)   
     }
  
  RECAP[[i]]=power_results
  
}
  write.table(RECAP,paste0("power RECAP, K=",K,"and alpha=",alpha,".txt"))
  return(RECAP)
}

RECAP_K2_alpha.01=summary(K=2,alpha=.01)
RECAP_K2_alpha.05=summary(K=2,alpha=.05)
RECAP_K2_alpha.1=summary(K=2,alpha=.1)
RECAP_K3_alpha.01=summary(K=3,alpha=.01)
RECAP_K3_alpha.05=summary(K=3,alpha=.05)
RECAP_K3_alpha.1=summary(K=3,alpha=.1)

############################## Graphics ############################## 

# Legend

par(mfrow=c(1,1),oma = c(0, 0, 0, 0))
plot(1:3,RECAP[[1]][1,4:6],bty="n",xaxt="n",,yaxt="n",ylim=c(.62,.67),main="",ylab="",xlab="",pch=19,type="o")
legend("center", legend=c("Chi-square and normal Left-skewed","Chi-square and normal Righ-skewed","Double exponential","Mixed normal","Normal","Normal Right-skewed and Normal Left-skewed","Normal right-skewed"),
       lty=1:7,pch=1:7,cex=1.1)

graphs=function(K,alpha){

if (K==2){
  if(alpha==.01){RECAP=RECAP_K2_alpha.01
  } else if (alpha==.05){RECAP=RECAP_K2_alpha.05
  } else if (alpha==.1){RECAP=RECAP_K2_alpha.1}
} else if (K==3){
  if(alpha==.01){RECAP=RECAP_K3_alpha.01
  } else if (alpha==.05){RECAP=RECAP_K3_alpha.05
  } else if (alpha==.1){RECAP=RECAP_K3_alpha.1}
}

subcategory=sapply(RECAP, '[[', 3)[,1]

index=c("a","b","c","d","e","f","g","h","i")
for (S in 1:length(subcategory)){
  
  par(xpd=FALSE,mar=c(3,4,4,1))  
  pow=NULL
    for (j in 1:length(RECAP)){
    pow=c(pow,RECAP[[j]][S,4:6])
      pow=unlist(pow)
    }
    MIN_Y=min(unlist(pow))-min(unlist(pow))%%.1
    if(MIN_Y>=.7){
      YMIN=.7
      YMAX=1} else{
        YMIN=MIN_Y
        YMAX=YMIN+.3  
      }

    setwd("C:/Users/Administrateur/Desktop/Plots W-test")
    png(file=paste0("Fig2",index[S],", K=",K,", alpha=",alpha,".png"),width=2000,height=1700, units = "px", res = 300)
    par(mfrow=c(1,2),oma = c(0, 0, 3, 0))
    plot(1:3,NULL,bty="n",ylim=c(YMIN,YMAX),xaxt="n",main="Power",cex.main=1.2,xlab="",ylab="averaged power",pch=19,type="o",col="white")
    axis(side=1,1:3,c("F-test","W-test","F*-test"))
    for (j in 1:length(RECAP)){ 
    lines(1:3,RECAP[[j]][S,4:6],bty="n",xaxt="n",pch=j,type="o",lty=j)}

    plot(1:3,NULL,bty="n",ylim=c(-1,2),xaxt="n",main="Consistency",cex.main=1.2,xlab="",ylab="averaged consistency",pch=19,type="o")
    axis(side=1,1:3,c("F-test","W-test","F*-test"))
    for (j in 1:length(RECAP)){ 
      lines(1:3,RECAP[[j]][S,7:9],bty="n",xaxt="n",pch=j,type="o",lty=j)
    }
    abline(h=0,lty=2,lwd=2,col="red")
    dev.off()
    
   } 
}   

graphs(K=2,alpha=.01)
graphs(K=2,alpha=.05)
graphs(K=2,alpha=.1)

graphs(K=3,alpha=.01)
graphs(K=3,alpha=.05)
graphs(K=3,alpha=.1)

getwd()