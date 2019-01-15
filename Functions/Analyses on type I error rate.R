files_path="G:/Welch's W ANOVA/ANOVA's Welch/Outputs of simulations/statistics_power and type 1 error rate/alpha/"
setwd(files_path)
files=list.files(files_path)

files_k2=files[grepl("k= 2",files)==TRUE]
files_k3=files[grepl("k= 3",files)==TRUE] 
files_k4=files[grepl("k= 4",files)==TRUE] 
files_k5=files[grepl("k= 5",files)==TRUE] 

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

summary=function(K){

  if (K==2){files=files_k2
  } else if (K==3){files=files_k3
  } else if (K==4){files=files_k4
  } else if (K==5){files=files_k5}
  
for (i in 1:length(files)){

  File=read.table(files[i],sep=";",dec=".")

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
  
  Hom_bal=apply(File[id_Hom_bal,c(4:6)],2,mean)
  Hom_nN=apply(File[id_Hom_nN,c(4:6)],2,mean)
  Hom_Nn=apply(File[id_Hom_Nn,c(4:6)],2,mean)
  sdSD_bal=apply(File[id_sdSD_bal,c(4:6)],2,mean)
  SDsd_bal=apply(File[id_SDsd_bal,c(4:6)],2,mean)
  sdSD_nN=apply(File[id_sdSD_nN,c(4:6)],2,mean)
  SDsd_Nn=apply(File[id_SDsd_Nn,c(4:6)],2,mean)
  SDsd_nN=apply(File[id_SDsd_nN,c(4:6)],2,mean)
  sdSD_Nn=apply(File[id_sdSD_Nn,c(4:6)],2,mean)
  
  subcateg=list(Hom_bal=Hom_bal,Hom_nN=Hom_nN,Hom_Nn=Hom_Nn,sdSD_bal=sdSD_bal,SDsd_bal=SDsd_bal,sdSD_nN=sdSD_nN,SDsd_Nn=SDsd_Nn,SDsd_nN=SDsd_nN,sdSD_Nn=sdSD_Nn)
  
  # Print results
  
  alpha_results=data.frame(matrix(0, ncol = 6, nrow = 9))
  colnames(alpha_results)=c("Distribution","K","subcategory","alpha_F","alpha_W","alpha_BF")

  alpha_results[,1]=File[1,1] # Distributions of simulations
  alpha_results[,2]=File[1,2] #  Number of groups
  for (j in 1:length(subcateg)){
    alpha_results[j,3]=ls(subcateg)[j]
    alpha_results[j,4:6]=round(subcateg[[j]],3)  
     }
  
  RECAP[[i]]=alpha_results
  
}
  return(RECAP)
  cat(capture.output(print(RECAP), file=paste0("alpha RECAP, K=",K,".txt")))
}

RECAP_K2=summary(K=2)
RECAP_K3=summary(K=3)
RECAP_K4=summary(K=4)
RECAP_K5=summary(K=5)

############################## Graphics ############################## 

# Legend

plot(1:3,RECAP[[1]][1,4:6],bty="n",xaxt="n",,yaxt="n",ylim=c(.62,.67),main="",xlab="",ylab="averaged alpha",pch=19,type="o")
legend("center", legend=c("Chi-square and normal Left-skewed","Chi-square and normal Righ-skewed","Double exponential","Mixed normal","Normal","Normal Right-skewed and Normal Left-skewed","Normal right-skewed"),
       lty=1:7,pch=1:7,cex=1.1)


graphs=function(K){

if (K==2){RECAP=RECAP_K2
} else if (K==3){RECAP=RECAP_K3
} else if (K==4){RECAP=RECAP_K4
} else if (K==5){RECAP=RECAP_K5}

subcategory=sapply(RECAP, '[[', 3)[,1]

##
## CHANGER LE NOM DES GRAPHIQUES
##

for (S in 1:length(subcategory)){
  
    if (grepl("Hom",subcategory[S])==TRUE){
      categ="equal sd across groups"
        if (grepl("bal",subcategory[S])==TRUE){n="equal n"
        } else if (grepl("nN",subcategory[S])==TRUE){n="the last group has the highest n"
        } else if (grepl("Nn",subcategory[S])==TRUE){n="the last group has the smallest n"}      
    
    } else if (grepl("sdSD",subcategory[S])==TRUE){
        if (grepl("bal",subcategory[S])==TRUE){
          categ="the last group has the highest SD"
          n="equal n"
        } else if (grepl("nN",subcategory[S])==TRUE){
          categ="positive correlation between n and sd"
          n="the last group has the biggest n and sd"
        } else if (grepl("Nn",subcategory[S])==TRUE){
          categ="negative correlation between n and sd"
          n="the last group has the smallest n and the biggest sd"}      
    
    } else if (grepl("SDsd",subcategory[S])==TRUE) {
        if (grepl("bal",subcategory[S])==TRUE){
          categ="the last group has the smallest SD"
          n="equal n"
        } else if (grepl("nN",subcategory[S])==TRUE){
          categ="negative correlation between n and sd"
          n="the last group has the biggest n and the smallest sd"
        } else if (grepl("Nn",subcategory[S])==TRUE){
          categ="positive correlation between n and sd"
          n="the last group has the smallest n and sd"}      
    }

  Title = paste0("Averaged type I error rate of 3 tests: ","\n",categ,"\n",n)

  par(xpd=FALSE,mar=c(3,3,4,1))  

    #png(file = paste0("Observed alpha, condition " ,subcategory[S]," when K=",K,".png"), width = 800, height = 700) 
    plot(1:3,RECAP[[1]][1,4:6],bty="n",ylim=c(0,.20),xaxt="n",main=Title,xlab="",ylab="averaged alpha",pch=7,type="o")
    axis(side=1,1:3,c("F-test","W-test","F*-test"))
    for (j in 1:length(RECAP)){ 
    lines(1:3,RECAP[[j]][S,4:6],bty="n",xaxt="n",main="Averaged alpha of 3 tests when n and sd are equal across groups",pch=j,type="o",lty=j)}
    abline(h=0.05,lty=2,lwd=2,col="red")
    
    #dev.off()
  } 

} 

### SEUL COUAC: PNG ET DEV.OFF, ça ne fonctionne pas. 
graphs(K=2)
graphs(K=3)
graphs(K=4)
graphs(K=5)
 