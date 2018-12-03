for (package in c("bda","moments","onewaytests","smoothmest","fGarch","stringr")) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}

theo_power=function(m,sd,n,alpha=.05){

    k=length(m)
    sd_pooled<-sqrt(sum((n-1)*sd^2)/(sum(n)-k)) #pooled standard deviation
    sd_BF<-sqrt(sum((1-n/sum(n))*sd^2)) #sd of B-F
    mu <- sum(n*m)/sum(n)
    mean_SS <- sum(n*(m-mu)^2)
    mean_sd <- sqrt(sum(n*(m-mu)^2)/sum(n))
    f <- mean_sd/sd_pooled
  
    #Power for ANOVA
    df1 <- k - 1 
    df2_anova <- sum(n) - k
    crit_f <- qf(1 - alpha, df1, df2_anova)
    ncp_anova <- mean_SS/(sd_pooled^2) #ncp_anova <- f^2 * sum(n)
    power_anova <- pf(crit_f, df1, df2_anova, ncp = ncp_anova, lower.tail = FALSE)

    #Power for Welch's ANOVA
    df1 <- k - 1 
    w_j <- n/sd^2 
    mu <- sum((w_j*m)/sum(w_j))
    df2_welch <- (k^2-1)/(3*sum((1-w_j/sum(w_j))^2/(n-1)))
    crit_w <- qf(1 - alpha, df1, df2_welch)
    ncp_w <- sum(w_j*(m-mu)^2)
    power_ANOVAwelch <- pf(crit_w, df1, df2_welch, ncp = ncp_w, lower.tail = FALSE)

    #Power for Brown-Forsythe's test of comparison of means
    df1 <- k - 1 
    df2_BF <- 1/sum((((1-n/sum(n))*sd^2)/sum((1-n/sum(n))*sd^2))^2/(n-1))
    crit_BF <- qf(1 - alpha, df1, df2_BF)
    ncp_BF <- mean_SS/(sd_BF^2)*(k-1) 
    power_BF <- pf(crit_BF, df1, df2_BF, ncp = ncp_BF, lower.tail = FALSE)

    return(c(round(power_anova,3),round(power_ANOVAwelch,3),round(power_BF,3)))

}

power_k2=rbind(theo_power(m=c(0,1),sd=c(2,1),n=c(20,20,10)),
theo_power(m=c(0,1),sd=c(2,2),n=c(20,20,10)),
theo_power(m=c(0,1),sd=c(2,4),n=c(20,20,10)),
theo_power(m=c(0,1),sd=c(2,8),n=c(20,20,10)),
theo_power(m=c(0,1),sd=c(2,1),n=c(20,20,20)),
theo_power(m=c(0,1),sd=c(2,2),n=c(20,20,20)),
theo_power(m=c(0,1),sd=c(2,4),n=c(20,20,20)),
theo_power(m=c(0,1),sd=c(2,8),n=c(20,20,20)),
theo_power(m=c(0,1),sd=c(2,1),n=c(20,20,30)),
theo_power(m=c(0,1),sd=c(2,2),n=c(20,20,30)),
theo_power(m=c(0,1),sd=c(2,4),n=c(20,20,30)),
theo_power(m=c(0,1),sd=c(2,8),n=c(20,20,30)),
theo_power(m=c(0,1),sd=c(2,1),n=c(20,20,40)),
theo_power(m=c(0,1),sd=c(2,2),n=c(20,20,40)),
theo_power(m=c(0,1),sd=c(2,4),n=c(20,20,40)),
theo_power(m=c(0,1),sd=c(2,8),n=c(20,20,40)),
theo_power(m=c(0,1),sd=c(2,1),n=c(30,30,15)),
theo_power(m=c(0,1),sd=c(2,2),n=c(30,30,15)),
theo_power(m=c(0,1),sd=c(2,4),n=c(30,30,15)),
theo_power(m=c(0,1),sd=c(2,8),n=c(30,30,15)),
theo_power(m=c(0,1),sd=c(2,1),n=c(30,30,30)),
theo_power(m=c(0,1),sd=c(2,2),n=c(30,30,30)),
theo_power(m=c(0,1),sd=c(2,4),n=c(30,30,30)),
theo_power(m=c(0,1),sd=c(2,8),n=c(30,30,30)),
theo_power(m=c(0,1),sd=c(2,1),n=c(30,30,45)),
theo_power(m=c(0,1),sd=c(2,2),n=c(30,30,45)),
theo_power(m=c(0,1),sd=c(2,4),n=c(30,30,45)),
theo_power(m=c(0,1),sd=c(2,8),n=c(30,30,45)),
theo_power(m=c(0,1),sd=c(2,1),n=c(30,30,60)),
theo_power(m=c(0,1),sd=c(2,2),n=c(30,30,60)),
theo_power(m=c(0,1),sd=c(2,4),n=c(30,30,60)),
theo_power(m=c(0,1),sd=c(2,8),n=c(30,30,60)),
theo_power(m=c(0,1),sd=c(2,1),n=c(40,40,20)),
theo_power(m=c(0,1),sd=c(2,2),n=c(40,40,20)),
theo_power(m=c(0,1),sd=c(2,4),n=c(40,40,20)),
theo_power(m=c(0,1),sd=c(2,8),n=c(40,40,20)),
theo_power(m=c(0,1),sd=c(2,1),n=c(40,40,40)),
theo_power(m=c(0,1),sd=c(2,2),n=c(40,40,40)),
theo_power(m=c(0,1),sd=c(2,4),n=c(40,40,40)),
theo_power(m=c(0,1),sd=c(2,8),n=c(40,40,40)),
theo_power(m=c(0,1),sd=c(2,1),n=c(40,40,60)),
theo_power(m=c(0,1),sd=c(2,2),n=c(40,40,60)),
theo_power(m=c(0,1),sd=c(2,4),n=c(40,40,60)),
theo_power(m=c(0,1),sd=c(2,8),n=c(40,40,60)),
theo_power(m=c(0,1),sd=c(2,1),n=c(40,40,80)),
theo_power(m=c(0,1),sd=c(2,2),n=c(40,40,80)),
theo_power(m=c(0,1),sd=c(2,4),n=c(40,40,80)),
theo_power(m=c(0,1),sd=c(2,8),n=c(40,40,80)),
theo_power(m=c(0,1),sd=c(2,1),n=c(50,50,25)),
theo_power(m=c(0,1),sd=c(2,2),n=c(50,50,25)),
theo_power(m=c(0,1),sd=c(2,4),n=c(50,50,25)),
theo_power(m=c(0,1),sd=c(2,8),n=c(50,50,25)),
theo_power(m=c(0,1),sd=c(2,1),n=c(50,50,50)),
theo_power(m=c(0,1),sd=c(2,2),n=c(50,50,50)),
theo_power(m=c(0,1),sd=c(2,4),n=c(50,50,50)),
theo_power(m=c(0,1),sd=c(2,8),n=c(50,50,50)),
theo_power(m=c(0,1),sd=c(2,1),n=c(50,50,75)),
theo_power(m=c(0,1),sd=c(2,2),n=c(50,50,75)),
theo_power(m=c(0,1),sd=c(2,4),n=c(50,50,75)),
theo_power(m=c(0,1),sd=c(2,8),n=c(50,50,75)),
theo_power(m=c(0,1),sd=c(2,1),n=c(50,50,100)),
theo_power(m=c(0,1),sd=c(2,2),n=c(50,50,100)),
theo_power(m=c(0,1),sd=c(2,4),n=c(50,50,100)),
theo_power(m=c(0,1),sd=c(2,8),n=c(50,50,100)),
theo_power(m=c(0,1),sd=c(2,1),n=c(100,50)),
theo_power(m=c(0,1),sd=c(2,2),n=c(100,50)),
theo_power(m=c(0,1),sd=c(2,4),n=c(100,50)),
theo_power(m=c(0,1),sd=c(2,8),n=c(100,50)),
theo_power(m=c(0,1),sd=c(2,1),n=c(100,100)),
theo_power(m=c(0,1),sd=c(2,2),n=c(100,100)),
theo_power(m=c(0,1),sd=c(2,4),n=c(100,100)),
theo_power(m=c(0,1),sd=c(2,8),n=c(100,100)),
theo_power(m=c(0,1),sd=c(2,1),n=c(100,150)),
theo_power(m=c(0,1),sd=c(2,2),n=c(100,150)),
theo_power(m=c(0,1),sd=c(2,4),n=c(100,150)),
theo_power(m=c(0,1),sd=c(2,8),n=c(100,150)),
theo_power(m=c(0,1),sd=c(2,1),n=c(100,200)),
theo_power(m=c(0,1),sd=c(2,2),n=c(100,200)),
theo_power(m=c(0,1),sd=c(2,4),n=c(100,200)),
theo_power(m=c(0,1),sd=c(2,8),n=c(100,200)))

write.table(power_k2,"power_k2.txt",sep=";",dec=",")


power_k3=rbind(theo_power(m=c(0,0,1),sd=c(2,2,1),n=c(20,20,10)),
theo_power(m=c(0,0,1),sd=c(2,2,2),n=c(20,20,10)),
theo_power(m=c(0,0,1),sd=c(2,2,4),n=c(20,20,10)),
theo_power(m=c(0,0,1),sd=c(2,2,8),n=c(20,20,10)),
theo_power(m=c(0,0,1),sd=c(2,2,1),n=c(20,20,20)),
theo_power(m=c(0,0,1),sd=c(2,2,2),n=c(20,20,20)),
theo_power(m=c(0,0,1),sd=c(2,2,4),n=c(20,20,20)),
theo_power(m=c(0,0,1),sd=c(2,2,8),n=c(20,20,20)),
theo_power(m=c(0,0,1),sd=c(2,2,1),n=c(20,20,30)),
theo_power(m=c(0,0,1),sd=c(2,2,2),n=c(20,20,30)),
theo_power(m=c(0,0,1),sd=c(2,2,4),n=c(20,20,30)),
theo_power(m=c(0,0,1),sd=c(2,2,8),n=c(20,20,30)),
theo_power(m=c(0,0,1),sd=c(2,2,1),n=c(20,20,40)),
theo_power(m=c(0,0,1),sd=c(2,2,2),n=c(20,20,40)),
theo_power(m=c(0,0,1),sd=c(2,2,4),n=c(20,20,40)),
theo_power(m=c(0,0,1),sd=c(2,2,8),n=c(20,20,40)),
theo_power(m=c(0,0,1),sd=c(2,2,1),n=c(30,30,15)),
theo_power(m=c(0,0,1),sd=c(2,2,2),n=c(30,30,15)),
theo_power(m=c(0,0,1),sd=c(2,2,4),n=c(30,30,15)),
theo_power(m=c(0,0,1),sd=c(2,2,8),n=c(30,30,15)),
theo_power(m=c(0,0,1),sd=c(2,2,1),n=c(30,30,30)),
theo_power(m=c(0,0,1),sd=c(2,2,2),n=c(30,30,30)),
theo_power(m=c(0,0,1),sd=c(2,2,4),n=c(30,30,30)),
theo_power(m=c(0,0,1),sd=c(2,2,8),n=c(30,30,30)),
theo_power(m=c(0,0,1),sd=c(2,2,1),n=c(30,30,45)),
theo_power(m=c(0,0,1),sd=c(2,2,2),n=c(30,30,45)),
theo_power(m=c(0,0,1),sd=c(2,2,4),n=c(30,30,45)),
theo_power(m=c(0,0,1),sd=c(2,2,8),n=c(30,30,45)),
theo_power(m=c(0,0,1),sd=c(2,2,1),n=c(30,30,60)),
theo_power(m=c(0,0,1),sd=c(2,2,2),n=c(30,30,60)),
theo_power(m=c(0,0,1),sd=c(2,2,4),n=c(30,30,60)),
theo_power(m=c(0,0,1),sd=c(2,2,8),n=c(30,30,60)),
theo_power(m=c(0,0,1),sd=c(2,2,1),n=c(40,40,20)),
theo_power(m=c(0,0,1),sd=c(2,2,2),n=c(40,40,20)),
theo_power(m=c(0,0,1),sd=c(2,2,4),n=c(40,40,20)),
theo_power(m=c(0,0,1),sd=c(2,2,8),n=c(40,40,20)),
theo_power(m=c(0,0,1),sd=c(2,2,1),n=c(40,40,40)),
theo_power(m=c(0,0,1),sd=c(2,2,2),n=c(40,40,40)),
theo_power(m=c(0,0,1),sd=c(2,2,4),n=c(40,40,40)),
theo_power(m=c(0,0,1),sd=c(2,2,8),n=c(40,40,40)),
theo_power(m=c(0,0,1),sd=c(2,2,1),n=c(40,40,60)),
theo_power(m=c(0,0,1),sd=c(2,2,2),n=c(40,40,60)),
theo_power(m=c(0,0,1),sd=c(2,2,4),n=c(40,40,60)),
theo_power(m=c(0,0,1),sd=c(2,2,8),n=c(40,40,60)),
theo_power(m=c(0,0,1),sd=c(2,2,1),n=c(40,40,80)),
theo_power(m=c(0,0,1),sd=c(2,2,2),n=c(40,40,80)),
theo_power(m=c(0,0,1),sd=c(2,2,4),n=c(40,40,80)),
theo_power(m=c(0,0,1),sd=c(2,2,8),n=c(40,40,80)),
theo_power(m=c(0,0,1),sd=c(2,2,1),n=c(50,50,25)),
theo_power(m=c(0,0,1),sd=c(2,2,2),n=c(50,50,25)),
theo_power(m=c(0,0,1),sd=c(2,2,4),n=c(50,50,25)),
theo_power(m=c(0,0,1),sd=c(2,2,8),n=c(50,50,25)),
theo_power(m=c(0,0,1),sd=c(2,2,1),n=c(50,50,50)),
theo_power(m=c(0,0,1),sd=c(2,2,2),n=c(50,50,50)),
theo_power(m=c(0,0,1),sd=c(2,2,4),n=c(50,50,50)),
theo_power(m=c(0,0,1),sd=c(2,2,8),n=c(50,50,50)),
theo_power(m=c(0,0,1),sd=c(2,2,1),n=c(50,50,75)),
theo_power(m=c(0,0,1),sd=c(2,2,2),n=c(50,50,75)),
theo_power(m=c(0,0,1),sd=c(2,2,4),n=c(50,50,75)),
theo_power(m=c(0,0,1),sd=c(2,2,8),n=c(50,50,75)),
theo_power(m=c(0,0,1),sd=c(2,2,1),n=c(50,50,100)),
theo_power(m=c(0,0,1),sd=c(2,2,2),n=c(50,50,100)),
theo_power(m=c(0,0,1),sd=c(2,2,4),n=c(50,50,100)),
theo_power(m=c(0,0,1),sd=c(2,2,8),n=c(50,50,100)),
theo_power(m=c(0,0,1),sd=c(2,2,1),n=c(100,100,50)),
theo_power(m=c(0,0,1),sd=c(2,2,2),n=c(100,100,50)),
theo_power(m=c(0,0,1),sd=c(2,2,4),n=c(100,100,50)),
theo_power(m=c(0,0,1),sd=c(2,2,8),n=c(100,100,50)),
theo_power(m=c(0,0,1),sd=c(2,2,1),n=c(100,100,100)),
theo_power(m=c(0,0,1),sd=c(2,2,2),n=c(100,100,100)),
theo_power(m=c(0,0,1),sd=c(2,2,4),n=c(100,100,100)),
theo_power(m=c(0,0,1),sd=c(2,2,8),n=c(100,100,100)),
theo_power(m=c(0,0,1),sd=c(2,2,1),n=c(100,100,150)),
theo_power(m=c(0,0,1),sd=c(2,2,2),n=c(100,100,150)),
theo_power(m=c(0,0,1),sd=c(2,2,4),n=c(100,100,150)),
theo_power(m=c(0,0,1),sd=c(2,2,8),n=c(100,100,150)),
theo_power(m=c(0,0,1),sd=c(2,2,1),n=c(100,100,200)),
theo_power(m=c(0,0,1),sd=c(2,2,2),n=c(100,100,200)),
theo_power(m=c(0,0,1),sd=c(2,2,4),n=c(100,100,200)),
theo_power(m=c(0,0,1),sd=c(2,2,8),n=c(100,100,200)))

write.table(power_k3,"power_k3.txt",sep=";",dec=",")





theoretical_power=function(k){

p_val<-data.frame(matrix(0,160,6,dimnames=list(1:160,c("one-way ANOVA","ANOVA's Welch test","ANOVA's B-F test","A-G test","Second James Order 5%","K-W test"))))
  # where k = number of groups
  # and p =  number of parameters
  #set up empty container for all estimated parameters

SDR=rep(c(0.5,1,2,4),40)

   for (i in 1:160){

      if(i%%2==0) {
		p_val[i,1]="obs"
		p_val[i,2]="obs"
		p_val[i,3]="obs"

          }

      if(i%%2==1) {
               if(i<=8) {
                     p_val[i,1]=theo_power(m=c(rep(0,k-1),1),n=c(rep(20,k-1),20*0.5),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[1]
                     p_val[i,2]=theo_power(m=c(rep(0,k-1),1),n=c(rep(20,k-1),20*0.5),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[2]
                     p_val[i,3]=theo_power(m=c(rep(0,k-1),1),n=c(rep(20,k-1),20*0.5),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[3]

                         }
               if(i>8 & i <=16) {
                     p_val[i,1]=theo_power(m=c(rep(0,k-1),1),n=c(rep(20,k-1),20*1),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[1]
                     p_val[i,2]=theo_power(m=c(rep(0,k-1),1),n=c(rep(20,k-1),20*1),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[2]
                     p_val[i,3]=theo_power(m=c(rep(0,k-1),1),n=c(rep(20,k-1),20*1),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[3]
                               }
               if(i>16 & i <=24) {
                     p_val[i,1]=theo_power(m=c(rep(0,k-1),1),n=c(rep(20,k-1),20*1.5),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[1]
                     p_val[i,2]=theo_power(m=c(rep(0,k-1),1),n=c(rep(20,k-1),20*1.5),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[2]
                     p_val[i,3]=theo_power(m=c(rep(0,k-1),1),n=c(rep(20,k-1),20*1.5),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[3]
                               }
               if(i>24 & i <=32) {
                     p_val[i,1]=theo_power(m=c(rep(0,k-1),1),n=c(rep(20,k-1),20*2),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[1]
                     p_val[i,2]=theo_power(m=c(rep(0,k-1),1),n=c(rep(20,k-1),20*2),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[2]
                     p_val[i,3]=theo_power(m=c(rep(0,k-1),1),n=c(rep(20,k-1),20*2),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[3]
                               }
               if(i>32 & i <=40) {
                   p_val[i,1]=theo_power(m=c(rep(0,k-1),1),n=c(rep(30,k-1),30*0.5),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[1]
                   p_val[i,2]=theo_power(m=c(rep(0,k-1),1),n=c(rep(30,k-1),30*0.5),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[2]
                   p_val[i,3]=theo_power(m=c(rep(0,k-1),1),n=c(rep(30,k-1),30*0.5),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[3]
                               }
               if(i>40 & i <=48) {
                     p_val[i,1]=theo_power(m=c(rep(0,k-1),1),n=c(rep(30,k-1),30*1),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[1]
                     p_val[i,2]=theo_power(m=c(rep(0,k-1),1),n=c(rep(30,k-1),30*1),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[2]
                     p_val[i,3]=theo_power(m=c(rep(0,k-1),1),n=c(rep(30,k-1),30*1),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[3]
                               }
               if(i>48 & i <=56) {
                     p_val[i,1]=theo_power(m=c(rep(0,k-1),1),n=c(rep(30,k-1),30*1.5),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[1]
                     p_val[i,2]=theo_power(m=c(rep(0,k-1),1),n=c(rep(30,k-1),30*1.5),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[2]
                     p_val[i,3]=theo_power(m=c(rep(0,k-1),1),n=c(rep(30,k-1),30*1.5),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[3]
                               }
               if(i>56 & i <=64) {
                     p_val[i,1]=theo_power(m=c(rep(0,k-1),1),n=c(rep(30,k-1),30*2),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[1]
                     p_val[i,2]=theo_power(m=c(rep(0,k-1),1),n=c(rep(30,k-1),30*2),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[2]
                     p_val[i,3]=theo_power(m=c(rep(0,k-1),1),n=c(rep(30,k-1),30*2),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[3]
                               }
               if(i>64 & i <=72) {
                     p_val[i,1]=theo_power(m=c(rep(0,k-1),1),n=c(rep(40,k-1),40*0.5),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[1]
                     p_val[i,2]=theo_power(m=c(rep(0,k-1),1),n=c(rep(40,k-1),40*0.5),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[2]
                     p_val[i,3]=theo_power(m=c(rep(0,k-1),1),n=c(rep(40,k-1),40*0.5),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[3]
                               }
               if(i>72 & i <=80) {
                     p_val[i,1]=theo_power(m=c(rep(0,k-1),1),n=c(rep(40,k-1),40*1),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[1]
                     p_val[i,2]=theo_power(m=c(rep(0,k-1),1),n=c(rep(40,k-1),40*1),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[2]
                     p_val[i,3]=theo_power(m=c(rep(0,k-1),1),n=c(rep(40,k-1),40*1),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[3]
                               }
               if(i>80 & i <=88) {
                     p_val[i,1]=theo_power(m=c(rep(0,k-1),1),n=c(rep(40,k-1),40*1.5),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[1]
                     p_val[i,2]=theo_power(m=c(rep(0,k-1),1),n=c(rep(40,k-1),40*1.5),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[2]
                     p_val[i,3]=theo_power(m=c(rep(0,k-1),1),n=c(rep(40,k-1),40*1.5),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[3]
                               }
               if(i>88 & i <=96) {
                     p_val[i,1]=theo_power(m=c(rep(0,k-1),1),n=c(rep(40,k-1),40*2),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[1]
                     p_val[i,2]=theo_power(m=c(rep(0,k-1),1),n=c(rep(40,k-1),40*2),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[2]
                     p_val[i,3]=theo_power(m=c(rep(0,k-1),1),n=c(rep(40,k-1),40*2),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[3]
                               }
               if(i>96 & i <=104) {
                     p_val[i,1]=theo_power(m=c(rep(0,k-1),1),n=c(rep(50,k-1),50*0.5),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[1]
                     p_val[i,2]=theo_power(m=c(rep(0,k-1),1),n=c(rep(50,k-1),50*0.5),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[2]
                     p_val[i,3]=theo_power(m=c(rep(0,k-1),1),n=c(rep(50,k-1),50*0.5),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[3]
                               }
               if(i>104 & i <=112) {
                     p_val[i,1]=theo_power(m=c(rep(0,k-1),1),n=c(rep(50,k-1),50*1),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[1]
                     p_val[i,2]=theo_power(m=c(rep(0,k-1),1),n=c(rep(50,k-1),50*1),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[2]
                     p_val[i,3]=theo_power(m=c(rep(0,k-1),1),n=c(rep(50,k-1),50*1),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[3]
                               }
               if(i>112 & i <=120) {
                     p_val[i,1]=theo_power(m=c(rep(0,k-1),1),n=c(rep(50,k-1),50*1.5),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[1]
                     p_val[i,2]=theo_power(m=c(rep(0,k-1),1),n=c(rep(50,k-1),50*1.5),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[2]
                     p_val[i,3]=theo_power(m=c(rep(0,k-1),1),n=c(rep(50,k-1),50*1.5),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[3]
                               }
               if(i>120 & i <=128) {
                     p_val[i,1]=theo_power(m=c(rep(0,k-1),1),n=c(rep(50,k-1),50*2),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[1]
                     p_val[i,2]=theo_power(m=c(rep(0,k-1),1),n=c(rep(50,k-1),50*2),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[2]
                     p_val[i,3]=theo_power(m=c(rep(0,k-1),1),n=c(rep(50,k-1),50*2),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[3]
                               }
               if(i>128 & i <=136) {
                     p_val[i,1]=theo_power(m=c(rep(0,k-1),1),n=c(rep(100,k-1),100*0.5),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[1]
                     p_val[i,2]=theo_power(m=c(rep(0,k-1),1),n=c(rep(100,k-1),100*0.5),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[2]
                     p_val[i,3]=theo_power(m=c(rep(0,k-1),1),n=c(rep(100,k-1),100*0.5),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[3]
                               }
               if(i>136 & i <=144) {
                     p_val[i,1]=theo_power(m=c(rep(0,k-1),1),n=c(rep(100,k-1),100*1),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[1]
                     p_val[i,2]=theo_power(m=c(rep(0,k-1),1),n=c(rep(100,k-1),100*1),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[2]
                     p_val[i,3]=theo_power(m=c(rep(0,k-1),1),n=c(rep(100,k-1),100*1),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[3]
                               }
               if(i>144 & i <=152) {
                     p_val[i,1]=theo_power(m=c(rep(0,k-1),1),n=c(rep(100,k-1),100*1.5),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[1]
                     p_val[i,2]=theo_power(m=c(rep(0,k-1),1),n=c(rep(100,k-1),100*1.5),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[2]
                     p_val[i,3]=theo_power(m=c(rep(0,k-1),1),n=c(rep(100,k-1),100*1.5),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[3]
                               }
               if(i>152 & i <=160) {
                     p_val[i,1]=theo_power(m=c(rep(0,k-1),1),n=c(rep(100,k-1),100*2),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[1]
                     p_val[i,2]=theo_power(m=c(rep(0,k-1),1),n=c(rep(100,k-1),100*2),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[2]
                     p_val[i,3]=theo_power(m=c(rep(0,k-1),1),n=c(rep(100,k-1),100*2),sd=c(rep(2,k-1),2*SDR[(i-1)/2+1]))[3]
                               }

                 }

        }
setwd(dir="C:/Users/mdela/Dropbox/ANOVA's Welch/theoretical power, copy paste for different k (make it more dynamic for final paper)") # destination file
write.table(p_val,str_c("theoretical power with k =",k,".txt"),sep=";",dec=",")
 }


theoretical_power(k=3)
theoretical_power(k=4)
theoretical_power(k=5)
theoretical_power(k=6)
theoretical_power(k=7)
theoretical_power(k=8)
theoretical_power(k=9)

View(p_val)


