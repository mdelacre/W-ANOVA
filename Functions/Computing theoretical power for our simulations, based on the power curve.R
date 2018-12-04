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

    return(c(power_anova,power_ANOVAwelch,power_BF))

}

