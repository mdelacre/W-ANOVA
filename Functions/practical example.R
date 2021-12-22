A1 <- c(19.5,31.5,14,27,25,26,9,24,24,36,9,26,24,38,9,22,19,33.5,24.5,27,23.5,21,12,39,15.5,39,10.5,21,14.5,34.5,27.5,32,16,22.5,16.5,35.5,10.5,35.5,27,38.5,24)
A2 <- c(27,19,24,24,24,27,27,22,23,20,17,21,18,19,27,22,27,25,25,21.5,23.5)
A3 <- c(31,27,39,27,15,19,21,25,25,21,35,35,29,39,35,29,27,29,37,21,23,31,25,23,21,23,29,29,25,17,25)

DV=c(A1,A2,A3)
IV=factor(c(rep(1,length(A1)),rep(2,length(A2)),rep(3,length(A3))))

database<-data.frame(IV = factor(IV), DV) 

#------------------------------------------------------------
                   Performing a Shapiro-Wilk test
#------------------------------------------------------------

shapiro.test(database$DV[database$IV==1]) # W=0.95365; p = .09 --> NRH0
shapiro.test(database$DV[database$IV==2]) # W=0.93191; p = .1503 --> NRH0
shapiro.test(database$DV[database$IV==3]) # W=0.96785; p = .462 --> NRH0

#------------------------------------------------------------
                  Performing a W-test
#------------------------------------------------------------

#oneway.test(dv.name ~ iv.name, data=data.name, var.equal=FALSE)

oneway.test(DV ~ IV, data=database, var.equal=FALSE)
