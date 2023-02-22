rm(list=ls(all.names=TRUE))
setwd("C:/Users/wang0075/Desktop/AHR/github")
source("cpDesign_fun.R")

## design parameters for PH trial
alpha <- 0.025  ## type I error rate to be controled
sided <- 1      ## 1=one sided test, 2=two sided test
beta <- 1-0.90  ## type II error rate or 1-power
r <- 1          ## allocation ratio bteween two treatment groups
b <- 20         ## accrual rate per time unit
B <- 34         ## accural length
lamc <- -log(0.5)/12    ## hazard rate for the control
HR.H0 <- 1              ## hazard ratio under H0
HR.H1 <- 0.75           ## hazard ratio under H1

theta.H1 <- log(HR.H1)
theta.H0 <- log(HR.H0)
AHR <- exp(theta.H1)
theta.cp <- theta.H1

## Create a reference PH design
## The reference design can be obtained using get.REFdesign() function
## or using the nSurv() function in gsDesign
S1 <- "FS"
s1 <- nSurv(sided=1, alpha=alpha, beta=beta, hr=HR.H1, lambdaC=lamc, eta=0, gamma=b, R=B)
N <- as.numeric(s1$R*s1$gamma)
nEvents <- as.numeric(s1$d)
B <- s1$R
L <- s1$R+s1$minfu
s1 <- data.frame(S=S1, N=N, B=B, L=L, nEvents=nEvents, AHRK.H1=exp(theta.H1), pow.nbind=s1$power)

## The reference design can be created using group sequential method to mornitor early stopping for efficacy
S2 <- "GSD"
s2 <- gsSurv(k=2, timing=1, sfu = sfLDOF, test.type=1, alpha=alpha, beta=beta, hr=HR.H1, lambdaC=lamc, eta=0, gamma=b, R=B)
N <- as.numeric(s2$R*s2$gamma)
nEvents <- as.numeric(s2$eDE[2]+s2$eDC[2])
B <- s2$R
L <- s2$R+s2$minfu
s2 <- data.frame(S=S2, N=N, B=B, L=L, nEvents=nEvents, AHRK.H1=exp(theta.H1), pow.nbind=1-s2$beta)

## Based on the FS and GSD design, we decide to use the following desing for the phase III trial
## Create futility stopping rules for the PH(0.75) design similar to Table 2
nEvents <- round(nEvents)
S3 <- paste0("PH(", HR.H1, ")")
lam1c <- lamc; HR1 <- HR2 <- HR.H1
L <- B+find.fu.PH(N, B, lamc, HR.H1, nEvents, alpha)
s3 <- get.bounds.summary(S3, N, B, L, 0, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha)

## Create futility stopping rules for the NPH(3,1,0.6938) design similar to Table 2
## This NPH has delayed treatment effect until td=3
## HR in the first time intervla is 1.0
## HR in the 2nd time interval is caculated so that the AHR=0.75
pow0 <- pow(r, nEvents, alpha, exp(theta.H1))
lam1c <- lamc; HR1 <- c(1); td <- 3
HR2 <- get.HR2(N, B, td, lam1c, HR1, AHR, pow0, alpha)
S4 <- paste0("NPH(3,1,",round(HR2,4),")")
L <- B+find.fu(N, B, td, lam1c, HR1, HR2, pow0, alpha)
nEvents <- find.ne.new(N, td, B, L, lam1c, HR1, HR2)$e
s4 <- get.bounds.summary(S4, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha)

ss <- bind_rows(s1,s2,s3,s4)
write.csv(ss, "application.csv")
