rm(list=ls(all.names = TRUE))
setwd("C:/Users/wang0075/Desktop/AHR/github")
source("cpDesign_fun.R")

######################################################################################
## design parameters for PH trial
alpha <- 0.025
sided <- 1
beta <- 0.1
r <- 1
B <- 34
N <- 680
lamc <- -log(0.5)/12
HR.H0 <- 1
HR.H1 <- 0.75
AHR <- HR.H1
theta.H0 <- log(1)
theta.cp <- log(AHR)

nEvents0 <- 508
pow0 <- pow(r, nEvents0, alpha, HR.H1)
L0 <- B+find.fu(N, B, 3, lamc, AHR, AHR, pow0, alpha)

S0 <- "PH(0.75)"
td <- 3
lam1c <- lamc
HR1 <- 0.75
HR2 <- 0.75
res0 <- get.fig1.data(S0, N, B, L0, td, nEvents0, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha)

lam1c <- lamc; HR1 <- c(1); td <- 3
HR2 <- get.HR2(N, B, td, lam1c, HR1, AHR, pow0, alpha)
S1 <- paste0("NPH(3,1.0,",round(HR2,4),")")
L <- B+find.fu(N, B, td, lam1c, HR1, HR2, pow0, alpha)
nEvents <- find.ne.new(N, td, B, L, lam1c, HR1, HR2)$e
res1 <- get.fig1.data(S1, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha)

lam1c <- lamc; HR1 <- c(1.3); td <- 3
HR2 <- get.HR2(N, B, td, lam1c, HR1, AHR, pow0, alpha)
S2 <- paste0("NPH(3,1.3,",round(HR2,4),")")
L <- B+find.fu(N, B, td, lam1c, HR1, HR2, pow0, alpha)
nEvents <- find.ne.new(N, td, B, L, lam1c, HR1, HR2)$e
res2 <- get.fig1.data(S2, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha)

lam1c <- lamc; HR1 <- c(1); td <- 6
HR2 <- get.HR2(N, B, td, lam1c, HR1, AHR, pow0, alpha)
S3 <- paste0("NPH(6,1,",round(HR2,4),")")
L <- B+find.fu(N, B, td, lam1c, HR1, HR2, pow0, alpha)
nEvents <- find.ne.new(N, td, B, L, lam1c, HR1, HR2)$e
res3 <- get.fig1.data(S3, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha)

res <- bind_rows(res0,res1,res2,res3)

plot.fig1(res, "Fig1.pdf")
