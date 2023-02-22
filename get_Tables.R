rm(list=ls(all.names = TRUE))
setwd("C:/Users/wang0075/Desktop/AHR/github")
source("cpDesign_fun.R")

################################################################################
################################################################################
################################################################################
## Create the reference design
ref <- get.REFdesign(alpha=0.025, power=0.90, sided=1, r=1, lamc=-log(0.5)/12, HR.H1=0.75, B=34, gamma=20)
print(ref)

################################################################################
tab2 <- get.table2(HR.H1=0.75, B=34, N=680, nEvents=508)
write.csv(tab2, "table2.csv")

################################################################################
tab3 <- get.table3(HR.H1=0.75, B=34, N=680, nEvents=508)
write.csv(tab3, "table3.csv")

tab4a <- get.table4(HR.H1=0.75, B=12, N=680, nEvents=508, look=1)
tab4b <- get.table4(HR.H1=0.75, B=34, N=680, nEvents=508, look=1)
tab4c <- get.table4(HR.H1=0.75, B=60, N=680, nEvents=508, look=1)
write.csv(rbind(tab4a, tab4b, tab4c), "table4.csv")
