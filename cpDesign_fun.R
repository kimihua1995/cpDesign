options(scipen=100, warn=-1)
library(dplyr)
library(gsDesign)
library(mvtnorm)
library(ggplot2)
library(grDevices)
library(gridExtra)
#options(scipen=0)

################################################################################
## N sample size
## B accrual time
## L study length
## td end time of delayed effect
## tia time of IA
## r allocation ratio
## alpha
## lambdac lambda for control
## first period HR1
## second period HR2
## method=num_int numberic integration method

################################################################################
ne <- function(r, alpha, beta, AHR) {
    ceiling((1+r)^2/r*(qnorm(1-alpha)+qnorm(1-beta))^2/log(AHR)^2)
}

################################################################################
pow <- function(r, nEvents, alpha, AHR) {
    pnorm(sqrt(nEvents)*log(1/AHR)/((1+r)/sqrt(r))-qnorm(1-alpha))
}

pow.new <- function(r, lamc, B, tia, N, alpha, HR) {
    r <- find.ne.PH.new(N, B, tia, lamc, HR)
    pow <- pnorm(log(1/HR)/sqrt(1/(r$Ntia/2*r$pe)+1/(r$Ntia/2*r$pc))-qnorm(1-alpha))
    return(list(pow=pow, nee=r$Ntia/2*r$pe, nec=r$Ntia/2*r$pc, ne=ceiling(r$Ntia/2*r$pe+r$Ntia/2*r$pc)))
}

################################################################################
peL1.num.PH <- function(B, tia, lam) {
    St <- function(x) exp(-lam*x)
    Sx <- function(x) St(tia-x)
    return(1-1/min(tia,B)*integrate(Sx, lower=0, upper=min(tia,B))$value)
}

###############################################################################
find.ne.PH <- function(N, B, tia, lamc, lame) {
    pc <- peL1.num.PH(B, tia, lamc)
    pe <- peL1.num.PH(B, tia, lame)
    Ntia <- N/B*min(tia,B)
    list(e=pc*Ntia/2+pe*Ntia/2, pc=pc, pe=pe, Ntia=Ntia)
}

find.ne.PH.new <- function(N, B, tia, lamc, HR) {
    lame <- lamc*HR
    pc <- peL1.num.PH(B, tia, lamc)
    pe <- peL1.num.PH(B, tia, lame)
    Ntia <- N/B*min(tia,B)
    list(e=pc*Ntia/2+pe*Ntia/2, pc=pc, pe=pe, Ntia=Ntia)
}

##############################################################################
find.fu.PH <- function(N, B, lamc, HR, nEvents, alpha) {
    f <- function(x) {
        nex <- find.ne.PH(N, B, B+x, lamc, lame)$e
        return(nex-nEvents)
    }
    lame <- lamc*HR
    uniroot(f, interval=c(0,100000))$root
}

################################################################################
get.time.PH <- function(pct.nEvents, N, B, L, nEvents, lamc, HR) {
 f <- function (x) {
    lame <- HR*lamc
    return(find.ne.PH(N, B, x, lamc, lame)$e-pct.nEvents*nEvents)
 }
 uniroot(f, interval=c(0.0001, 1.5*L))$root
}

###############################################################################
### conditional power
cp <- function(Zk, Ik, IK, theta, alpha){
    pnorm((-Zk*sqrt(Ik)-qnorm(1-alpha)*sqrt(IK)-theta*(IK-Ik))/sqrt(IK-Ik))
}

###############################################################################
fck <- function(Ik, IK, thetak, thetaK, alpha, method="Wieand", etak=0, theta.cp=log(0.75), gk=NULL) {
    if (method=="Wieand" || method=="KF") ck <- thetak*sqrt(Ik)-etak
    else if (method=="cp") { etak <- -(qnorm(1-alpha)*sqrt(IK)+theta.cp*(IK-Ik)+qnorm(gk)*sqrt(IK-Ik))/sqrt(Ik); ck <- thetak*sqrt(Ik)-etak }
    return(ck)
}


##############################################################################
PF1.new <- function(Ik, IK, thetak, thetaK, alpha, etak=0, theta.cp, gk=NULL) {
    if (!is.null(gk)) etak <- (-qnorm(gk)*sqrt(IK-Ik)-qnorm(1-alpha)*sqrt(IK)-theta.cp*(IK-Ik))/sqrt(Ik)
    zk <- etak
    cpk <- pnorm((-zk*sqrt(Ik)-qnorm(1-alpha)*sqrt(IK)-theta.cp*(IK-Ik))/sqrt(IK-Ik))
    ck <- -etak+thetak*sqrt(Ik)
    pr.S <- pnorm(ck)
    cK <- qnorm(1-alpha)+thetaK*sqrt(IK)
    pow.nbind <- 1-pnorm(cK)
    sigma <- diag(2)
    sigma[2,1] <- sigma[1,2] <- sqrt(Ik/IK)
    pow.bind <- pmvnorm(lower=c(ck, cK), upper=c(Inf, Inf), mean=c(0,0), sigma=sigma)
    return(list(pr.S=pr.S, cpk=cpk, pow.nbind=pow.nbind, pow.bind=pow.bind, powloss=pow.nbind-pow.bind))
}

###############################################################################
PF2.new <- function(Ij, Ik, IK, thetaj, thetak, thetaK, alpha, etaj=0, etak=0, theta.cp, gj=NULL, gk=NULL) {
    if (!is.null(gk)) {
        etaj <- (-qnorm(gj)*sqrt(IK-Ij)-qnorm(1-alpha)*sqrt(IK)-theta.cp*(IK-Ij))/sqrt(Ij)
        etak <- (-qnorm(gk)*sqrt(IK-Ik)-qnorm(1-alpha)*sqrt(IK)-theta.cp*(IK-Ik))/sqrt(Ik)
    }
    zj <- etaj; zk <- etak
    cpj <- pnorm((-zj*sqrt(Ij)-qnorm(1-alpha)*sqrt(IK)-theta.cp*(IK-Ij))/sqrt(IK-Ij))
    cpk <- pnorm((-zk*sqrt(Ik)-qnorm(1-alpha)*sqrt(IK)-theta.cp*(IK-Ik))/sqrt(IK-Ik))
    cj <- -etaj+thetaj*sqrt(Ij)
    ck <- -etak+thetak*sqrt(Ik)
    xj <- pnorm(cj)
    xk <- pnorm(ck)
    sigma <- diag(2)
    sigma[2,1] <- sigma[1,2] <- sqrt(Ij/Ik)
    xjk <- pmvnorm(lower=c(-Inf, -Inf), upper=c(cj, ck), mean=c(0,0), sigma=sigma)
    cK <- qnorm(1-alpha)+thetaK*sqrt(IK)
    pow.nbind <- 1-pnorm(cK)
    sigma <- diag(3)
    sigma[2,1] <- sigma[1,2] <- sqrt(Ij/Ik)
    sigma[3,1] <- sigma[1,3] <- sqrt(Ij/IK)
    sigma[3,2] <- sigma[2,3] <- sqrt(Ik/IK)
    pow.bind <- pmvnorm(lower=c(cj, ck, cK), upper=c(Inf, Inf, Inf), mean=c(0,0,0), sigma=sigma)
    return(list(pr.S1=xj, pr.S2cS1=xk-xjk, pr.S=xj+xk-xjk, cpj=cpj, cpk=cpk, pow.nbind=pow.nbind, pow.bind=pow.bind, powloss=pow.nbind-pow.bind))
    }

## Functions for obtaining average HR
################################################################################
peL1.num <- function(td, B, tia, lam) {
    St <- function(x) {I(x<td)*exp(-lam*x)+I(x>=td)*exp(-lam*td)}
    Sx <- function(x) St(tia-x)
    return(1-1/min(tia,B)*integrate(Sx, lower=0, upper=min(tia,B))$value)
}

################################################################################
peL.num <- function(td, B, tia, lam1, lam2) {
    St <- function(x) {I(x<td)*exp(-lam1*x)+I(x>=td)*exp(-lam1*td-lam2*(x-td))}
    Sx <- function(x) St(tia-x)
    return(1-1/min(tia,B)*integrate(Sx, lower=0, upper=min(tia,B))$value)
}


################################################################################
## find num of events occuring before td or the end of delayed treatment
find1.ne <- function(N, td, B, tia, lam1c, lam1e) {
    pc <- peL1.num(td, B, tia, lam1c)
    pe <- peL1.num(td, B, tia, lam1e)
    Ntia <- N/B*min(tia,B)
    list(e1=pc*Ntia/2+pe*Ntia/2, pc1=pc, pe1=pe, Ntia=Ntia)
}

################################################################################
find.ne <- function(N, td, B, tia, lam1c, lam2c, lam1e, lam2e) {
    pc <- peL.num(td, B, tia, lam1c, lam2c)
    pe <- peL.num(td, B, tia, lam1e, lam2e)
    Ntia <- N/B*min(tia,B)
    list(e=pc*Ntia/2+pe*Ntia/2, pc=pc, pe=pe, Ntia=Ntia)
}

################################################################################
find.ne.new <- function(N, td, B, tia, lam1c, HR1, HR2) {
    lam2c <- lam1c
    lam1e <- HR1*lam1c
    lam2e <- HR2*lam2c
    pc <- peL.num(td, B, tia, lam1c, lam2c)
    pe <- peL.num(td, B, tia, lam1e, lam2e)
    Ntia <- N/B*min(tia,B)
    list(e=pc*Ntia/2+pe*Ntia/2, pc=pc, pe=pe, Ntia=Ntia)
}

## get2.time <- function(pct.nEvents, pct.fu3, N, B, L, td, nEvents, lam1c, HR1, HR2)
## tmin <- get.time(pct.nEvents, N, B, L, td, nEvents, lam1c, HR1, HR2)
## Ntmin <- N/B*min(tmin,B)
## 2/3*Ntmin
################################################################################
get.time <- function(pct.nEvents, N, B, L, td, nEvents, lam1c, HR1, HR2) {
 f <- function (x) {
    lam2c <- lam1c
    lam1e <- HR1*lam1c
    lam2e <- HR2*lam2c
    return(find.ne(N, td, B, x, lam1c, lam2c, lam1e, lam2e)$e-pct.nEvents*nEvents)
 }
 uniroot(f, interval=c(0.0001, 1.5*L))$root
}

get.St <- function(td, tx, lam1c, HR1, HR2) {
    lam2c <- lam1c
    lam1e <- HR1*lam1c
    lam2e <- HR2*lam2c
    return(I(tx<td)*exp(-lam1e*tx)+I(tx>=td)*exp(-lam1e*td-lam2e*(tx-td)))
}

################################################################################
get.AHR <- function(N, B, L, td, tia, lam1c, HR1, HR2) {
    lam2c <- lam1c
    lam1e <- HR1*lam1c
    lam2e <- HR2*lam2c
    e1 <- find1.ne(N, td, B, tia, lam1c, lam1e)
    e <- find.ne(N, td, B, tia, lam1c, lam2c, lam1e, lam2e)
    e2 <- e$e-e1$e1
    p1 <- e1$e1/e$e
    p2 <- e2/e$e
    AHR <- exp(p1*log(HR1)+p2*log(HR2))
}

################################################################################
get.Ek <- function(N, B, L, td, tia, lam1c, HR1, HR2) {
    lam2c <- lam1c
    lam1e <- HR1*lam1c
    lam2e <- HR2*lam2c
    e <- find.ne(N, td, B, tia, lam1c, lam2c, lam1e, lam2e)
}

################################################################################
find.fu <- function(N, B, td, lam1c, HR1, HR2, pow0, alpha) {
    f <- function(x) {
        thetax <- log(get.AHR(N, B, B+x, td, B+x, lam1c, HR1, HR2))
        Ex <- find.ne.new(N, td, B, B+x, lam1c, HR1, HR2)$e
        powx <- 1-pnorm(qnorm(1-alpha)+thetax*sqrt(Ex/4))
        return(powx-pow0)
    }
    thetaB <- log(get.AHR(N, B, B, td, B, lam1c, HR1, HR2))
    EB <- find.ne.new(N, td, B, B, lam1c, HR1, HR2)$e
    powB <- 1-pnorm(qnorm(1-alpha)+thetaB*sqrt(EB/4))
    if(powB>pow0 || abs(powB-pow0)<0.000001) return(0)
    else uniroot(f, interval=c(0,100000))$root
}

get.HR2 <- function(N, B, td, lam1c, HR1, AHR, pow0, alpha) {
    f <- function(x) {
        L <- B+find.fu(N, B, td, lam1c, HR1, x, pow0, alpha)
        AHRx <- get.AHR(N, B, L, td, L, lam1c, HR1, x)
        return(AHRx-AHR)
    }
 ##the upper limit 0.70 needs further work to be correct
 if(HR1>1) uniroot(f, interval=c(0.01, 0.65))$root
 else uniroot(f, interval=c(0.01, 0.70))$root
 ##else uniroot(f, interval=c(0.01, 0.71))$root
 ##if(HR1>1) uniroot(f, interval=c(0.01, 0.5))$root
 ##else uniroot(f, interval=c(0.01, 0.55))$root
}


################################################################################
find.pctfu2td <- function(td, B, tia) {
    St <- function(x) I(x<td)*1+I(x>=td)*0
    Sx <- function(x) St(tia-x)
    pct.fu2td <- (1-1/min(tia,B)*integrate(Sx, lower=0, upper=min(tia,B))$value)
}

################################################################################
## get the time that at least 2/3 patients have been followed more than 3 months
get.EFtime <- function(qk, N, B, L, td, EK, lam1c, HR1, HR2) {
    tk0 <- get.time(qk, N, B, L, td, EK, lam1c, HR1, HR2)
    pctk0 <- find.pctfu2td(td, B, tk0)
    if (pctk0>2/3) tEFk <- tk0
    else {
        f <- function(x) find.pctfu2td(td, B, tk0+x)-2/3
        x <- uniroot(f, c(0, 1000))$root
        tEFk <- tk0+x
    }
    pctEFk <- find.pctfu2td(td, B, tEFk)
    Nk <- N/B*min(tk0, B)
    NEFk <- N/B*min(tEFk, B)
    res <- cbind(qk, tk0, tEFk, pctk0, pctEFk, Nk, NEFk)
    res
}

################################################################################
## get the time that at least 2/3 events occurred after 3 months from randomization.
## this is ad hoc rule and after 3 months why not 6 months for scenario 6
## this is a big limitation for EF rule, as for some scenarios, the time to observe >2/3 events occurs after 3 months will
## exceeds the study druation.
find.pctnetd <- function(N, td, B, tia, lam1c, HR1, HR2) {
    lam2c <- lam1c
    lam1e <- HR1*lam1c
    lam2e <- HR2*lam2c
    ## e1 <- find1.ne(N, td, B, tia, lam1c, lam1e)
    e1 <- find1.ne(N, 3, B, tia, lam1c, lam1e)
    e <- find.ne(N, td, B, tia, lam1c, lam2c, lam1e, lam2e)
    (e$e-e1$e1)/e$e
}

get.KF2time <- function(qk, N, B, L, td, EK, lam1c, HR1, HR2) {
    tk0 <- get.time(qk, N, B, L, td, EK, lam1c, HR1, HR2)
    pctnetdk0 <- find.pctnetd(N, td, B, tk0, lam1c, HR1, HR2)
    if (pctnetdk0>2/3) tEFk <- tk0
    else {
        f <- function(x) find.pctnetd(N, td, B, tk0+x, lam1c, HR1, HR2)-2/3
        x <- uniroot(f, c(0, 1000))$root
        tEFk <- tk0+x
    }
    pctnetdEFk <- find.pctnetd(N, td, B, tEFk, lam1c, HR1, HR2)
    Nk <- N/B*min(tk0, B)
    NEFk <- N/B*min(tEFk, B)
    tEFk
}

################################################################################
Wieand.look1 <- function(qk, etak, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha, gk=NULL) {
    P1 <- r/(1+r)
    qK <- 1
    Ek <- nEvents*qk
    EK <- nEvents*qK
    Ik <- Ek*P1*(1-P1)
    IK <- EK*P1*(1-P1)
    if (!is.null(gk)) etak <- (-qnorm(gk)*sqrt(IK-Ik)-qnorm(1-alpha)*sqrt(IK)-theta.cp*(IK-Ik))/sqrt(Ik)

    tk <- get.time(qk, N, B, L, td, EK, lam1c, HR1, HR2)
    thetak.H1 <- log(get.AHR(N, B, L, td, tk, lam1c, HR1, HR2))
    tK <- get.time(qK, N, B, L, td, EK, lam1c, HR1, HR2)
    Nk <- N/B*min(tk,B)
    NK <- N/B*min(tK,B)
    thetaK.H1 <- log(get.AHR(N, B, L, td, tK, lam1c, HR1, HR2))
    H1k <- PF1.new(Ik, IK, thetak.H1, thetaK.H1, alpha, etak, theta.cp)
    H1.AE <- H1k$pr.S*Ek +(1-H1k$pr.S)*EK
    H1.AN <- H1k$pr.S*Nk +(1-H1k$pr.S)*NK
    H1.AL <- H1k$pr.S*tk +(1-H1k$pr.S)*tK

    HR.H0 <- exp(theta.H0)
    Ek0 <- Ek; EK0 <- EK
    Ik0 <- Ek0*P1*(1-P1); IK0 <- EK0*P1*(1-P1)
    tk0 <- get.time(qk, N, B, L, td, EK, lam1c, HR.H0, HR.H0)
    tK0 <- get.time(qK, N, B, L, td, EK, lam1c, HR.H0, HR.H0)
    Nk0 <- N/B*min(tk0,B); NK0 <- N/B*min(tK0,B)
    H0k <- PF1.new(Ik0, IK0, theta.H0, theta.H0, alpha, etak, theta.cp)

    H0.AE <- H0k$pr.S*Ek0+(1-H0k$pr.S)*EK0
    H0.AN <- H0k$pr.S*Nk0+(1-H0k$pr.S)*NK0
    H0.AL <- H0k$pr.S*tk0+(1-H0k$pr.S)*tK0
    if(etak>100) rule="nostop" else if(etak==0) rule="Wieand1" else if (etak==-0.011) rule="OF1" else if (etak==-0.490 || !is.null(gk)) rule="Xi1"
    res <- data.frame(rule, N, B, L, nEvents, etak, qk, qK, tk, tK, AHRk.H1=exp(thetak.H1), AHRK.H1=exp(thetaK.H1),
    H1.pr.S=H1k$pr.S, H1.cpk=H1k$cpk, H1.AE, H1.AN, H1.AL, H1.pow.nbind=H1k$pow.nbind, H1.pow.bind=H1k$pow.bind, H1.powloss=H1k$pow.nbind-H1k$pow.bind,
    H0.pr.S=H0k$pr.S, H0.cpk=H0k$cpk, H0.AE, H0.AN, H0.AL)
}


################################################################################
Wieand.look2 <- function(qj, qk, etaj, etak, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha, gj=NULL, gk=NULL) {
    P1 <- r/(1+r)
    qK <- 1
    Ej <- nEvents*qj
    Ek <- nEvents*qk
    EK <- nEvents*qK
    Ij <- Ej*P1*(1-P1)
    Ik <- Ek*P1*(1-P1)
    IK <- EK*P1*(1-P1)

    if (!is.null(gk)) {
        etaj <- (-qnorm(gj)*sqrt(IK-Ij)-qnorm(1-alpha)*sqrt(IK)-theta.cp*(IK-Ij))/sqrt(Ij)
        etak <- (-qnorm(gk)*sqrt(IK-Ik)-qnorm(1-alpha)*sqrt(IK)-theta.cp*(IK-Ik))/sqrt(Ik)
    }

    tj <- get.time(qj, N, B, L, td, EK, lam1c, HR1, HR2)
    thetaj.H1 <- log(get.AHR(N, B, L, td, tj, lam1c, HR1, HR2))
    tk <- get.time(qk, N, B, L, td, EK, lam1c, HR1, HR2)
    thetak.H1 <- log(get.AHR(N, B, L, td, tk, lam1c, HR1, HR2))
    tK <- get.time(qK, N, B, L, td, EK, lam1c, HR1, HR2)
    Nj <- N/B*min(tj,B)
    Nk <- N/B*min(tk,B)
    NK <- N/B*min(tK,B)
    thetaK.H1 <- log(get.AHR(N, B, L, td, tK, lam1c, HR1, HR2))
    H1jk <- PF2.new(Ij, Ik, IK, thetaj.H1, thetak.H1, thetaK.H1, alpha, etaj, etak, theta.cp)

    H1.AE <- H1jk$pr.S1*Ej +H1jk$pr.S2cS1*Ek +(1-H1jk$pr.S1-H1jk$pr.S2cS1)*EK
    H1.AN <- H1jk$pr.S1*Nj +H1jk$pr.S2cS1*Nk +(1-H1jk$pr.S1-H1jk$pr.S2cS1)*NK
    H1.AL <- H1jk$pr.S1*tj +H1jk$pr.S2cS1*tk +(1-H1jk$pr.S1-H1jk$pr.S2cS1)*tK

    HR.H0 <- exp(theta.H0)
    Ej0 <- Ej; Ek0 <- Ek; EK0 <- EK
    tj0 <- get.time(qj, N, B, L, td, EK, lam1c, HR.H0, HR.H0)
    tk0 <- get.time(qk, N, B, L, td, EK, lam1c, HR.H0, HR.H0)
    tK0 <- get.time(qK, N, B, L, td, EK, lam1c, HR.H0, HR.H0)
    Nj0 <- N/B*min(tj0,B); Nk0 <- N/B*min(tk0,B); NK0 <- N/B*min(tK0,B)
    Ij0 <- Ej0*P1*(1-P1); Ik0 <- Ek0*P1*(1-P1); IK0 <- EK0*P1*(1-P1)
    H0jk <- PF2.new(Ij0, Ik0, IK0, theta.H0,  theta.H0,  theta.H0,  alpha, etaj, etak, theta.cp)

    H0.AE <- H0jk$pr.S1*Ej0+H0jk$pr.S2cS1*Ek0+(1-H0jk$pr.S1-H0jk$pr.S2cS1)*EK0
    H0.AN <- H0jk$pr.S1*Nj0+H0jk$pr.S2cS1*Nk0+(1-H0jk$pr.S1-H0jk$pr.S2cS1)*NK0
    H0.AL <- H0jk$pr.S1*tj0+H0jk$pr.S2cS1*tk0+(1-H0jk$pr.S1-H0jk$pr.S2cS1)*tK0

    if(etaj>100 && etak>100) rule="nostop" else if(etaj==0 && etak==0) rule="Wieand2" else if (etaj==-0.011 && etak==-0.864) rule="OF2" else if ((etaj==-0.02 && etak==-0.78) || !is.null(gk)) rule="Xi2"

    res <- data.frame(rule, N, B, L, nEvents, etaj, etak, qj, qk, qK, tj, tk, tK, AHRj.H1=exp(thetaj.H1), AHRk.H1=exp(thetak.H1), AHRK.H1=exp(thetaK.H1),
    H1.pow.nbind=H1jk$pow.nbind, H1.pow.bind=H1jk$pow.bind, H1.powloss=H1jk$pow.nbind-H1jk$pow.bind,
    H1.pr.S1=H1jk$pr.S1, H1.pr.S=H1jk$pr.S, H1.cpj=H1jk$cpj, H1.cpk=H1jk$cpk, H1.AE, H1.AN, H1.AL,
    H0.pr.S1=H0jk$pr.S1, H0.pr.S=H0jk$pr.S, H0.cpj=H0jk$cpj, H0.cpk=H0jk$cpk, H0.AE, H0.AN, H0.AL)
}


################################################################################
KF.look1 <- function(qk, etak, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha) {
    rule="KF1"
    P1 <- r/(1+r)
    qK <- 1
    Ek <- nEvents*qk
    EK <- nEvents*qK
    Ik <- Ek*P1*(1-P1)
    IK <- EK*P1*(1-P1)

    tk <- get.KF2time(qk, N, B, L, td, EK, lam1c, HR1, HR2)
    thetak.H1 <- log(get.AHR(N, B, L, td, tk, lam1c, HR1, HR2))
    tK <- get.time(qK, N, B, L, td, EK, lam1c, HR1, HR2)
    Nk <- N/B*min(tk,B); NK <- N/B*min(tK,B)
    thetaK.H1 <- log(get.AHR(N, B, L, td, tK, lam1c, HR1, HR2))

    Ek <- find.ne.new(N, td, B, tk, lam1c, HR1, HR2)$e
    EK <- find.ne.new(N, td, B, tK, lam1c, HR1, HR2)$e
    Ik <- Ek*P1*(1-P1); IK <- EK*P1*(1-P1)
    H1k <- PF1.new(Ik, IK, thetak.H1, thetaK.H1, alpha, etak, theta.cp)
    H1.AE <- H1k$pr.S*Ek +(1-H1k$pr.S)*EK
    H1.AN <- H1k$pr.S*Nk +(1-H1k$pr.S)*NK
    H1.AL <- H1k$pr.S*tk +(1-H1k$pr.S)*tK

    HR.H0 <- exp(theta.H0)
    tk0 <- get.KF2time(qk, N, B, L, td, EK, lam1c, HR.H0, HR.H0)
    tK0 <- get.time(qK, N, B, L, td, EK, lam1c, HR.H0, HR.H0)
    Nk0 <- N/B*min(tk0,B); NK0 <- N/B*min(tK0,B)
    Ek0 <- find.ne.new(N, td, B, tk0, lam1c, HR.H0, HR.H0)$e
    EK0 <- find.ne.new(N, td, B, tK0, lam1c, HR.H0, HR.H0)$e
    Ik0 <- Ek0*P1*(1-P1); IK0 <- EK0*P1*(1-P1)
    H0k <- PF1.new(Ik0, IK0, theta.H0, theta.H0, alpha, etak, theta.cp)

    H0.AE <- H0k$pr.S*Ek0+(1-H0k$pr.S)*EK0
    H0.AN <- H0k$pr.S*Nk0+(1-H0k$pr.S)*NK0
    H0.AL <- H0k$pr.S*tk0+(1-H0k$pr.S)*tK0

    res <- data.frame(rule, N, B, L, nEvents, etak, qk=Ik0/IK0, qK, tk, tK, AHRk.H1=exp(thetak.H1), AHRK.H1=exp(thetaK.H1),
    H1.pr.S=H1k$pr.S, H1.cpk=H1k$cpk, H1.AE, H1.AN, H1.AL, H1.pow.nbind=H1k$pow.nbind, H1.pow.bind=H1k$pow.bind, H1.powloss=H1k$powloss,
    H0.pr.S=H0k$pr.S, H0.cpk=H0k$cpk, H0.AE, H0.AN, H0.AL)
}


################################################################################
KF.look2 <- function(qj, qk, etaj, etak, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha) {
    rule <- "KF2"
    P1 <- r/(1+r)
    qK <- 1
    Ej <- nEvents*qj
    Ek <- nEvents*qk
    EK <- nEvents*qK
    Ij <- Ej*P1*(1-P1)
    Ik <- Ek*P1*(1-P1)
    IK <- EK*P1*(1-P1)
    tj <- get.KF2time(qj, N, B, L, td, EK, lam1c, HR1, HR2)
    thetaj.H1 <- log(get.AHR(N, B, L, td, tj, lam1c, HR1, HR2))
    tk <- get.KF2time(qk, N, B, L, td, EK, lam1c, HR1, HR2)
    thetak.H1 <- log(get.AHR(N, B, L, td, tk, lam1c, HR1, HR2))
    tK <- get.time(qK, N, B, L, td, EK, lam1c, HR1, HR2)
    Nj <- N/B*min(tj,B); Nk <- N/B*min(tk,B); NK <- N/B*min(tK,B)
    Ej <- find.ne.new(N, td, B, tj, lam1c, HR1, HR2)$e
    Ek <- find.ne.new(N, td, B, tk, lam1c, HR1, HR2)$e
    EK <- find.ne.new(N, td, B, tK, lam1c, HR1, HR2)$e
    Ij <- Ej*P1*(1-P1); Ik <- Ek*P1*(1-P1); IK <- EK*P1*(1-P1)
    thetaK.H1 <- log(get.AHR(N, B, L, td, tK, lam1c, HR1, HR2))
    H1jk <- PF2.new(Ij, Ik, IK, thetaj.H1, thetak.H1, thetaK.H1, alpha, etaj, etak, theta.cp)
    H1.AE <- H1jk$pr.S1*Ej +H1jk$pr.S2cS1*Ek +(1-H1jk$pr.S1-H1jk$pr.S2cS1)*EK
    H1.AN <- H1jk$pr.S1*Nj +H1jk$pr.S2cS1*Nk +(1-H1jk$pr.S1-H1jk$pr.S2cS1)*NK
    H1.AL <- H1jk$pr.S1*tj +H1jk$pr.S2cS1*tk +(1-H1jk$pr.S1-H1jk$pr.S2cS1)*tK

    HR.H0 <- exp(theta.H0)
    tj0 <- get.KF2time(qj, N, B, L, td, EK, lam1c, HR.H0, HR.H0)
    tk0 <- get.KF2time(qk, N, B, L, td, EK, lam1c, HR.H0, HR.H0)
    tK0 <- get.time(qK, N, B, L, td, EK, lam1c, HR.H0, HR.H0)
    Nj0 <- N/B*min(tj0,B); Nk0 <- N/B*min(tk0,B); NK0 <- N/B*min(tK0,B)
    Ej0 <- find.ne.new(N, td, B, tj0, lam1c, HR.H0, HR.H0)$e
    Ek0 <- find.ne.new(N, td, B, tk0, lam1c, HR.H0, HR.H0)$e
    EK0 <- find.ne.new(N, td, B, tK0, lam1c, HR.H0, HR.H0)$e
    Ij0 <- Ej0*P1*(1-P1); Ik0 <- Ek0*P1*(1-P1); IK0 <- EK0*P1*(1-P1)
    H0jk <- PF2.new(Ij0, Ik0, IK0, theta.H0,  theta.H0,  theta.H0,  alpha, etaj, etak, theta.cp)

    H0.AE <- H0jk$pr.S1*Ej0+H0jk$pr.S2cS1*Ek0+(1-H0jk$pr.S1-H0jk$pr.S2cS1)*EK0
    H0.AN <- H0jk$pr.S1*Nj0+H0jk$pr.S2cS1*Nk0+(1-H0jk$pr.S1-H0jk$pr.S2cS1)*NK0
    H0.AL <- H0jk$pr.S1*tj0+H0jk$pr.S2cS1*tk0+(1-H0jk$pr.S1-H0jk$pr.S2cS1)*tK0

    res <- data.frame(rule, N, B, L, nEvents, etaj, etak, qj=Ij0/IK0, qk=Ik0/IK0, qK, tj, tk, tK, AHRj.H1=exp(thetaj.H1), AHRk.H1=exp(thetak.H1), AHRK.H1=exp(thetaK.H1),
    H1.pow.nbind=H1jk$pow.nbind, H1.pow.bind=H1jk$pow.bind, H1.powloss=H1jk$pow.nbind-H1jk$pow.bind,
    H1.pr.S1=H1jk$pr.S1, H1.pr.S=H1jk$pr.S, H1.cpj=H1jk$cpj, H1.cpk=H1jk$cpk, H1.AE, H1.AN, H1.AL,
    H0.pr.S1=H0jk$pr.S1, H0.pr.S=H0jk$pr.S, H0.cpj=H0jk$cpj, H0.cpk=H0jk$cpk, H0.AE, H0.AN, H0.AL)
}


################################################################################
opt.look1 <- function(pks, gammas, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha) {
    cat("1 stop at pks !!\n")
    res <- NULL
    P1 <- r/(1+r)
    qK <- 1
    rule.id <- 1
    for (ttk in 1:length(pks)) {
            Ek <- nEvents*pks[ttk]
            EK <- nEvents*qK
            Ik <- Ek*P1*(1-P1)
            IK <- EK*P1*(1-P1)
            tk <- get.time(pks[ttk], N, B, L, td, EK, lam1c, HR1, HR2)
            thetak.H1 <- log(get.AHR(N, B, L, td, tk, lam1c, HR1, HR2))
            tK <- get.time(qK, N, B, L, td, EK, lam1c, HR1, HR2)
            thetaK.H1 <- log(get.AHR(N, B, L, td, tK, lam1c, HR1, HR2))
            Nk <- N/B*min(tk,B)
            NK <- N/B*min(tK,B)

            HR.H0 <- exp(theta.H0)
            Ek0 <- Ek; EK0 <- EK
            tk0 <- get.time(pks[ttk], N, B, L, td, EK, lam1c, HR.H0, HR.H0)
            tK0 <- get.time(qK, N, B, L, td, EK, lam1c, HR.H0, HR.H0)
            Nk0 <- N/B*min(tk0,B); NK0 <- N/B*min(tK0,B)

            for (gk in c(1:length(gammas))) {
                    H1k <- PF1.new(Ik, IK, thetak.H1, thetaK.H1, alpha, 0, theta.cp, gammas[gk])
                    H1.AE <- H1k$pr.S*Ek +(1-H1k$pr.S)*EK
                    H1.AN <- H1k$pr.S*Nk +(1-H1k$pr.S)*NK
                    H1.AL <- H1k$pr.S*tk +(1-H1k$pr.S)*tK

                    H0k <- PF1.new(Ik, IK, theta.H0, theta.H0, alpha, 0, theta.cp, gammas[gk])
                    H0.AE <- H0k$pr.S*Ek0+(1-H0k$pr.S)*EK0
                    H0.AN <- H0k$pr.S*Nk0+(1-H0k$pr.S)*NK0
                    H0.AL <- H0k$pr.S*tk0+(1-H0k$pr.S)*tK0

                    ##admissible rules
                    rule.id <- rule.id+1
                    a.rule <- ifelse(H0k$pr.S>=0.50 && ((H1k$pow.nbind-H1k$pow.bind)<=0.03), 1, 0)
                    ##a.rule <- ifelse(H0k$pr.S>=0.50 && (H1k$pow.nbind-H1k$pow.bind)<=0.031, 1, 0)
                    ##practical rules
                    ##p.rule <- ifelse(a.rule==1 && pks[ttk]>=0.20 && pks[ttk]<=0.80, 1, 0)
                    p.rule <- a.rule
                    tmp <- c(rule.id, a.rule, p.rule, pks[ttk], qK, tk, tK, exp(thetak.H1), exp(thetaK.H1), gammas[gk],
                    H0k$pr.S, H1k$pr.S, H0k$cpk, H1k$cpk, H1k$pow.nbind, H1k$pow.bind, H1k$pow.nbind-H1k$pow.bind, H0.AE, H1.AE, H0.AN, H1.AN, H0.AL, H1.AL)
                    res <- rbind(res, tmp)
            }}
    colnames(res) <- c("rule.id", "a.rule", "p.rule", "qk", "qK", "tk", "tK", "AHRk.H1", "AHRK.H1", "gk", "H0.pr.S", "H1.pr.S", "H0.cpk", "H1.cpk", "H1.pow.nbind", "H1.pow.bind", "H1.powloss", "H0.AE", "H1.AE", "H0.AN", "H1.AN",  "H0.AL", "H1.AL")
    res.a.rules <- data.frame(res)
    write.csv(res.a.rules,'look1_admissible_rules.csv')
    res.p.rules <- res.a.rules[res.a.rules$p.rule==1,]
    write.csv(res.p.rules,'look1_practical_rules.csv')
    r1 <- res.p.rules[which.min(res.p.rules$H0.AE),]
    r2 <- res.p.rules[which.min(res.p.rules$H0.AN),]
    r3 <- res.p.rules[which.min(res.p.rules$H0.AL),]
    res.o.rules <- rbind(r1,r2,r3)
    rownames(res.o.rules) <- c("CP1AE", "CP1AN", "CP1AL")
    write.csv(res.o.rules,'look1_optimal_rules.csv')
    o <- res.o.rules; o$rule <- rownames(o)
    res <- data.frame(o$rule, N, B, L, nEvents, o$gk, o$qk, o$qK, o$tk, o$tK, o$AHRk.H1, o$AHRK.H1,
    o$H1.pow.nbind, o$H1.pow.bind, o$H1.pow.nbind-o$H1.pow.bind,
    o$H1.pr.S, H1.cpk=o$H1.cpk, o$H1.AE, o$H1.AN, o$H1.AL,
    o$H0.pr.S, H0.cpk=o$H0.cpk, o$H0.AE, o$H0.AN, o$H0.AL
)
    colnames(res) <- c("rule", "N", "B", "L", "nEvents", "gk", "qk", "qK", "tk", "tK", "AHRk.H1", "AHRK.H1",
    "H1.pow.nbind", "H1.pow.bind", "H1.powloss",
    "H1.pr.S", "H1.cpk", "H1.AE", "H1.AN", "H1.AL",
    "H0.pr.S", "H0.cpk", "H0.AE", "H0.AN", "H0.AL")
    res
    }

################################################################################
opt.look2 <- function(pjs, pks, gammas, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha) {
    cat("two stops at pjs and pks !!\n")
    res <- NULL
    P1 <- r/(1+r)
    qK <- 1
    rule.id <- 1
    for (ttj in 1:length(pjs)) {
        for (ttk in c(ttj:length(pks))) {
            Ej <- nEvents*pjs[ttj]
            Ek <- nEvents*pks[ttk]
            EK <- nEvents*qK
            Ij <- Ej*P1*(1-P1)
            Ik <- Ek*P1*(1-P1)
            IK <- EK*P1*(1-P1)
            tj <- get.time(pjs[ttj], N, B, L, td, EK, lam1c, HR1, HR2)
            thetaj.H1 <- log(get.AHR(N, B, L, td, tj, lam1c, HR1, HR2))
            tk <- get.time(pks[ttk], N, B, L, td, EK, lam1c, HR1, HR2)
            thetak.H1 <- log(get.AHR(N, B, L, td, tk, lam1c, HR1, HR2))
            tK <- get.time(qK, N, B, L, td, EK, lam1c, HR1, HR2)
            thetaK.H1 <- log(get.AHR(N, B, L, td, tK, lam1c, HR1, HR2))
            Nj <- N/B*min(tj,B)
            Nk <- N/B*min(tk,B)
            NK <- N/B*min(tK,B)

            HR.H0 <- exp(theta.H0)
            Ej0 <- Ej; Ek0 <- Ek; EK0 <- EK
            Ij0 <- Ej0*P1*(1-P1); Ik0 <- Ek0*P1*(1-P1); IK0 <- EK0*P1*(1-P1)
            tj0 <- get.time(pjs[ttj], N, B, L, td, EK, lam1c, HR.H0, HR.H0)
            tk0 <- get.time(pks[ttk], N, B, L, td, EK, lam1c, HR.H0, HR.H0)
            tK0 <- get.time(qK, N, B, L, td, EK, lam1c, HR.H0, HR.H0)
            Nj0 <- N/B*min(tj0,B); Nk0 <- N/B*min(tk0,B); NK0 <- N/B*min(tK0,B)

            for (tgj in c(1:length(gammas))) {
                for (tgk in c(1:length(gammas))) {
                    H1jk <- PF2.new(Ij, Ik, IK, thetaj.H1, thetak.H1, thetaK.H1, alpha, 0, 0, theta.cp, gammas[tgj], gammas[tgk])
                    H1.AE <- H1jk$pr.S1*Ej +H1jk$pr.S2cS1*Ek +(1-H1jk$pr.S1-H1jk$pr.S2cS1)*EK
                    H1.AN <- H1jk$pr.S1*Nj +H1jk$pr.S2cS1*Nk +(1-H1jk$pr.S1-H1jk$pr.S2cS1)*NK
                    H1.AL <- H1jk$pr.S1*tj +H1jk$pr.S2cS1*tk +(1-H1jk$pr.S1-H1jk$pr.S2cS1)*tK

                    H0jk <- PF2.new(Ij0, Ik0, IK0, theta.H0,  theta.H0,  theta.H0,  alpha, 0, 0, theta.cp, gammas[tgj], gammas[tgk])
                    H0.AE <- H0jk$pr.S1*Ej0+H0jk$pr.S2cS1*Ek0+(1-H0jk$pr.S1-H0jk$pr.S2cS1)*EK0
                    H0.AN <- H0jk$pr.S1*Nj0+H0jk$pr.S2cS1*Nk0+(1-H0jk$pr.S1-H0jk$pr.S2cS1)*NK0
                    H0.AL <- H0jk$pr.S1*tj0+H0jk$pr.S2cS1*tk0+(1-H0jk$pr.S1-H0jk$pr.S2cS1)*tK0

                    ##admissible rules
                    rule.id <- rule.id+1
                    ##a.rule <- ifelse((H0jk$pr.S>=0.50 && (H1jk$pow.nbind-H1jk$pow.bind)<=0.03), 1, 0)
                    a.rule <- ifelse((H0jk$pr.S>=0.50 && (H1jk$pow.nbind-H1jk$pow.bind)<=0.03), 1, 0)
                    ##practical rules
                    ##p.rule <- ifelse((a.rule==1 && pjs[ttj]>=0.20 && pks[ttk]<=0.80 && pks[ttk]-pjs[ttj]>=0.20 && H0jk$pr.S1>=0.10), 1, 0)
                    p.rule <- a.rule
                    if (a.rule==1) res <- rbind(res, c(rule.id, a.rule, p.rule, pjs[ttj], pks[ttk], qK, tj, tk, tK, exp(thetaj.H1), exp(thetak.H1), exp(thetaK.H1),
                    gammas[tgj], gammas[tgk], H0jk$pr.S1, H0jk$pr.S, H0jk$cpj, H0jk$cpk, H1jk$pr.S1, H1jk$pr.S, H0jk$cpj, H0jk$cpk, H1jk$pow.nbind, H1jk$pow.bind, H1jk$pow.nbind-H1jk$pow.bind, H0.AE, H1.AE, H0.AN, H1.AN, H0.AL, H1.AL))
            }}
            }}
    colnames(res) <- c("rule.id", "a.rule", "p.rule", "qj", "qk", "qK", "tj", "tk", "tK", "AHRj.H1", "AHRk.H1", "AHRK.H1", "gj", "gk",
    "H0.pr.S1", "H0.pr.S", "H0.cpj", "H0.cpk",
    "H1.pr.S1", "H1.pr.S", "H1.cpj", "H1.cpk",
    "H1.pow.nbind", "H1.pow.bind", "H1.powloss", "H0.AE", "H1.AE", "H0.AN", "H1.AN",  "H0.AL", "H1.AL")
    res.a.rules <- data.frame(res)
    write.csv(res.a.rules,'look2_admissible_rules.csv')
    res.p.rules <- res.a.rules[res.a.rules$p.rule==1,]
    write.csv(res.p.rules,'look2_practical_rules.csv')
    r1 <- res.p.rules[which.min(res.p.rules$H0.AE),]
    r2 <- res.p.rules[which.min(res.p.rules$H0.AN),]
    r3 <- res.p.rules[which.min(res.p.rules$H0.AL),]
    res.o.rules <- rbind(r1,r2,r3)
    rownames(res.o.rules) <- c("CP2AE", "CP2AN", "CP2AL")
    write.csv(res.o.rules,'look2_optimal_rules.csv')
    o <- res.o.rules; o$rule <- rownames(o)
    res <- data.frame(o$rule, N, B, L, nEvents, o$gj, o$gk, o$qj, o$qk, o$qK, o$tj, o$tk, o$tK,
    o$AHRj.H1, o$AHRk.H1, o$AHRK.H1,
    o$H1.pow.nbind, o$H1.pow.bind, o$H1.pow.nbind-o$H1.pow.bind,
    o$H1.pr.S1, o$H1.pr.S, o$H1.cpj, o$H1.cpk, o$H1.AE, o$H1.AN, o$H1.AL,
    o$H0.pr.S1, o$H0.pr.S, o$H0.cpj, o$H0.cpk, o$H0.AE, o$H0.AN, o$H0.AL)
    colnames(res) <- c("rule", "N", "B", "L", "nEvents", "gj", "gk", "qj", "qk", "qK", "tj", "tk", "tK",
    "AHRj.H1", "AHRk.H1", "AHRK.H1",
    "H1.pow.nbind", "H1.pow.bind", "H1.powloss",
    "H1.pr.S1", "H1.pr.S", "H1.cpj", "H1.cpk", "H1.AE", "H1.AN", "H1.AL",
    "H0.pr.S1", "H0.pr.S", "H0.cpj", "H0.cpk", "H0.AE", "H0.AN", "H0.AL")
    res
}



###############################################################################
bounds <- function(Ik, IK, thetak, thetaK, alpha, etak, theta.cp, gk=NULL) {
    if (!is.null(gk)) etak <- (-qnorm(gk)*sqrt(IK-Ik)-qnorm(1-alpha)*sqrt(IK)-theta.cp*(IK-Ik))/sqrt(Ik)
    zk <- etak
    cpk <- pnorm((-zk*sqrt(Ik)-qnorm(1-alpha)*sqrt(IK)-theta.cp*(IK-Ik))/sqrt(IK-Ik))
    ppk <- pnorm((-zk*sqrt(Ik)-qnorm(1-alpha)*sqrt(Ik))/sqrt(IK-Ik))
    ck <- -etak+thetak*sqrt(Ik)
    pr.S <- pnorm(ck)
    ahrk <- exp(-ck/sqrt(Ik))
    cK <- qnorm(1-alpha)+thetaK*sqrt(IK)
    pow.nbind <- 1-pnorm(cK)
    sigma <- diag(2)
    sigma[2,1] <- sigma[1,2] <- sqrt(Ik/IK)
    pow.bind <- pmvnorm(lower=c(ck, cK), upper=c(Inf, Inf), mean=c(0,0), sigma=sigma)
    return(list(qk=Ik/IK, cpk=cpk, ahrk=ahrk, etak=etak, ppk=ppk, pr.S=pr.S, pow.nbind=pow.nbind, pow.bind=pow.bind, powloss=pow.nbind-pow.bind))
}

###############################################################################
get.bounds <- function(rule, qk, etak, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha, gk=NULL) {
    P1 <- r/(1+r)
    qK <- 1

    Ek <- nEvents*qk
    EK <- nEvents*qK
    Ik <- Ek*P1*(1-P1)
    IK <- EK*P1*(1-P1)
    tk <- get.time(qk, N, B, L, td, EK, lam1c, HR1, HR2)
    thetak.H1 <- log(get.AHR(N, B, L, td, tk, lam1c, HR1, HR2))
    tK <- get.time(qK, N, B, L, td, EK, lam1c, HR1, HR2)
    Nk <- N/B*min(tk,B)
    NK <- N/B*min(tK,B)
    thetaK.H1 <- log(get.AHR(N, B, L, td, tK, lam1c, HR1, HR2))
    H1k <- bounds(Ik, IK, thetak.H1, thetaK.H1, alpha, etak, theta.cp, gk)

    HR.H0 <- exp(theta.H0)
    Ek0 <- Ek; EK0 <- EK
    Ik0 <- Ek0*P1*(1-P1); IK0 <- EK0*P1*(1-P1)
    tk0 <- get.time(qk, N, B, L, td, EK, lam1c, HR.H0, HR.H0)
    tK0 <- get.time(qK, N, B, L, td, EK, lam1c, HR.H0, HR.H0)
    Nk0 <- N/B*min(tk0,B); NK0 <- N/B*min(tK0,B)
    H0k <- bounds(Ik0, IK0, theta.H0, theta.H0, alpha, etak, theta.cp, gk)

    H0.AE <- H0k$pr.S*Ek0+(1-H0k$pr.S)*EK0
    H0.AN <- H0k$pr.S*Nk0+(1-H0k$pr.S)*NK0
    H0.AL <- H0k$pr.S*tk0+(1-H0k$pr.S)*tK0
    res <- data.frame(rule, N, B, L, nEvents, tk, tK, AHRk.H1=exp(thetak.H1), AHRK.H1=exp(thetaK.H1), H0.AE=H0.AE, H0.AN=H0.AN, H0.AL=H0.AL, H0AER=H0.AE/nEvents, H0.ANR=H0.AN/N, H0ALR=H0.AL/L, qk=H0k$qk, cpk=H0k$cpk, ahrk=H0k$ahrk, etak=H0k$etak, ppk=H0k$ppk, pr.S=H0k$pr.S, pow.nbind=H1k$pow.nbind, pow.bind=H1k$pow.bind, powloss=H1k$powloss)
}

################################################################################
get.bounds.summary <- function(S, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha) {
r1 <- get.bounds(rule="Wieand1", 0.50, 0, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha, gk=NULL)
r2 <- get.bounds(rule="KF1",     0.50, 0, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha, gk=NULL)
r3 <- get.bounds(rule="OF1",     1/3, -0.011, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha, gk=NULL)
r4 <- get.bounds(rule="Xi1",     0.40, -0.49, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha, gk=NULL)
r5 <- get.bounds(rule="CP1AE.PH(0.75)", 0.42, 0, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha, gk=0.64)
r6 <- get.bounds(rule="CP1AE.NPH(3,1.0,0.6938)", 0.47, 0, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha, gk=0.50)
r7 <- get.bounds(rule="CP1AE.NPH(3,1.3,0.6295)", 0.55, 0, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha, gk=0.31)
r8 <- get.bounds(rule="CP1AE.NPH(6,1.0,0.623)", 0.54, 0, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha, gk=0.34)
res <- bind_rows(r1, r2, r3, r4, r5, r6, r7, r8)
res <- cbind(S=S, res)
}

################################################################################
table3 <- function(design_name, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha, look=1) {
    ##res0 <- Wieand.look1(qk=0.5, etak=Inf, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha)
    res1a <- Wieand.look1(qk=0.5, etak=0, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha)
    res1b <- Wieand.look1(qk=1/3, etak=-0.011, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha)
    res1c <- KF.look1(qk=0.5, etak=0, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha)
    res1d <- Wieand.look1(qk=0.40, etak=-0.49, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha)
    gammas <- seq(0.05, 0.85, by=0.05)
    pks <- seq(0.1, 0.9, by=0.05)
    res1e <- opt.look1(pks, gammas, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha)
    res <- bind_rows(res1e, res1d, res1c, res1b, res1a)

if (look==2) {
    res2a <- Wieand.look2(qj=0.5, qk=0.75, etaj=0, etak=0, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha)
    res2b <- Wieand.look2(qj=1/3, qk=2/3, etaj=-0.011, etak=-0.864, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha)
    res2c <- KF.look2(qj=0.5, qk=0.75, etaj=0, etak=0, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha)
    res2d <- Wieand.look2(qj=0.30, qk=0.60, etaj=-0.02, etak=-0.78, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha)
    gammas <- seq(0.05, 0.85, by=0.05)
    pjs <- seq(0.1, 0.9, by=0.05)
    pks <- pjs+0.05
    res2e <- opt.look2(pjs, pks, gammas, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha)
    res <- bind_rows(res2e, res2d, res2c, res2b, res2a)
}
res
}

################################################################################
## Create the reference design
get.REFdesign <- function(alpha=0.025, power=0.90, sided=1, r=1, lamc=-log(0.5)/12, HR.H1=0.75, B=34, gamma=20) {
    theta.H1 <- log(HR.H1)
    theta.H0 <- log(1)
    nEvents <- ceiling(ne(r, alpha, 1-power, exp(theta.H1)))
    nEvents(hr=exp(theta.H1), alpha=alpha, beta=1-power, ratio=r, sided=sided, hr0=exp(theta.H0))
    AHR <- exp(theta.H1)
    theta.cp <- theta.H1
    N <- B*gamma
    r0 <- nSurv(lambdaC=lamc, hr=exp(theta.H1), hr0=exp(theta.H0), R=B, gamma=N/B, alpha=alpha, beta=1-power, sided=sided)
    res <- cbind(r0$n, r0$d, r0$B, r0$T, 1-r0$beta, r0$hr)
}

################################################################################
## Get Table 2 in the manuscript
get.table2 <- function(HR.H1=0.75, B=34, N=680, nEvents=508) {
## fixed parameters
    alpha <- 0.025
    beta <- 0.10
    sided <- 1
    r <- 1
    lamc <- -log(0.5)/12
    theta.H1 <- log(HR.H1)
    theta.H0 <- log(1)
    nEvents0 <- ceiling(ne(r, alpha, beta, exp(theta.H1)))
    pow0 <- pow(r, nEvents, alpha, exp(theta.H1))

    AHR <- exp(theta.H1)
    theta.cp <- theta.H1

    ## We will run a PH trial with N, nEvents, etc
    print(S3 <- paste0("PH(", HR.H1, ")"))
    lam1c <- lamc; HR <- HR.H1
    L <- B+find.fu.PH(N, B, lamc, HR, nEvents, alpha)
    HR1 <- HR2 <- HR
    s3 <- get.bounds.summary(S3, N, B, L, 3, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha)

    ## We will run a NPH trial with N, nEvents, etc
    lam1c <- lamc; HR1 <- 1; td <- 3
    HR2 <- get.HR2(N, B, td, lam1c, HR1, AHR, pow0, alpha)
    print(S4 <- paste0("NPH(3,1,",round(HR2,4),")"))
    L <- B+find.fu(N, B, td, lam1c, HR1, HR2, pow0, alpha)
    nEvents_check <- find.ne.new(N, td, B, L, lam1c, HR1, HR2)$e
    s4 <- get.bounds.summary(S4, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha)

    ## We will run a NPH trial with N, nEvents, and theta.A
    lam1c <- lamc; HR1 <- 1.3; td <- 3
    HR2 <- get.HR2(N, B, td, lam1c, HR1, AHR, pow0, alpha)
    print(S5 <- paste0("NPH(3,1.3,",round(HR2,4),")"))
    L <- B+find.fu(N, B, td, lam1c, HR1, HR2, pow0, alpha)
    nEvents_check <- find.ne.new(N, td, B, L, lam1c, HR1, HR2)$e
    s5 <- get.bounds.summary(S5, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha)

    ## We will run a NPH trial with N, nEvents, and theta.A
    lam1c <- lamc; HR1 <- 1; td <- 6
    HR2 <- get.HR2(N, B, td, lam1c, HR1, AHR, pow0, alpha)
    print(S6 <- paste0("NPH(6,1,",round(HR2,4),")"))
    L <- B+find.fu(N, B, td, lam1c, HR1, HR2, pow0, alpha)
    nEvents_check <- find.ne.new(N, td, B, L, lam1c, HR1, HR2)$e
    s6 <- get.bounds.summary(S6, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha)
    res <- bind_rows(s3, s4, s5, s6)
}

################################################################################
## Get Table 3 in the manuscript
get.table3 <- function(HR.H1=0.75, B=34, N=680, nEvents=508) {
## fixed parameters
    alpha <- 0.025
    beta <- 0.10
    sided <- 1
    r <- 1
    lamc <- -log(0.5)/12
    theta.H1 <- log(HR.H1)
    theta.H0 <- log(1)
    nEvents0 <- ceiling(ne(r, alpha, beta, exp(theta.H1)))
    pow0 <- pow(r, nEvents, alpha, exp(theta.H1))

    AHR <- exp(theta.H1)
    theta.cp <- theta.H1

    ## We will run a PH trial with N, nEvents, etc
    print(S3 <- paste0("PH(", HR.H1, ")"))
    lam1c <- lamc; HR <- HR.H1
    L <- B+find.fu.PH(N, B, lamc, HR, nEvents, alpha)
    HR1 <- HR2 <- HR
    s3 <- table3(S3, N, B, L, 3, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha)

    ## We will run a NPH trial with N, nEvents, etc
    lam1c <- lamc; HR1 <- 1; td <- 3
    HR2 <- get.HR2(N, B, td, lam1c, HR1, AHR, pow0, alpha)
    print(S4 <- paste0("NPH(3,1,",round(HR2,4),")"))
    L <- B+find.fu(N, B, td, lam1c, HR1, HR2, pow0, alpha)
    nEvents_check <- find.ne.new(N, td, B, L, lam1c, HR1, HR2)$e
    s4 <- table3(S4, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha)

    ## We will run a NPH trial with N, nEvents, and theta.A
    lam1c <- lamc; HR1 <- 1.3; td <- 3
    HR2 <- get.HR2(N, B, td, lam1c, HR1, AHR, pow0, alpha)
    print(S5 <- paste0("NPH(3,1.3,",round(HR2,4),")"))
    L <- B+find.fu(N, B, td, lam1c, HR1, HR2, pow0, alpha)
    nEvents_check <- find.ne.new(N, td, B, L, lam1c, HR1, HR2)$e
    s5 <- table3(S5, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha)


    ## We will run a NPH trial with N, nEvents, and theta.A
    lam1c <- lamc; HR1 <- 1; td <- 6
    HR2 <- get.HR2(N, B, td, lam1c, HR1, AHR, pow0, alpha)
    print(S6 <- paste0("NPH(6,1,",round(HR2,4),")"))
    L <- B+find.fu(N, B, td, lam1c, HR1, HR2, pow0, alpha)
    nEvents_check <- find.ne.new(N, td, B, L, lam1c, HR1, HR2)$e
    s6 <- table3(S6, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha)

    res <- bind_rows(s3, s4, s5, s6)
}


################################################################################
## Get Table 4 in the manuscript
get.table4 <- function(HR.H1=0.75, B=34, N=608, nEvents=508, look=2) {
## fixed parameters
    alpha <- 0.025
    beta <- 0.10
    sided <- 1
    r <- 1
    lamc <- -log(0.5)/12
    theta.H1 <- log(HR.H1)
    theta.H0 <- log(1)
    nEvents0 <- ceiling(ne(r, alpha, beta, exp(theta.H1)))
    pow0 <- pow(r, nEvents, alpha, exp(theta.H1))

    AHR <- exp(theta.H1)
    theta.cp <- theta.H1

    ## We will run a NPH trial with N, nEvents, etc
    lam1c <- lamc; HR1 <- 1; td <- 3
    HR2 <- get.HR2(N, B, td, lam1c, HR1, AHR, pow0, alpha)
    print(S4 <- paste0("NPH(3,1,",round(HR2,4),")"))
    L <- B+find.fu(N, B, td, lam1c, HR1, HR2, pow0, alpha)
    nEvents_check <- find.ne.new(N, td, B, L, lam1c, HR1, HR2)$e
    s4 <- table3(S4, N, B, L, 3, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha, look)
}

###############################################################################################
plot.fig1 <- function(res, file) {
    g_legend<-function(a.gplot){
      tmp <- ggplot_gtable(ggplot_build(a.gplot))
      leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
      legend <- tmp$grobs[[leg]]
      return(legend)
    }

        fa=ggplot(res,aes(x=pk,y=AHRk,linetype=factor(design), group=factor(design)))+
          xlab("Information")+
          ylab("AHR")+
          ggtitle("(a) AHR vs. Information")+
          geom_line()+
          scale_x_continuous(limits=c(0,1))+
          scale_y_continuous(limits=c(0.5,1.3))+
          scale_color_discrete(name="design")+
          theme(legend.position="bottom")+
          theme(legend.text=element_text(size=9.5, face="bold"))+
          theme(legend.title=element_blank())

        fb=ggplot(res,aes(x=pk,y=powk,linetype=factor(design), group=factor(design)))+
          xlab("Information")+
          ylab("Power")+
          ggtitle("(b) Power vs. Information")+
          geom_line()+
          scale_x_continuous(limits=c(0,1))+
          scale_y_continuous(limits=c(0,1))+
          scale_color_discrete(name="design")+
          theme(legend.position="bottom")+
          theme(legend.text=element_text(size=9.5, face="bold"))+
          theme(legend.title=element_blank())

        fc=ggplot(res,aes(x=tk,y=AHRk,linetype=factor(design), group=factor(design)))+
          xlab("Time (months)")+
          ylab("AHR")+
          ggtitle("(c) AHR vs. Calendar Time")+
          geom_line()+
          scale_x_continuous(limits=c(0,50))+
          scale_y_continuous(limits=c(0.5,1.3))+
          scale_color_discrete(name="design")+
          theme(legend.position="bottom")+
          theme(legend.text=element_text(size=9.5, face="bold"))+
          theme(legend.title=element_blank())

        fd=ggplot(res,aes(x=tk,y=powk,linetype=factor(design), group=factor(design)))+
          xlab("Time (months)")+
          ylab("Power")+
          ggtitle("(d) Power vs. Calendar Time")+
          geom_line()+
          scale_x_continuous(limits=c(0,50))+
          scale_y_continuous(limits=c(0,1))+
          scale_color_discrete(name="design")+
          theme(legend.position="bottom")+
          theme(legend.text=element_text(size=9.5, face="bold"))+
          theme(legend.title=element_blank())

          ml <- g_legend(fa)

        pdf(file, width=13, height=11)
        grid.arrange(fa, fb, fc, fd, ncol=2, nrow=2)
        dev.off()
}

###############################################################################################
get.fig1.data <- function(design, N, B, L, td, nEvents, lam1c, HR1, HR2, r, theta.H0, theta.cp, alpha) {
pks <- seq(0.0001, 1, by=0.01)
res <- data.frame()
for (pk in pks) {
    tk <- get.time(pk, N, B, L, td, nEvents, lam1c, HR1, HR2)
    tx <- pk*L
    nek <- find.ne.new(N, td, B, tk, lam1c, HR1, HR2)$e
    Sk <- get.St(td, tx, lam1c, HR1, HR2)
    AHRk <- get.AHR(N, B, L, td, tk, lam1c, HR1, HR2)
    powk <- pow(r, nek, alpha, AHRk)
    resk <- data.frame(design, pk, tk, tx, Sk, nek, AHRk, powk, L=ceiling(L))
    res <- rbind(res, resk)
}
return(res)
}

###############################################################################################
###############################################################################################
plot.figS1 <- function(res0, res1, res2, res3, file) {
    g_legend<-function(a.gplot){
      tmp <- ggplot_gtable(ggplot_build(a.gplot))
      leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
      legend <- tmp$grobs[[leg]]
      return(legend)
    }

        f0=ggplot(res0,aes(x=tx,y=Sk,linetype=factor(A), group=factor(A)))+
          xlab("Follow-up Time (months)")+
          ylab("Survival Prob.")+
          ggtitle(paste0("(a) ", res0$design[1]))+
          geom_line()+
          scale_x_continuous(limits=c(0,50))+
          scale_y_continuous(limits=c(0,1))+
          scale_linetype_discrete(name = "Treatment")+
          theme(legend.position="bottom")

        f1=ggplot(res1,aes(x=tx,y=Sk,linetype=factor(A), group=factor(A)))+
          xlab("Follow-up Time (months)")+
          ylab("Survival Prob.")+
          ggtitle(paste0("(b) ", res1$design[1]))+
          geom_line()+
          scale_x_continuous(limits=c(0,50))+
          scale_y_continuous(limits=c(0,1))+
          scale_linetype_discrete(name = "Treatment")+
          theme(legend.position="bottom")

        f2=ggplot(res2,aes(x=tx,y=Sk,linetype=factor(A), group=factor(A)))+
          xlab("Follow-up Time (months)")+
          ylab("Survival Prob.")+
          ggtitle(paste0("(c) ", res2$design[1]))+
          geom_line()+
          scale_x_continuous(limits=c(0,50))+
          scale_y_continuous(limits=c(0,1))+
          scale_linetype_discrete(name = "Treatment")+
          theme(legend.position="bottom")

        f3=ggplot(res3,aes(x=tx,y=Sk,linetype=factor(A), group=factor(A)))+
          xlab("Follow-up Time (months)")+
          ylab("Survival Prob.")+
          ggtitle(paste0("(d) ", res3$design[1]))+
          geom_line()+
          scale_x_continuous(limits=c(0,50))+
          scale_y_continuous(limits=c(0,1))+
          scale_linetype_discrete(name = "Treatment")+
          theme(legend.position="bottom")

          ml <- g_legend(f0)

        pdf(file, width=11, height=11)
        grid.arrange(f0, f1, f2, f3, ncol=2, nrow=2)
        dev.off()
}

###############################################################################################
get.figS1.data <- function(design, L, td, lam1c, HR1, HR2) {
    pks <- seq(0.01, 1, by=0.01)
    res <- data.frame()
    for (pk in pks) {
        tx <- pk*L
        S0 <- get.St(td, tx, lam1c, 1, 1)
        S1 <- get.St(td, tx, lam1c, HR1, HR2)
        resk0 <- data.frame(design=design, A=0, tx=tx, Sk=S0)
        resk1 <- data.frame(design=design, A=1, tx=tx, Sk=S1)
        res <- rbind(res, resk0, resk1)
    }
    res
}
