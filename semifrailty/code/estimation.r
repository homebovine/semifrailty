library(MASS)
library(BB)
library(pracma)
library(parallel)
##########General functions#########
G <- function(mu, bg, br){
    - log(sum(exp(bg * mu) * br))
}
dG <- function(mu, bg, br){
    sum(bg * exp(- bg * mu) * br) /sum(exp(bg * mu) * br)
}
ddG <- function(mu, bg, br){
    Q <- sum(exp(bg * mu) * br)
    dQ <- sum(-bg * exp(-bg * mu) * br)
    ddQ <- sum(bg^2 * exp(-bg * mu) * br)
    - (ddQ * Q - dQ^2)/(dQ^2)
}
iniA <- function(t){
    0.5 * t ^2
}
inia <- function(t){
     t
}
A <- iniA
a <- inia
########estimation for alpha############
########function for e2##############
#A, a functions for hazards, bg, chosen gamma values, br, the probability.

fb2 <- function(A, a, bb, bg, x, t){
    expbbx <- exp(x %*% bb)
    integrand1 <- function(s, i){
        As <- A(s)
        as <- a(s)
        mu <-  As *   expbbx
        ddG(mu, bg, br)/dG(mu, bg, br) * expbbx * 1 * exp(bg[i] * mu) * (bg[i] * as * expbbx- dG(mu, bg, br)* as *  expbbx ) 
     }

    integrand2 <- function(s, i){
        As <- A(s)
        as <- a(s)
        mu <-  As *   expbbx
         (ddG(mu, bg, br)/dG(mu, bg, br) * expbbx * (s- t) + 1/as) * exp(bg[i] * mu) * (bg[i] * as * expbbx- dG(mu, bg, br)* as *  expbbx ) 
  }
    integrand1 <- Vectorize(integrand1)
    integrand2 <- Vectorize(integrand2)
    b2 <- function(i){
        c(integrate(integrand1, 0, tau, i)$value, integrate(integrand2, 0, tau, i)$value)
    }
    t(do.call(rbind, lapply(1 : m, b2)))
}
fA2 <- function(A, a, bb, bg, br, x){
    mA2 <- matrix(NA, m, m)
    expbbx <- exp(x %*% bb)
    integrand <- function(s, i, j){
        As <- A(s)
        as <- a(s)
        mu <-  As *  expbbx
        expbgbbx <- exp(bg * As * expbbx)
        sumexpc <-  sum(exp(bg * As * expbbx) *br )
        sumexprc <-  sum(exp(bg * As * expbbx) * bg *br )
        res <- (expbgbbx[j] * bg[j] * br[j]/sumexprc -  expbgbbx[j]  * br[j]/sumexpc)  * expbgbbx[i] * (bg[i] * as * expbbx- dG(mu, bg, br)*as*  expbbx )
        if(is.na(res)|abs(res) == Inf|is.nan(res)){
            browser()
        }
        res
    }
    vintegrand <- Vectorize(integrand)
    A2 <- function(ij){
        i <- ij[1]
        j <- ij[2]
        res <- try(integral(vintegrand, 0, tau,method = c("Kronrod"), vectorized = TRUE, arrayValued = FALSE, reltol = 1e-08, abstol = 0, i, j))
        if(class(res) == "try-error"){
            browser()
        }
        res
    }
    for(itr in 1 : (m^2)){
        ij <- mij[itr, ]
        mA2[ij[1], ij[2]] <- A2(ij)
    }
    return(mA2)
     
}

gete2 <- function(i, Afun, afun, bb, bg, br, t, data){
    x <- data[i, 3:ncol(data)]
    A2 <- fA2(Afun, afun, bb, bg, br, x)
    b2 <- fb2(Afun, afun, bb, bg, x, t)
    e2 <- t(ginv(t(A2)%*% A2) %*% (t(A2) %*% t(b2)))
}


scorealpha <- function(i, data,  ba, bb, Afun, afun, bg, br, t, me2){
    d <- data[i, 1]
    y <- data[i, 2]
    x <- data[i, 3: ncol(data)]
    expbbx <- exp(x %*% bb)
    Aexpbbx <- ba[1] *expbbx
    a <- ba[2]
    A <- ( ba[1] + ba[2]* (y - t)) 
    mu <- ( ba[1] + ba[2]* (y - t)) *   expbbx
    expbgbbx <- exp(bg * A *  expbbx)
    sumexpc <-  sum(exp(bg * A * expbbx) *br )
    sumexprc <-  sum(exp(bg  * A *  expbbx) *bg * br )
    e2 <- me2[[i]]
    numf <- apply(cbind(e2[1, ] * bg *  expbgbbx * br, e2[2, ] * bg * expbgbbx * br), 2, sum)
    nums <- apply(cbind(e2[1, ] *   expbgbbx * br, e2[2, ] *  expbgbbx * br), 2, sum)
    if(d == 1){
       res <-  ddG(mu, bg, br)/dG(mu, bg, br)  * expbbx * c(1, y - t) + c(0, 1)/ a- ((numf/ sumexprc) - nums/(sumexpc)) - dG(mu, bg, br) * expbbx * c(1, y - t) -  nums/(sumexpc)
    }else{
        res <- - dG(mu, bg, br) * expbbx * c(1, y - t) -  nums/(sumexpc)
    }
    t(res) %*% diag(dnorm((y - t)/h), 2, 2)%*% res
}
scorealphasum <- function(ba, bb, data, A, a, bg, br, t, me2){
 #   print(ba)
    ba <- c((ba[1]), ba[2])
   # apply(do.call(rbind, lapply(1 :n, data, 1, scorealpha, ba, bb, A, a, bg, c, t)), 2, mean, na.rm = T)
    mean(do.call(rbind, lapply(1 :n, scorealpha, data,   ba, bb, A, a, bg, br, t, me2)))
}

tryscorealphasum <- function(ba, bb, data, A, a, bg, br, t, me2){
    res <- try(scorealphasum(ba, bb, data, A, a, bg, br, t, me2))
    if(class(res) == "try-error")
        browser()
    res
}
##get ba function
fmb <- function(i, mt, bb, A, a, bg, br){
    t <- mt[i]
    me2 <- mclapply(1 : n, gete2, A, a, bb, bg, br, t, data, mc.cores = 2)
    if(i != 1){
    spg(c(A(t), a(t)), tryscorealphasum, gr=NULL, method=3, lower=c(A(mt[i-1]), a(mt[i-1])), upper=Inf, project=NULL, projectArgs=NULL, control=list(), quiet=FALSE,  bb, data, A, a, bg, br, t, me2)$par
}else{
    spg(c(A(t), a(t)), tryscorealphasum, gr=NULL, method=3, lower=c(0, 0), upper=Inf, project=NULL, projectArgs=NULL, control=list(), quiet=FALSE,  bb, data, A, a, bg, br, t, me2)$par
}
}


######function for bb###############
fb1 <- function(A, a, bb, bg, x){
    expbbx <- exp(x %*% bb)
    integrand <- function(s, i, x){
        mu <-  A(s) *   expbbx
        (ddG(mu, bg, br)/dG(mu, bg, c) * mu * x + x) * exp(bg[i] * mu) * (bg[i] * a(s) * expbbx- dG(mu, bg, br)* a(s)*  expbbx ) 
     }

   
    integrand <- Vectorize(integrand)
    b1 <- function(i){
        res <- 0
        for(q in 1 : p){
        res <- c(res, integrate(integrand, 0, tau, i, x[q])$value)
    }
        res <- res[-1]
    }
    t(do.call(rbind, lapply(1 : m, b1)))
}
fA1 <- function(A, a, bb, bg, br, x){
    mA1 <- matrix(NA, m, m)
    expbbx <- exp(x %*% bb)
    integrand <- function(s, i, j){
        mu <-  A(s) *  expbbx
        expbgbbx <- exp(bg * A(s) * expbbx)
        sumexpc <-  sum(exp(bg * A(s) * expbbx) *c )
        sumexprc <-  sum(exp(bg * A(s) * expbbx) * bg *c )
        res <- (expbgbbx[j] * bg[j] * br[j]/sumexprc -  expbgbbx[j]  * br[j]/sumexpc)  * expbgbbx[i] * (bg[i] * a(s) * expbbx- dG(mu, bg, br)* a(s)*  expbbx )
        if(is.an(res)|abs(res) == Inf|is.nan(res)){
            browser()
        }
    }
    vintegrand <- Vectorize(integrand)
    A1 <- function(ij){
        i <- ij[1]
        j <- ij[2]
        res <- try(integrate(vintegrand, 0, tau, i, j)$value)
        if(class(res) == "try-error"){
            browser()
        }
        res
    }
    for(itr in 1 : (m^2)){
        ij <- mij[itr, ]
        mA1[ij[1], ij[2]] <- A1(ij)
    }
    return(mA1)
     
}

gete1 <- function(i, Afun, afun, bb, bg, br,  data){
    x <- data[i, 3:ncol(data)]
    A1 <- fA1(Afun, afun, bb, bg, br, x)
    b1 <- fb1(Afun, afun, bb, bg, x)
    e2 <- t(ginv(t(A1)%*% A1) %*% (t(A1) %*% t(b1)))
}

scorebeta <- function(i, data,   bb,   Afun, afun, bg, br,  me1){
    d <- data[i, 1]
    y <- data[i, 2]
    x <- data[i, 3: ncol(data)]
    expbbx <- exp(x %*% bb)
    a <- afun(y)#ba[2]
    A <- Afun(y)#( ba[1] + ba[2]* (y - t)) 
    mu <- A *   expbbx
    expbgbbx <- exp(bg * A *  expbbx)
    sumexpc <-  sum(exp(bg * A * expbbx) *br )
    sumexprc <-  sum(exp(bg  * A *  expbbx) *bg * br)
    e1 <- me1[[i]]
   
    numf <-  apply(matrix(apply(e1, 1, "*", bg *  expbgbbx * br), ncol = p), 2, sum)
    nums <-  apply(matrix(apply(e1, 1, "*",  expbgbbx * br), ncol = p), 2, sum)
    if(d == 1){
       res <-  (ddG(mu, bg, br)/dG(mu, bg, br)  * expbbx * A * x + x) - ((numf/ sumexprc) - nums/(sumexpc)) - dG(mu, bg, br) * A * expbbx * x-  nums/(sumexpc)
    }else{
        res <- - dG(mu, bg, br) * A * expbbx * x-  nums/(sumexpc)
    }
    t(res) %*% res
}
scorebetasum <- function( bb,  data, A, a, bg, br,  me1){
   
   # apply(do.call(rbind, lapply(1 :n, data, 1, scorealpha, ba, bb, A, a, bg, c, t)), 2, mean, na.rm = T)
    mean(do.call(rbind, lapply(1 :n, scorebeta, data,    bb,   A, a, bg, br,  me1)))
}

tryscorebetasum <- function( bb,  data, A, a, bg, br, me1){
    res <- try(scorebetasum( bb,  data, A, a, bg, br, me1))
    if(class(res) == "try-error")
        browser()
    res
}

fmbeta <- function( data, bb, A, a, bg, br){
    me1 <- lapply(1 : n, gete1, A, a, bb, bg, br,  data)

    spg(bb, tryscorebetasum, gr=NULL, method=3, lower= -Inf, upper=Inf, project=NULL, projectArgs=NULL, control=list(), quiet=FALSE,   data, A, a, bg, br,  me1)$par
}
#########simulation###############
simu <- function(n, v, bb,lambda,  alpha, tau){
   u <- runif(n)
   x <- matrix(rbeta(n, 0.5, 1), nrow = 100)
   rg <- rgamma(n, 1/v, 1/v)
   t <- (- log(u)/ (rg * lambda * exp(x%*% bb)))^(1/alpha)
   d <- t <= tau
   cbind(d, pmin(t, tau), x)
   
}
####test code###########
v <- 0.5
tau <- 1.6
data <- simu(100, v, bb, 1, 2, tau)
rg <- rgamma(10000, 1/v, 1/v)
m <- 10
bg <- as.vector(quantile(rg, seq(0.20, 0.80, length.out = m)))
br <- dgamma(bg, 1/v, 1/v)/ sum(dgamma(bg, 1/v, 1/v))
p <- ncol(x)
n <- nrow(data)
h <- bw.nrd(data[, 2]) * n^{-1/15}
mij <- as.matrix(expand.grid(1:m, 1:m))
lmt <- 22
mt <- seq(min(data[, 2]), max(data[data[, 1] == 1, 2]), length.out = lmt)
#mt <- data[, 2]
omba <- mba

for(itr in 1:100){
    mba<- do.call(rbind, lapply(1 :lmt, fmb, mt, bb, A, a, bg, br))
  #  A <- function(s) predict(smooth.spline(mt, mba[, 1]), s)$y
   # a <- function(s) predict(smooth.spline(mt, mba[, 2]), s)$y
    A <- approxfun(mt, mba[, 1], rule = 2)
     a <- approxfun(mt, mba[, 2], rule = 2)
    if(itr != 1 & (sum(abs(mba[, 1] - omba[, 1])) / lmt <= 1e-3) ){
        break
    }
    omba <- mba
}

#dfsane(ba,tryscorealphasum, method = 2, control = list(), quiet = FALSE, bb, data, A, a, bg, c, t)$par
for(itr  in 1:100){
mbeta <- fmbeta(data, bb, A, a, bg, br)
if(sum(abs(bb - mbeta))/p <= 1e-5){
    break
}
bb <- mbeta

}
