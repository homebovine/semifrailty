

library(MASS)
library(BB)
library(pracma)
library(parallel)
library(caTools)
##########General functions#########
G <- function(mu, bg, br){
    - log(sum(exp(-bg * mu) * br))
}
dG <- function(mu, bg, br){
    sum(bg * exp(- bg * mu) * br) /sum(exp(-bg * mu) * br)
}
ddG <- function(mu, bg, br){
    Q <- sum(exp(-bg * mu) * br)
    dQ <- sum(-bg * exp(-bg * mu) * br)
    ddQ <- sum(bg^2 * exp(-bg * mu) * br)
    res <- - (ddQ * Q - dQ^2)/(Q^2)
    ## if(is.nan(res)){
    ##     browser()
    ## }
    res
}
iniA <- function(t){
     0.5 * t^2 
}
inia <- function(t){
    t
}
A <- iniA
a <- inia
########estimation for alpha############
########function for e2##############
#A, a functions for hazards, bg, chosen gamma values, br, the probability.

fb2 <- function(A, a,  bb, bg, x, t){
    expbbx <- exp(x %*% bb)
    integrand1 <- function(s, i){
        As <- A(s)#A(t) + a(t) * (s - t)#A(s)
        as <- a(s)#a(t)#a(s)
        mu <-  As *   expbbx
        ddG(mu, bg, br)/dG(mu, bg, br) * expbbx * 1 * exp(-bg[i] * mu) * (bg[i] * as * expbbx- dG(mu, bg, br)* as *  expbbx ) 
     }

    integrand2 <- function(s, i){
        As <- A(t) + a(t)  * (s - t)#A(s)
        as <- a(t)#a(s)
        mu <-  As *   expbbx
         (ddG(mu, bg, br)/dG(mu, bg, br) * expbbx * (s- t) + 1/as) * exp(-bg[i] * mu) * (bg[i] * as * expbbx- dG(mu, bg, br)* as *  expbbx ) 
  }
    integrand1 <- Vectorize(integrand1)
    integrand2 <- Vectorize(integrand2)
    b2 <- function(i){
        c(myintegral(integrand1, 0.001, tau, i), myintegral(integrand2, 0.001, tau, i))
        #c(integral(integrand1, 0, tau,  method = "Kronrod", vectorized = T, arrayValued = F, waypoints = NULL, reltol = 1e-05, abstol = 1e-5, i), integral(integrand2, 0, tau,  method = "Kronrod", vectorized = T, arrayValued = F, waypoints = NULL, reltol = 1e-05, abstol = 1e-5, i))
    }
    t(do.call(rbind, lapply(1 : m, b2)))
}
fA2 <- function(A, a, bb, bg, br, x, t){
    mA2 <- matrix(NA, m, m)
    expbbx <- exp(x %*% bb)
    integrand <- function(s, i, j){
        As <- A(s)#A(t) + a(t)  * (s - t)#A(s)#ba[1] + ba[2](s - t)
        as <- a(s)#a(t)#a(s)#ba[2]
        mu <-  As *  expbbx
        expbgbbx <- exp(-bg * As * expbbx)
        sumexpc <-  sum(exp(-bg * As * expbbx) *br )
        sumexprc <-  sum(exp(-bg * As * expbbx) * bg *br )
        res <- (expbgbbx[j] * bg[j] * br[j]/sumexprc -  expbgbbx[j]  * br[j]/sumexpc)  * expbgbbx[i] * (bg[i] * as * expbbx- dG(mu, bg, br)*as*  expbbx )
        if(is.na(res)|abs(res) == Inf|is.nan(res)){
            #browser()
        }
        res
    }
    vintegrand <- Vectorize(integrand)
    A2 <- function(ij){
        i <- ij[1]
        j <- ij[2]
        res <- myintegral(vintegrand, 0, tau, i, j)
        #res <- try(integral(vintegrand, 0, tau,method = c("Kronrod"), vectorized = TRUE, arrayValued = FALSE, waypoints = NULL, reltol = 1e-05, abstol = 1e-5, i, j))
        ## if(class(res) == "try-error"){
        ##     browser()
        ## }
        res
    }
    for(itr in 1 : (m^2)){
        ij <- mij[itr, ]
        mA2[ij[1], ij[2]] <- A2(ij)
    }
    return(mA2)
     
}



fb21 <- function(ba,  bb, Afun, afun,  bg, x, t){
    expbbx <- exp(x %*% bb)
    integrand1 <- function(s, i){
        As <- ba[1] + ba[2] * (s - t)#A(s)
        as <- ba[2]#a(s)
        if(t== s){
        oAs <- ba[1] + ba[2] * (s - t)#A(s)
        oas <- ba[2]#a(s)
        }else{
            oAs <- Afun(s)
            oas <- afun(s)
        }
    #    As <- oAs
     #   as <- oas
        
        mu <-  As *   expbbx
        omu <-  oAs *   expbbx
        ddG(omu, bg, br)/dG(omu, bg, br) * expbbx * 1 * exp(-bg[i] * omu) * (bg[i] * oas * expbbx- dG(omu, bg, br)* oas *  expbbx ) 
     }

    integrand2 <- function(s, i){
        As <- ba[1] + ba[2] * (s - t)#A(s)
        as <- ba[2]#a(s)
        if(t== s){
        oAs <- ba[1] + ba[2] * (s - t)#A(s)
        oas <- ba[2]#a(s)
        }else{
            oAs <- Afun(s)
            oas <- afun(s)
        }
    #    As <- oAs
     #   as <- oas
        
        mu <-  As *   expbbx
        omu <-  oAs *   expbbx
         (ddG(omu, bg, br)/dG(omu, bg, br) * expbbx * c( s- t) + c(1/as))  * exp(-bg[i] * omu) * (bg[i] * oas * expbbx- dG(omu, bg, br)* oas *  expbbx ) 
  }
    integrand1 <- Vectorize(integrand1)
    integrand2 <- Vectorize(integrand2)
    b2 <- function(i){
        c(myintegral(integrand1, 0, mtau, i), myintegral(integrand2, 0, mtau, i))
        #c(integral(integrand2, 0.0001, mtau,  method = "Simpson", vectorized = T, arrayValued = T, waypoints = NULL, reltol = 1e-08, abstol = 0, i))
    }
    t(do.call(rbind, lapply(1 : m, b2)))
}
fA21 <- function(ba, bb, Afun, afun, bg, br, x, t){
    mA2 <- matrix(NA, m, m)
    expbbx <- exp(x %*% bb)
    integrand <- function(s, i, j){
        As <- ba[1] + ba[2] * (s - t)#A(s)
        as <- ba[2]#a(s)
        if(t== s){
        oAs <- ba[1] + ba[2] * (s - t)#A(s)
        oas <- ba[2]#a(s)
        }else{
            oAs <- Afun(s)
            oas <- afun(s)
        }
        As <- oAs
        as <- oas
        
        mu <-  As *  expbbx
        omu <-  oAs *   expbbx
        expbgbbx <- exp(-bg * As * expbbx)
        oexpbgbbx <- exp(-bg * oAs * expbbx)
        sumexpc <-  sum(exp(-bg * As * expbbx) *br )
        sumexprc <-  sum(exp(-bg * As * expbbx) * bg *br )
        res <- (expbgbbx[j] * bg[j] * br[j]/sumexprc -  expbgbbx[j]  * br[j]/sumexpc)  * oexpbgbbx[i] * (bg[i] * oas * expbbx- dG(omu, bg, br)*oas*  expbbx )
        if(is.na(res)|abs(res) == Inf|is.nan(res)){
           # browser()
        }
        res
    }
    vintegrand <- Vectorize(integrand)
    A2 <- function(ij){
        i <- ij[1]
        j <- ij[2]
        res <- myintegral(vintegrand, 0, mtau, i, j)
    #    res <- try(integral(vintegrand, 0.0001, mtau,method = c("Simpson"), vectorized = TRUE, arrayValued = FALSE, waypoints = NULL, reltol = 1e-08, abstol = 0, i, j))
        ## if(class(res) == "try-error"){
        ##     browser()
        ## }
        res
    }
    for(itr in 1 : (m^2)){
        ij <- mij[itr, ]
        mA2[ij[1], ij[2]] <- A2(ij)
    }
    return(mA2)
     
}
fb22 <- function(ba,  bb, Afun, afun,  bg, x, t){
    expbbx <- exp(x %*% bb)
    integrand1 <- function(s, i){
        if(TRUE){
        As <- ba[1] + ba[2] * (s - t)#A(s)
        as <- ba[2]#a(s)
        }else{
            As <- Afun(s)
            as <- afun(s)
        }
        mu <-  As *   expbbx
        A <- Afun(s)
        a <- afun(s)
        ((ddG(mu, bg, br)/dG(mu, bg, br) * expbbx * c(1) + c(0)) * (s <= tau)  - dG(mu, bg, br)* c(1)*expbbx) * (bg[i] * a * expbbx) * exp(-A* expbbx * bg[i])
     }
    
    integrand2 <- function(s, i){
        if(TRUE){
        As <- ba[1] + ba[2] * (s - t)#A(s)
        as <- ba[2]#a(s)
        }else{
            As <- Afun(s)
            as <- afun(s)
        }
        mu <-  As *   expbbx
        A <- Afun(s)
        a <- afun(s)
         ((ddG(mu, bg, br)/dG(mu, bg, br) * expbbx * (s- t) + 1/as) * (s <= tau)  - dG(mu, bg, br) *expbbx * (s- t)) * (bg[i] * a * expbbx) * exp(-A * expbbx * bg[i])
  }
    integrand1 <- Vectorize(integrand1)
    integrand2 <- Vectorize(integrand2)
    b2 <- function(i){
        c(integrate(integrand1, 0.001, mtau, i)$value, integrate(integrand2, 0.001, mtau, i)$value)
        #c(integral(integrand1, 0.001, mtau,  method = "Simpson", vectorized = T, arrayValued = T, reltol = 1e-05, abstol = 1e-5, i))
    }
    t(do.call(rbind, lapply(1 : m, b2)))
}
fA22 <- function(ba, bb, Afun, afun, bg, br, x, t){
    mA2 <- matrix(NA, m, m)
    expbbx <- exp(x %*% bb)
    integrand <- function(s, i, j){
        if(TRUE){
        As <- ba[1] + ba[2] * (s - t)#A(s)
        as <- ba[2]#a(s)
        }else{
            As <- Afun(s)
            as <- afun(s)
        }
        A <- Afun(s)
        a <- afun(s)
        mu <-  As *  expbbx
        expbgbbx <- exp(-bg * As * expbbx)
        sumexpc <-  sum(exp(-bg * As * expbbx) *br )
        sumexprc <-  sum(exp(-bg * As * expbbx) * bg *br )
        res <- ((expbgbbx[j] * bg[j] * br[j]/sumexprc -  expbgbbx[j]  * br[j]/sumexpc) * (s <= tau) +expbgbbx[j]  * br[j]/sumexpc )  * (bg[i] * a * expbbx) * exp(-A * expbbx * bg[i])
        if(is.na(res)|abs(res) == Inf|is.nan(res)){
            #browser()
        }
        res
    }
    vintegrand <- Vectorize(integrand)
    A2 <- function(ij){
        i <- ij[1]
        j <- ij[2]
        res <- integrate(vintegrand, 0, mtau, i, j)$value


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
    A2 <- fA2(Afun, afun, bb, bg, br, x, t)
    b2 <- fb2(Afun, afun, bb, bg, x, t)
    e2 <- t(ginv(t(A2)%*% A2) %*% (t(A2) %*% t(b2)))
}


scorealpha <- function(i, data,  ba, bb, Afun, afun, bg, br, t, me2){
    d <- data[i, 1]
    y <- data[i, 2]
    x <- data[i, 3: ncol(data)]
    expbbx <- exp(x %*% bb)
    a <- ba[2]
    A <- ( ba[1] + ba[2]* (y - t)) 
    mu <- ( ba[1] + ba[2]* (y - t)) *   expbbx
    expbgbbx <- exp(-bg * A *  expbbx)
    sumexpc <-  sum(exp(-bg * A * expbbx) *br )
    sumexprc <-  sum(exp(-bg  * A *  expbbx) *bg * br )
    
    e2 <- me2[[x + 1.5]]
    numf <- apply(cbind(e2[1, ] * bg *  expbgbbx * br, e2[2, ] * bg * expbgbbx * br), 2, sum)
    nums <- apply(cbind(e2[1, ] *   expbgbbx * br, e2[2, ] *  expbgbbx * br), 2, sum)
    compens <- function(s, flg){
        A <- ba[1] + ba[2] * (s - t)#A(s)
        a <- ba[2]#a(s)
        if(t == s){
        oA <- ba[1] + ba[2] * (s - t)#A(s)
        oa <- ba[2]#a(s)
        }else{
            oA <- Afun(s)
            oa <- afun(s)
        }
      # A <- oA
      # a <- oa
#        oA<- A
 #       oa <- a
        mu <- oA *   expbbx
        omu <- oA *   expbbx
       expbgbbx <- exp(-bg * oA *  expbbx)
        oexpbgbbx <- exp(-bg * oA *  expbbx)
    sumexpc <-  sum(exp(-bg * oA * expbbx) *br )
       sumexprc <-  sum(exp(-bg  * oA *  expbbx) *bg * br )
       numf <- apply(cbind(e2[1, ] * bg *  expbgbbx * br, e2[2, ] * bg * expbgbbx * br), 2, sum)
       nums <- apply(cbind(e2[1, ] *   expbgbbx * br, e2[2, ] *  expbgbbx * br), 2, sum)
        if(flg == 1){
       (((ddG(mu, bg, br)/dG(mu, bg, br)  * expbbx * c(1, s- t) + c(0, 1)/ a) - ((numf/ sumexprc) - nums/(sumexpc))) * dG(omu, bg, br) * expbbx * oa)[1]
       }else{
           (((ddG(mu, bg, br)/dG(mu, bg, br)  * expbbx * c(1, s- t) + c(0, 1)/ a) - ((numf/ sumexprc) - nums/(sumexpc))) * dG(omu, bg, br) * expbbx * oa)[2]
           }
        
        }
    compens <- Vectorize(compens)
    if(d == 1){
        temp <- c(myintegral(compens, 0, y, 1), myintegral(compens, 0, y, 2))
#        res <-  ddG(mu, bg, br)/dG(mu, bg, br)  * expbbx * c(1, y - t) + c(0, 1)/ a- ((numf/ sumexprc) - nums/(sumexpc)) - dG(mu, bg, br) * expbbx * c(1, y- t) -  nums/(sumexpc)
        res <- (ddG(mu, bg, br)/dG(mu, bg, br)  * expbbx * c(1, y - t) + c(0, 1)/ a- ((numf/ sumexprc) - nums/(sumexpc))) - temp#integral(compens, 0.0001, y, method = "Simpson", vectorized = TRUE, arrayValued = TRUE, waypoints = NULL,reltol = 1e-08, abstol = 0)
        

    }else{
      temp <- c(myintegral(compens, 0, y, 1), myintegral(compens, 0, y, 2))
  #   res <- - dG(mu, bg, br) * expbbx * c(1, y - t) -  nums/(sumexpc)
        res <- - temp # -integral(compens, 0.0001, y, method = "Simpson", vectorized = TRUE, arrayValued = TRUE, waypoints = NULL,reltol = 1e-08, abstol = 0)
      
    }
    (dnorm((y - t)/h)) * t(res) %*% res
}

dscorealpha <- function(i, data,  ba, bb, Afun, afun, bg, br, t, me2){
    d <- data[i, 1]
    y <- data[i, 2]
    x <- data[i, 3: ncol(data)]
    expbbx <- exp(x %*% bb)
    a <- ba[2]
    A <- ( ba[1] + ba[2]* (y - t)) 
    mu <- ( ba[1] + ba[2]* (y - t)) *   expbbx
    expbgbbx <- exp(-bg * A *  expbbx)
    sumexpc <-  sum(exp(-bg * A * expbbx) *br )
    sumexprc <-  sum(exp(-bg  * A *  expbbx) *bg * br )
    
    e2 <- me2[[x + 1.5]]
    numf <- apply(cbind(e2[1, ] * bg *  expbgbbx * br, e2[2, ] * bg * expbgbbx * br), 2, sum)
    nums <- apply(cbind(e2[1, ] *   expbgbbx * br, e2[2, ] *  expbgbbx * br), 2, sum)
    compens <- function(s){
        A <- ba[1] + ba[2] * (s - t)#A(s)
        a <- ba[2]#a(s)
        if(t == s){
        oA <- ba[1] + ba[2] * (s - t)#A(s)
        oa <- ba[2]#a(s)
        }else{
            oA <- Afun(s)
            oa <- afun(s)
        }
      #  A <- oA
      #  a <- oa
        mu <- A *   expbbx
        omu <- oA *   expbbx
       expbgbbx <- exp(-bg * A *  expbbx)
        oexpbgbbx <- exp(-bg * oA *  expbbx)
    sumexpc <-  sum(exp(-bg * A * expbbx) *br )
       sumexprc <-  sum(exp(-bg  * A *  expbbx) *bg * br )
       numf <- apply(cbind(e2[1, ] * bg *  expbgbbx * br, e2[2, ] * bg * expbgbbx * br), 2, sum)
       nums <- apply(cbind(e2[1, ] *   expbgbbx * br, e2[2, ] *  expbgbbx * br), 2, sum)
       ((ddG(mu, bg, br)/dG(mu, bg, br)  * expbbx * c(1, s - t) + c(0, 1)/ a) - ((numf/ sumexprc) - nums/(sumexpc))) * dG(omu, bg, br) * expbbx * oa 
        }
    compens <- Vectorize(compens)
    if(d == 1){

     #res <-  ddG(mu, bg, br)/dG(mu, bg, br)  * expbbx * c(1, y - t) + c(0, 1)/ a- ((numf/ sumexprc) - nums/(sumexpc)) - dG(mu, bg, br) * expbbx * c(1, y- t) -  nums/(sumexpc)
     res <- (ddG(mu, bg, br)/dG(mu, bg, br)  * expbbx * c(1, y - t) + c(0, 1)/ a- ((numf/ sumexprc) - nums/(sumexpc))) - integral(compens, 0.0001, y, method = "Simpson", vectorized = TRUE, arrayValued = TRUE,  reltol = 1e-08, abstol = 0)
  #      res2nd <- t(ba) %*% (((numf/ sumexprc) - nums/(sumexpc)) +  nums/sumexpc) 
   #     res <- -(log(dG(mu, bg, br)) +  log(a * expbbx) - G(mu, bg, br) - (res2nd)  )

    }else{
      
     # res <- - dG(mu, bg, br) * expbbx * c(1, y - t) -  nums/(sumexpc)
      res <- -integral(compens, 0.0001, y, method = "Simpson", vectorized = TRUE, arrayValued = TRUE,  reltol = 1e-08, abstol = 0)
#        res2nd <- t(ba) %*% (nums/sumexpc) 
 #       res <- G(mu, bg, br) + res2nd#t(res2nd) %*% ba
    }
    (dnorm((y - t)/h)) *  res
}

scorealphasum <- function(ba, bb, data, A, a, bg, br, t, me2){
    #print(ba)

    A2 <- fA21(ba, bb, A, a, bg, br, -0.5, t)
    b2 <- fb21(ba, bb, A, a, bg, -0.5, t)
    e2 <- t(ginv(t(A2)%*% A2) %*% (t(A2) %*% t(b2)))
    me2[[1]] <- e2
    A2 <- fA21(ba, bb, A, a, bg, br, 0.5, t)
    b2 <- fb21(ba, bb, A, a, bg, 0.5, t)
    e2 <- t(ginv(t(A2)%*% A2) %*% (t(A2) %*% t(b2)))
    me2[[2]] <- e2
    me2 <<- me2
    vsum <- apply(do.call(cbind, mclapply(1 :n,  scorealpha, data, ba, bb, A, a, bg, br, t, me2, mc.cores = 20)), 1, sum,  na.rm = T) #/n^(1/2)#* (oba)
    return(vsum)
#    1/2 * sqrt(t(vsum) %*% vsum)
    #mean(do.call(rbind, lapply(1 :n, scorealpha, data,   ba, bb, A, a, bg, br, t, me2)))
}

dscorealphasum <- function(ba, bb, data, A, a, bg, br, t, me2){
    #print(ba)
    ## if(i > 1){
    ##     oba <- exp(ba)
    ##     ba <- exp(ba) + c(exp((mba[i - 1, 1])), 0)
    ##     ba <- oba
    ## }else{
    ##     ba <- exp(ba)
    ##     oba <- ba
    ## }
    A2 <- fA22(ba, bb, A, a, bg, br, -0.5, t)
    b2 <- fb22(ba, bb, A, a, bg, -0.5, t)
    e2 <- t(ginv(t(A2)%*% A2) %*% (t(A2) %*% t(b2)))
    me2[[1]] <- e2
    A2 <- fA22(ba, bb, A, a, bg, br, 1.5, t)
    b2 <- fb22(ba, bb, A, a, bg, 1.5, t)
    e2 <- t(ginv(t(A2)%*% A2) %*% (t(A2) %*% t(b2)))
    me2[[2]] <- e2
    vsum <- apply(do.call(cbind, mclapply(1 :n,  dscorealpha, data, ba, bb, A, a, bg, br, t, me2, mc.cores = 20)), 1, mean,  na.rm = T) #* (oba)
    return(vsum)
#    sqrt(t(vsum) %*% vsum)
    #mean(do.call(rbind, lapply(1 :n, scorealpha, data,   ba, bb, A, a, bg, br, t, me2)))
}

tryscorealphasum <- function(ba, bb, data, A, a, bg, br, t, me2){
    res <- try(scorealphasum(ba, bb, data, A, a, bg, br, t, me2))

        #browser()
    res
}
##get ba function
fmb <- function(i, mt, bb, A, a, bg, br){
    t <- mt[i]
    me2 <- vector("list")# mclapply(1 : n, gete2, A, a, bb, bg, br, t, data, mc.cores = 10)
    if(t <= 2){
         temp <- spg(c(A(t), a(t)), tryscorealphasum, gr= NULL, method=1, lower=c(0, 0.001), upper=c(2, 2), project=NULL, projectArgs=NULL, control=list(maxit = 30, ftol = 1e-8, gtol = 1e-5), quiet=FALSE,  bb, data, A, a, bg, br, t, me2)#dfsane((c(A(t) , a(t))), dscorealphasum, method = 2, control = list(tol = 1e-5), quiet = FALSE,  bb, data, A, a, bg, br, t, me2)##dfsane((c(A(t), a(t))),tryscorealphasum, method = 2, control = list(tol = 1e-5), quiet = FALSE,  bb, data, A, a, bg, br, t, me2)#####m## #### optim(c(A(t), a(t)) , tryscorealphasum, gr = NULL, bb, data, A, a, bg, br, t, me2, method = "L-BFGS-B", lower = c(0, 0.001), upper = c(2, 2), control = list(), hessian = FALSE)# ####
      #  print(temp$convergence)
      #  print(temp$value)
       # print(temp$gradient)
        #print(temp$message)
    temp$par
}else{
        temp <- spg(c(A(t), a(t)), tryscorealphasum, gr= NULL, method=1, lower=c(0, 0.001), upper=c(30, 8), project=NULL, projectArgs=NULL, control=list(maxit = 30, ftol = 1e-8, gtol = 1e-5), quiet=FALSE,  bb, data, A, a, bg, br, t, me2)#dfsane((c(A(t) , a(t))), dscorealphasum, method = 2, control = list(tol = 1e-5), quiet = FALSE,  bb, data, A, a, bg, br, t, me2)##dfsane((c(A(t), a(t))),tryscorealphasum, method = 2, control = list(tol = 1e-5), quiet = FALSE,  bb, data, A, a, bg, br, t, me2)#####m## #### optim(c(A(t), a(t)) , tryscorealphasum, gr = NULL, bb, data, A, a, bg, br, t, me2, method = "L-BFGS-B", lower = c(0, 0.001), upper = c(2, 2), control = list(), hessian = FALSE)# ####
       # print(temp$convergence)
       # print(temp$value)
       # print(temp$gradient)
       # print(temp$message)
    temp$par
    
}
}
tryfmb <- function(i, mt, bb, A, a, bg, br){
    res <- try(fmb(i, mt, bb, A, a, bg, br))
    if(class(res) == "try-error"){
        return(c(NA, NA))
    }
    res
}

######function for bb###############
fb1 <- function(A, a, bb, bg, x){
    expbbx <- exp(x %*% bb)
    integrand <- function(s, i, x){
        mu <-  A(s) *   expbbx
        (ddG(mu, bg, br)/dG(mu, bg, br) * mu * x + x) * exp(-bg[i] * mu) * (bg[i] * a(s) * expbbx- dG(mu, bg, br)* a(s)*  expbbx ) 
     }

   
    integrand <- Vectorize(integrand)
    b1 <- function(i){
        res <- 0
        for(q in 1 : p){
        res <- c(res, myintegral(integrand, 0, mtau, i, x[q]))
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
        expbgbbx <- exp(-bg * A(s) * expbbx)
        sumexpc <-  sum(exp(-bg * A(s) * expbbx) *br )
        sumexprc <-  sum(exp(-bg * A(s) * expbbx) * bg *br )
        res <- (expbgbbx[j] * bg[j] * br[j]/sumexprc -  expbgbbx[j]  * br[j]/sumexpc)  * expbgbbx[i] * (bg[i] * a(s) * expbbx- dG(mu, bg, br)* a(s)*  expbbx )
        if(is.na(res)|abs(res) == Inf|is.nan(res)){
            #browser()
        }
        res
    }
    vintegrand <- Vectorize(integrand)
    A1 <- function(ij){
        i <- ij[1]
        j <- ij[2]
        res <- try(myintegral(vintegrand, 0, mtau, i, j))
        if(class(res) == "try-error"){
            #browser()
        }
        res
    }
    for(itr in 1 : (m^2)){
        ij <- mij[itr, ]
        mA1[ij[1], ij[2]] <- A1(ij)
    }
    return(mA1)
     
}

fb12 <- function(A, a, bb, bg, x){
    expbbx <- exp(x %*% bb)
    integrand <- function(s, i, x){
        mu <-  A(s) *   expbbx
        ((ddG(mu, bg, br)/dG(mu, bg, br) * mu * x + x) * (s<= tau) - dG(mu, bg, br)*  mu * x) * (bg[i] * a(s)* expbbx) * exp(-A(s)* expbbx * bg[i])
     }

   
    integrand <- Vectorize(integrand)
    b1 <- function(i){
        res <- 0
        for(q in 1 : p){
        res <- c(res, myintegral(integrand, 0, mtau, i, x[q]))
    }
        res <- res[-1]
    }
    t(do.call(rbind, lapply(1 : m, b1)))
}

fA12 <- function(A, a, bb, bg, br, x){
    mA1 <- matrix(NA, m, m)
    expbbx <- exp(x %*% bb)
    integrand <- function(s, i, j){
        mu <-  A(s) *  expbbx
        expbgbbx <- exp(-bg * A(s) * expbbx)
        sumexpc <-  sum(exp(-bg * A(s) * expbbx) *br )
        sumexprc <-  sum(exp(-bg * A(s) * expbbx) * bg *br )
        res <- ((expbgbbx[j] * bg[j] * br[j]/sumexprc -  expbgbbx[j]  * br[j]/sumexpc)* (s<=tau ) + expbgbbx[j]  * br[j]/sumexpc ) * (bg[i] * a(s) * expbbx) * exp(-A(s)* expbbx * bg[i])
        if(is.na(res)|abs(res) == Inf|is.nan(res)){
            #browser()
        }
        res
    }
    vintegrand <- Vectorize(integrand)
    A1 <- function(ij){
        i <- ij[1]
        j <- ij[2]
        res <- try(myintegral(vintegrand, 0, mtau, i, j))
        if(class(res) == "try-error"){
            #browser()
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
    e1 <- t(ginv(t(A1)%*% A1) %*% (t(A1) %*% t(b1)))
}

scorebeta <- function(i, data,   bb,   Afun, afun, bg, br,  me1){
    d <- data[i, 1]
    y <- data[i, 2]
    x <- data[i, 3: ncol(data)]
    expbbx <- exp(x %*% bb)
    a <- afun(y)#ba[2]
    A <- Afun(y)#( ba[1] + ba[2]* (y - t)) 
    mu <- A *   expbbx
    expbgbbx <- exp(-bg * A *  expbbx)
    sumexpc <-  sum(exp(-bg * A * expbbx) *br )
    sumexprc <-  sum(exp(-bg  * A *  expbbx) *bg * br)
    #A1 <- fA1(Afun, afun, bb, bg, br, x)
    #b1 <- fb1(Afun, afun, bb, bg, x)
    #e1 <- t(ginv(t(A1)%*% A1) %*% (t(A1) %*% t(b1)))
   e1 <- me1[[x+ 1.5]]
    numf <-  apply(matrix(apply(e1, 1, "*", bg *  expbgbbx * br), ncol = p), 2, sum)
    nums <-  apply(matrix(apply(e1, 1, "*",  expbgbbx * br), ncol = p), 2, sum)
    compens <- function(s){
       A <- Afun(s)
       a <- afun(s)
       mu <- A *   expbbx
       expbgbbx <- exp(-bg * A *  expbbx)
       sumexpc <-  sum(exp(-bg * A * expbbx) *br )
       sumexprc <-  sum(exp(-bg  * A *  expbbx) *bg * br)
       numf <-  apply(matrix(apply(e1, 1, "*", bg *  expbgbbx * br), ncol = p), 2, sum)
    nums <-  apply(matrix(apply(e1, 1, "*",  expbgbbx * br), ncol = p), 2, sum)
       ((ddG(mu, bg, br)/dG(mu, bg, br)  * expbbx * A * x + x) - ((numf/ sumexprc) - nums/(sumexpc))) * dG(mu, bg, br) * expbbx * a
        }

    
    compens <- Vectorize(compens)
    if(d == 1){
     res <-  (ddG(mu, bg, br)/dG(mu, bg, br)  * expbbx * A * x + x) - ((numf/ sumexprc) - nums/(sumexpc)) - dG(mu, bg, br) * A * expbbx * x-  nums/(sumexpc)
#      res <-  (ddG(mu, bg, br)/dG(mu, bg, br)  * expbbx * A * x + x) - ((numf/ sumexprc) - nums/(sumexpc)) - integrate(compens, 0, y)$value
    }else{
        #res <- - integrate(compens, 0, y)$value##- dG(mu, bg, br) * A * expbbx * x-  nums/(sumexpc)#
        res <- - dG(mu, bg, br) * A * expbbx * x-  nums/(sumexpc)
    }
    1/2 * (res) ^2
}

dscorebeta <- function(i, data,   bb,   Afun, afun, bg, br,  me1){
    d <- data[i, 1]
    y <- data[i, 2]
    x <- data[i, 3: ncol(data)]
    expbbx <- exp(x %*% bb)
    a <- afun(y)#ba[2]
    A <- Afun(y)#( ba[1] + ba[2]* (y - t)) 
    mu <- A *   expbbx
    expbgbbx <- exp(-bg * A *  expbbx)
    sumexpc <-  sum(exp(-bg * A * expbbx) *br )
    sumexprc <-  sum(exp(-bg  * A *  expbbx) *bg * br)
    #A1 <- fA1(Afun, afun, bb, bg, br, x)
    #b1 <- fb1(Afun, afun, bb, bg, x)
    #e1 <- t(ginv(t(A1)%*% A1) %*% (t(A1) %*% t(b1)))
   e1 <- me1[[x+ 1.5]]
    numf <-  apply(matrix(apply(e1, 1, "*", bg *  expbgbbx * br), ncol = p), 2, sum)
    nums <-  apply(matrix(apply(e1, 1, "*",  expbgbbx * br), ncol = p), 2, sum)
    compens <- function(s){
       A <- Afun(s)
       a <- afun(s)
       mu <- A *   expbbx
       expbgbbx <- exp(-bg * A *  expbbx)
       sumexpc <-  sum(exp(-bg * A * expbbx) *br )
       sumexprc <-  sum(exp(-bg  * A *  expbbx) *bg * br)
       numf <-  apply(matrix(apply(e1, 1, "*", bg *  expbgbbx * br), ncol = p), 2, sum)
    nums <-  apply(matrix(apply(e1, 1, "*",  expbgbbx * br), ncol = p), 2, sum)
      c1 <-  ((ddG(mu, bg, br)/dG(mu, bg, br)  * expbbx * A * x + x) - ((numf/ sumexprc) - nums/(sumexpc))) * dG(mu, bg, br) * expbbx * a
#       c2 <- ((ddG(mu, bg, br)/dG(mu, bg, br)  * expbbx * A * x + x) - ((numf/ sumexprc) - nums/(sumexpc))) * (ddG(mu, bg, br) * A * expbbx ^ 2 * x * a + dG(mu, bg, br) * expbbx * x * a)
       return(c(c1))
        }
    compens <- Vectorize(compens)
#    temp <- integral(compens, 0, y, method = "Simpson", vectorized = T, arrayValued = T, waypoints = NULL, reltol = 1e-05, abstol = 1e-5)
    if(d == 1){
     res <-  (ddG(mu, bg, br)/dG(mu, bg, br)  * expbbx * A * x + x) - ((numf/ sumexprc) - nums/(sumexpc)) - dG(mu, bg, br) * A * expbbx * x-  nums/(sumexpc)
        
#      res <-  ((ddG(mu, bg, br)/dG(mu, bg, br)  * expbbx * A * x + x) - ((numf/ sumexprc) - nums/(sumexpc)) - temp[1]) #* temp[2]#integrate(compens, 0, y)$value
    }else{
 # res <- - temp[1] #* temp[2]#- integrate(compens, 0, y)$value
       res <- - dG(mu, bg, br) * A * expbbx * x-  nums/(sumexpc)
    }
    (res) 
}

scorebetasum <- function( bb,  data, A, a, bg, br){
    me1 <- vector("list")

   #bb <- exp(bb)
   #print(bb)
   A1 <- (fA12(A, a, bb, bg, br, -0.5))
   b1 <- fb12(A, a, bb, bg, -0.5)
   me1[[1]] <- t(ginv(t(A1)%*% A1) %*% (t(A1) %*% t(b1)))
   A1 <- (fA12(A, a, bb, bg, br, 0.5))
   b1 <- fb12(A, a, bb, bg, 0.5)
   me1[[2]] <- t(ginv(t(A1)%*% A1) %*% (t(A1) %*% t(b1)))
   me1 <<- me1
  (apply(do.call(cbind, mclapply(1 :n, scorebeta, data,  bb, A, a, bg, br, me1, mc.cores = 20)), 1, sum, na.rm = T)) 
    #mean(do.call(rbind, lapply(1 :n, scorebeta, data,    bb,   A, a, bg, br,  me1)))
}

dscorebetasum <- function( bb,  data, A, a, bg, br){
   ##  obb <- bb
   ##  expbb<-  exp(bb)/(1 + exp(bb))
   ##  if(is.nan(expbb)){
   ##      expbb <- 1
   ##      }
   ## bb <- 4 * expbb- 2
   #print(bb)
   me1 <- vector("list")
   A1 <- (fA12(A, a, bb, bg, br, -0.5))
   b1 <- fb12(A, a, bb, bg, -0.5)
   me1[[1]] <- t(ginv(t(A1)%*% A1) %*% (t(A1) %*% t(b1)))
   A1 <- (fA12(A, a, bb, bg, br, 0.5))
   b1 <- fb12(A, a, bb, bg, 0.5)
   me1[[2]] <- t(ginv(t(A1)%*% A1) %*% (t(A1) %*% t(b1)))
  (apply(do.call(cbind, mclapply(1 :n, dscorebeta, data,  bb, A, a, bg, br, me1, mc.cores = 20)), 1, sum, na.rm = T))#* 4 * exp(obb)/(1 + exp(obb))^2 
    #mean(do.call(rbind, lapply(1 :n, scorebeta, data,    bb,   A, a, bg, br,  me1)))
}

tryscorebetasum <- function( bb,  data, A, a, bg, br){
    res <- try(scorebetasum( bb,  data, A, a, bg, br))

    res
}

fmbeta <- function( data, bb, A, a, bg, br){
 #   me1 <- vector("list")#mclapply(1 : n, gete1, A, a, bb, bg, br,  data, mc.cores = 10)
 #temp <- dfsane(bb, dscorebetasum, method = 2, control = list(maxit = 30), quiet = FALSE, data, A, a, bg, br)#
   temp <- spg(bb, tryscorebetasum, gr= NULL, method=3, lower= -1, upper=2, project=NULL, projectArgs=NULL, control=list(maxit = 15, ftol = 1e-6, gtol = 0.0001), quiet=FALSE,   data, A, a, bg, br)#optim(0.5, tryscorebetasum, gr = NULL, data, A, a, bg, br,   method = "L-BFGS-B", lower = -1, upper = 2, control = list(maxit = 30, REPORT = 1, pgtol = 1e-5, trace = 0), hessian = FALSE)#### 
 #   temp <- try(uniroot(dscorebetasum, c(-1, 2),  data, A, a, bg, br))
   # print(temp$convergence)
 #   temp[[1]]
    if(class(temp) == "try-error"){
        return(NA)
        }
    return(c(temp[[1]], temp$convergence))
}
#########simulation###############
simu <- function(n, v, bb,lambda,  alpha, tau){
   u <- runif(n)
   x <- matrix(rbinom(n, 1, 0.3), nrow = n) - 0.5
   rg <<- runif(n, 0.1, 5)#rgamma(n, 1/v, 1/v)
   t <<- (- log(1 - u )/ (rg * lambda * exp(x%*% bb)))^(1/alpha)
   d <- t <= tau
   cbind(d, pmin(t, tau), x)
   
}

myintegral <- function(fun, low, upper, ...){
    x <-  seq(low, upper, length.out = 200)
    y <- fun(x, ...)
    idx = 2:length(x)
    return (as.double( (x[idx] - x[idx-1]) %*% (y[idx] + y[idx-1])) / 2)
    
}

####test code###########
sta <- 1
mbeta <- matrix(NA, 100, 2)
lmba <- vector("list")
for(itr  in sta:(sta + 19)){
    print(itr)
    set.seed(itr + 2015)
    m <- 15
    mtau <- 25
    tau  <- 2.5
    bb <- 0.5
    obb <- 0.5
    data <- simu(500, v, bb, 0.5, 2, tau)
    bg <- seq(0.1, 5, length.out = m)# seq(qgamma(0.002, 1/v, 1/v), qgamma(0.998, 1/v, 1/v), length.out = m)#as.vector(quantile(rg, seq(0, 1, length.out = m)))
    br <- rep(1, m) / m #dgamma(bg, 1/v, 1/v)/ sum(dgamma(bg, 1/v, 1/v))#dnorm(bg, mean(rg), sd(rg))/sum(dnorm(bg, mean(rg), sd(rg)))##rep(1, 30)/30#
    n <- nrow(data)
    h <- bw.nrd(data[data[, 1] == 1, 2]) * n^{-1/15}
    mij <- as.matrix(expand.grid(1:m, 1:m))
    lmt <- 10
    #mt <- seq(quantile(data[data[, 1] == 1, 2], 0.1), quantile(data[data[, 1] == 1, 2], 0.8), length.out = lmt)
    mt <- c(seq(quantile(t, 0), 2, length.out = lmt - 1), max(t))
    mba <- matrix(0, lmt, 2)
    omba <- mba
    mtau <- max(t)
    A <- iniA
    a <- inia
for(j in 1: 2){
for(i in 1 : lmt){
    ## mba <- (do.call(rbind, mclapply(1:lmt, tryfmb, mt, bb, A, a, bg, br, mc.cores = 10)))
    mba[i, ] <- fmb(i, mt, bb, A, a, bg, br)
}
        A <- approxfun(mt, mba[, 1], method = "linear", yleft = 0, rule = 2)
        a <- approxfun(mt, mba[, 2], method = "linear",  rule = 2)    

mbb <- try(fmbeta(data, obb, A, a, bg, br))
    if(class(mbb) == "try-error"& j == 1){
        bb <- obb
    }else if(class(mbb) == "try-error"& j == 2){
    bb <- bb
}else{
    bb <- mbb[1]
}
    
}
 
save(mba, mbb, file = paste(paste("../res/", itr, sep = ""),  obb, "res", sep = "_"))
}

