##==============================================================#
##              MPA smulation Juan fernandez                    #
## Creator  :  Javier Porobic                                   #
## Date     :  6-Aug-2015                                       #
## Last mod :  6-Aug-2015                                       #
##==============================================================#

trans.matrix <- function(size.vec, k.g = NULL, linf = NULL, beta.g = NULL){
##################################################################################
    ## Normalize Gamma shaped Transitin matrix                                      ##
    ## size.vec:   vector of legths
    ## beta.g  :   Beta from the gamma distribution
    ## k.g     :   Kappa from the gamma disgribution
    ## linf    :  maximum length in the polulation
##################################################################################
    n.sizes <- length(size.vec)
    intervals <- (size.vec + diff(size.vec/2))
    intervals[length(intervals)] <- size.vec[length(size.vec)] # limits
    ## mean growth
    mean.growth <-(linf - intervals) * (1 - exp( - k.g))
    ## Transition matrix
    trans.matrix <- matrix(NA, n.sizes, n.sizes)
    for( j in 1 : n.sizes){
        for( i in 1 : n.sizes){
            trans.matrix[j, i] <- ifelse(size.vec[i] > size.vec[j] & mean.growth[j] < 0.00001, 0,
                                         dgamma(intervals[i] - size.vec[j], mean.growth[j] / beta.g, beta.g))
        }
    }
    normalize <- apply(trans.matrix, 2, function(x) x / rowSums(trans.matrix))
    normalize[which(is.na(normalize))] <- 0
    return(normalize)
}

read.dat <- function(file){
    options(warn=-1)
    ##-----------------------------------------------------------#
    ##  Function to read data file in a ".dat" files             #
    ##        Creador:  Steve Martell  Mod: Javier Porobic       #
    ##-----------------------------------------------------------#
    ifile <- scan(file, what = "character", flush = TRUE, blank.lines.skip = FALSE, quiet = TRUE)
    idx   <- sapply(as.double(ifile), is.na)
    vnam  <- ifile[idx]	                    #list names
    nv    <- length(vnam)	                    #number of objects
    A     <- list()
    ir    <- 0
    for(i in 1 : nv)
    {
        ir <- match(vnam[i], ifile)
        if(i != nv) irr <- match(vnam[i + 1], ifile) else irr = length(ifile) + 1 #next row
        dum <- NA
        if(irr-ir==2) dum <- as.double(scan(file, skip = ir, nlines = 1, quiet = TRUE, what = ""))
        if(irr-ir>2)  dum <- as.matrix(read.table(file, skip = ir, nrow = irr-ir-1, fill = TRUE))
        if(is.numeric(dum))#Logical test to ensure dealing with numbers
        {
            A[[ vnam[i ]]] <- dum
        }
    }
    return(A)
}

fec <- function(a, b, size){
    ## Calculation of the fecundity at age
    vect <- vector()
    size <- c(size, tail(size,1) + diff(size)[1])
    for( i in 1 : (length(size) - 1)){
        vect[i] <- a * mean(size[i : (i+ 1)]) ^ b
    }
    return(vect)
}

matu <- function(a, b, size){
    ## Calculation of Maturity vat age
    vect <- vector()
    size <- c(size, tail(size,1) + diff(size)[1])
    for( i in 1 : (length(size) - 1)){
        vect[i] <- 1 / (1 + exp(a + b * mean(size[i : (i + 1)])))
    }
    return(vect)
}

alometry.f <- function(a, b, size){
    ## Alometric function
    vect <- vector()
    size <- c(size, tail(size,1) + diff(size)[1])
    for( i in 1 : (length(size) - 1)){
        vect[i] <- exp(a) * (mean(size[i : (i + 1)]) / 10) ^ b
    }
    return(vect)
}

selectivity.f <- function(l50, l95, size){
    vect <- 1 / (1 + exp(( - log(19) * (size - l50)) / (l95 - l50)))
    return(vect)
}

## Calculation of recruitment
recruitment <- function(rec.a, rec.b, spawners, vt){
    ## this function calculate the recruitmen equation
    rec.out <- spawners / (rec.a + rec.b * spawners) * vt
    return(rec.out)
}

## Recruitment pattern
rec.pat <- function(size, beta, alfa, normalize = TRUE){
    s.l     <- length(size)
    size <- c(size, tail(size,1) + diff(size)[1])
    pattern <- vector()
    for( i in 1 : s.l){
        pattern[i] <- dgamma(x = mean(size[i : (i + 1)]), shape = alfa, scale = beta)
    }
    if(isTRUE(normalize)) pattern <- pattern / sum(pattern)
    return(pattern)
}

## Calculation of recruitment
rec.shp <- function(rec.areas, conect_matrix){
    ## Matrix multiplication to obtain the new shape in rec
    shape.rec.out <- rec.areas %*% conect_matrix
    return(shape.rec.out)
}

## Abundance
abund <- function(alive, trans.matrix, patterns, recruit, n.zones = NULL){
    if(is.null(n.zones)){
        out <- (alive %*% trans.matrix) + patterns * recruit
    } else {
        out <- alive * NA
        for(z in 1 : n.zones){
            out[z, ] <- (alive[z, ] %*% trans.matrix) + patterns * recruit[z]
        }
    }
    return(out)
}

## Spawning stock
spawn <- function(fecundity, maturity, abundance){
    output <- sum(fecundity * maturity * abundance)
    return(output)
}

## Abundance after Benthic migration
benth.mig <- function(abundance, connect.mat, m.pattern, zone){
    ## Inmigrants
    inmig <- (colSums(t(t(abundance * connect.mat[, zone]) * mig.pat)))
    ## Migrants
    mig   <- (colSums(t(t(abundance[rep(zone, 8),  ] * connect.mat[zone,]) * mig.pat)))
                                        #Total migration
    migration <- abundance[zone, ] + inmig - mig
    return(migration)
}

## Catches
catch.f <- function(abundance.y, selec.fish, F, effort, M, zone){
    ## Baranov Catch Equation
    ## Using the total effort by area
    total.catch <- abundance.y[zone, ]  * ((selec.fish * F * effort[zone]) / (M + F * effort[zone] * selec.fish)) *
        (1 - exp( - (M + F * effort[zone] * selec.fish)))
    return(total.catch)
}

alive.f <- function(abun, M, Sel.F, F, effort, zone){
    out <- abun[zone, ] * exp( - (M + Sel.F * F * effort[zone]))
    return(out)
}


simulation <- function(M, R0, stpnss.h, n.years, maturity, w.l, f.selec, trap.s, f.cur, projection = FALSE, fecundity, pela.mat, bent.mat, n.zones, mig.pat,
                       n.rec.pat, normal.t.matrix, res.Rec, t.rec.a = NULL, t.rec.b = NULL, abun.sim = NULL, f.end = NULL, res.rec.end=NULL, qCPUE){
    abundance <- array(0, dim = c(n.zones, n.years, length(fecundity)))
    aamig     <- catch <- lfd <- alive <- abundance
    rec       <- matrix(NA, ncol= n.years, nrow = n.zones)
    b.catch   <- v.biomass <- spawners <- rec
    if(!isTRUE(projection)){
        ## |||||~~~~ At equilibrum ~~~~||||
        ## ~ initial recruitment
        pela.avg    <- rowMeans(pela.mat, dim = 2)
        rec[, 1]      <- R0 / n.zones
        rec.shape     <- rec.shp(rec[, 1], pela.avg)
        for(z in 1 : n.zones){ # loop for the zones
            ##~ Abundance
            abundance[z, 1, ] <- (n.rec.pat %*% solve(diag(1, dim(normal.t.matrix)) - normal.t.matrix * exp(-M))) * rec.shape[z]
            ##~ Biomass
            v.biomass[z, 1] <- sum(abundance[z, 1, ] *  w.l * f.selec) / 1000000
            ## Catch
            catch[z, 1, ]   <- catch.f(abundance[, 1, ], f.selec, exp(f.cur[1]), h.effort, M, z)
            ## ldf
            lfd[z, 1, ]     <- catch.f(abundance[, 1, ], trap.s, exp(f.cur[1]), h.effort, M, z)
            ## Catch biomass
            b.catch[z, 1]   <- sum(catch[z, 1, ] * w.l) / 1000000
        }
        ##~ Spawning biomass at equilibrium
        S0            <- spawn(fecundity, maturity, colSums(abundance[1 : 4, 1, ]))
        t.rec.a       <-  (S0 / R0 )  * (1 - (stpnss.h - 0.2) / (0.8 * stpnss.h))
        t.rec.b       <-  (stpnss.h - 0.2) / (0.8 * stpnss.h * R0)
        spawners[, 1] <- S0 / n.zones
    } else {
        if(any(c(is.null(t.rec.a), is.null(t.rec.b)))){
            cat('\n\tfor projection you need to provide the values for alfa and beta in the spawning equation\n')
        }
        if(any(c(is.null(alive), is.null(abun.sim)))){
            cat('\n\t You need to provide the initial abundance')
        }
        abundance[1 : n.zones, 1, ] <- abun.sim
        for(z in 1 : n.zones){ # loop for the zones
            v.biomass[z, 1]   <- sum(abundance[z, 1, ] *  w.l * f.selec) / 1000000
            catch[z, 1, ]     <- catch.f(abundance[, 1, ], f.selec, f.end, h.effort, M, z)
            b.catch[z, 1]     <- sum(catch[z, 1, ] * w.l) / 1000000
            lfd[z, 1, ]       <- catch.f(abundance[, 1, ], trap.s, exp(f.cur[1]), h.effort, M, z)
            spawners[z, 1]    <- spawn(fecundity, maturity, abundance[z, 1, ])
        }
        rec[, 1]    <- recruitment(t.rec.a, t.rec.b, spawners[, 1], exp(res.rec.end))
        rec.shape   <- rec.shp(rec[, 1], pela.mat)
    }
    for(year in 2 : n.years){
        ##~ Recruitment
        rec[, year]  <-  recruitment(t.rec.a, t.rec.b, spawners[, (year - 1)], exp(res.Rec[year]))
        if(year < 51 || year > (50 + dim(pela.mat)[3])){
            pela.con <- pela.avg
        } else {
            pela.con <- pela.mat[, , year - 50]
        }
        rec.shape  <- rec.shp(rec[, year], pela.con)
        ##~ Recruitment by zone after connectivity
        for(z in 1 : n.zones){
            ## Abundance after Benthic migration
            aamig[z, year, ]     <- benth.mig(abundance[, year-1, ], bent.mat, mig.pat, z)
            ## Alive after catch, Natural moratality and migration
            alive[z, year, ]     <- alive.f(aamig[, year, ], M, f.selec, f.cur[year - 1], h.effort, z)
            ## Abundance
            abundance[z, year, ] <- abund(alive[z, year, ], normal.t.matrix, n.rec.pat, rec.shape[z])
            ## Catch
            catch[z, year, ]     <- catch.f(abundance[, year, ], f.selec, f.cur[year], h.effort, M, z)
            ## lfd
            lfd[z, year, ]       <- catch.f(abundance[, year, ], trap.s, exp(f.cur[year]), h.effort, M, z)
            ## Biomass catch
            b.catch[z, year]     <- sum(catch[z, year, ] * w.l) / 1000000
            ## Spawners
            spawners[z, year]    <- spawn(fecundity, maturity, abundance[z, year, ])
            ##~ Biomass
            v.biomass[z, year]   <- sum(abundance[z, year, ] *  w.l * f.selec) / 1000000
        }
    }
    lfd.t <- colSums(lfd[, pcll, ], na.rm = TRUE, 1)
    lfd.t <- t(apply(lfd.t, 1, function(x) x / sum(x, na.rm = TRUE)))
    cpue.est <- colSums(v.biomass[, pcpue]) * exp(qCPUE)
        output <- list(Abundance  = abundance,
                   Catch      = catch,
                   B.catch    = b.catch,
                   Spawners   = spawners,
                   V.biomass  = v.biomass,
                   Trap.lfd   = lfd.t,
                   CPUE.cpp   = cpue.est,
                   Rec        = rec)
    if(!isTRUE(projection)){
        output$Parameters <- list(S0 = S0, t.rec.a = t.rec.a, t.rec.b = t.rec.b)
    }
    return(output)
}

aut.rd <- function(years, sigmaz, phi, seed = NULL){
    if(!is.null(seed)) set.seed(seed)
    rd     <- rnorm(years, sd = 3)
    sigmae <- sigmaz * (1 - phi ^ 2)
    zt     <- vector(mode = 'numeric',  length = years)
    zt[1]  <- sqrt(sigmaz) * rd[1]
    for(dt in 2 : years){
        zt[dt] <- zt[dt - 1] * phi + rd[dt] * sqrt(sigmae)
    }
    return(zt)
}


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Maximum Likelihood Estimator of parameters of multinomial distribution
##' @param obs Observed values
##' @param est estimated values
##' @param ssize sample size
##' @return likelihood
##' @author Demiurgo
multinomial <- function(obs, est, ssize){
    offset <- sum( - 1 *  (ssize * ((obs + 0.001) * log(obs + 0.001))), na.rm = TRUE)
    v      <- sum(ssize * ((obs + 0.001) * log(est + 0.001)), na.rm = TRUE)
    output <- sum(-1 * v - offset, na.rm = TRUE)
    return(output)
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Maximum likelihood for CPUE
##' @param obs obser4ved values
##' @param est Estimated values
##' @param cv coeficient of variation
##' @return likelihood
##' @author Demiurgo
like.normal <- function(obs, est, cv){
    LL <- 0.5 * sum(log(cv ^ 2 * 2 * pi), na.rm = TRUE) + sum(((log(obs) - log(est)) ^ 2 / (2 * (cv) ^ 2)), na.rm = TRUE)
    return(LL)
}

like.rec <- function(rec, cv){
    LL <- 0.5 * cv * sum(rec ^ 2)
    return(LL)
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Function for the estimation of paramters
##' @param Par Parameter to stimated,  for this model correspond to R0,  qCPUE, F and recdev
##' @param M Natural mortality
##' @param stpnss.h Steepness
##' @param n.years number of years
##' @param maturity maturity vector
##' @param w.l weight at length
##' @param f.selec fishery selectivity (knife - edge)
##' @param f.fish Selectivity of the trap
##' @param fecundity fecundity
##' @param pela.mat matrix of pelagic connectivity
##' @param bent.mat matrix of benthic connectivity
##' @param n.zones number of zones
##' @param mig.pat pattern of migration
##' @param n.rec.pat recruitment pattern
##' @param normal.t.matrix Normlaized size transition matrix
##' @param t.rec.a recruitment parameters (alpha)
##' @param t.rec.b recruitment parameter (beta)
##' @param abun.sim (simulated abundance)
##' @param f.end final fishing mortality
##' @param res.rec.end residual recruitments
##' @return estimated parameters
##' @author Demiurgo
hindcast <- function(Par, M, stpnss.h, n.years, maturity, w.l, f.selec, trap.s, fecundity, pela.mat, bent.mat, n.zones, mig.pat,
                     n.rec.pat, normal.t.matrix, phase = 1, cvs){
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    ## ~        At equilibrium    ~ ##
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    if(length(dim(pela.mat)) > 2)
    abundance <- array(0, dim = c(n.zones, n.years, length(fecundity)))
    aamig     <- catch <- alive <- lfd <- abundance
    rec       <- matrix(NA, ncol= n.years, nrow = n.zones)
    b.catch   <- v.biomass <- spawners <- rec
    ## Paramter for stimation. this part is hardwired
    R0      <- 13.3416
    f.cur   <- c(rep(1, 29), rep(1, 51), rep(1, 19), rep(1, 12), rep(1, 5))
    res.Rec <- rep(0, length(Par[8 : 123]))
    ## updating values
    if(phase >= 1) qCPUE   <- Par[2]
    if(phase >= 2) R0      <- Par[1]
    if(phase >= 3) f.cur   <- c(rep(Par[3], 29), rep(Par[4], 51), rep(Par[5], 19), rep(Par[6], 12), rep(Par[7], 5))
    if(phase >= 4) res.Rec <- Par[8 : 123]
    ## |||||~~~~ At equilibrum ~~~~||||
    ## ~ initial recruitment
    rec[, 1]      <- exp(R0) / n.zones
    pela.avg    <- rowMeans(pela.mat, dim=2)
    rec.shape     <- rec.shp(rec[, 1], pela.avg)
    for(z in 1 : n.zones){ # loop for the zones
        ##~ Abundance
        abundance[z, 1, ] <- (n.rec.pat %*% solve(diag(1, dim(normal.t.matrix)) - normal.t.matrix * exp(-M))) * rec.shape[z]
        ##~ Biomass
        v.biomass[z, 1] <- sum(abundance[z, 1, ] *  w.l * f.selec) / 1000000
        ## Catch
        catch[z, 1, ]   <- catch.f(abundance[, 1, ], f.selec, exp(f.cur[1]), h.effort, M, z)
        ## ldf
        lfd[z, 1, ]     <- catch.f(abundance[, 1, ], trap.s, exp(f.cur[1]), h.effort, M, z)
        ## Catch biomass
        b.catch[z, 1]   <- sum(catch[z, 1, ] * w.l) / 1000000
    }
    ##~ Spawning biomass at equilibrium
    S0            <- spawn(fecundity, maturity, colSums(abundance[1 : n.zones, 1, ]))
    t.rec.a       <-  (S0 / exp(R0) )  * (1 - (stpnss.h - 0.2) / (0.8 * stpnss.h))
    t.rec.b       <-  (stpnss.h - 0.2) / (0.8 * stpnss.h * exp(R0))
    spawners[, 1] <- S0 / n.zones
    for(year in 2 : n.years){
        ##~ Recruitment
        rec[, year]  <-  recruitment(t.rec.a, t.rec.b, spawners[, (year - 1)], exp(res.Rec[year]))
        ##~ Recruitment by zone after connectivity
        if(year < 51 || year > (50 + dim(pela.mat)[3])){
            pela.con <- pela.avg
        } else {
            pela.con <- pela.mat[, , year - 50]
        }
        rec.shape  <- rec.shp(rec[, year], pela.con)
        for(z in 1 : n.zones){
            ## Abundance after Benthic migration
            aamig[z, year, ]     <- benth.mig(abundance[, year-1, ], bent.mat, mig.pat, z)
            ## Alive after catch, Natural moratality and migration
            alive[z, year, ]     <- alive.f(aamig[, year, ], M, f.selec, exp(f.cur[year - 1]), h.effort, z)
            ## Abundance
            abundance[z, year, ] <- abund(alive[z, year, ], normal.t.matrix, n.rec.pat, rec.shape[z])
            ## Catch
            catch[z, year, ]     <- catch.f(abundance[, year, ], f.selec, exp(f.cur[year]), h.effort, M, z)
            ## lfd
            lfd[z, year, ]       <- catch.f(abundance[, year, ], trap.s, exp(f.cur[year]), h.effort, M, z)
            ## Biomass catch
            b.catch[z, year]     <- sum(catch[z, year, ] * w.l) / 1000000
            ## Spawners
            spawners[z, year]    <- spawn(fecundity, maturity, abundance[z, year, ])
            ##~ Biomass
            v.biomass[z, year]   <- sum(abundance[z, year, ] *  w.l * f.selec) / 1000000
        }
    }
    ##----------------------------
    ## Likelihoods
    ##---------------------------
    ## transformation LFD to proportions
    #browser()
    lfd.t <- colSums(lfd[, pcll, ], na.rm = TRUE, 1)
    ldf.t <- t(apply(lfd.t, 1, function(x) x / sum(x, na.rm = TRUE)))
    ## LLCPUE
    cpue.est <- colSums(v.biomass[, pcpue]) * exp(qCPUE)
    cpue.ll  <- like.normal(obs = datos$cpue, est = cpue.est, cv = rep(cvs[1], length(pcpue)))
    ## LLCATCH
    catch.ll <- like.normal(obs = datos$total_annual_catch, est = colSums(b.catch[, pcatch]), cv = rep(cvs[2], length(pcatch)))
    ## LLSIZE
    lfd.ll   <- multinomial(obs = datos$cll, est = ldf.t, ssize = cvs[3])
    ## recruitment
    rec.ll <- like.rec(res.Rec, cvs[4])
    ## Total
    llike    <- sum(cpue.ll, catch.ll, lfd.ll,  rec.ll)
    ## likeout
    return(llike)
}
