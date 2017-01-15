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
                                        #-----------------------------------------------------------#
                                        #  Function to read data file in a ".dat" files             #
                                        #        Creador:  Steve Martell  Mod: Javier Porobic       #
                                        #-----------------------------------------------------------#
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


simulation <- function(M, R0, stpnss.h, n.years, maturity, w.l, f.selec, f.cur, projection = FALSE, fecundity, pela.mat, bent.mat, n.zones, mig.pat,
                       n.rec.pat, normal.t.matrix, res.Rec, t.rec.a = NULL, t.rec.b = NULL, abun.sim = NULL, f.end = NULL, res.rec.end=NULL){
    abundance <- array(0, dim = c(n.zones, n.years, length(fecundity)))
    aamig     <- catch <- alive <- abundance
    rec       <- matrix(NA, ncol= n.years, nrow = n.zones)
    b.catch   <- v.biomass <- spawners <- rec
    if(!isTRUE(projection)){
        ## |||||~~~~ In equilibrum ~~~~||||
        ## ~ initial recruitment
        rec[, 1]      <- R0
        rec.shape     <- rec.shp(rec[, 1], pela.mat)
        for(z in 1 : n.zones){ # loop for the zones
            ##~ Abundance
            abundance[z, 1, ] <- (n.rec.pat %*% solve(diag(1, dim(normal.t.matrix)) - normal.t.matrix * exp(-M))) * rec.shape[z]
            ##~ Biomass
            v.biomass[z, 1] <- sum(abundance[z, 1, ] *  w.l * f.selec) / 1000000
        }
        ##~ Spawning biomass at equilibrium
        S0            <- spawn(fecundity, maturity, colSums(abundance[1 : 4, 1, ])) / 4
        t.rec.a       <-  (S0 / R0 )  * (1 - (stpnss.h - 0.2) / (0.8 * stpnss.h))
        t.rec.b       <-  (stpnss.h - 0.2) / (0.8 * stpnss.h * R0)
        spawners[, 1] <- S0
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
            spawners[z, 1]    <- spawn(fecundity, maturity, abundance[z, 1, ])
        }
        rec[, 1]    <- recruitment(t.rec.a, t.rec.b, spawners[, 1], exp(res.rec.end))
        rec.shape   <- rec.shp(rec[, 1], pela.mat)
    }
    for(year in 2 : n.years){
        ##~ Recruitment
        rec[, year]  <-  recruitment(t.rec.a, t.rec.b, spawners[, (year - 1)], exp(res.Rec[year]))
        ##~ Recruitment by zone after connectivity
        rec.shape  <- rec.shp(rec[, year], pela.mat)
        for(z in 1 : n.zones){
            ## Abundance after Benthic migration
            aamig[z, year, ]     <- benth.mig(abundance[, year-1, ], bent.mat, mig.pat, z)
            ## Alive after catch, Natural moratality and migration
            alive[z, year, ]     <- alive.f(aamig[, year, ], M, f.selec, f.cur[year - 1], h.effort, z)
            ## Abundance
            abundance[z, year, ] <- abund(alive[z, year, ], normal.t.matrix, n.rec.pat, rec.shape[z])
            ## Catch
            catch[z, year, ]     <- catch.f(abundance[, year, ], f.selec, f.cur[year], h.effort, M, z)
            ## Biomass catch
            b.catch[z, year]     <- sum(catch[z, year, ] * w.l) / 1000000
            ## Spawners
            spawners[z, year]    <- spawn(fecundity, maturity, abundance[z, year, ])
            ##~ Biomass
            v.biomass[z, year]   <- sum(abundance[z, year, ] *  w.l * f.selec) / 1000000
        }
    }
    output <- list(Abundance  = abundance,
                   Catch      = catch,
                   B.catch    = b.catch,
                   Spawners   = spawners,
                   V.biomass  = v.biomass,
                   Rec        = rec)
    if(!isTRUE(projection)){
        output$Parameters <- list(S0 = S0, t.rec.a = t.rec.a, t.rec.b = t.rec.b)
    }
                   ## Parameters = ifelse(!isTRUE(projection), list(S0      = S0,
                   ##                                               t.rec.a = t.rec.a,
                   ##                                               t.rec.b = t.rec.b), NA))

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
