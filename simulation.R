#==============================================================#
#              MPA smulation Juan fernandez                    #
# Creator  :  Javier Porobic                                   #
# Date     :  6-Aug-2015                                       #
# Last mod :  6-Aug-2015                                       #
#==============================================================#

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

aut.rd <- function(years, sigmaz, phi, seed = NULL){
    if(!is.null(seed)) set.seed(seed)
    rd     <- rnorm(50, sd = 3)
    sigmae <- sigmaz * (1 - phi ^ 2)
    zt     <- vector(mode = 'numeric',  length = 50)
    zt[1]  <- sqrt(sigmaz) * rd[1]
    for(dt in 2 : years){
        zt[dt] <- zt[dt - 1] * phi + rd[dt] * sqrt(sigmae)
    }
    return(zt)
}
