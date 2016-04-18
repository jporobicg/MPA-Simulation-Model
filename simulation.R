#==============================================================#
#              MPA smulation Juan fernandez                    #
# Creator  :  Javier Porobic                                   #
# Date     :  6-Aug-2015                                       #
# Last mod :  6-Aug-2015                                       #
#==============================================================#



trans.matrix <- function(size.vec, k.g = NULL, linf = NULL, beta.g = NULL){
    ##################################################################################
    ## Normalize Transitin matrix                                                   ##
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

    return(normalize)
}

read.dat <- function(file)
    {
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
    ## Calculation of the fecundity vector=a_fec*(AVERAGE(H38:I38)^b_fec)
    vect <- rep(NA,length(size) - 1)
    for( i in 1 : length(vect)){
        vect[i] <- a * mean(size[i : (i+ 1)]) ^ b
    }
    return(vect)
}


## Calculation of recruitment
recruitment <- function(rec.a, rec.b, spawners, vt){
    ## this function calculate the recruitmen equation
    rec.out <- spawner / (log(rec.a) + log(rec.b) * spawners) * vt
    return(rec.out)
}

## Calculation of recruitment
rec.shp <- function(rec.areas, conect_matrix){
    ## Matrix multiplication to obtain the new shape in rec
    shape.rec.out <- rec.areas %*% conect_matrix
    return(shape.rec.out)
}

## Spawning stock
spawners <- function(fecundity, maturity, abundance){
    spawn <- sum(fecundity * maturity * abundance)
    return(spawn)
}
