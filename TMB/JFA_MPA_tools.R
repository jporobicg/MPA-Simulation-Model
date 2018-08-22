##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Read data files
##' @param file Dat file
##' @return a list with all the information contained on the dat file
##' @author Demiurgo
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

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Normalize Gamma shaped Transition matrix
##' @param size.vec Vector of lengths
##' @param k.g Kappa paratmer from the gamma distribution
##' @param linf Average maximum length in the popularion
##' @param beta.g Beta from the gamma distribution
##' @return A transition matrix
##' @author Demiurgo
trans.matrix <- function(size.vec, k.g = NULL, linf = NULL, beta.g = NULL){
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

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Patern of recruitment at size
##' @param size Vector of sizes
##' @param beta Beta paramter for the gamma function
##' @param alpha alpha paramter for the gamma function
##' @param normalized default TRUE. if True the value of the vector are normalized to 1 if not, keep the original values
##' @return A verctor with the redistribution of the recruitment
##' @author Demiurgo
rec.pat <- function(size, beta, alpha, normalized = TRUE){
    s.l     <- length(size)
    size <- c(size, tail(size,1) + diff(size)[1])
    pattern <- vector()
    for( i in 1 : s.l){
        pattern[i] <- dgamma(x = mean(size[i : (i + 1)]), shape = alpha, scale = beta)
    }
    if(isTRUE(normalized)) pattern <- pattern / sum(pattern)
    return(pattern)
}
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Calculating the maturity at size
##' @param alpha Apha paramter
##' @param beta Beta paramter
##' @param size Size class
##' @return a verctor with the maturity al size
##' @author Demiurgo
mature <- function(alpha, beta, size){
    ## Calculation of Maturity vat age
    vect <- vector()
    size <- c(size, tail(size, 1) + diff(size)[1])
    for( i in 1 : (length(size) - 1)){
        vect[i] <- 1 / (1 + exp(alpha + beta * mean(size[i : (i + 1)])))
    }
    return(vect)
}
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title allometric function
##' @param alpha Alpha paramter
##' @param beta Beta parameter
##' @param size Vector of sizes
##' @return A vector with the maturity at size
##' @author Demiurgo
alometry.f <- function(alpha, beta, size){
    ## Alometric function
    vect <- vector()
    size <- c(size, tail(size, 1) + diff(size)[1])
    for( i in 1 : (length(size) - 1)){
        vect[i] <- exp(alpha) * (mean(size[i : (i + 1)]) / 10) ^ beta
    }
    return(vect)
}
