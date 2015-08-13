#==============================================================#
#              MPA smulation Juan fernandez                    #
# Creator  :  Javier Porobic                                   #
# Date     :  6-Aug-2015                                       #
# Last mod :  6-Aug-2015                                       #
#==============================================================#



trans.matrix <- function(sequence = NULL, L.ini = NULL, L.end = NULL, lag = NULL, k.g = 0.06985, linf = 213.495, beta.g = 1){

    if(is.null(sequence)){
        vec <- seq(from  = 10,  to = 214, by = 2)
    } else {
        if(!isTRUE(any(L.ini, L.end, lag))) stop ('one of the parameters have it is null')
        vec <- seq(from = L.ini, to = L.end, by = lag)
    }

    ## Matrices
    mat2 <- mat3 <- mat1 <- matrix(NA, length(vec), length(vec))
    mean.growth <- vec * NA
    ## Fill the matrix with the the gamma - Growth
    for(i in 1 : (length(vec) - 1)){    # Rows
        mean.growth[i] <- (linf - mean(vec[i : (i + 1)])) * (1 - exp( - k.g))
        if(i == (length(vec) - 1)){
                    mean.growth[i + 1] <- (linf - vec[i + 1]) * (1 - exp( - k.g))
                }
        for(j in 1 : (length(vec) - 1)){
            mat1[i, j] <- ifelse(vec[j] > vec[i] && mean.growth[i] <= 0.0001, 0, dgamma(mean(vec[j : (j + 1)]) - vec[i], mean.growth[i] / beta.g, 1))
            if(j ==                     #Fix the at value on the column
        }
    }

    plot(vec, mat1[1, ])
    mean.growth
