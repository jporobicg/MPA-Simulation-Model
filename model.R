## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~                        Model to simulate the effect of an MPA in JFRE                    ~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
rm(list = ls())
## Run model
source('simulation.R')
## Read data file
datos <- read.dat('data/jf_f3.dat')
## Stimated parameters
param <- read.dat('data/jf_f.par')
## Read estimation from ADMB
est   <- read.dat('data/jf_f.rep')
## Just as reminder the order of the zones are
zones <- c('RC-SC NO', 'RC-SC NE', 'RC-SC SO', 'RC-SC SE', 'AS-NE', 'AS-SE', 'AS-SO', 'AS-NO')


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##          Conectivity Matrix            ##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## Pelagic connectivity matrix
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#pela.mat <- read.dat('data/conect_mat.dat')[[3]]
pela.mat <- diag(1, 1, 1) ## Matrix to test the model, with one population

## Benthic connectivity matrix
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#bent.mat <- read.dat('data/conect_mat.dat')[[5]]
bent.mat <- diag(0, 1, 1) ## Matrix to test the model, with one population

##~~~~~~~~~~~~~~~~~~~~##
##  Known parameters  ##
##~~~~~~~~~~~~~~~~~~~~##
## Fecundity
a.fec    <- 0.00047
b.fec    <- 4.4005
## Maturity
a.mat    <- 26.03
b.mat    <- -0.37
## Length - weight (Alometric relation)
a.wl     <-  0.0031
b.wl     <-  2.6672
# Selectivity of the fishery (almost knife - edge)
f.l50    <-  114.5
f.l95    <-  115
# Selectivity of the trap
p.l50    <- 95.2134056031
p.l95    <- 106.000000240
# Natural mortality
M        <- 0.18
# Steepness
stpnss.h <- 0.7

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##     Estimated in ADMB                  ##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

R0         <- est$Rec[1] # / 4  # Virginal Recruitment
size.vec   <- with(datos, seq(from = lowersize, to = bigsize, by = sizedelta))
res.Rec    <- log(est$Recruitment_Residuals)             # Recruitment residual
res.Rec2   <- 0           # Recruitment residual
#Rec        <- est$Rec
Fs         <- with(param, exp(c(logF_1, logF_2, logF_3, logF_4)))




##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##      Equation and models               ##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
start.year <- 1901
end.year   <- 2011

fecundity       <- fec(a.fec, b.fec, size.vec)                                       # Fecundity
maturity        <- matu(a.mat, b.mat, size.vec)                                      # Maturity
w.l             <- alometry.f(a.wl, b.wl, size.vec)                                  # Alometric relationship
f.selec         <- selectivity.f(f.l50, f.l95, size.vec)                             # Fishery Selectivity
p.selec         <- selectivity.f(p.l50, p.l95, size.vec)                             # Trap Selectivity
#p.selec <- est$Selectivity_pot

normal.t.matrix <- trans.matrix(size.vec, k.g = 0.06985, linf = 213.495, beta.g = 1) # Normalize transition matrix
mig.pat         <- c(rep(0.3, 28), rep(0.6, 8), rep(1, 67))                          # Migration Pattern (Puerulus - juvenil - adult)
n.rec.pat       <- rec.pat(size.vec, beta = 0.3, alfa = 45)                          # Normalize recruitment pattern
#h.effort <- c(0.84, 0.24, 2.32, 0.60, 0.72, 0.84, 1.08, 1.40) # Historical effort by zone
h.effort <- c(1, 1, 1, 1, 1, 1, 1, 1) # Historical effort by zone


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##          Simulation                    ##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
n.zones   <- 1
years.sim <- start.year : end.year
reclut    <- matrix(NA, ncol= length(years.sim), nrow = n.zones)
rec       <- matrix(NA, ncol= length(years.sim), nrow = n.zones)
rec2      <- matrix(NA, ncol= length(years.sim), nrow = n.zones)
rec3      <- matrix(NA, ncol= length(years.sim), nrow = n.zones)
b.catch   <- v.biomass <- spawners  <- rec
## Array for the abundance and the biomass
abundance2 <- abundance <- array(0, dim = c(length(zones), length(years.sim), length(size.vec)))
aamig2    <- aamig <- catch <- alive <- abundance

f.cur <- c(rep(Fs[1], 29), rep(Fs[2], 51), rep(Fs[3], 19), rep(Fs[4], 12))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~        At equilibrium    ~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

##~ initial recruitment
rec[, 1]  <- R0
rec2[, 1] <- R0
rec3[, 1] <- R0

rec.shape     <- rec.shp(rec[, 1], pela.mat)
for(z in 1 : n.zones){ # loop for the zones
    ##~ Abundance
    abundance[z, 1, ] <- (n.rec.pat %*% solve(diag(1, dim(normal.t.matrix)) - normal.t.matrix * exp(-M))) * rec.shape[z]
    abundance2[z, 1, ] <- (n.rec.pat %*% solve(diag(1, dim(normal.t.matrix)) - normal.t.matrix * exp(-M))) * rec.shape[z]
    #abundance[z, 1, ] <- est$Number[1,]
    #abundance2[z, 1, ] <- est$Number[1,]
    ##~ Biomass
    v.biomass[z, 1] <- sum(abundance[z, 1, ] *  w.l * f.selec) / 1000000
}
##~ Spawning biomass at equilibrium
S0       <- spawn(fecundity, maturity, abundance[1, 1, ])#  / 4
## S0M       <- spawn(fecundity, est$Maturity, sum(abundance[1 : 4, 1, ]))#  / 4
## S0F       <- spawn(est$Fecundity, maturity, sum(abundance[1 : 4, 1, ]))#  / 4
## S0N       <- spawn(fecundity, maturity, sum(est$N[1, ]))#  / 4
## SOA <- spawn(est$Fecundity, est$Maturity, sum(est$N[1, ]))#  / 4

## sum(abundance[1, 1, ] * fecundity * maturity)
## sum(abundance[1, 1, ] * est$Reprod)

## plot(aundance(fecundity * maturity)))

## S0
## S0M
## S0F
## S0N
## SOA


#S0      <- est$Spawners[1]
t.rec.a <-  (S0 / R0)  * (1 - (stpnss.h - 0.2) / (0.8 * stpnss.h))
t.rec.b <-  (stpnss.h - 0.2) / (0.8 * stpnss.h * R0)
spawners[z, 1] <- S0 / n.zones

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~           population with catch        ~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
for(year in 2 : length(years.sim)){ # years with fishery
    ##~ Recruitment
    rec[, year]  <-  recruitment(t.rec.a, t.rec.b, spawners[, (year - 1)], exp(res.Rec[year]))
    ##~ Recruitment by zone after connectivity
    rec.shape  <- rec.shp(rec[, year], pela.mat)
    for(z in n.zones){
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

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##         Projections 50 years            ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

year.proj      <- 50
seq.years      <- seq(from=2012, length=50)
rec.proj       <- matrix(NA, ncol= year.proj, nrow = n.zones)
b.catch.proj   <- v.biomass.proj <- spawners.proj  <- rec.proj
## Array for the abundance and the biomass
abundance.proj <- array(0, dim = c(length(zones), year.proj, length(size.vec)))
aamig.proj     <- catch.proj <- alive.proj <- abundance.proj
## Autocorrelation recruitment
res.Rec.proj  <- aut.rd(years = 50, sigmaz = abs(mean(tail(res.Rec, 10))), phi = 0.95, seed = 800)
cat(t(res.Rec.proj))
plot(exp(res.Rec), type = 'l')
lines(exp(res.Rec.proj), col = 2, lty = 2)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~             Projection                       ~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##~ First Year
rec.proj[, 1] <- recruitment(t.rec.a, t.rec.b, spawners[, length(years.sim)], exp(res.Rec[length(years.sim)]))
rec.shape     <- rec.shp(rec.proj[, 1], pela.mat)

for(z in 1 : n.zones){ # loop for the zones
        abundance.proj[z, 1, ] <- abund(alive[z, length(years.sim), ], normal.t.matrix, n.rec.pat, rec.shape[z],)
        v.biomass.proj[z, 1]   <- sum(abundance.proj[z, 1, ] *  w.l * f.selec) / 1000000
        spawners.proj[z, 1]    <- spawn(fecundity, maturity, abundance.proj[z, 1, ])
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~           population with catch        ~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
for(proj in 2 : year.proj){ # years with fishery
    ##~ Recruitment
    rec.proj[, proj]  <- recruitment(t.rec.a, t.rec.b, spawners.proj[, (proj - 1)], exp(res.Rec.proj[proj]))
    ##~ Recruitment by zone after connectivity
    rec.shape        <- rec.shp(rec.proj[, proj], pela.mat)
        for(z in n.zones){
        ## Abundance after Benthic migration
        aamig.proj[z, proj, ]     <- benth.mig(abundance[, proj-1, ], bent.mat, mig.pat, z)
        ## Alive after catch, Natural moratality and migration
        alive.proj[z, proj, ]     <- alive.f(aamig[, proj, ], M, f.selec, f.cur[proj - 1], h.effort, z)
        ## Abundance
        abundance.proj[z, proj, ] <- abund(alive.proj[z, proj, ], normal.t.matrix, n.rec.pat, rec.shape[z])
        ## Catch
        catch.proj[z, proj, ]     <- catch.f(abundance[, proj, ], f.selec, f.cur[proj], h.effort, M, z)
        ## Biomass catch
        b.catch.proj[z, proj]     <- sum(catch.proj[z, proj, ] * w.l) / 1000000
        ## Spawners
        spawners.proj[z, proj]    <- spawn(fecundity, maturity, abundance.proj[z, proj, ])
        ##~ Biomass
        v.biomass.proj[z, proj]   <- sum(abundance.proj[z, proj, ] *  w.l * f.selec) / 1000000
    }
}



par(mar = c(1, 1, 1, 1), oma = c(1, 1, 1, 1), mfrow = c(2, 2))
plot(years.sim, b.catch[1 , ], type = 'o', xlim = c(1900, 2060))
lines(years.sim, est$Catch, col = 2)
lines(seq.years, b.catch.proj[, ], col = 'blue')
plot(years.sim, v.biomass[1, ], type = 'o', xlim = c(1900, 2060))
lines(years.sim, est$Biomass[-112], col = 2)
lines(seq.years, v.biomass.proj[1, ], col = 'blue')

plot(years.sim, log(rec[1, ]), type = 'o', xlim = c(1900, 2060))
rec.est <- log(est$Rec * c(1,exp(est$Recruitment_Residuals)[-111]))
lines(years.sim, rec.est[-112], col = 2)
lines(seq.years, log(rec.proj))

plot(years.sim, log(spawners[1 , ]), type = 'o', xlim = c(1900, 2060))
lines(years.sim, log(est$Spawners[-112]), col = 2)
lines(seq.years, log(spawners.proj))
