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
pela.mat <- read.dat('data/conect_mat.dat')[[2]]
#pela.mat <- diag(1, 8) ## Matrix to test the model, with one population

## Benthic connectivity matrix
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bent.mat <- read.dat('data/conect_mat.dat')[[5]]
#bent.mat <- diag(0, 8) ## Matrix to test the model, with one population

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

R0         <- est$Rec[1] / 8            # Virginal Recruitment
size.vec   <- with(datos, seq(from = lowersize, to = bigsize, by = sizedelta))
res.Rec    <- log(est$Recruitment_Residuals)             # Recruitment residual
Fs         <- with(param, exp(c(logF_1, logF_2, logF_3, logF_4)))


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##      Equation and models               ##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
start.year      <- 1901
end.year        <- 2016
fecundity       <- fec(a.fec, b.fec, size.vec)                                       # Fecundity
maturity        <- matu(a.mat, b.mat, size.vec)                                      # Maturity
w.l             <- alometry.f(a.wl, b.wl, size.vec)                                  # Alometric relationship
f.selec         <- selectivity.f(f.l50, f.l95, size.vec)                             # Fishery Selectivity
p.selec         <- selectivity.f(p.l50, p.l95, size.vec)                             # Trap Selectivity
normal.t.matrix <- trans.matrix(size.vec, k.g = 0.06985, linf = 213.495, beta.g = 1) # Normalize transition matrix
mig.pat         <- c(rep(0.3, 28), rep(0.6, 8), rep(1, 67))                          # Migration Pattern (Puerulus - juvenil - adult)
n.rec.pat       <- rec.pat(size.vec, beta = 0.3, alfa = 45)                          # Normalize recruitment pattern
h.effort        <- c(0.84, 0.24, 2.32, 0.6, 0.72, 0.8, 1.08, 1.4)
f.cur           <- c(rep(Fs[1], 29), rep(Fs[2], 51), rep(Fs[3], 19), rep(Fs[4], 12))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##          Simulation                    ##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
n.zones   <- 8
years.sim <- start.year : end.year
pcll   <- which(years.sim %in% datos$yearCLL)
pcpue  <- which(years.sim %in% datos$yearCPUE)
pcatch <- which(years.sim %in% datos$year)

n.years   <- length(years.sim)

## reclut    <- matrix(NA, ncol= length(years.sim), nrow = n.zones)
## rec       <- matrix(NA, ncol= length(years.sim), nrow = n.zones)
## rec2      <- matrix(NA, ncol= length(years.sim), nrow = n.zones)
## rec3      <- matrix(NA, ncol= length(years.sim), nrow = n.zones)
## b.catch   <- v.biomass <- spawners  <- rec
## ## Array for the abundance and the biomass
## abundance <- array(0, dim = c(length(zones), length(years.sim), length(size.vec)))
## aamig2    <- aamig <- catch <- alive <- abundance


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~          Estimation function       ~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Run model
source('simulation.R')

                                        # I need to deffine the initial condition
qCPUE   <- log(1e-7)
log.Rec <- c(res.Rec, rep(1, 5))
Par     <- c(log(R0), qCPUE, log(c(Fs, 1)), log.Rec)
## debug(hindcast)
## hindcast(Par, M = M, stpnss.h = stpnss.h, n.years = n.years, maturity = maturity, w.l = w.l,
##                           f.selec = f.selec, fecundity = fecundity, pela.mat = pela.mat, bent.mat = bent.mat, n.zones = 8,
##                           mig.pat = mig.pat, n.rec.pat = n.rec.pat, normal.t.matrix = normal.t.matrix)#
## library(compiler)
## compile(hindcast)
result <- optim(par = Par, fn = hindcast, M = M, stpnss.h = stpnss.h, n.years = n.years, maturity = maturity, w.l = w.l,
                          f.selec = f.selec, trap.s = p.selec, fecundity = fecundity, pela.mat = pela.mat, bent.mat = bent.mat, n.zones = 8,
                          mig.pat = mig.pat, n.rec.pat = n.rec.pat, normal.t.matrix = normal.t.matrix, method = "BFGS")

#save(result, file = 'Estimation.RData')
load(file = 'Estimation.RData')
R0      <- exp(result$par[1])
qCPUE   <- exp(result$par[2])
f.cur   <- c(rep(exp(result$par[3]), 29), rep(exp(result$par[4]), 51), rep(exp(result$par[5]), 19), rep(exp(result$par[6]), 12), rep(exp(result$par[7]), 5))
res.Rec2 <- result$par[8 : 123]

cvs <- (0.2, 0.05, 100, 0.6) ## csv for cpue, catch, lfd and recdev
output.simu <- simulation(M = M, R0 = R0, qCPUE = log(qCPUE), stpnss.h = stpnss.h, n.years = n.years, maturity = maturity, f.cur = f.cur, w.l = w.l,
                          f.selec = f.selec, trap.s = p.selec, fecundity = fecundity, pela.mat = pela.mat, bent.mat = bent.mat, n.zones = 8,
                          mig.pat = mig.pat, n.rec.pat = n.rec.pat, normal.t.matrix = normal.t.matrix, res.Rec = res.Rec2,
                          projection = FALSE)

x = data.frame(rec = unlist(res.Rec2))
type(x)
write(t(unlist(x)) , 'out.csv')
dim(output.simu$Catch)
par(mfrow = c(2, 2))
# Total Catch
plot(years.sim,datos$total_annual_catch)
lines(years.sim, colSums(output.simu$B.catch[1 : 8, ]), col = 2)
## CPP
plot(datos$cpue)
lines(output.simu$CPUE.cpp)
## Recruitment
plot(colSums(output.simu$V.biomass))
## Recruitment deviations
plot(log.Rec)
lines(exp(res.Rec2))

datos$cpue

names(output.simu)

par(mfrow = c(3, 4))
for(cll in 1 :  length(pcll)){
    plot(size.vec, datos$cll[cll, ])
    lines(size.vec, output.simu$Trap.lfd[cll, ], col = 2)
    abline(v = 115, lty = 2)
}

sum(output.simu$Trap.lfd[1, ])

dim(output.simu$Trap.lfd)

plot(years.sim, colSums(output.simu$B.catch[1 : 8, ]))


##-------------------
## Projections
##------------------
summary(output.simu)
#length(f.cur)
pars <- output.simu$Parameters
## Abundance of the last year of simulation
abun.sim     <- output.simu$Abundance[1 : 8, 116, ]
## Last f for the simulation
f.end        <- f.cur[116]
## last recruitment
res.rec.end  <- res.Rec[116]
## Years of projection
n.years      <- 100
seq.years    <- seq(from = 2012, length = n.years)
## F-projection
f.proj <- rep(f.end, n.years)
## Residual deviations
res.Rec.proj <- aut.rd(years = n.years, sigmaz = abs(mean(tail(res.Rec, 10))), phi = 0.95, seed = 800)
## projection
output.proj  <- simulation(M = M, R0 = R0, stpnss.h = stpnss.h, n.years = n.years, maturity = maturity, f.cur = f.proj, w.l = w.l,
                          f.selec = f.selec, fecundity = fecundity, pela.mat = pela.mat, bent.mat = bent.mat, n.zones = 8,
                          mig.pat = mig.pat, n.rec.pat = n.rec.pat, normal.t.matrix = normal.t.matrix, res.Rec = res.Rec.proj,
                          projection = TRUE, t.rec.a = pars$t.rec.a, t.rec.b = pars$t.rec.b, abun.sim = abun.sim, f.end = f.end,
                          res.rec.end = res.rec.end)





## spgred   <- rgb(t(col2rgb("red")), alpha=10, maxColorValue=255)
## plot(log(res.Rec  -  min(res.Rec)), type = 'l', ylim = c( - 4,  4))
## for(i in 1 : 1000){
##     res.Rec.proj <- aut.rd(years = n.years, sigmaz = abs(mean(tail(res.Rec, 10))), phi = 0.95, seed = i)
##     lines(res.Rec.proj, col = spgred)
## }


## summary(output.proj)

## output.proj$Abundance
datos

## plot proportion by Age


## Plots
par(mar = c(1, 1, 1, 1), oma = c(1, 1, 1, 1), mfrow = c(2, 2))
with(output.simu, plot(years.sim, colSums(B.catch[1 : 4 , ]), type = 'o', xlim = c(1900, 2060)))

with(output.proj, lines(seq.years, colSums(B.catch[1 : 4, ]), col = 'blue'))

with(output.simu, plot(years.sim, colSums(V.biomass[1 : 4, ]), type = 'o', xlim = c(1900, 2060)))
with(output.proj, lines(seq.years, colSums(V.biomass[1 : 4, ]), col = 'blue'))

with(output.simu, plot(years.sim, log(colSums(Rec[1:4, ])), type = 'o', xlim = c(1900, 2060)))
with(output.proj, lines(seq.years, log(colSums(Rec[1 : 4, ])), col = 4))

with(output.simu, plot(years.sim, log(colSums(Spawners[1 : 4 , ])), type = 'o', xlim = c(1900, 2060)))
with(output.proj, lines(seq.years, log(colSums(Spawners[1 : 4, ])), col = 4))




















##~ initial recruitment
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

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~           population with catch        ~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
for(year in 2 : length(years.sim)){ # years with fishery
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
        B.catch[z, year]     <- sum(catch[z, year, ] * w.l) / 1000000
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
B.catch.proj   <- v.biomass.proj <- spawners.proj  <- rec.proj
## Array for the abundance and the biomass
abundance.proj <- array(0, dim = c(length(zones), year.proj, length(size.vec)))
aamig.proj     <- catch.proj <- alive.proj <- abundance.proj
## Autocorrelation recruitment
res.Rec.proj  <- aut.rd(years = 50, sigmaz = abs(mean(tail(res.Rec, 10))), phi = 0.95, seed = 800)
## Fishing mortality
f.proj <- rep(tail(f.cur, 1), year.proj)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~             Projection                       ~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##~ First Year
rec[, 1] <- recruitment(t.rec.a, t.rec.b, spawners[, length(years.sim)], exp(res.Rec[length(years.sim)]))
rec.shape     <- rec.shp(rec.proj[, 1], pela.mat)

for(z in 1 : n.zones){ # loop for the zones
        abundance.proj[z, 1, ] <- abund(alive[z, length(years.sim), ], normal.t.matrix, n.rec.pat, rec.shape[z])
        v.biomass.proj[z, 1]   <- sum(abundance.proj[z, 1, ] *  w.l * f.selec) / 1000000
        catch.proj[z, 1, ]     <- catch.f(abundance.proj[, 1, ], f.selec, f.proj[1], h.effort, M, z)
        B.catch.proj[z, 1]     <- sum(catch.proj[z, 1, ] * w.l) / 1000000
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

        for(z in 1 : n.zones){
        ## Abundance after Benthic migration
        aamig.proj[z, proj, ]     <- benth.mig(abundance[, proj - 1, ], bent.mat, mig.pat, z)
        ## Alive after catch, Natural moratality and migration
        alive.proj[z, proj, ]     <- alive.f(aamig[, proj, ], M, f.selec, f.proj[proj - 1], h.effort, z)
        ## Abundance
        abundance.proj[z, proj, ] <- abund(alive.proj[z, proj, ], normal.t.matrix, n.rec.pat, rec.shape[z])
        ## Catch
        catch.proj[z, proj, ]     <- catch.f(abundance.proj[, proj, ], f.selec, f.proj[proj], h.effort, M, z)
        ## Biomass catch
        B.catch.proj[z, proj]     <- sum(catch.proj[z, proj, ] * w.l) / 1000000
        ## Spawners
        spawners.proj[z, proj]    <- spawn(fecundity, maturity, abundance.proj[z, proj, ])
        ##~ Biomass
        v.biomass.proj[z, proj]   <- sum(abundance.proj[z, proj, ] *  w.l * f.selec) / 1000000
    }
}



x11()
par(mar = c(1, 1, 1, 1), oma = c(1, 1, 1, 1), mfrow = c(2, 2))
plot(years.sim, colSums(B.catch[1 : 4 , ]), type = 'o', xlim = c(1900, 2060))
lines(years.sim, est$Catch, col = 2)
lines(seq.years, colSums(B.catch.proj[1 : 4, ]), col = 'blue')
##lines(seq.years, B.catch.proj[1, ], col = 'blue')
plot(years.sim, colSums(v.biomass[1 : 4, ]), type = 'o', xlim = c(1900, 2060))
lines(years.sim, est$Biomass[-112], col = 2)
lines(seq.years, colSums(v.biomass.proj[1 : 4, ]), col = 'blue')

plot(years.sim, log(colSums(rec[1:4, ])), type = 'o', xlim = c(1900, 2060))
rec.est <- log(est$Rec * c(1, exp(est$Recruitment_Residuals)[-111]))
lines(years.sim, rec.est[-112], col = 2)
lines(seq.years, log(colSums(rec.proj[1 : 4, ])), col = 4)

plot(years.sim, log(colSums(spawners[1 : 4 , ])), type = 'o', xlim = c(1900, 2060))
lines(years.sim, log(est$Spawners[-112]), col = 2)
lines(seq.years, log(colSums(spawners.proj[1 : 4, ])), col = 4)

64
116
166
214

breaks <- c(0, 116.4591393657, 202.6449897483, 214.4175236044,500)
prop <- c(sum(est$Number[1, which(size.vec < 64)]),
          sum(est$Number[1, which(size.vec >= 64  & size.vec < 116)]),
          sum(est$Number[1, which(size.vec >= 116 & size.vec < 166)]),
          sum(est$Number[1, which(size.vec >= 166 & size.vec <=  214)]))
dec <- c(0.2, 0.02, 0.02, 0.002, 0.0002, 0.00002)
prop <- c(prop, tail(prop, 1) * dec)
prop <- prop / sum(prop)

/ sum(est$Number[1, ])



hist(est$Number[1,], breaks=breaks)
ls()




n.years=length(years.sim)
abun.sim <- abundance[1 : 8, 111, ]



f.end <- f.cur[111]
res.rec.end <- res.Rec[111]
projection=T
n.years=50



salida <- simulation(M = M, R0 = R0, stpnss.h = stpnss.h, n.years = n.years, maturity = maturity, w.l = w.l,
           f.selec = f.selec, fecundity = fecundity, pela.mat = pela.mat, bent.mat = bent.mat, n.zones = 8,
           n.rec.pat = n.rec.pat, normal.t.matrix = normal.t.matrix, res.Rec = res.Rec,
           projection = FALSE, t.rec.a = NULL, t.rec.b = NULL, abun.sim = NULL, f.end = NULL,
           res.rec.end = NULL)
