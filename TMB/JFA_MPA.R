require(TMB)
source('JFA_MPA_tools.R')
dat <- read.dat('JFA_MPA.dat')
##~~~~~~~~~~~~~~~~~~~~##
##  Known parameters  ##
##~~~~~~~~~~~~~~~~~~~~##

fec      <- dat$fec                  # Fecundity
mat      <- dat$maturity             # Maturity
wl       <- dat$weight              # Length - weight (Alometric relation)
sel.f    <- c(dat$l50s, dat$l95s)   # Selectivity of the fishery (almost knife - edge)
M        <- 0.18                     # Natural mortality
stpnss.h <- 0.7                      # Steepness
Lvec <- with(dat, seq(from = lowersize, to = bigsize, by = sizedelta))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##      Equation and models               ##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
start.year      <- dat$startyr
end.year        <- dat$end_year
Nyears          <- start.year : end.year
Nclass          <- length(Lvec)
maturity        <- mature(mat[1], mat[2], Lvec)                                      # Maturity
w.l             <- alometry.f(wl[1], wl[2], Lvec)                                  # Alometric relationship
n.rec.pat       <- rec.pat(Lvec, beta = 0.3, alpha = 45)                          # Normalize recruitment pattern
pcll            <- which(Nyears %in% dat$yearCLL)
pcpue           <- which(Nyears %in% dat$yearCPUE)
pcatch          <- which(Nyears %in% dat$year)
names(dat)


compile('JFA_MPA.cpp')
