
library(ggplot2)
library(shiny)
library(tidyr)
library(dplyr)
library(deSolve)
library(adaptivetau)

################################### STOCHASTIC BURDEN SIMULATOR #####################################
#####   Stochastically plot fraction of engineered cells in a population over generation time  ######
#####   according to the burden and the mutation rate of the construct                         ######
#####################################################################################################

pars<- c(
  b  = 0.2,  # burden
  u  = 1e-7  # mutation rate
)


burdenmod <- function(t, y, p) {
  with(as.list(c(y, p)), {
    dEt <- ((1-b) * Et) - (u * (1-b) * Et)  # rate at which engineered cells (Et) decrease while growing
    dFt <- Ft + (u * (1-b) * Et)  # rate at which failed cells (Ft) accumulate
    return(list(c(dEt, dFt)))
  })
}

yini  <- c(Et = 1, Ft = 0)  # initial values of Et and Ft. Begin with one engineered cell.
times <- seq(0,200)  # arbitrary number of time steps

## Run model
out   <- ode(yini, times, burdenmod, pars)
  # the model simulates engineered cell growth, mutation, and depletion 

out <- data.frame(out)
out$fraction <- (out$Et/(out$Et+out$Ft))  # fraction of engineered cells
out$generation <- log2(out$Et+out$Ft)  # number of cell doublings
summary(out)

## Plot model
plot.mod= ggplot(out, aes(x=generation, y=fraction)) +
  geom_line() +
# ylim(c(0, 1)) +
#  xlim(c(0, 200)) +
  ggtitle("Modeling Burden")+
  NULL
plot.mod

######################################################################################
#one way to use adaptivetau below
stoburd <- function(x, params, t) {
  with(params, {
    return(c(((1-b) * x[1]) - (u * (1-b) * x[1]),  # engineered cells
             x[2] + (u * (1-b) * x[2])))  # failed cells
  })
}
par=list(b=0.2, u=10^-6, t = 200);
sbs=ssa.adaptivetau(c(1,0),
                  matrix(c(0,200, 1,0, 0,-1), nrow=2), stoburd, par, tf=200)
matplot(r[,"t"], sbs[,c("x1","x2")], type='l', xlab='generations', ylab='Counts')

#####################################################################################
#more explicit code here
b=0.1
u=10^-6
transitions = list(c(e = +1), # trans 1: e cells grow
                    c(e = -1, f = +1), # trans 2: e cells mutate into f cells
                    c(f = +1)) # trans 3: f cells grow

lvRateF <- function(x, params, t) {
    return(c(params$r * x["e"], # rate of e cells growing
               params$beta * x["e"],  # rate of e cells turning into f cells
               params$delta * x["f"]))  # rate of f cells growing
    }

set.seed(4) # set random number generator seed to be reproducible
simResults = ssa.adaptivetau(init.values = c(e = 1, f = 0),
                                transitions, lvRateF,
                                params = list(r=(1-b), beta=(u*(1-b)), delta=1),
                                tf=100)
matplot(simResults[,"time"], simResults[,c("e","f")], type='l',
             xlab='Time', ylab='Counts (log scale)', log='y')
legend("bottomleft", legend=c("e", "f"), lty=1:2, col=1:2)
