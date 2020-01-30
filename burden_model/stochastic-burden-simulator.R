
################################### STOCHASTIC BURDEN SIMULATOR #####################################
#####  Stochastically plot fraction of engineered cells in a population versus cell doublings  ######
#####   depending only on the mutation rate and burden                                         ######
#####################################################################################################

library(deSolve)
library(adaptivetau)
library(tidyverse)

################ Build ODE model #################

pars<- c(
  b = 0.4,  # burden value
  u = 1e-5  # mutation rate
)

burdenode <- function(t, y, p) {
  with(as.list(c(y, p)), {
    dEt <- ((1-b) * Et) - (u * (1-b) * Et)  # rate at which engineered cells (Et) grow/mutate
    dFt <- Ft + (u * (1-b) * Et)            # rate at which failed cells (Ft) accumulate/grow
    return(list(c(dEt, dFt)))
  })
}

yini  <- c(Et = 1, Ft = 0)  # initial values of Et and Ft. Begin with one engineered cell.
times <- seq(0,200, 0.1)    # time sequence for ODE solver

## Run model
out   <- ode(yini, times, burdenode, pars)
out <- data.frame(out)                    # save output into dataframe
out$fraction <- (out$Et/(out$Et+out$Ft))  # fraction of engineered cells remaining in population
out$generation <- log2(out$Et+out$Ft)     # number of cell doublings
summary(out)

## Plot model
ggplot(out, aes(x=generation, y=fraction)) +
  geom_line()

################# Build stochastic simulator  #################

transition.list = list(
  c(e = +1),           # trans 1: e cells grow
  c(f = +1),           # trans 3: f cells grow
  c(e = -1, f = +1)    # trans 2: e cells mutate into f cells
)

rates.function <- function(x, params, t) {
    return(c(
              (1 - params$b) * x["e"],             # rate of e cells growing
              x["f"],                              # rate of f cells growing
              (1 - params$b) * x["e"] * params$u   # rate of e cells mutating into f cells
          )) 
}


## Store all results in a dataframe for easier graphing
sim.results.df = data.frame()

## Loops over seeds and parameters to generate all possible graphs (takes ~5 mins) -GAM
for (this.u in 8:5) {
for (this.b in 1:4) {
for (this.seed in 1:20) {
  
  set.seed(this.seed) # set random number generator seed to be reproducible
  par=list(b=this.b/10, u=1*10^-this.u);
  sim.results = ssa.adaptivetau(
    init.values = c(e = 1, f = 0),
    transition.list,
    rates.function,
    par,
    tf=200
  )
  
  # reformat as a data frame and calculate generations
  this.sim.results.df =   data.frame (
    time = sim.results[,c("time")],
    generation = log2(sim.results[,c("e")] + sim.results[,c("f")]),
    e = sim.results[,c("e")],
    f = sim.results[,c("f")],
    fr.e = sim.results[,c("e")] / (sim.results[,c("e")] + sim.results[,c("f")]),
    seed = this.seed,
    b = this.b/10,
    u = 1*10^-this.u,
    method = "stochastic"
  )
  sim.results.df = rbind(sim.results.df, this.sim.results.df)
}
}
}

sim.results.df$seed = as.factor(sim.results.df$seed)

################### Generate graphs ################

## Compare ODE v stochastic methods
# filter results which match ODE parameters
b4.u5.sim = filter(sim.results.df, b == 0.4 & u == 1e-5)
comp.plot= ggplot() +
  geom_line(data = b4.u5.sim, aes(x=generation, y=fr.e, group=seed, color = "black"))+
  geom_line(data = out, aes(x=generation, y= fraction, color = "red"))+
  labs(title = "Comparison of Models", x = "Cell Doublings", y = "Fraction of Engineered Cells")+
  scale_colour_manual(values = c('black', 'red'), labels = c("Stochastic", "ODE"), name = "Type")+
  scale_x_continuous(limits = c(0,100))
comp.plot

## Show distribution of failure generations created in simulator
f5.sim.results = data.frame()
# filter fraction results rounding to 0.5
f5.sim.results = filter(sim.results.df, fr.e < 0.64 & fr.e > 0.45)
f5.sim.results$b = f5.sim.results$b*100
f5.sim.results$b = as.factor(f5.sim.results$b)

f5.plot= ggplot(f5.sim.results, aes(b, generation, group = b)) +
  geom_boxplot(aes(color = b), outlier.shape =1, outlier.alpha =0.1) + facet_grid(.~u)+ theme(legend.position = "none")+
  labs(title = "Stochastic Distribution Faceted by Mutation Rate", x = "Burden (%)", y = "Cell Doublings to 50% Failure")
f5.plot

