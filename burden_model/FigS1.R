
####  Code separated from main script to reproduce Fig. S1

################################### STOCHASTIC BURDEN SIMULATOR #####################################
#####  Stochastically plot fraction of engineered cells in a population versus cell doublings  ######
#####   depending only on the mutation rate and burden                                         ######
#####################################################################################################

library(deSolve)
library(adaptivetau)
library(tidyverse)

################ Build ODE model ####################

burdenode <- function(t, y, p) {
  with(as.list(c(y, p)), {
    dEt <- ((1-b) * Et) - (u * (1-b) * Et)  # rate at which engineered cells (Et) grow/mutate
    dFt <- Ft + (u * (1-b) * Et)            # rate at which failed cells (Ft) accumulate/grow
    return(list(c(dEt, dFt)))
  })
}

yini  <- c(Et = 1, Ft = 0)    # initial values of Et and  Ft. Begin with one engineered cell.
times <- seq(0,200)   # time sequence for ODE solver
log10uvec <- seq(-8,-4,by=1)
burdenvec <- seq(0.1,0.5,by=0.1)                     #vector of burden values, loops only by whole number
results<- vector(length(burdenvec), mode="list")  #store results
out<-data.frame()                                 #store in dataframe when combining parameters

for(u.i in seq_along(log10uvec)){
  for (burd.i in seq_along(burdenvec)){
    this.log10u = log10uvec[u.i]
    this.u = 10^this.log10u
    this.b = burdenvec[burd.i]
    results<-ode(yini, times, burdenode,
                 parms=c(b=this.b, u=this.u))
    df = as.data.frame(results)
    df$burden = this.b
    df$mutation = this.u
    out=rbind(out, df)
  }
}
out$burden<-as.numeric(out$burden)
out$fraction <- (out$Et/(out$Et+out$Ft))  #fraction of engineered cells remaining in population
out$generation <- log2(out$Et+out$Ft)     #number of cell doublings


################# Build stochastic simulator  ####################

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

for(u.i in seq_along(log10uvec)){
  for (burd.i in seq_along(burdenvec)){
    
    this.log10u = log10uvec[u.i]
    this.u = 10^this.log10u
    this.b = burdenvec[burd.i]
    cat(this.u, ", ", this.b, "\n")
    
    for (this.seed in 1:20) {
      
      set.seed(this.seed) # set random number generator seed to be reproducible
      par=list(b=this.b, u=this.u);
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
        b = this.b,
        u = this.u,
        method = "stochastic"
      )
      sim.results.df = rbind(sim.results.df, this.sim.results.df)
    }
  }
}

sim.results.df$seed = as.factor(sim.results.df$seed)

#tidy more data
out = out %>% rename(b = burden)
out = out %>% rename(u =mutation)
sim.results.df = sim.results.df %>% rename(fraction = fr.e)
out = filter(out, b != 0)

## Plot showing ODE v stochastic curves for all whole number parameters
fig1s <- ggplot()+
  geom_line(data = sim.results.df, aes(generation, fraction, group = seed))+
  geom_line(data = out, aes(generation, fraction), color = "red") +
  facet_grid(u ~ b)+
  scale_x_continuous(limits = c(0,100))+
  labs(title = "Curves Produced with Varying Parameters", x = "Cell Doublings", y = "Fraction of Engineered Cells") +
  theme_classic() +
  NULL
fig1s
