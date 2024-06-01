####  Code separated from main script to reproduce Fig. S1

library(deSolve)
library(adaptivetau)
library(tidyverse)

################################### ODE MODEL ###########################################
# Plot fraction of unmutated cells over time.
# Parameters:
#   burden (b), expressed as fractional decrease in growth rate
#   failure mutation rate (u), in units of mutations per cell doubling

simulation.time = 400
stochastic.replicates = 10000

ODE_failure_model_function <- function(t, y, p) {
  with(as.list(c(y, p)), {
    
    # dE(t) is the change in engineered (functional) cells over time.
    #   It is equal to the rate at which they increase due to replication
    #   minus the rate at which some offspring are converted to mutated cells.
    dEt <- ((1-b) * Et) - (u * (1-b) * Et)  
    
    # dF(t) is the change in failed (mutated) cells over time.
    #   It is equal to the rate at which they increase due to replication
    #   plus the rate at which some engineered cells are converted to mutated cells.
    dFt <- Ft + (u * (1-b) * Et)
    return(list(c(dEt, dFt)))
  })
}
yini  <- c(Et = 1, Ft = 0)    # initial values of Et and  Ft. Begin with one engineered cell.
times <- seq(0,simulation.time, by=0.1)   # time sequence for ODE solver
log10uvec <- seq(-8,-4,by=1)
burdenvec <- seq(0.1,0.5,by=0.1)                     #vector of burden values
results<- vector(length(burdenvec), mode="list")  #store results
ode.model.results<-data.frame()                                 #store in dataframe when combining parameters

for(u.i in seq_along(log10uvec)){
  for (burd.i in seq_along(burdenvec)){
    this.log10u = log10uvec[u.i]
    this.u = 10^this.log10u
    this.b = burdenvec[burd.i]
    results<-ode(yini, times, ODE_failure_model_function,
                 parms=c(b=this.b, u=this.u))
    df = as.data.frame(results)
    df$b = this.b
    df$u = this.u
    df$fraction <- (df$Et/(df$Et+df$Ft))
    df$total.cell.doublings <- log2(df$Et+df$Ft)
    ode.model.results=rbind(ode.model.results, df)
  }
}

## Plot Distributions of 50% cell divisions in each model
ode.model.cell.doublings.to.50 = ode.model.results %>% filter(fraction<=0.5) %>% group_by(b, u) %>% slice_head(n=1)
ode.model.cell.doublings.to.50$line.color = factor(ode.model.cell.doublings.to.50$b)

################################ STOCHASTIC MODEL #######################################

stochastic.failure.model.transitions = list(
  c(e = +1),           # trans 1: e cells grow
  c(f = +1),           # trans 3: f cells grow
  c(e = -1, f = +1)    # trans 2: e cells mutate into f cells
)

stochastic_failure_model_function <- function(x, params, t) {
  return(c(
    (1 - params$b) * x["e"],             # rate of e cells growing
    x["f"],                              # rate of f cells growing
    (1 - params$b) * x["e"] * params$u   # rate of e cells mutating into f cells
  )) 
}


## Store all results in a dataframe for easier graphing

stochastic.model.cell.doublings.to.50 = data.frame()
for(u.i in seq_along(log10uvec)){
  for (burd.i in seq_along(burdenvec)){
    
    this.log10u = log10uvec[u.i]
    this.u = 10^this.log10u
    this.b = burdenvec[burd.i]
    cat(this.u, ", ", this.b, "\n")
    
    for (this.seed in 1:stochastic.replicates) {
      
      if (this.seed %% 100 == 0) {
        cat("  replicate:", this.seed, "\n")
      }
      
      set.seed(this.seed) # set random number generator seed to be reproducible
      par=list(b=this.b, u=this.u);
      sim.results = ssa.adaptivetau(
        init.values = c(e = 1, f = 0),
        stochastic.failure.model.transitions,
        stochastic_failure_model_function,
        par,
        tf=simulation.time
      )
      
      # reformat as a data frame and calculate generations
      this.stochastic.model.results =   data.frame (
        time = sim.results[,c("time")],
        total.cell.doublings = log2(sim.results[,c("e")] + sim.results[,c("f")]),
        e = sim.results[,c("e")],
        f = sim.results[,c("f")],
        fraction = sim.results[,c("e")] / (sim.results[,c("e")] + sim.results[,c("f")]),
        seed = this.seed,
        b = this.b,
        u = this.u,
        method = "stochastic"
      )
      
      this.stochastic.model.cell.doublings.to.50 = this.stochastic.model.results %>% filter(fraction<=0.5) %>% slice_head(n=1)
      
      stochastic.model.cell.doublings.to.50 = rbind(stochastic.model.cell.doublings.to.50, this.stochastic.model.cell.doublings.to.50)
      
    }
  }
}
stochastic.model.cell.doublings.to.50$line.color = factor(stochastic.model.cell.doublings.to.50$b)

## Plot showing ODE v stochastic curves for all whole number parameters
fig1s <- ggplot(ode.model.cell.doublings.to.50) +
  geom_vline(xintercept = 22.93156857, color = "#6BAED5", size = 0.8) +
  geom_vline(xintercept = 34.21928095, color = "#4292C5", size = 0.8) +
  geom_vline(xintercept = 41.12617154, color = "#2271B5", size = 0.8) +
  geom_vline(xintercept = 56.15084952, color = "#194891", size = 0.8) +
 # geom_vline(aes(xintercept=total.cell.doublings, color=line.color), linetype="dashed", size = 0.8) +
  scale_color_manual(values = c("#ffcccc", "#ff9999", "#ff3232", "#cc0000", "#7f0000")) +
  stat_ecdf(data=stochastic.model.cell.doublings.to.50, aes(x=total.cell.doublings, group=b, color=line.color)) +
  facet_grid(rows=vars(-u)) +
  coord_cartesian(xlim=c(0,100)) +
  scale_x_continuous(breaks = seq(0, 100, 5), labels=c("0", "", "", "", "", "25", "", "", "", "", "50", "", "", "", "", "75", "", "", "", "", "100")) +
  
  scale_y_continuous() +

  labs(x = "Cell Divisions", y = "Fraction Engineered Cells") +
  theme_linedraw() +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  NULL
  
fig1s
ggsave("Fig3.pdf")

write_csv(stochastic.model.cell.doublings.to.50, "stochastic.model.cell.doublings.to.50.csv")
write_csv(ode.model.cell.doublings.to.50, "ode.model.cell.doublings.to.50.csv")
