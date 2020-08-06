
library(tidyverse)
library(deSolve)
library(cowplot)

################################### ODE MODEL ###########################################
# Plot fraction of unmutated cells over time.
# Parameters:
#   burden (b), expressed as fractional decrease in growth rate
#   failure mutation rate (u), in units of mutations per cell doubling


ODE_failure_model <- function(t, y, p) {
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

# Create a data frame with the results for several different values of burden and a mutation rate of 1E-6
u = 1e-6
all.output = data.frame()
for (b in seq(from=0, to=0.5, by=0.1)) {
  pars<- c(
    b  = b,  # burden
    u  = u  # mutation rate
  )

  yini  <- c(Et = 1, Ft = 0)  # initial values of Et and Ft. Begin with one engineered cell.
  times <- seq(from=0, to=200, by=0.1)  # arbitrary number of time steps
  
  ## Run model
  out   <- ode(yini, times, ODE_failure_model, pars)
  
  out <- data.frame(out)
  out$fraction.engineered.cells <- (out$Et/(out$Et+out$Ft))  # fraction of engineered cells
  out$total.cell.doublings <- log2(out$Et+out$Ft)  # calculate the actual number of cell doublings
  summary(out)
  
  out$b = b
  out$u = u
  
  all.output = all.output %>% bind_rows(out)
}


all.output$line.color = factor(all.output$b)
  
## Plot model
plot.mod= ggplot(all.output, aes(x=total.cell.doublings, y=fraction.engineered.cells, group=b, color=line.color), clip="off") +
  geom_hline(aes(yintercept=0.5), linetype="dashed") +
  geom_vline(xintercept = 22.93156857, color = "#6BAED5", size = 0.8) +
  geom_vline(xintercept = 34.21928095, color = "#4292C5", size = 0.8) +
  geom_vline(xintercept = 41.12617154, color = "#2271B5", size = 0.8) +
  geom_vline(xintercept = 56.15084952, color = "#194891", size = 0.8) +
  geom_line(size=1.5) +
  scale_color_manual(values = c("black", "#ffcccc", "#ff9999", "#ff3232", "#cc0000", "#7f0000")) +
  scale_x_continuous(limits=c(0,100), expand = expansion(add = c(0, 0))) +
  scale_y_continuous(limits=c(0,1), expand = expansion(add = c(0, 0)), breaks = seq(0, 1, 0.25)) +
  theme_linedraw() +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  NULL
plot.mod

#This used to be Fig. 1C in a prior version of the paper
ggsave("Fig. 1C.pdf", plot = plot.mod)