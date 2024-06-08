#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(gridExtra))
suppressMessages(library(cowplot))
suppressMessages(library(optparse))

## Run this from the main directory

## Code what you want plotted
biobricks =   c("BB205", "BB296", "BB272", "BB220",  "BB288", "BB298", "BB328")
experiments = c("exp009","exp028","exp020","exp050","exp020","exp022","exp011")
              
max_time=420
time_point_delta=40
OD_offset_for_window_graph = 0.75
OD_separation_for_window_graph = 0.02

plot_array = list()

for (i in 1:length(biobricks)) {
 
  this.biobrick = biobricks[i]
  this.experiment = experiments[i]
  
 # this_biobrick_curves_path = file.path("01-output-plate-fits", this.experiment, paste0(this.experiment,"-plots"),paste0(this.biobrick,"__NA.rates.csv"))
  
  this_experiment_all_rates_path =file.path("01-output-plate-fits", this.experiment, paste0(this.experiment,".rates.all.csv"))
  
  Y = read_csv(this_experiment_all_rates_path) 
  Y = Y %>% filter(strain==this.biobrick)
  Y$well=factor(Y$well)
  Y$plot_y= OD_offset_for_window_graph - OD_separation_for_window_graph * as.numeric(Y$well)

  this_experiment_add_values_path = file.path("01-output-plate-fits", this.experiment, paste0(this.experiment,".tidy.background.subtracted.measurements.csv"))

  X = read_csv(this_experiment_add_values_path) 
  
  X = X %>% filter(strain==this.biobrick)
  
  
  okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  p = ggplot(X, aes(x=time.min, y=OD, color=well)) + 
    #geom_point() + 
    geom_line(size=2) +
    scale_colour_manual(values=okabe) + 
    scale_x_continuous(limits=c(-0.01,max_time), breaks=c(0:7)*(60)) +
    scale_y_continuous(limits=c(-0.05,0.75), breaks=c(0:3)*(0.25)) +
    geom_segment(data=Y, aes(x=max.growth.rate.time-time_point_delta, y=plot_y, xend=max.growth.rate.time+time_point_delta, yend=plot_y, color=well)) +
    theme_cowplot(12) + 
    panel_border(color = "black") +
    theme(legend.position="none")
  
  plot_array[[i]] = p
}


q = gridExtra::grid.arrange(plot_array[[1]], plot_array[[2]], plot_array[[3]], 
                        plot_array[[4]], plot_array[[5]], plot_array[[6]],
                        ncol=2)


ggsave(file.path("12-OD-curves", "supplemental_OD_curves.pdf"), q, width=8, height=8)
