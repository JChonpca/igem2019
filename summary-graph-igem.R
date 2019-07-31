#!/usr/bin/env Rscript

##############################################################
#
# Using RStudio to run this code? Follow these steps:
#
# 1) Define the variable 'input.file.string' in your Rstudio shell 
#    and then all of the command-line option code will be 
#    skipped and you can run this script within RStudio!
#
#    Example command
#    input.file.string = "input1.csv,input2.csv"
#
# 2) Set the working directory to a folder on your computer
#    that contains files 'input1.csv' and 'input2.csv'
#
#   Expects these files to exist:
#      input1.csv
#      input2.csv
#
##############################################################

library(tidyverse)
library(plotly)
library(gridExtra)
library(cowplot)
library(optparse)
library(xtable)
library(ggrepel)

input.file.string = "exp006.rates.summary.csv,exp007.rates.summary.csv,exp008.rates.summary.csv,exp009.rates.summary.csv,exp010.rates.summay.csv,exp011.rates.summary.csv"

# display instructions at the command line to use this script 
if (!exists("input.file.string")) {
  
  option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL, 
                help="Input file. You may separate values by commas to analyze multiple input files together", metavar="input.csv")
  )
  
  usage_string = paste(
    "summary-graph.R -i input1,input2,input3\n\n",
    sep = ""
  ) 
  
  opt_parser = OptionParser(usage=usage_string, option_list=option_list)
  opt = parse_args(opt_parser)
  
  if (is.null(opt$input)) {
    print_help(opt_parser)
    stop("You must supply the -i|--input argument for the input", call.=FALSE)
  }
  
  input.file.string = opt$input
}

# Allow loading of multiple different results

input.file.names = strsplit(input.file.string,",")

readings = c("growth", "GFP") 

all.data = data.frame()
for (this.file.name in input.file.names[[1]]) {
# this loop iterates over each input file name and reads them into
# data frame this.data, then combines them into one data frame all.data

# allow multiple input file types (comma or tab separated files)  
  if (tools::file_ext(this.file.name) == "csv") {
    this.data <- read_csv(this.file.name)
    output.prefix= sub(".csv", "", input.file.names[[1]][1])
  } else if (tools::file_ext(this.file.name) == "tsv") {
    this.data <- read_tsv(this.file.name)
    output.prefix= sub(".tsv", "", input.file.names[[1]][1])
  } else {
    stop(paste0("Unrecognized file extension (must be csv or tsv) for file: ", this.file.name))
  }
  
  if("graph" %in% colnames(this.data)) {
    this.data = this.data %>% filter(graph==1)
  }
  all.data = rbind(all.data, this.data)
}


#convert isolate to factor, so its not treated as numeric and mess with plotting
all.data$isolate=as.factor(all.data$isolate)

# fit line and show entire range

fit.data = all.data
if("fit" %in% colnames(fit.data)) {
  fit.data = this.data %>% filter(fit==1)
}

#Find 90% confidence interval and subset data
calibrations <- subset(fit.data, strain == "JEB1204"| strain == "JEB1205"|
                        strain =="JEB1206"| strain == "JEB1207"| strain =="JEB1208")
fit_fixed_zero = lm(GFP.rate~growth.rate + 0, calibrations)
slope_fixed_zero = coef(fit_fixed_zero)
fitin = lm(GFP.rate~growth.rate, calibrations)
confident<- confint(fitin, 'growth.rate', level = 0.90)
oddballs <- subset(fit.data, GFP.rate > confident[2] | GFP.rate < confident[1])

##Need to also account for (x,y) uncertainty, then create p-value correction for multiple testing##
#  sse <- sum((all.data$GFP.rate - fit_fixed_zero$fitted.values)^2)
#  mse <- sse / (n - 2)
#  t.val <- qt(0.975, n - 2) 
#  x_new <- 1:max(all.data$growth.rate)
#  y.fit <- x_new + slope_fixed_zero
#  se <- sqrt(sum((all.data$GFP.rate - y.fit)^2) / (n - 2)) * sqrt(1 / n + (all.data$growth.rate - mean(all.data$growth.rate))^2 / sum((all.data$growth.rate - mean(all.data$growth.rate))^2))
#  x_new2 <- 1:max(all.data$growth.rate)
#  y.fit2 <- x_new2 + slope_fixed_zero
#  slope.upper <- suppressWarnings(y.fit2 + t.val * se)
#  slope.lower <- suppressWarnings(y.fit2 - t.val * se)
#  bands <- data.frame(cbind(slope.lower, slope.upper))
#  outliers<- subset(all.data, GFP.rate < bands$slope.lower | GFP.rate > bands$slope.upper, select = growth.rate:strain)

##Residuals to determine %translational v. %transl.+other burden
# resid.data<-residuals(fit.data)
# resid.plot = ggplot(resid.data, aes(x= as.name(paste0(readings[1], ".rate)), y=as.name(paste0(readings[2], ".rate)),...
# resid.plot
##Shiny app/plotly to hover over points lying outside corrected confint that gives %other+%transl. burden ???

# plot showing all strains burden vs growth 
burdenVsGrowthPlot = ggplot(all.data, aes_(x=as.name(paste0(readings[1], ".rate")), y=as.name(paste0(readings[2], ".rate")), color=as.name("strain"), shape=as.name("isolate")))  +
  geom_errorbarh(aes(xmin=growth.rate-growth.rate.sd, xmax=growth.rate+growth.rate.sd), height=0) +
  geom_errorbar(aes(ymin=GFP.rate-GFP.rate.sd, ymax=GFP.rate+GFP.rate.sd), width=0) + 
  geom_point(size=2)  +
  scale_x_continuous(limits = c(0, max(all.data$growth.rate+all.data$growth.rate.sd))) + 
  scale_y_continuous(limits = c(0, max(all.data$GFP.rate+all.data$GFP.rate.sd))) + 
  geom_abline(intercept=0, slope = slope_fixed_zero) +
  geom_label_repel(data = oddballs, aes(label = strain), label.size = 0.10, nudge_x = -2, direction = "y", hjust = 1, force = 10, segment.alpha = 0.5, segment.colour = "grey" ) + theme(legend.position = "none") +
  NULL

ggsave(paste0(output.prefix, ".burden_vs_growth_rates.pdf"))

burdenVsGrowthPlotly <- ggplotly(burdenVsGrowthPlot)
htmlwidgets::saveWidget(as_widget(burdenVsGrowthPlotly), paste0(output.prefix, ".burden_vs_growth_rates.html"))

# We need to convert NA isolates to a factor number.
all.data$isolate=factor(all.data$isolate, levels=c(1,levels(all.data$isolate)))
all.data$isolate[is.na(all.data$isolate)] = 1

all.data$isolate=as.factor(all.data$isolate)

# Plot growth rate for every strain
growthRatePlot = ggplot(all.data, aes_(x=as.name("strain"), y=as.name(paste0(readings[1], ".rate")), fill=as.name("isolate")))  +  
	geom_bar(size=3, stat="identity", position=position_dodge()) +
	geom_errorbar(aes(ymin=growth.rate-growth.rate.sd, ymax=growth.rate+growth.rate.sd), position=position_dodge()) + 
	scale_y_continuous(limits = c(0, max(all.data$growth.rate+all.data$growth.rate.sd))) +
	NULL	

ggsave(paste0(output.prefix, ".growth_rates.pdf"))

growthRatePlotly <- ggplotly(growthRatePlot)
htmlwidgets::saveWidget(as_widget(growthRatePlotly), paste0(output.prefix, ".growth_rates.html"))

write_csv(oddballs, paste0(output.prefix, ".outliers.csv"))
