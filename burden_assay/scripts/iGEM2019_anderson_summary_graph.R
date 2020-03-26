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

suppressMessages(library(tidyverse))
suppressMessages(library(plotly))
suppressMessages(library(gridExtra))
suppressMessages(library(cowplot))
suppressMessages(library(optparse))
suppressMessages(library(xtable))
suppressMessages(library(ggrepel))

input.file.string = "exp_057.rates.summary.csv"

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

readings = c("growth", "GFP", "RFP") 

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
if (!("isolate" %in% colnames(all.data))) {
  all.data$isolate=rep("1", nrow(all.data))
}
all.data$isolate=as.factor(all.data$isolate)

# fit line and show entire range

fit.data = all.data
if("fit" %in% colnames(fit.data)) {
  fit.data = this.data %>% filter(fit==1)
}

fit_fixed_zero = lm(GFP.rate~growth.rate + 0, fit.data)
slope_fixed_zero = coef(fit_fixed_zero)


# plot showing all strains burden vs growth 
burdenVsgrowthPlot = ggplot(all.data, aes_(x=as.name(paste0(readings[1], ".rate")), y=as.name(paste0(readings[2], ".rate")), color=all.data$promoter.strength, shape=as.name("isolate")))  +
  geom_errorbarh(aes(xmin=growth.rate-growth.rate.sd, xmax=growth.rate+growth.rate.sd), height=0) +
  geom_errorbar(aes(ymin=GFP.rate-GFP.rate.sd, ymax=GFP.rate+GFP.rate.sd), width=0) + 
  geom_point(size=5, shape=16)  +
  scale_x_continuous(limits = c(0, max(all.data$growth.rate+all.data$growth.rate.sd))) + 
  scale_y_continuous(limits = c(0, max(all.data$GFP.rate+all.data$GFP.rate.sd))) + 
  geom_abline(intercept=0, slope = slope_fixed_zero) +
  scale_color_gradient(low = "lightpink", high = "red3", na.value = "grey") +
  labs(color = "Promoter Strength") + 
  geom_label_repel(aes(label = all.data$strain),size = 2.5,box.padding   = 0.35,point.padding = 0.5,segment.color = 'grey50') +
  ggtitle("Anderson Series Burden Plot") + xlab("Growth Rate") + ylab("GFP Rate")
  NULL
burdenVsgrowthPlot

burdenVsgrowthPlotly <- ggplotly(burdenVsgrowthPlot)
htmlwidgets::saveWidget(as_widget(burdenVsgrowthPlotly), paste0(output.prefix, ".burden_vs_growth_rates.html"))
ggsave(paste0(output.prefix, ".burden_vs_growth_rates.pdf"))

# We need to convert NA isolates to a factor number.
all.data$isolate=factor(all.data$isolate)
all.data$isolate[is.na(all.data$isolate)] = levels(all.data$isolate)[1]

all.data$isolate=as.factor(all.data$isolate)


# growth BAR PLOT 
growthRatePlot = ggplot(all.data, aes(x = reorder(all.data$strain, -(all.data$growth.rate)), y=all.data$growth.rate, fill=all.data$promoter.strength)) +  
  geom_bar(size=3, stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=growth.rate-growth.rate.sd, ymax=growth.rate+growth.rate.sd), position=position_dodge()) + 
  scale_y_continuous(limits = c(0, max(all.data$growth.rate+all.data$growth.rate.sd))) +   
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 12)) +
  scale_fill_gradient(low = "lightpink", high = "red3", na.value = "gray") +
  #scale_color_viridis_c(direction=-1) +
  labs(fill = "Promoter Strength") + 
  ggtitle("Anderson Series Growth Rates") + xlab("Strain") + ylab("Growth Rate") 
  NULL	
growthRatePlot 
ggsave(paste0(output.prefix, ".growth_rates.pdf"))
growthRatePlotly <- ggplotly(growthRatePlot)
htmlwidgets::saveWidget(as_widget(growthRatePlotly), paste0(output.prefix, ".growth_rates.html"))


#RFP Plot
rfpRatePlot = ggplot(all.data, aes(x = reorder(all.data$strain, -(all.data$other.rate)), y=all.data$other.rate, fill=all.data$promoter.strength)) +  
  geom_bar(size=3, stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=other.rate-other.rate.sd, ymax=other.rate+other.rate.sd), position=position_dodge()) + 
  scale_y_continuous(limits = c(0, max(all.data$other.rate+all.data$other.rate.sd))) +   
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 12)) +
  scale_fill_gradient(low = "lightpink", high = "red3", na.value = "gray") +
  #scale_color_viridis_c(direction=-1) +
  labs(fill = "Promoter Strength") + 
  ggtitle("RFP Rates") + xlab("Strain") + ylab("RFP Rate")
NULL	
rfpRatePlot 



