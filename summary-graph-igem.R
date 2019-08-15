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
#      strain.list.csv
#
##############################################################

library(tidyverse)
library(plotly)
library(gridExtra)
library(cowplot)
library(optparse)
library(xtable)
library(ggpubr)
library(forcats)
library(ggridges)
library(plyr)

input.file.string =
"exp005.rates.summary.csv,exp006.rates.summary.csv,exp007.rates.summary.csv,
exp008.rates.summary.csv,exp009.rates.summary.csv,exp010.rates.summary.csv,
exp011.rates.summary.csv,exp012.rates.summary.csv,exp014.rates.summary.csv,
exp015.rates.summary.csv,exp017.rates.summary.csv,exp018.rates.summary.csv,
exp020.rates.summary.csv,exp021.rates.summary.csv,exp022.rates.summary.csv,
exp024.rates.summary.csv,exp027.rates.summary.csv,exp028.rates.summary.csv"

strain.file.string = "strain.list.csv"
#TODO: allow for this to be called from command line
  
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

#Clean and subset data
calibrations <- subset(fit.data, strain == "JEB1204"| strain == "JEB1205"|
                        strain =="JEB1206"| strain == "JEB1207"| strain =="JEB1208")
other_strains <- subset(fit.data, strain != "JEB1204"& strain != "JEB1205"&
                          strain !="JEB1206"& strain != "JEB1207"& strain !="JEB1208")
strain.list = read.csv(strain.file.string)
strain.data <- inner_join(other_strains, strain.list, by = "strain")
control.data <- mutate(calibrations, part.type = "Control")
total.data <- bind_rows(strain.data, control.data)
total.data$part.type = as.factor(total.data$part.type)
total.data$part.type<- recode_factor(total.data$part.type, DNA = "Part", Tag="Part", Regulatory="Part", Terminator="Part", RBS="Part",
                                     Protein_Domain="Part", Coding="Part", RNA="Part", Translational_Unit="Part", Generator = "Device",
                                     Signaling = "Device", Reporter = "Device", Measurement ="Device", 
                                     Device = "Device", Signalling = "Device", Intermediate = "Intermediate", Composite = "Intermediate", .default = "Other")
levels(total.data$part.type)
strain.data <- semi_join(total.data, strain.data, by = "strain")
clean.control.data <- subset(calibrations, growth.rate.sd <0.2)
clean.total.data <- subset(fit.data, growth.rate.sd < 0.2)
clean.data <- subset(strain.data, growth.rate.sd < 0.2)
dirty.data <- subset(total.data, growth.rate.sd > 0.2)
write_csv(dirty.data, paste0(output.prefix, ".outliers.csv"))

#Create regression
fit_fixed_zero= lm(GFP.rate~growth.rate + 0, clean.control.data)
slope_fixed_zero = coef(fit_fixed_zero)

#TODO: need to create test to assess duplicate difference for retention or removal of data

# Plot showing all strains burden vs growth 
burdenVsGrowthPlot = ggplot(clean.data, aes_(x=as.name(paste0(readings[1], ".rate")), y=as.name(paste0(readings[2], ".rate")), color = as.name("part.type")))  +
  geom_errorbarh(aes(xmin=growth.rate-growth.rate.sd, xmax=growth.rate+growth.rate.sd), height=0) +
  geom_errorbar(aes(ymin=GFP.rate-GFP.rate.sd, ymax=GFP.rate+GFP.rate.sd), width=0) + 
  geom_point(size=2)  +
  scale_x_continuous(limits = c(0, max(clean.data$growth.rate+clean.data$growth.rate.sd))) + 
  scale_y_continuous(limits = c(0, max(clean.data$GFP.rate+clean.data$GFP.rate.sd))) + 
  geom_abline(intercept=0, slope = slope_fixed_zero) +
  geom_vline(xintercept = mean(clean.total.data$growth.rate), linetype = "dashed") +
  NULL

burdenVsGrowthPlot

ggsave(paste0(output.prefix, ".burden_vs_growth_rates.pdf"))

burdenVsGrowthPlotly <- ggplotly(burdenVsGrowthPlot)
htmlwidgets::saveWidget(as_widget(burdenVsGrowthPlotly), paste0(output.prefix, ".burden_vs_growth_rates.html"))

# We need to convert NA isolates to a factor number.
clean.data$isolate=factor(clean.data$isolate, levels=c(1,levels(clean.data$isolate)))
clean.data$isolate[is.na(clean.data$isolate)] = 1
clean.data$isolate=as.factor(clean.data$isolate)

# Plot growth rate for every strain
growthRatePlot = ggplot(clean.data, aes_(x= reorder(clean.data$strain, -clean.data$growth.rate), y=as.name(paste0(readings[1], ".rate")), fill=as.name("part.type")))  +  
  geom_bar(size=3, stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=growth.rate-growth.rate.sd, ymax=growth.rate+growth.rate.sd), position=position_dodge()) + 
  scale_y_continuous(limits = c(0, max(clean.data$growth.rate+clean.data$growth.rate.sd))) +
  geom_hline(yintercept = mean(clean.total.data$growth.rate)) +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  NULL

growthRatePlot

ggsave(paste0(output.prefix, ".growth_rates.pdf"))

growthRatePlotly <- ggplotly(growthRatePlot)
htmlwidgets::saveWidget(as_widget(growthRatePlotly), paste0(output.prefix, ".growth_rates.html"))

#Plot distribution of growth rates
growthRateHist = ggplot(clean.data, aes(x=growth.rate, fill = part.type)) +
  stat_bin(bins = 45)
growthRateHist

ggsave(paste0(output.prefix, ".growth_rate_distribution.pdf"))

growthRateHistPlotly <- ggplotly(growthRateHist)
htmlwidgets::saveWidget(as_widget(growthRateHistPlotly), paste0(output.prefix, ".growth_rate_distribution.html"))

#Plot growth rate density per part type
growthRateDensity <- ggplot(data = clean.data, aes(x = growth.rate, y = stat(count), group = part.type, fill = part.type)) +
                              geom_density_ridges2(alpha=0.7, scale=1.5) +
                              geom_vline(xintercept = mean(clean.total.data$growth.rate), linetype = "dashed")+
                              geom_rug(sides = "b")
growthRateDensity

ggsave(paste0(output.prefix, ".growth_rate_density.pdf"))

growthRateDensityPlotly <- ggplotly(growthRateDensity)
htmlwidgets::saveWidget(as_widget(growthRateDensityPlotly), paste0(output.prefix, ".growth_rate_density.html"))

#Assess normality
growthRateTotalDensity <- ggdensity(data = clean.data, x = "growth.rate", add = "mean", rug = TRUE)
growthRateTotalDensity

ggsave(paste0(output.prefix, ".growth_rate_total_density.pdf"))

growthRateTotalDensityPlotly <- ggplotly(growthRateTotalDensity)
htmlwidgets::saveWidget(as_widget(growthRateTotalDensityPlotly), paste0(output.prefix, ".growth_rate_total_density.html"))

shapiro.test(strain.data$growth.rate)
