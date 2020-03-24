#!/usr/bin/env Rscript

# Load many different summary files to determine with fitting method gives lowest variation

library(tidyverse)
library(ggplot2)
library(cowplot)

if (!exists("input.path")) {
   suppressMessages(library(optparse))
   
   option_list = list(
      make_option(
         c("-i", "--input"), type="character", default=NULL, 
         help="Input path of merged summary CSV file to be analyzed.", metavar="merged.csv"
      ),
      make_option(
         c("-o", "--output"), type="character", default=NULL, 
         help="Base name of output. ", metavar="merged.csv"
      )
   )
   
   usage_string = paste(
      "burden_merge.R -i path_to_results -o merged.csv\n\n",
      sep = ""
   ) 
   
   opt_parser = OptionParser(usage=usage_string, option_list=option_list)
   opt = parse_args(opt_parser)
   
   if (is.null(opt$input) || is.null(opt$output)) {
      print_help(opt_parser)
      stop("You must supply the -i, -o arguments", call.=FALSE)
   }
   
   input.path = opt$input
   output.path = opt$output 
}

##############  Read in the input files

input.summary.file.names  = list.files(path = input.path)
all.data = data.frame()
for (this.summary.file.name in input.summary.file.names) {
   this.file.path = file.path(input.path,this.summary.file.name)
   this.data <- read.csv(this.file.path)
   
   
   this.data$setting=this.summary.file.name
   this.data$setting = sub(".csv", "", this.data$setting)
   
   #Let's move the setting column to the leftmost!
   this.data <- this.data[c("setting", colnames(this.data)[1:(length(colnames(this.data))-1)])]
   
   all.data = rbind(all.data, this.data)
}

all.data = all.data %>% filter(replicates>1)
all.data = all.data %>% filter(!is.na(all.data$growth.rate))
all.data = all.data %>% filter(!is.na(all.data$GFP.rate))

all.data$growth.rate.cv = all.data$growth.rate.sd / all.data$growth.rate
all.data$GFP.rate.cv = all.data$GFP.rate.sd / all.data$GFP.rate

### No graph

dir.create(output.path)

setting.data = all.data %>% group_by(setting) %>% summarize(n=n(), mean.growth.rate.cv = mean(growth.rate.cv), mean.GFP.rate.cv = mean(GFP.rate.cv))


write.csv(setting.data, file.path(output.path, paste0("setting_comparison.csv")))
print(setting.data)


p = ggplot(setting.data, aes(x=mean.growth.rate.cv, y=mean.GFP.rate.cv, color=setting)) + geom_point()
ggsave(file.path(output.path, paste0("variation.pdf")), p)


p = ggplot(setting.data, aes(x=setting, y=n)) + geom_bar(stat="identity") 
ggsave(file.path(output.path, paste0("fits.pdf")), p)