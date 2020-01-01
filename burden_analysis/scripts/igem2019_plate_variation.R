#!/usr/bin/env Rscript

# Examines variation across different plates

library(tidyverse)
library(ggplot2)
library(cowplot)

## Uncomment to use in interactive mode
#input.file.path = "minOD_0.1_maxmethod_2.csv"
#output.path = "../03-plate-variation/minOD_0.1_maxmethod_2"

if (!exists("input.file.path")) {
   suppressMessages(library(optparse))
   
   option_list = list(
      make_option(
         c("-i", "--input"), type="character", default=NULL, 
         help="Input path of merged summary CSV file to be analyzed.", metavar="merged.csv"
      ),
      make_option(
         c("-o", "--output"), type="character", default="output", 
         help="Base name of output. ", metavar="output"
      )
   )
   
   usage_string = paste(
      "burden_merge.R -i summary.csv -o output_path\n\n",
      sep = ""
   ) 
   
   opt_parser = OptionParser(usage=usage_string, option_list=option_list)
   opt = parse_args(opt_parser)
   
   if (is.null(opt$input) || is.null(opt$output)) {
      print_help(opt_parser)
      stop("You must supply the -i, -o arguments", call.=FALSE)
   }
   
   input.file.path = opt$input
   output.path = opt$output
}


cat("Graphing plate variation statistics...\n")
cat("Input file of merged plates:  ", input.file.path,"\n")
cat("Output file path: ", output.path,"\n")

##############  Read in the input files

all.data = read.csv(input.file.path)
all.data$growth.rate.cv = all.data$growth.rate.sd / all.data$growth.rate

all.data$plate.strain = paste0(all.data$plate, "_", all.data$strain)


## Sort based on CV
dir.create(output.path)



##############  Make a general QC graph that highlights the controls

control.data=all.data %>% filter(strain %in% c("JEB1204", "JEB1205", "JEB1206", "JEB1207", "JEB1208"))

## Collect top ones
sorted.data = all.data %>% arrange(-growth.rate.cv)
top.data=sorted.data[1:20,]
top.data = top.data %>% arrange(plate)

p = ggplot(top.data , aes(x=plate.strain, y=growth.rate.cv)) + geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(file.path(output.path, paste0("most_variation.pdf")), p)


## Graph them all
for (this.plate in levels(all.data$plate)) {
   p = ggplot(all.data %>% filter(plate == this.plate), aes(x=strain, y=growth.rate.cv)) + 
      geom_bar(stat="identity") +
      ylim(0.4,1.6)
   ggsave(file.path(output.path, paste0(this.plate, ".cv.pdf")), p) 
}




##############  Make a  QC graph that connects the controls

control.data.means = control.data %>% group_by(plate, strain) %>% summarize(mean.growth.rate=mean(growth.rate))

p = ggplot(control.data, aes(x=plate, y=growth.rate, color=strain)) + 
  geom_point(data=control.data, size=1, aes(x=plate, y=growth.rate, color=strain)) +
  geom_line(data=control.data.means, aes(x=plate, y=mean.growth.rate, color=strain, group=strain)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(file.path(output.path, paste0("control.growth.rates.pdf")), p) 

