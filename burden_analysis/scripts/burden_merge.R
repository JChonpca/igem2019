#!/usr/bin/env Rscript

library(tidyverse)


##############  Merge burden values from multiple plates

# Then, run...
metadata.file.path = NA

#debug
#input.path="02-output-all-merged"
#output.file.path = "test.csv"
#metadata.file.path = "igem2019_strain_metadata.csv"

if (!exists("input.path")) {
  suppressMessages(library(optparse))
  
  option_list = list(
    make_option(
      c("-i", "--input"), type="character", default=NULL, 
      help="Input path containing CSV files to be merged.", metavar="merged.csv"
    ),
    make_option(
      c("-o", "--output"), type="character", default=NULL, 
      help="Output merged CSV file.", metavar="merged.csv"
    ),
    make_option(
      c("-m", "--metadata"), type="character", default=NA, 
      help="Optional CSV file containing one line per strain analyzed. Must contain a column called 'GPF.interference.metadata'. All strains with this tag will have theif GFP rate set to NA.", metavar="strain_metadata.csv"
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
  output.file.path = opt$output
  metadata.file.path = opt$metadata
}

input.rates.file.names  = list.files(path = input.path)

all.data = data.frame()
for (this.file.name in input.rates.file.names) {
  #cat(this.file.name, "\n")
  # this loop iterates over each input file name and reads them into
  # data frame this.data, then combines them into one data frame all.data
  
  this.file.path = file.path(input.path,this.file.name)
  this.data <- read.csv(this.file.path)
  this.data$strain=as.character(this.data$strain)
  this.data$isolate=as.character(this.data$isolate)
  
  if ("wells" %in% colnames(this.data)) {
    this.data$wells=as.character(this.data$wells)
  }
  
  if ("well" %in% colnames(this.data)) {
    this.data$well=as.character(this.data$well)
  }
  
  if("graph" %in% colnames(this.data)) {
    this.data = this.data %>% filter(graph==1)
  }
  
  ### Two modes for "all" and "summary" rates files
  
  # We know we are using a summary if the "replicates" column exists:
  
  if ("replicates" %in% colnames(this.data)) {
    
    #we have to add other.rate if it is missing to join all data together
    if(!("other.rate" %in% colnames(this.data))) {
      this.data$other.rate=c(NA)
    }
    
    if(!("other.rate.sd" %in% colnames(this.data))) {
      this.data$other.rate.sd=c(NA)
    }
    
  } else {
    
    #we have to add other.rate if it is missing to join all data together
    if(!("other.rate" %in% colnames(this.data))) {
      this.data$other.rate=c(NA)
    }
    
    if(!("max.other.rate.time" %in% colnames(this.data))) {
      this.data$max.other.rate.time=c(NA)
    }
  }
  

  
  this.data$plate=this.file.name
  this.data$plate = sub(".rates.all.csv", "", this.data$plate)
  this.data$plate = sub(".rates.summary.csv", "", this.data$plate)
  #this.data$plate = sub("exp", "", this.data$plate)

  #Let's move the run column to the leftmost!
  this.data <- this.data[c("plate", colnames(this.data)[1:(length(colnames(this.data))-1)])]
  
  all.data = rbind(all.data, this.data)
}

############# Remove GFP measurements from contains.GFP parts

if (!is.na(metadata.file.path)) {
  cat("Setting GFP.rate to NA for strains marked with 'GFP.inteference' in metadata file:", metadata.file.path, "\n")
  metadata = read.csv(metadata.file.path)
  metadata$GFP.interference = as.logical(metadata$GFP.inteference)
  GPF.interference.metadata = metadata %>% filter(GFP.interference==T)
  all.data$GFP.rate[all.data$strain %in% GPF.interference.metadata$strain] = NA
}


############## Write out the merged file

write.csv(all.data, output.file.path, row.names=F)

