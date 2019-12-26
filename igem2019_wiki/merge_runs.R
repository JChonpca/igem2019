##############  Read in the input files
# How to use this script:
#
# Set the input.path variable to the path containing folders of runs
#    input.path = "_rates_all_0.0_"
#
# Optionally, override the output.file.path assigned below
output.file.path = "rates_all_0.0_merged.csv"
# Then, run...

input.rates.file.names  = list.files(path = input.path)

all.data = data.frame()
for (this.file.name in input.rates.file.names) {
  # this loop iterates over each input file name and reads them into
  # data frame this.data, then combines them into one data frame all.data
  
  this.file.path = file.path(input.path,this.file.name)
  this.data <- read.csv(this.file.path)
  
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
  

  
  this.data$experiment=this.file.name
  this.data$experiment = sub(".rates.all.csv", "", this.data$experiment)
  this.data$experiment = sub(".rates.summary.csv", "", this.data$experiment)
  #this.data$experiment = sub("exp", "", this.data$experiment)

  #Let's move the run column to the leftmost!
  this.data <- this.data[c("experiment", colnames(this.data)[1:(length(colnames(this.data))-1)])]
  
  all.data = rbind(all.data, this.data)
}

############## Write out the merged file

write.csv(all.data, output.file.path, row.names=F)

