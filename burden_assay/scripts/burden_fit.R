#!/usr/bin/env Rscript

##############################################################
#
# Using RStudio to run this code? Follow these steps:
#
# 1) Define the variable 'input.prefix' in your Rstudio shell 
#    and then all of the command-line option code will be 
#    skipped and you can run this script within RStudio!
#
#    Example command
#      input.prefix = "my_exp"
#
# 2) Set the working directory to a folder on your computer
#    that contains comma separated value (CSV) input files. 
#    It should have a "measurements file saved from the 
#    platereader and metadata file that you created that
#    describes the samples (see full docs for format)
#    beginning with your prefix:
#    
#       'input.prefix.measurements.csv'
#       'input.prefix.metadata.csv'
#
#     You can also use tab-separated value (TSV) files that end in *.tsv
#
# 3) Set any of these variables ahead of time to override the defaults.
#
#    Example commands:
#      max.method = 1
#      output.prefix = "my_output"
#
##############################################################

suppressMessages(library(tidyverse))
suppressMessages(library(gridExtra))
suppressMessages(library(cowplot))


#debug
#input.prefix="input-plate-data/exp011/exp011"
#output.prefix="01-output-plate-fits/exp011"
#time.point.delta = 4
#FP.per.OD.at.time = T
#max.method = 1
#two.point.fit = F
#minimum.OD = 0.03

##############################################################
#### Error Handler
##############################################################
# This function will pause executionin Rstudio by entering an
# infinite loop so that you can see the actual fatal error

fatal_error <- function(msg) {
  
  message("======= FATAL ERROR=======")
  if (interactive()) {
    message(msg)
    message("Click the STOP sign to end after this fatal error")
    while(T) {}
  } else {
    message(msg)
    options("show.error.messages" = F)
    stop()
  }
  
}


##############################################################
#### Load command line options and set global parameters
##############################################################


## Defaults:
default.two.point.fit = F
default.minimum.OD = 0.0
default.maximum.time = 600
default.offset.readings.to.average.time.span.end = 60
default.growth.rate.method ="exponential"
default.fluorescence.rate.method = "exponential-with-production"

if (!exists("input.prefix")) {
  suppressMessages(library(optparse))
  option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL, 
                help="Input file prefix. Expects to find the files <input>.metadata.tsv and <input>.measurement.tsv", metavar="input.csv"),
    make_option(c("-o", "--output"), type="character", default=NULL, 
                help="Output file prefix. Output files of the form <output>.* will be created. If this option is not provided the input file prefix will be used.", metavar="output_prefix"),
    make_option(c("-m", "--max-method"), type="numeric", default=2, 
                help="Maximum picking method. 1=pick max growth rate time and use other rates at this same time. 2=pick maximum value of each curve at whatever time it occurs (may be different for growth rate and for fluorescence).", metavar="1/2"),
    make_option(c("--two-point-fit"), type="logical", action="store", default=default.two.point.fit, 
                help="Fit growth and fluorescence rates using just the first and last points in the intervals. The default without this flag is to use all time points in the interval."),
    make_option(c("-g", "--growth-rate-method"), type="character", action="store", default=NA, 
                help="Method used for fitting fluorescence production rates. Options are log-linear or exponential. Only works when --two-point-fit is FALSE. Fits to the exponential model may fail to converge. (Default = exponential)"),
    make_option(c("-f", "--fluorescence-rate-method"), type="character", action="store", default=NA, 
                help="Method used for fitting fluorescence production rates. Options are log-linear-per-OD or exponential-with-production. Only works when --two-point-fit is FALSE. Fits to the exp-production model may fail to converge. (Default = exponential-production)"),
    make_option(c("-t", "--time-point-delta"), type="numeric", default=3, 
                help="Time point delta used in calculating rates. This many points before and after the current point will be used for fitting rates. So, if you recorded data every 10 minutes and choose a value of 2 (the default), then 50 minute long windows spanning five data point are used for fitting rates.", metavar="int"),
    make_option(c("-d", "--min-OD"), type="numeric", default=default.minimum.OD, 
                help="Only use points with at least this minimum OD measurement", metavar="value"),
    make_option(c("-x", "--max-time"), type="numeric", default=default.maximum.time, 
                help="Only consider data from the beginning of growth up to this time", metavar="value"),
    make_option(c("-e", "--avg-end-time"), type="numeric", default=default.offset.readings.to.average.time.span.end, 
                help="Offset the readings of every well of a sample to the average value of all readings in all wells of a  sample over the interval from time zero up through and including this final time. Set to -1 to turn off this initial timepoint averaging.", metavar="value"),
    make_option(c("-b", "--background"), type="character", default=NA, 
                            help="Constant background measurements to subtract from measurements in place. Must be a comma-separated list with as many background values as there are measurement types. If this is not provided (the default), then the time course from BLANK wells is subtracted.", metavar="value1,value2,[value3]"),
     make_option(c("-p", "--no-plots"), type="logical", action="store_true", 
                help="Don't create plots for quicker execution.")
  )
  
  usage_string = paste(
    "burden.R -i input -o output",
    sep = ""
  ) 
  
  opt_parser = OptionParser(usage=usage_string, option_list=option_list);
  opt = parse_args(opt_parser);
  
  if (is.null(opt$input)) {
    print_help(opt_parser)
    stop("You must supply the -i|--input argument for the input", call.=FALSE)
  }
  input.prefix = opt$input
  
  if (!is.null(opt$output)) {
    output.prefix = opt$output
  }
  
  if (!is.null(opt$'max-method')) {
    max.method = opt$'max-method' 
  }
  
  if (!is.null(opt$'two-point-fit')) {
    two.point.fit = opt$'two-point-fit' 
  }
  
  if (!is.null(opt$'growth-rate-method' )) {
    growth.rate.method = opt$'growth-rate-method' 
  }
  
  if (!is.null(opt$'fluorescence-rate-method' )) {
    fluorescence.rate.method = opt$'fluorescence-rate-method' 
  }
  
  if (!is.null(opt$'time-point-delta')) {
    time.point.delta = opt$'time-point-delta' 
  }
  
  if (!is.null(opt$'max-time')) {
    maximum.time = opt$'max-time' 
  }
  
  if (!is.null(opt$'no-plots')) {
    no.plots = opt$'no-plots'
  }
  
  if (!is.null(opt$'min-OD')) {
    minimum.OD = opt$'min-OD'
  }
  
  if (!is.null(opt$'avg-end-time')) {
    offset.readings.to.average.time.span.end = opt$'avg-end-time'
  }
  
  if (!is.null(opt$'background')) {
    set.background.string = opt$'background'
  }
  
} else {
  # RStudio mode 
  # You have manually assigned input.prefix to reach this block
  # Nothing else is necessary...
}


# Set values of these variables to command-line values if provided.
# "option.____" variable names are the ones that should be used throughout
# the code below in case the settings are changed within RStudio
# between runs of the code

option.input.prefix = input.prefix

# Default to the same output prefix as input prefix. 
if (exists("output.prefix")) {
  option.output.prefix = output.prefix
} else {
  option.output.prefix = input.prefix
}

#Default to max.method 2 (for now)
if (exists("max.method")) {
  option.max.method = max.method
} else {
  option.max.method = 2
}


if (exists("two.point.fit")) {
  option.two.point.fit = two.point.fit
} else {
  option.two.point.fit = default.two.point.fit
}


if (option.two.point.fit & exists("growth.rate.method")) {
  fatal_error("Cannot have both --growth-rate-method and --two-point-fit set to TRUE")
}

if (exists("growth.rate.method") && !is.na(growth.rate.method) ) {
  option.growth.rate.method = growth.rate.method
} else {
  option.growth.rate.method = default.growth.rate.method
}
option.growth.rate.method = toupper(option.growth.rate.method)

if (!(option.growth.rate.method %in% c("LOG-LINEAR", "EXPONENTIAL"))) {
  fatal_error("--growth-rate-method must be \"log-linear\" or \"exponential\"")
}

if (option.two.point.fit & exists("fluorescence.rate.method")) {
  fatal_error("Cannot have both --fluorescence-rate-method and --two-point-fit set to TRUE")
}

if (exists("fluorescence.rate.method") && !is.na(fluorescence.rate.method)) {
  option.fluorescence.rate.method = fluorescence.rate.method
} else {
  option.fluorescence.rate.method = default.fluorescence.rate.method
}
option.fluorescence.rate.method = toupper(option.fluorescence.rate.method)

if (!(option.fluorescence.rate.method %in% c("LINEAR", "LOG-LINEAR-PER-OD", "EXPONENTIAL-WITH-PRODUCTION"))) {
  fatal_error("--fluorescence-rate-method must be \"log-linear-per-OD\" or \"exponential-with-production\"")
}
  

if (exists("time.point.delta")) {
  option.time.point.delta = time.point.delta
} else {
  option.time.point.delta = 2
}

if (exists("no.plots")) {
  option.no.plots = no.plots
} else {
  option.no.plots = F
}

if (exists("minimum.OD")) {
  option.minimum.OD = minimum.OD 
} else {
  option.minimum.OD = default.minimum.OD
}

if (exists("maximum.time")) {
  option.maximum.time = maximum.time
} else {
  option.maximum.time = default.maximum.time
}

if (exists("offset.readings.to.average.time.span.end")) {
  options.offset.readings.to.average.time.span.end = offset.readings.to.average.time.span.end
} else {
  options.offset.readings.to.average.time.span.end = default.offset.readings.to.average.time.span.end
}

option.set.background.values = c()
if (exists("set.background.string") && !is.na(set.background.string)) {
  option.set.background.values = as.vector(strsplit(set.background.string, ",")[[1]])
}

##############################################################
#### Advanced options that are currently hard-coded
##############################################################


#Gets set automatically
num.readings = NA

# Fix all initial measurements to be equal across all non-blank samples in a category

if (options.offset.readings.to.average.time.span.end >= 0) {
  offset.readings.to.average = T
} else {
  offset.readings.to.average = F
}

offset.readings.to.average.time.span = c(0,options.offset.readings.to.average.time.span.end)


##############################################################
#### Print the Settings
##############################################################

#Create a dataframe of settings so that we can write them out
# For sdf (settings data frame), the columns are ["name", "value"]

add_setting <- function(sdf, n, v) {
  sdf = suppressWarnings(bind_rows(sdf, data.frame(name=as.character(n), value=as.character(v))))
  return(sdf)
}

print_settings_df <- function(sdf) {
  for (i in 1:nrow(sdf)) {
    cat(paste0(sdf$name[i], " = ", sdf$value[i],"\n"))
  }
}

# If possible (script has not been moved), get the repository/version of the code
# and output this as a setting so that we can know which version was used!
#
# @JEB: This code finds the path of the current script so that we can offset to the 
# Git file containing the currently checked out version. This is a pretty complex
# thing to do and I have not tested this across platforms. It may need more work
# so that we can allow it to gracefully fail w/o erroring out of the script.

this_file_path <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    #cat("Rscript\n")
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else if (!is.null(sys.frames()[[1]]$ofile)) {
    # 'source'd via R console
    #cat("Sourced\n")
    return(normalizePath(sys.frames()[[1]]$ofile))
  } else {
    # Rstudio
    #cat("RStudio\n")
    return(rstudioapi::getActiveDocumentContext()$path)
  }
}

#Returns data frame with 'version' and 'repository' fields
git_repo_version <- function() {
  
  grv = c()
  git.fetch.head.path = file.path(dirname(this_file_path()), "..", "..", ".git", "FETCH_HEAD")
  if(file.exists(git.fetch.head.path)) {
    mystring <- read_file(git.fetch.head.path)
    mystring = str_remove_all(mystring, "[\r\n]")
    mystringlist = strsplit(mystring, "\t")
    grv$repository = mystringlist[[1]][3]
  } else {
    grv$repository = "unknown"
  }
  
  git.log.head.path = file.path(dirname(this_file_path()), "..", "..", ".git", "logs", "HEAD")
  if(file.exists(git.log.head.path)) {
    mystring <- read_file(git.log.head.path)
    linelist = strsplit(mystring, "\n")
    lastline = linelist[[1]][length(linelist[[1]])]
    lastline = str_remove_all(lastline, "[\r\n]")
    lastlinelist = strsplit(lastline, " ")
    grv$version = lastlinelist[[1]][2]
  } else {
    grv$version = "unknown"
  }
  
  return(grv)
}
git.repo.version.info = git_repo_version()

sdf = data.frame()
sdf = add_setting(sdf, "Run At", Sys.time())
sdf = add_setting(sdf, "Repository", git.repo.version.info$repository)
sdf = add_setting(sdf, "Version", git.repo.version.info$version)
sdf = add_setting(sdf, "Input Prefix", option.input.prefix)
sdf = add_setting(sdf, "Output Prefix", option.output.prefix)

# Translate the method to something meaningful
option.max.method.string = "unknown"

if (option.max.method == 1) {
  option.max.method.string = "1 (Find time point with maximum growth rate. Report all rates at this time point.) [one time = max growth rate]"
} else if (option.max.method == 2) {
  option.max.method.string = "2 (Find time point with maximum rate and report that for each curve separately.) [multiple times = max of each curve]"  
}
sdf = add_setting(sdf, "Maximum Picking Method", option.max.method.string)

sdf = add_setting(sdf, "Offset Readings to Initial Average", offset.readings.to.average)
sdf = add_setting(sdf, "Offset Readings to Initial Average Time Span", 
                  paste(offset.readings.to.average.time.span, collapse="-"))
sdf = add_setting(sdf, "Minimum OD", option.minimum.OD)
sdf = add_setting(sdf, "Maximum Time", option.maximum.time)

if (option.two.point.fit) {
  option.fit.method.string = "using two points (first and last in time window)"
  sdf = add_setting(sdf, "Fitting method", option.fit.method.string)
} else {
  option.fit.method.string = "using all points in time window"
  sdf = add_setting(sdf, "Fitting method", option.fit.method.string)
  
  sdf = add_setting(sdf, "  Growth rate method", option.growth.rate.method)
  sdf = add_setting(sdf, "  Fluorescence rate method", option.fluorescence.rate.method)
}      


sdf = add_setting(sdf, "Time Point Span", 2 * option.time.point.delta + 1)
sdf = add_setting(sdf, "Time Point Delta", option.time.point.delta)



option.use.blanks.for.background = T
if (length(option.set.background.values) > 0) {
  option.use.blanks.for.background = F
}
sdf = add_setting(sdf, "Use Blank Wells", option.use.blanks.for.background)
sdf = add_setting(sdf, "Fixed Blank Values", paste0(option.set.background.values, "", collapse=","))

sdf = add_setting(sdf, "No plots", option.no.plots)

cat("=== SETTINGS ==\n")
print_settings_df(sdf)

##############################################################
#### Check that measurement and metadata files exist
##############################################################

if (   !file.exists(paste0(option.input.prefix, ".measurements.csv")) 
       && !file.exists(paste0(option.input.prefix, ".measurements.tsv")) ) {
  fatal_error(paste0("Could not find a valid measurements file. Tried:\n", paste0(option.input.prefix, ".measurements.csv"), "\n", paste0(option.input.prefix, ".measurements.tsv")))
}

if (   !file.exists(paste0(option.input.prefix, ".metadata.csv")) 
       && !file.exists(paste0(option.input.prefix, ".metadata.tsv")) ) {
  fatal_error(paste0("Could not find a valid measurements file. Tried:\n", paste0(option.input.prefix, ".measurements.csv"), "\n", paste0(option.input.prefix, ".measurements.tsv")))
}

##############################################################
#### Input files existed, now create the plot directory and write  settings
##############################################################

write.csv(sdf, paste0(option.output.prefix, ".settings.csv"), row.names=F)

if (!option.no.plots) {
  plot.directory = paste0(option.output.prefix, "-plots")
  dir.create(file.path(plot.directory), showWarnings = FALSE)
}

##############################################################
#### Load and tidy measurement file
##############################################################

#SPL Note: had to run "dos2unix" to clean up the file from our platereader before it could be read properly. Not sure how to fix this on the windows end.
#TODO: Check for incorrect line endings and warn user

all_data = data.frame()

# We know one of these exists b/c we checked above!
if (file.exists(paste0(option.input.prefix, ".measurements.csv"))) {
  all_data <- read_csv(paste0(option.input.prefix, ".measurements.csv"), col_names=F, comment = "#" )
} else if (file.exists(paste0(option.input.prefix, ".measurements.tsv"))) {
  all_data <- read_tsv(paste0(option.input.prefix, ".measurements.tsv"), col_names=F, comment = "#" )
}

# If there is no header row of wells, assume standard order by row from top of plate.
# Otherwise, take the first row as the well labels.
if ((all_data$X1[1] == "0s" || as.character(all_data$X1[1]) == "0")) {
  names(all_data) = c("Time",
                      "A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11", "A12",
                      "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B11", "B12",
                      "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12",
                      "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "D11", "D12",
                      "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11", "E12",
                      "F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", "F11", "F12",
                      "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9", "G10", "G11", "G12",
                      "H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11", "H12"
  )
  
  
  
} else {
  names(all_data) = all_data[1,]
  all_data = all_data[-1,]
}

#Remove row thats are NA in time (this is the trialing trash data). 
#Remove "s" and convert to integer, remove row that dont convert properly and leave NA

all_data$Time <- str_replace(all_data$Time, "s", "")
all_data$Time <- as.integer(all_data$Time)
all_data <- all_data %>% filter(!is.na(Time))

#Separate data into data frames to make tidy

#Find the zero times? (they separate the different types of readings)
zeros <- which(all_data$Time == "0")
zeros = c(zeros, nrow(all_data)+1)

stopifnot(length(zeros) <= 4) 
num.readings = length(zeros) - 1

#generic names for the readings
reading.names = c("OD", "GFP", "other")

tidy_all = data.frame()

## Need to keep all times in lock step despite some rounding
## so we use the first chunk times
definitive_times = c()

for(i in 1:(length(zeros)-1)) {
  this.data.chunk = all_data[zeros[i]:(zeros[i+1]-1),]
  
  if (i==1) {
    definitive_times = floor(this.data.chunk$Time/60)
  }
  
  ## If the file was truncated - readings stopped in the middle of a row by the user,
  ## Then later chunks may have fewer values. In this case, we need to shorten the earlier
  ## ones and the number of definitive times.
  if (nrow(this.data.chunk) < length(definitive_times)) {
    definitive_times = definitive_times[1:nrow(this.data.chunk)]
    tidy_all = tidy_all %>% filter(time.min < definitive_times[length(definitive_times)])
  }
  
  tidy.data.chunk <- gather(this.data.chunk, key = "well", value = "value", -Time)
  tidy.data.chunk$reading <- reading.names[i]
  
  tidy.data.chunk$time.min = rep(definitive_times, nrow(tidy.data.chunk)/length(definitive_times))
  tidy.data.chunk = tidy.data.chunk %>% select (-Time)
  
  tidy_all = tidy_all %>% bind_rows(tidy.data.chunk)
}

# Remove overflow readings
tidy_all$value=as.numeric(tidy_all$value)


#Write the tidy data
write_csv( tidy_all, paste0(option.output.prefix, ".tidy.measurements.csv"))

##############################################################
#### Load and tidy metadata file
##############################################################

metadata <- data.frame()

# We know one of these exists b/c we checked above!
if (file.exists(paste0(option.input.prefix, ".metadata.csv"))) {
  metadata <- read_csv(paste0(option.input.prefix, ".metadata.csv"), col_names=T, comment = "#" )
} else if (file.exists(paste0(option.input.prefix, ".metadata.tsv"))) {
  metadata <- read_tsv(paste0(option.input.prefix, ".metadata.tsv"), col_names=T, comment = "#" )
}

## Check for required columns
if(! ("well" %in% colnames(metadata)) ) {
  stop("Metadata does not have 'well' column (case-sensitive)", call.=FALSE)
}
if(! ("strain" %in% colnames(metadata)) ) {
  stop("Metadata does not have 'strain' column (case-sensitive)", call.=FALSE)
}

#Add optional columns if missing
if(! ("include" %in% colnames(metadata)) ) {
  metadata$include = T
}
if(! ("isolate" %in% colnames(metadata)) ) {
  metadata$isolate = NA
}
if(! ("description" %in% colnames(metadata)) ) {
  metadata$description = ""
}

#convert boolean fields
metadata$include = as.logical(metadata$include)


## Determine what wells we are ignoring and remove from metadata
ignore.wells = metadata$well[metadata$include==F]
metadata = metadata %>% filter(include==T)

# Fix the names of strains that are synonyms for blank (case-insensitive)
metadata$strain.upper = toupper(metadata$strain)
metadata$strain[grepl("^(B|BLANK)$", metadata$strain.upper, perl=T)] = "blank"
metadata = metadata %>% select(-strain.upper)

# Blanks can't have isolate numbers
metadata$isolate[metadata$strain=="blank"] = NA

metadata$strain.isolate = paste0(metadata$strain, "__", metadata$isolate)

#Fix wells that have blank back to blank
metadata$strain.isolate = sub("blank__NA", "blank", metadata$strain.isolate)

write_csv( metadata, paste0(option.output.prefix, ".tidy.metadata.csv"))

##############################################################
#### Combine tidy measurement and metadata files
##############################################################

#Now,  load the data file and wrangle it
X = tidy_all
X = X %>% filter(time.min < option.maximum.time)
X = X %>% filter(!(well %in% ignore.wells))
X$value=as.numeric(X$value)
Y = X %>% spread(reading, value)

#join sample data to master tidy data
Y <- full_join(metadata, Y, by = "well")
Y$well = factor(Y$well)
Y = Y %>% filter(time.min != 0)


# Analysis starts here ------->
##############################################################
#### Background Correction
#### * Subtracts the average of ALL blank wells at a given time from readings
##############################################################

if (!option.use.blanks.for.background) {
  
  if (length(option.set.background.values) !=  num.readings) {
    fatal_error("WARNING: The number of background values provided does not match the number of measurements.\n")
  }
  
  Y$OD = as.numeric(Y$OD) - as.numeric(option.set.background.values[1])
  Y$GFP = as.numeric(Y$GFP) - as.numeric(option.set.background.values[2])
  
  if (num.readings == 3) {
    Y$other = as.numeric(Y$other) - as.numeric(option.set.background.values[3])
  }
  
  Z = Y
  
} else {
  
  if ( nrow(metadata %>% filter(strain.isolate == "blank")) == 0) {
    fatal_error("No blank wells are included. Check your metadata file or add the --background option.\n")
  }
  
  Y$GFP = as.numeric(Y$GFP)
  Y$OD = as.numeric(Y$OD)
  
  BG = Y %>% filter(strain=="blank")
  
  BG$well = droplevels(BG$well)
  
  if (!option.no.plots) {
    plot1 = ggplot(BG, aes(time.min, OD, color=well)) + geom_point() 
    plot2 = ggplot(BG, aes(time.min, GFP, color=well)) + geom_point() 
    
    
    if (num.readings == 2) {
      p = grid.arrange(plot1, plot2, nrow=2)
    } else if (num.readings == 3) {
      plot3 = ggplot(BG, aes(time.min, other, color=well)) + geom_point() 
      p = grid.arrange(plot1, plot2, plot3, nrow=3)
    }
    
    ggsave(file.path(plot.directory, "background.pdf"), plot=p)
  }
  
  Z = Y
  if (num.readings== 3) {
    
    BG.means = BG %>% group_by(time.min) %>% summarize(other.bg = mean(other), GFP.bg = mean(GFP), OD.bg = mean(OD))
    
    Z = Z %>% left_join(BG.means, by="time.min")
    
    Z$OD = Z$OD - Z$OD.bg
    Z$GFP = Z$GFP - Z$GFP.bg
    Z$other = Z$other - Z$other.bg
    
  } else if (num.readings == 2) {
    
    BG.means = BG %>% group_by(time.min) %>% summarize(GFP.bg = mean(GFP), OD.bg = mean(OD))
    
    Z = Z %>% left_join(BG.means, by="time.min")
    
    Z$OD = Z$OD - Z$OD.bg
    Z$GFP = Z$GFP - Z$GFP.bg
    
  }
  
  ## There is an error in group_by with some versions of the tidyverse and it does
  ## not average per time point and only gives an overall average
  
  if (nrow(BG.means)==1) {
    fatal_error("Averaging the background per time point faied. Please update your tidyverse packages.\n   install.packages(\"tidyverse\")")
  }
}

#Ditch the blanks here
Z = Z %>% filter(strain != "blank")

#For debugging
#write_tsv(Z, "temp.tsv")

##############################################################
#### Average Per-Sample Offset Correction
#### * Equalizes the beginning of each curve for a given strain
####   This is a kludge to correct for well-to-well background differences
####   TODO:  figure out a way to experimentally measure these differences!
##############################################################
if (offset.readings.to.average) {
  
  cat("Offsetting to average measured values\n")
  
  #calculate the average over all times for a sample
  if (num.readings == 2) {
    
    ZZ = Z %>% filter((time.min >= offset.readings.to.average.time.span[1]) & (time.min <= offset.readings.to.average.time.span[2]))
    
    strain.means = ZZ %>% group_by(strain.isolate) %>% summarize(mean.strain.OD = mean(OD), mean.strain.GFP = mean(GFP))
    well.means = ZZ %>% group_by(well, strain.isolate) %>% summarize(mean.well.OD = mean(OD), mean.well.GFP = mean(GFP))
    
    well.means = well.means %>% left_join(strain.means, by="strain.isolate")
    
    well.means$well.OD.offset = well.means$mean.strain.OD - well.means$mean.well.OD
    well.means$well.GFP.offset = well.means$mean.strain.GFP - well.means$mean.well.GFP
    
    well.means = well.means %>% select(well, well.OD.offset, well.GFP.offset)
    
    Z = Z %>% left_join(well.means, by="well")
    
    Z$OD = Z$OD + Z$well.OD.offset
    Z$GFP = Z$GFP + Z$well.GFP.offset
    
    
  } else if (num.readings == 3) {
    
    ZZ = Z %>% filter((time.min >= offset.readings.to.average.time.span[1]) & (time.min <= offset.readings.to.average.time.span[2]))
    
    strain.means = ZZ %>% group_by(strain.isolate) %>% summarize(mean.strain.OD = mean(OD), mean.strain.GFP = mean(GFP), mean.strain.other = mean(other))
    well.means = ZZ %>% group_by(well, strain.isolate) %>% summarize(mean.well.OD = mean(OD), mean.well.GFP = mean(GFP), mean.well.other = mean(other))
    
    well.means = well.means %>% left_join(strain.means, by="strain.isolate")
    
    well.means$well.OD.offset = well.means$mean.strain.OD - well.means$mean.well.OD
    well.means$well.GFP.offset = well.means$mean.strain.GFP - well.means$mean.well.GFP
    well.means$well.other.offset = well.means$mean.strain.other - well.means$mean.well.other
    
    well.means = well.means %>% select(well, well.OD.offset, well.GFP.offset, well.other.offset)
    
    Z = Z %>% left_join(well.means, by="well")
    
    Z$OD = Z$OD + Z$well.OD.offset
    Z$GFP = Z$GFP + Z$well.GFP.offset
    Z$other = Z$other + Z$well.other.offset
    
  }
}


##############################################################
#### Calculate and graph rates per strain
##############################################################

## Some notes on the fitting equations and procedures

## C = cells (assumed ~ OD)
## F = fluorescence
## t = time
## r = specific growth rate for cells (units of 1/t)
## P = production rate of fluorescence (units of AU/t)

## Growth equation
## dC/dt = r * C(t)
## log[ C(t) ] = r * t + k
## => C(t) = k * exp(r * t)
## C(0) = k * exp (0) = k
## => C(t) = C0 * exp(r * t)
## => log[ C(t) ] = r * t + log[ C(0) ]
## slope = r; y-intercept = log[ C(0) ]
##
## Fits are to linearized data - this probably has little effect versus fitting
## directly to exponential equations given the short time windows

## Fluorescence equation
## dF/dt = P * C(t) = P * C(0) * exp(r * t)
## => F(t) = P * C(0) * 1/r * exp(r * t) + k
## => F(t) = P * 1/r * C(t) + k

## F(0) = P * C(0) * 1/r * exp (0) + k
## k = F(0) - P * C(0) * 1/r 
## => F(t) = P * C(0) * 1/r * exp(rt) + F(0) - P * C(0) * 1/r


Z = Z %>% filter (OD >= option.minimum.OD)
final.table.summary = data.frame()
final.table.all.wells = data.frame()
for (strain.of.interest in unique(Z$strain.isolate) )
{
  cat("STRAIN:", strain.of.interest, "\n") 
  strain.data = Z %>% filter(strain.isolate==strain.of.interest)
  head(strain.data)
  
  ## graph
  
  if (!option.no.plots) {
    plot1 = ggplot(strain.data, aes(time.min, OD, color=well)) + geom_point() 
    plot2 = ggplot(strain.data, aes(time.min, GFP, color=well)) + geom_point() 
    
    if (num.readings == 2) {
      p = grid.arrange(plot1, plot2, nrow=2)
    } else if (num.readings == 3) {
      plot3 = ggplot(strain.data, aes(time.min, other, color=well)) + geom_point() 
      p = grid.arrange(plot1, plot2, plot3, nrow=3)
    }
    
    ggsave(file.path(plot.directory, paste0(strain.of.interest, ".pdf")), plot=p)
  }
  
  strain.data$logOD = log(strain.data$OD)
  #ggplot(strain.data, aes(time.min, logOD, color=well)) + scale_x_continuous(limits = c(0, 400))+ geom_point() 
  
  
  #Calculate the growth rate over a given time period and graph
  growth.rate.data = strain.data
  growth.rate.data$well = droplevels(growth.rate.data$well)
  
  
  if (option.two.point.fit) {
    
    growth.rate.data = growth.rate.data %>% 
      group_by(well) %>% 
      mutate(t1 = lag(time.min, option.time.point.delta, order_by=well)) %>% 
      mutate( t2 = lead(time.min, option.time.point.delta, order_by=well)) %>% 
      mutate( logOD.t1 = lag(logOD, option.time.point.delta, order_by=well)) %>% 
      mutate( logOD.t2 = lead(logOD, option.time.point.delta, order_by=well))
    
    growth.rate.data = growth.rate.data %>% filter(!is.na(t1) & !is.na(t2))
    
    growth.rate.data = growth.rate.data %>% mutate(time.min = (t1+t2)/2)
    growth.rate.data = growth.rate.data %>% mutate(growth.rate = 60 * (logOD.t2 - logOD.t1) /  (t2 - t1))
    
  } else {
    
    new.growth.rate.data = data.frame()
    
    for (this.well in levels(growth.rate.data$well)) {
      #cat(this.well, "\n")
      replicate.growth.rate.data = growth.rate.data %>% filter(well == this.well)
      replicate.growth.rate.data$growth.rate = NA
      replicate.growth.rate.data$initial.OD = NA
      replicate.growth.rate.data$growth.rate.adj.r.squared = NA
      replicate.growth.rate.data$fit.OD = NA
      replicate.growth.rate.data$growth.fit.ss.residuals = NA
      
      
      if (nrow(replicate.growth.rate.data) >= 2 * option.time.point.delta + 1) {
        
        for(i in (1+option.time.point.delta):(nrow(replicate.growth.rate.data)-option.time.point.delta)) {
          
          if (option.growth.rate.method == "LOG-LINEAR") {
            
            times = replicate.growth.rate.data$time.min[(i-option.time.point.delta):(i+option.time.point.delta)]
            logODs = replicate.growth.rate.data$logOD[(i-option.time.point.delta):(i+option.time.point.delta)]
            
            fit = lm(logOD~time, data.frame(time=times, logOD=logODs))
            replicate.growth.rate.data$growth.rate[i] = 60 * coef(fit)[["time"]]
            replicate.growth.rate.data$initial.OD[i] = exp(coef(fit)[["(Intercept)"]])
            replicate.growth.rate.data$growth.rate.adj.r.squared[i] = summary(fit)$adj.r.squared
            replicate.growth.rate.data$fit.OD[i] = exp(fit$fitted.values[option.time.point.delta+1])
            
            replicate.growth.rate.data$growth.fit.ss.residuals[i] = sum(residuals(nls.fit)^2)
            
          } else if (option.growth.rate.method == "EXPONENTIAL") {
            
            nls.df = data.frame(
              t = replicate.growth.rate.data$time.min[(i-option.time.point.delta):(i+option.time.point.delta)] - replicate.growth.rate.data$time.min[i-option.time.point.delta],
              OD = replicate.growth.rate.data$OD[(i-option.time.point.delta):(i+option.time.point.delta)]
            )
            
            ## F(t) = P * 1/r * C(t) + k
            r_guess = (log(nls.df$OD[nrow(nls.df)]) - log(nls.df$OD[1]))/(nls.df$t[nrow(nls.df)]-nls.df$t[1])
            C0_guess = replicate.growth.rate.data$OD[i-option.time.point.delta]
            
            try(
              nls.fit <- nls(OD ~ C0 * exp(r * t),
                             data = nls.df, 
                             start = list(
                               r=r_guess,
                               C0=C0_guess
                             ), 
                             control=nls.control(maxiter=1000, tol = 1e-5, minFactor = 1/100000, warnOnly=TRUE),
                             trace = FALSE)
            )
            
            this.converged = nls.fit$convInfo$isConv
            
            if (this.converged) {
              replicate.growth.rate.data$growth.rate[i] = coef(nls.fit)[["r"]] * 60
              replicate.growth.rate.data$initial.OD[i] = coef(nls.fit)[["C0"]]
              replicate.growth.rate.data$fit.OD[i] = predict(nls.fit)[option.time.point.delta+1]
              replicate.growth.rate.data$growth.fit.ss.residuals[i] = sum(residuals(nls.fit)^2)
            }
          }
          
        }
      }
      
      new.growth.rate.data = rbind(new.growth.rate.data, replicate.growth.rate.data)
    }
    
    # Join values back to main data frame for graphing and use
    growth.rate.data = new.growth.rate.data
  }
  
  if (!option.no.plots) {
    plot1 = ggplot(growth.rate.data , aes(time.min, growth.rate, color=well)) +  geom_point() 
  }
  
  #Calculate the growth rate over a given time period and graph
  fluorescence.data = growth.rate.data
  
  if (option.two.point.fit) {
    
    fluorescence.data = fluorescence.data %>% 
      group_by(well) %>% 
      mutate(t1 = lag(time.min, option.time.point.delta, order_by=well)) %>% 
      mutate( t2 = lead(time.min, option.time.point.delta, order_by=well)) %>% 
      mutate( GFP.t1 = lag(GFP, option.time.point.delta, order_by=well)) %>% 
      mutate( GFP.t2 = lead(GFP, option.time.point.delta, order_by=well)) 
    
    if (num.readings == 3) {
      fluorescence.data = fluorescence.data %>% 
        mutate( other.t1 = lag(other, option.time.point.delta, order_by=well)) %>% 
        mutate( other.t2 = lead(other, option.time.point.delta, order_by=well))
    }
    
    fluorescence.data = fluorescence.data %>% filter(!is.na(t1) & !is.na(t2))
    fluorescence.data = fluorescence.data %>% mutate(time.min = (t1+t2)/2)
    fluorescence.data = fluorescence.data %>% mutate(GFP.rate = (GFP.t2 - GFP.t1) / (t2 - t1) / OD * 60)
    
    if (num.readings == 3) {
      fluorescence.data = fluorescence.data %>% mutate(other.rate = (other.t2 - other.t1) / (t2 - t1) / OD * 60)
    }
    
  } else {
    new.fluorescence.data = data.frame()
    
    for (this.well in levels(fluorescence.data$well)) {
      #cat(this.well, "\n")
      replicate.fluorescence.data = fluorescence.data %>% filter(well == this.well)
      replicate.fluorescence.data$GFP.rate = NA
      replicate.fluorescence.data$GFP.fit.adj.r.squared = NA
      replicate.fluorescence.data$GFP.fit.ss.residuals = NA
      
      replicate.fluorescence.data$other.rate = NA
      replicate.fluorescence.data$other.fit.adj.r.squared = NA
      replicate.fluorescence.data$other.fit.ss.residuals = NA
      
      if (nrow(replicate.fluorescence.data) >= 2 * option.time.point.delta + 1) {
        
        for(i in (1+option.time.point.delta):(nrow(replicate.fluorescence.data)-option.time.point.delta)) {
          
          #cat(i, "\n")
          #create a matrix
          
          if (option.fluorescence.rate.method == "LINEAR") {
            
            times = replicate.fluorescence.data$time.min[(i-option.time.point.delta):(i+option.time.point.delta)]
            fluorescences = replicate.fluorescence.data$GFP[(i-option.time.point.delta):(i+option.time.point.delta)]
            #don't do this with NAs => caused by overflow of fluorescent challenges
            if (!is.na(mean(fluorescences)) && !is.infinite(mean(fluorescences))) {
              fit = lm(fluorescence~time, data.frame(time=times, fluorescence=fluorescences))
              replicate.fluorescence.data$GFP.rate[i] = coef(fit)[["time"]] / replicate.fluorescence.data$fit.OD[i] * 60
              replicate.fluorescence.data$GFP.fit.adj.r.squared[i] = summary(fit)$adj.r.squared
            }
          }
          
          else  if (option.fluorescence.rate.method == "LOG-LINEAR-PER-OD") {
            
            times = replicate.fluorescence.data$time.min[(i-option.time.point.delta):(i+option.time.point.delta)]
            log_fluorescence = log(replicate.fluorescence.data$GFP[(i-option.time.point.delta):(i+option.time.point.delta)])
            
            ## Remove negative/missing values and fit if half are remaining...
            times = times[!is.na(log_fluorescence) & !is.nan(log_fluorescence) & !is.infinite(log_fluorescence)]
            log_fluorescence = log_fluorescence[!is.na(log_fluorescence) & !is.nan(log_fluorescence) & !is.infinite(log_fluorescence)]
            
            if (length(log_fluorescence) > option.time.point.delta) {
              fit = lm(fluorescence~time, data.frame(time=times, fluorescence=log_fluorescence))
              replicate.fluorescence.data$GFP.rate[i] = coef(fit)[["time"]] * exp(fit$fitted.values[option.time.point.delta+1]) / replicate.fluorescence.data$fit.OD[i] * 60
              
              replicate.fluorescence.data$GFP.fit.adj.r.squared[i] = summary(fit)$adj.r.squared
              replicate.fluorescence.data$GFP.fit.ss.residuals[i] = sum((exp(fit$fitted.value) - exp(log_fluorescence))^2)
              
            }
          }
          
          else if (option.fluorescence.rate.method == "EXPONENTIAL-WITH-PRODUCTION") {
            
            times = replicate.fluorescence.data$time.min[(i-option.time.point.delta):(i+option.time.point.delta)] - replicate.fluorescence.data$time.min[i-option.time.point.delta]
            
            fluorescences = replicate.fluorescence.data$GFP[(i-option.time.point.delta):(i+option.time.point.delta)]
            
            ## Remove negative/missing values and fit if half are remaining...
            times = times[!is.na(fluorescences) & !is.nan(fluorescences) & !is.infinite(fluorescences)]
            fluorescences = fluorescences[!is.na(fluorescences) & !is.nan(fluorescences) & !is.infinite(fluorescences)]
            
            if (
              (length(fluorescences) > option.time.point.delta) 
              && (!is.na(replicate.fluorescence.data$growth.rate[i]))
              && (!is.na(replicate.fluorescence.data$initial.OD[i]))
            ) {
              
              # specific growth rate (units of 1/min)
              r = replicate.fluorescence.data$growth.rate[i] / 60
              C0 = replicate.fluorescence.data$initial.OD[i]
              
              nls.df = data.frame(t=times, fl = fluorescences)
              nls.df$C = C0 * exp(r * nls.df$t)
              
              ## F(t) = P * 1/r * C(t) + k
              
              # Production rate (units of 1/min)
              P_guess = (nls.df$fl[nrow(nls.df)] - nls.df$fl[1]) / ((nls.df$C[nrow(nls.df)] + nls.df$C[1])/2) / 60 
              k_guess = 0
              #k_guess = nls.df$fl[1] - P_guess * 1/r * C0
              
              #Try a first fit
              try(
                nls.fit <- nls(
                  fl ~ P * (1/r) * C + k,
                  data = nls.df, 
                  start = list(P=P_guess, k=0), 
                  control=nls.control(maxiter=1000, tol = 1e-5, minFactor = 1/100000, warnOnly=TRUE),
                  trace = FALSE)
              )
              
              this.converged = nls.fit$convInfo$isConv
              
              # Try a second fit
              if (!this.converged) {
                P_guess = 2*(nls.df$fl[nrow(nls.df)] - nls.df$fl[1]) / ((nls.df$C[nrow(nls.df)] + nls.df$C[1])/2) / 60 
                k_guess = -100
                
                try(
                  nls.fit <- nls(
                    fl ~ P * (1/r) * C + k,
                    data = nls.df, 
                    start = list(P=P_guess, k=k_guess), 
                    control=nls.control(maxiter=1000, tol = 1e-5, minFactor = 1/100000, warnOnly=TRUE),
                    trace = FALSE)
                )
                this.converged = nls.fit$convInfo$isConv
              }
              
              if (this.converged) {
                replicate.fluorescence.data$GFP.rate[i] = coef(nls.fit)[["P"]] * 60
                replicate.fluorescence.data$GFP.fit.ss.residuals[i] = sum(resid(nls.fit)^2)
              }
              
            }
            
          }
          
          if (num.readings == 3) {
            
            if (option.fluorescence.rate.method == "LINEAR") {
              
              times = replicate.fluorescence.data$time.min[(i-option.time.point.delta):(i+option.time.point.delta)]
              fluorescences = replicate.fluorescence.data$other[(i-option.time.point.delta):(i+option.time.point.delta)]
              #don't do this with NAs => caused by overflow of fluorescent challenges
              if (!is.na(mean(fluorescences)) && !is.infinite(mean(fluorescences))) {
                fit = lm(fluorescence~time, data.frame(time=times, fluorescence=fluorescences))
                replicate.fluorescence.data$other.rate[i] = coef(fit)[["time"]] / replicate.fluorescence.data$fit.OD[i] * 60
                replicate.fluorescence.data$other.fit.adj.r.squared[i] = summary(fit)$adj.r.squared
              }
            }
            
            else  if (option.fluorescence.rate.method == "LOG-LINEAR-PER-OD") {
              
              times = replicate.fluorescence.data$time.min[(i-option.time.point.delta):(i+option.time.point.delta)]
              log_fluorescence = log(replicate.fluorescence.data$other[(i-option.time.point.delta):(i+option.time.point.delta)])
              
              ## Remove negative/missing values and fit if half are remaining...
              times = times[!is.na(log_fluorescence) & !is.nan(log_fluorescence) & !is.infinite(log_fluorescence)]
              log_fluorescence = log_fluorescence[!is.na(log_fluorescence) & !is.nan(log_fluorescence) & !is.infinite(log_fluorescence)]
              
              if (length(log_fluorescence) > option.time.point.delta) {
                fit = lm(fluorescence~time, data.frame(time=times, fluorescence=log_fluorescence))
                replicate.fluorescence.data$other.rate[i] = coef(fit)[["time"]] * exp(fit$fitted.values[option.time.point.delta+1]) / replicate.fluorescence.data$fit.OD[i] * 60
                
                replicate.fluorescence.data$other.fit.adj.r.squared[i] = summary(fit)$adj.r.squared
                replicate.fluorescence.data$other.fit.ss.residuals[i] = sum((exp(fit$fitted.value) - exp(log_fluorescence))^2)
                
              }
            }
            
            else if (option.fluorescence.rate.method == "EXPONENTIAL-WITH-PRODUCTION") {
              
              times = replicate.fluorescence.data$time.min[(i-option.time.point.delta):(i+option.time.point.delta)] - replicate.fluorescence.data$time.min[i-option.time.point.delta]
              
              fluorescences = replicate.fluorescence.data$other[(i-option.time.point.delta):(i+option.time.point.delta)]
              
              ## Remove negative/missing values and fit if half are remaining...
              times = times[!is.na(fluorescences) & !is.nan(fluorescences) & !is.infinite(fluorescences)]
              fluorescences = fluorescences[!is.na(fluorescences) & !is.nan(fluorescences) & !is.infinite(fluorescences)]
              
              if (
                (length(fluorescences) > option.time.point.delta) 
                && (!is.na(replicate.fluorescence.data$growth.rate[i]))
                && (!is.na(replicate.fluorescence.data$initial.OD[i]))
              ) {
                
                # specific growth rate (units of 1/min)
                r = replicate.fluorescence.data$growth.rate[i] / 60
                C0 = replicate.fluorescence.data$initial.OD[i]
                
                nls.df = data.frame(t=times, fl = fluorescences)
                nls.df$C = C0 * exp(r * nls.df$t)
                
                ## F(t) = P * 1/r * C(t) + k
                
                # Production rate (units of 1/min)
                P_guess = (nls.df$fl[nrow(nls.df)] - nls.df$fl[1]) / ((nls.df$C[nrow(nls.df)] + nls.df$C[1])/2) / 60 
                k_guess = 0
                #k_guess = nls.df$fl[1] - P_guess * 1/r * C0
                
                #Try a first fit
                try(
                  nls.fit <- nls(
                    fl ~ P * (1/r) * C + k,
                    data = nls.df, 
                    start = list(P=P_guess, k=0), 
                    control=nls.control(maxiter=1000, tol = 1e-5, minFactor = 1/100000, warnOnly=TRUE),
                    trace = FALSE)
                )
                
                this.converged = nls.fit$convInfo$isConv
                
                # Try a second fit
                if (!this.converged) {
                  P_guess = 2*(nls.df$fl[nrow(nls.df)] - nls.df$fl[1]) / ((nls.df$C[nrow(nls.df)] + nls.df$C[1])/2) / 60 
                  k_guess = -100
                  
                  try(
                    nls.fit <- nls(
                      fl ~ P * (1/r) * C + k,
                      data = nls.df, 
                      start = list(P=P_guess, k=k_guess), 
                      control=nls.control(maxiter=1000, tol = 1e-5, minFactor = 1/100000, warnOnly=TRUE),
                      trace = FALSE)
                  )
                  this.converged = nls.fit$convInfo$isConv
                }
                
                if (this.converged) {
                  replicate.fluorescence.data$other.rate[i] = coef(nls.fit)[["P"]] * 60
                  replicate.fluorescence.data$other.fit.ss.residuals[i] = sum(resid(nls.fit)^2)
                }
                
              }
            }
          }
        }
      }
      
      new.fluorescence.data = rbind(new.fluorescence.data, replicate.fluorescence.data)
    }
    
    # Join values back to main data frame for graphing and use
    fluorescence.data = new.fluorescence.data
    
    if (!option.no.plots) {
      write.csv(fluorescence.data, file.path(plot.directory, paste0(strain.of.interest, ".rates.csv")))
    }
    
  }
  
  if (!option.no.plots) {
    
    plot2 = ggplot(fluorescence.data , aes(time.min, GFP.rate, color=well)) +  geom_point() 
    if (num.readings == 2) {
      p = grid.arrange(plot1, plot2, nrow=2)
    } else if (num.readings == 3) {
      plot3 = ggplot(fluorescence.data , aes(time.min, other.rate, color=well)) +  geom_point() 
      p = grid.arrange(plot1, plot2, plot3, nrow=3)
    }
    
    ggsave(file.path(plot.directory, paste0(strain.of.interest, ".rates.pdf")), plot=p)
  }

  growth.rate.data$well = droplevels(growth.rate.data$well)
  
  strain.max.values = data.frame()
  
  
  cat("STRAIN__ISOLATE:", strain.of.interest, "\n\n")
  
  for (this.well in levels(growth.rate.data$well)) {
    
    #cat("Well:", this.well, "\n")
    
    replicate.growth.rate.data = growth.rate.data %>% filter(well == this.well)
    
    #identify the best growth rate row/time
    max.growth.rate.row = replicate.growth.rate.data[which.max(replicate.growth.rate.data$growth.rate),]
    
    #Case for failing to find a max
    if (nrow(max.growth.rate.row) == 0) {
      max.growth.rate.row  = replicate.growth.rate.data[1,]
      max.growth.rate.row$time.min[1] = NA
    }
    
    cat("Time of maximum growth rate:", max.growth.rate.row$time.min[1], "min\n")
    cat("Max growth rate:", max.growth.rate.row$growth.rate, "per hour\n")
    
    if (option.max.method==1) {
      max.GFP.fluorescence.data.row = fluorescence.data %>% filter(well == this.well) %>% filter(time.min == max.growth.rate.row$time.min[1])
    } else if (option.max.method==2) {
      
      replicate.fluorescence.data = fluorescence.data %>% filter(well == this.well)
      max.GFP.fluorescence.data.row = replicate.fluorescence.data[which.max(replicate.fluorescence.data$GFP.rate),]
    }
    
    #Case for failing to find a max
    if (nrow(max.GFP.fluorescence.data.row) == 0) {
      max.GFP.fluorescence.data.row  = fluorescence.data[1,]
      max.GFP.fluorescence.data.row$time.min[1] = NA
    }
    
    cat("Time of GFP production rate:", max.GFP.fluorescence.data.row$time.min[1], "min\n")
    cat("Max GFP production rate:", max.GFP.fluorescence.data.row$GFP.rate, "per hour\n")
    
    if (num.readings == 3) {
      
      if (option.max.method==1) {
        max.other.fluorescence.data.row = fluorescence.data %>% filter(well == this.well) %>% filter(time.min == max.growth.rate.row$time.min[1]) 
      } else if (option.max.method==2) {
        max.other.fluorescence.data.row = replicate.fluorescence.data[which.max(replicate.fluorescence.data$other.rate),]
      }
      
      #Case for failing to find a max
      if (nrow(max.other.fluorescence.data.row) == 0) {
        max.other.fluorescence.data.row  = fluorescence.data[1,]
        max.other.fluorescence.data.row$time.min[1] = NA
      }
      
      cat("Time of other production rate:", max.other.fluorescence.data.row$time.min[1], "min\n")
      cat("Max other production:", max.other.fluorescence.data.row$other.rate, "per hour\n")
    }
    
    
    
    if (num.readings == 2) {
      
      if(!is.na(max.growth.rate.row$growth.rate) && !is.na(max.GFP.fluorescence.data.row$GFP.rate)) {
        strain.max.values = bind_rows(strain.max.values, 
                                      data.frame(
                                        well = this.well,
                                        growth.rate = max.growth.rate.row$growth.rate,
                                        GFP.rate = max.GFP.fluorescence.data.row$GFP.rate,
                                        max.growth.rate.time = max.growth.rate.row$time.min[1],
                                        max.GFP.rate.time = max.GFP.fluorescence.data.row$time.min[1],
                                        growth.fit.ss.residuals = max.growth.rate.row$growth.fit.ss.residuals[1],
                                        GFP.fit.ss.residuals = max.GFP.fluorescence.data.row$GFP.fit.ss.residuals[1]
                                      )
        )
      }
      
    } else  if (num.readings == 3) {
      
      if(!is.na(max.growth.rate.row$growth.rate) && !is.na(max.GFP.fluorescence.data.row$GFP.rate)) {
        
        strain.max.values = bind_rows(strain.max.values, 
                                      data.frame(
                                        well = this.well,
                                        growth.rate = max.growth.rate.row$growth.rate,
                                        GFP.rate = max.GFP.fluorescence.data.row$GFP.rate,
                                        other.rate = max.other.fluorescence.data.row$other.rate,
                                        max.growth.rate.time = max.growth.rate.row$time.min[1],
                                        max.GFP.rate.time = max.GFP.fluorescence.data.row$time.min[1],
                                        max.other.rate.time = max.other.fluorescence.data.row$time.min[1],
                                        growth.fit.ss.residuals = max.growth.rate.row$growth.fit.ss.residuals[1],
                                        GFP.fit.ss.residuals = max.GFP.fluorescence.data.row$GFP.fit.ss.residuals[1],
                                        other.fit.ss.residuals = max.other.fluorescence.data.row$other.fit.ss.residuals[1]
                                      )
        )
      }
    }
  }
  
  cat("\n\n")
  cat("Growth rate:", mean(strain.max.values$growth.rate), "±", sd(strain.max.values$growth.rate), "\n")
  cat("GFP rate:", mean(strain.max.values$GFP.rate), "±", sd(strain.max.values$GFP.rate), "\n")
  
  if (num.readings == 3) {
    cat("Other rate:", mean(strain.max.values$other.rate), "±", sd(strain.max.values$other.rate), "\n")
  }
  cat("\n\n")
  
  if (!is.na(mean(strain.max.values$growth.rate)) && length(strain.max.values$growth.rate) > 1) {
    ttest.growth.rate = t.test(strain.max.values$growth.rate)
  } else {
    ttest.growth.rate = c()
    ttest.growth.rate$conf.int = c(NA, NA)
  }
  
  if (!is.na(mean(strain.max.values$GFP.rate)) && length(strain.max.values$GFP.rate) > 1) {
    ttest.GFP.rate = t.test(strain.max.values$GFP.rate)
  } else {
    ttest.growth.rate = c()
    ttest.growth.rate$conf.int = c(NA, NA)
  }
  
  if (num.readings == 2) {
    

    final.table.summary = bind_rows(final.table.summary, 
                                    data.frame(
                                      strain.isolate = strain.of.interest,
                                      replicates = nrow(strain.max.values),
                                      wells = paste(strain.max.values$well, collapse=","),
                                      growth.rate = mean(strain.max.values$growth.rate),
                                      growth.rate.sd = sd(strain.max.values$growth.rate),
                                      growth.rate.sem = sd(strain.max.values$growth.rate)/sqrt(length(strain.max.values$growth.rate)),
                                      growth.rate.95L = ttest.growth.rate$conf.int[1],
                                      growth.rate.95U = ttest.growth.rate$conf.int[2],
                                      GFP.rate = mean(strain.max.values$GFP.rate),
                                      GFP.rate.sd = sd(strain.max.values$GFP.rate),
                                      GFP.rate.sem = sd(strain.max.values$GFP.rate)/sqrt(length(strain.max.values$GFP.rate)),
                                      GFP.rate.95L = ttest.GFP.rate$conf.int[1],
                                      GFP.rate.95U = ttest.GFP.rate$conf.int[2]
                                    )
    )
    
    final.table.all.wells = bind_rows(final.table.all.wells, 
                                      data.frame(
                                        strain.isolate = rep(strain.of.interest, nrow(strain.max.values)),
                                        well = strain.max.values$well,
                                        growth.rate = strain.max.values$growth.rate,
                                        GFP.rate = strain.max.values$GFP.rate,
                                        max.growth.rate.time = strain.max.values$max.growth.rate.time,
                                        max.GFP.rate.time = strain.max.values$max.GFP.rate.time,
                                        growth.fit.ss.residuals = strain.max.values$growth.fit.ss.residuals,
                                        GFP.fit.ss.residuals = strain.max.values$GFP.fit.ss.residuals
                                      )
    )
    
  } else  if (num.readings == 3) {
    
    if (!is.na(mean(strain.max.values$other.rate)) && length(strain.max.values$other.rate) > 1) {
      ttest.other.rate = t.test(strain.max.values$other.rate)
    } else {
      ttest.other.rate = c()
      ttest.other.rate$conf.int = c(NA, NA)
    }
    
    final.table.summary = bind_rows(final.table.summary, 
                                    data.frame(
                                      strain.isolate = strain.of.interest,
                                      replicates = nrow(strain.max.values),
                                      wells = paste(strain.max.values$well, collapse=","),
                                      growth.rate = mean(strain.max.values$growth.rate),
                                      growth.rate.sd = sd(strain.max.values$growth.rate),
                                      growth.rate.sem = sd(strain.max.values$growth.rate)/sqrt(length(strain.max.values$growth.rate)),
                                      growth.rate.95L = ttest.growth.rate$conf.int[1],
                                      growth.rate.95U = ttest.growth.rate$conf.int[2],
                                      GFP.rate = mean(strain.max.values$GFP.rate),
                                      GFP.rate.sd = sd(strain.max.values$GFP.rate),
                                      GFP.rate.sem = sd(strain.max.values$GFP.rate)/sqrt(length(strain.max.values$GFP.rate)),
                                      GFP.rate.95L = ttest.GFP.rate$conf.int[1],
                                      GFP.rate.95U = ttest.GFP.rate$conf.int[2],
                                      other.rate = mean(strain.max.values$other.rate),
                                      other.rate.sd = sd(strain.max.values$other.rate),
                                      other.rate.sem = sd(strain.max.values$other.rate)/sqrt(length(strain.max.values$other.rate)),
                                      other.rate.95L = ttest.other.rate$conf.int[1],
                                      other.rate.95U = ttest.other.rate$conf.int[2]
                                    )
    )
    
    final.table.all.wells = bind_rows(final.table.all.wells, 
                                      data.frame(
                                        strain.isolate = rep(strain.of.interest, nrow(strain.max.values)),
                                        well = strain.max.values$well,
                                        growth.rate = strain.max.values$growth.rate,
                                        GFP.rate = strain.max.values$GFP.rate,
                                        other.rate = strain.max.values$other.rate,                              
                                        max.growth.rate.time = strain.max.values$max.growth.rate.time,
                                        max.GFP.rate.time = strain.max.values$max.GFP.rate.time,
                                        max.other.rate.time = strain.max.values$max.other.rate.time,
                                        growth.fit.ss.residuals = strain.max.values$growth.fit.ss.residuals,
                                        GFP.fit.ss.residuals = strain.max.values$GFP.fit.ss.residuals,
                                        other.fit.ss.residuals = strain.max.values$other.fit.ss.residuals
                                      )
    )
  }
  
  
}

##Split out strain and replicates, fix column order, sort rows

#### Summmry table
final.table.summary$isolate= sub("^.+__", "", final.table.summary$strain.isolate, perl = T)
final.table.summary$strain = sub("__.+$", "", final.table.summary$strain.isolate, perl = T)
final.table.summary = final.table.summary %>% select(-strain.isolate)

final.column.order = c("strain", "isolate", "replicates", "wells", "growth.rate", "growth.rate.sd", "growth.rate.sem", "growth.rate.95L", "growth.rate.95U", "GFP.rate", "GFP.rate.sd", "GFP.rate.sem", "GFP.rate.95L", "GFP.rate.95U")
if (num.readings == 3) {
  final.column.order = c(final.column.order, "other.rate", "other.rate.sd", "other.rate.sem", "other.rate.95L", "other.rate.95U")
}
final.table.summary = final.table.summary %>% select(all_of(final.column.order)) %>% arrange(strain, isolate)

#### All table
final.table.all.wells$isolate= sub("^.+__", "", final.table.all.wells$strain.isolate, perl = T)
final.table.all.wells$strain = sub("__.+$", "", final.table.all.wells$strain.isolate, perl = T)
final.table.all.wells = final.table.all.wells %>% select(-strain.isolate)

final.column.order = c("strain", "isolate", "well", "growth.rate", "GFP.rate")
if (num.readings== 3) {
  final.column.order = c(final.column.order, "other.rate")
}

final.column.order = c(final.column.order, "max.growth.rate.time", "max.GFP.rate.time")
if (num.readings== 3) {
  final.column.order = c(final.column.order, "max.other.rate.time")
}

final.column.order = c(final.column.order, "growth.fit.ss.residuals", "GFP.fit.ss.residuals")
if (num.readings== 3) {
  final.column.order = c(final.column.order, "other.fit.ss.residuals")
}

final.table.all.wells = final.table.all.wells %>% select(all_of(final.column.order)) %>% arrange(strain, isolate, well)

### Write both final tables
write_csv(final.table.summary, paste0(option.output.prefix, ".rates.summary.csv"))
write_csv(final.table.all.wells, paste0(option.output.prefix, ".rates.all.csv"))

