#!/usr/bin/env Rscript

library(tidyverse)
library(ggplot2)
library(cowplot)

input.file.string = "05-burden-final-output/output.final.csv"
sequencing.results.file.string = "igem2019_sequencing_results.csv"
output.prefix = "05-burden-final-output/output"

final.output = read_csv(input.file.string)
sequencing.results = read_csv(sequencing.results.file.string, comment = "#")

final.output = final.output %>% select(accession)
final.output$burden.measured = TRUE

sequencing.results = sequencing.results %>% select(-burden.measured) %>% left_join(final.output, by="accession")

write_csv(sequencing.results, paste0(output.prefix, ".sequencing.results.with.burden.measured.csv"), na = "")

# Now analyze a bit

distinct.sequencing.results = sequencing.results

#Remove genomes
nrow(distinct.sequencing.results)
distinct.sequencing.results = distinct.sequencing.results %>% filter(accession !="E. coli genome")
nrow(distinct.sequencing.results)

distinct.sequencing.results = distinct.sequencing.results %>% filter(!(accession %in% c("K3174002", "K3174003", "K3174004", "K3174006", "K3174007")))
nrow(distinct.sequencing.results)

distinct.sequencing.results = distinct.sequencing.results %>% filter(!(accession %in%
c(
"J23100",
"J23101",
"J23102",
"J23103",
"J23104",
"J23105",
"J23106",
"J23107",
"J23110",
"J23112",
"J23113",
"J23114",
"J23115",
"J23116",
"J23117",
"J23118"
)
))
nrow(distinct.sequencing.results)

distinct.sequencing.results = distinct.sequencing.results %>% group_by(accession, strain, backbone, host, sequencing.update, biobrick.prefix,	biobrick.suffix,	backbone.sequence,	backbone.discrepancies,	biobrick.sequence,	biobrick.discrepancies,	burden.measured) %>% summarize()

nrow(distinct.sequencing.results)

length(unique(distinct.sequencing.results$strain))

write_csv(distinct.sequencing.results, paste0(output.prefix, ".distinct.sequencing.results.with.burden.measured.csv"), na = "")

#check for duplicates by strain
distinct.sequencing.results %>% group_by(strain) %>% summarize(n=n()) %>% filter(n>1)

#check for duplicate BioBricks
distinct.sequencing.results.summary = distinct.sequencing.results %>% group_by(accession) %>% summarize(n=n())

distinct.sequencing.results = distinct.sequencing.results %>% left_join(distinct.sequencing.results.summary, by="accession")

write_csv(distinct.sequencing.results, paste0(output.prefix, ".distinct.strain.sequencing.results.with.burden.measured.csv"), na = "")

distinct.sequencing.results = distinct.sequencing.results %>% group_by(accession, backbone, host, sequencing.update, biobrick.prefix,	biobrick.suffix,	backbone.sequence,	backbone.discrepancies,	biobrick.sequence,	biobrick.discrepancies,	burden.measured) %>% summarize()

write_csv(distinct.sequencing.results, paste0(output.prefix, ".distinct.accession.sequencing.results.with.burden.measured.csv"), na = "")

# Only one that is included did not have its burden measured
distinct.sequencing.results %>% group_by(burden.measured) %>% summarize(n=n())


distinct.sequencing.results %>% group_by(biobrick.sequence) %>% summarize(n=n())

 