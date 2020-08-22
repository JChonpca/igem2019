#!/usr/bin/env Rscript

#igem2019 specific analyses

library(tidyverse)
library(ggplot2)
library(cowplot)

#debug
input.all.wells.file.string="04-normalization/output.no.burden.and.control.normalized.all.wells.csv"
input.part.file.string="04-normalization/output.part.burden.csv"

output.base.name="05-burden-final-output/output"



##############  Read in the input files

all.wells.data = read_csv(input.all.wells.file.string)
part.data = read_csv(input.part.file.string)


############## Filter out the ones with 2 replicates only!

part.data = part.data %>% filter(replicates>=2)
remaining.parts = unique(part.data$accession )
all.wells.data = all.wells.data %>% filter (accession %in% remaining.parts)


# Recalculate burden in all.wells.data
all.wells.data$burden = 1 - all.wells.data$normalized.growth.rate

# reorder!
part.data = part.data %>% arrange(-burden.mean)
part.data$accession = factor(part.data$accession, levels=part.data$accession)

all.wells.data$accession = factor(all.wells.data$accession, levels=part.data$accession)

accession.order = part.data

############## Huge bar graphs

p = ggplot(all.wells.data, aes(x=accession, y = burden)) +
   geom_bar(data=part.data, aes(x=accession, y = burden.mean, fill=burden.category), stat="identity", color=NA) +
   geom_point(alpha=0.5, shape = 16, size=0.5) +
   scale_y_continuous(breaks = seq(-0.1, 0.6, by=0.1), labels= c("-10%", "0%", "10%", "20%", "30%", "40%", "50%", "60%")) +
   scale_fill_manual(values=c("#4147F7", "#E1E1E1", "#FF6C6C")) +
   theme_linedraw() +
   theme_classic() +
   theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
   NULL

p

ggsave(paste0(output.base.name, ".burden.bar.graphs.pdf"), plot=p, width=14, height=7)

