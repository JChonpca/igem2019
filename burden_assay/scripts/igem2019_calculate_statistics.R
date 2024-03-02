#!/usr/bin/env Rscript

library(tidyverse)
library(ggplot2)
library(cowplot)

input.prefix = "04-normalization/output"
output.prefix = "05-burden-final-output/output"
strain.metadata.file.string = "igem2019_strain_metadata.csv"
part.metadata.file.string = "igem2019_part_metadata.csv"
all.data.file.string = paste0(input.prefix, ".no.burden.and.control.normalized.all.wells.csv")

part.input.file.string = paste0(output.prefix, ".burden_confidence.csv")
parts.means = read.csv(part.input.file.string)
parts.means = parts.means %>% filter( (replicates>=2) )

strain.metadata = read.csv(strain.metadata.file.string)
part.metadata = read.csv(part.metadata.file.string)

all.data = read.csv(all.data.file.string)

########################################################################
### Write final output files (missing the 1-off replicate ones)
parts.means = parts.means %>% filter( (replicates>=2) )
annotated.parts.means = parts.means %>% left_join(part.metadata, by="accession")
write_csv(annotated.parts.means, paste0(output.prefix, ".final.part.measurements.csv"))


#Clean up ones with only one replicate and remove unused columns
final.all.measurements.n = all.data %>% group_by(vector, accession) %>% 
  summarize(replicates=n())

final.all.measurements = all.data %>% left_join(final.all.measurements.n, by=c("vector", "accession"))
final.all.measurements = final.all.measurements %>% filter(replicates>=2)

final.all.measurements = final.all.measurements %>% select(-replicates, -plate.strain, -other.rate, -max.other.rate.time, -growth.fit.ss.residuals, -GFP.fit.ss.residuals, -isolate)

write_csv(final.all.measurements, paste0(output.prefix, ".final.all.measurements.csv"))
########################################################################
### Calculate KS test on different vectors

#we need to combine these to be BioBrick means!
strain.means = all.data %>% 
  filter(!strain %in% c("JEB1204", "JEB1205", "JEB1206", "JEB1207", "JEB1208")) %>% 
  group_by(vector, accession) %>% 
  summarize(normalized.growth.rate.mean=mean(normalized.growth.rate), replicates=n()) %>% 
   filter(replicates>=2)

ks.test(
  (strain.means %>% filter(vector=="pSB1A2"))$normalized.growth.rate.mean,
  (strain.means %>% filter(vector=="pSB1C3"))$normalized.growth.rate.mean
)

#Graph

p = ggplot(data=strain.means, aes(x=normalized.growth.rate.mean, group=vector, color=vector)) +
  stat_ecdf(geom = "step") +
  theme_cowplot(12) +
  panel_border(color = "black") +
  scale_x_continuous(limits=c(0.4, 1.2), breaks=seq(0.4,1.2,by=0.1)) +
  scale_color_manual(values=c("#FEC013", "#000000")) +
  NULL
p
ggsave(paste0(output.prefix, ".growth.rate.CDF.by.vector.pdf"), plot=p)

########################################################################
### Compare strains measured in both vectors

measured.twice = parts.means %>% filter(vectors=="pSB1A2,pSB1C3")
measured.twice.strains = (strain.metadata %>% filter(accession %in% measured.twice$accession))$strain

measured.twice.data = all.data %>% filter(strain %in% measured.twice.strains)

measured.twice.data = measured.twice.data %>%left_join(strain.metadata, by="strain")

fit = lm(normalized.growth.rate~accession+vector, measured.twice.data )
fit2 = lm(normalized.growth.rate~accession, measured.twice.data )
anova(fit,fit2)

measured.twice.data.mean = measured.twice.data %>% group_by(accession, vector) %>% summarize(normalized.growth.rate.mean=mean(normalized.growth.rate))

p = ggplot(data=measured.twice.data, aes(x=accession, y=normalized.growth.rate, color=vector, shape=vector)) +
  scale_color_manual(values=c("#FEC013", "#000000"))+
  geom_jitter(width=0.4, size=2) +
  geom_crossbar(data=measured.twice.data.mean, aes(y = normalized.growth.rate.mean, ymin = normalized.growth.rate.mean, ymax = normalized.growth.rate.mean),
                size=0.6, width = .5) +
  theme_cowplot(12) +
  panel_border(color = "black") +
  scale_y_continuous(limits=c(0.7, 1.2), breaks=seq(0.7,1.2,by=0.1)) +
  NULL
p
ggsave(paste0(output.prefix, ".growth.rate.parts.in.both.vectors.pdf"), plot=p, height=6, width=12)


########################################################################
### Calculate one-tailed tests of more than certain burden cutoffs

not.control.data = parts.means %>% filter(burden.category != "control")

cutoff.df = data.frame()

cat("burden cutoff", "# sig> than", "% sig > than", "\n")
for (b in seq(0, 0.5, by=0.05)) {
  
  p.values = pt((b-not.control.data$burden.mean)/not.control.data$burden.sem, not.control.data$replicates-1)
  #Uncomment for multiple testing correction
  #p.values = p.adjust(p.values, method="BH")
  num.significant = sum(p.values <= 0.05, na.rm=T)
  
  cutoff.df = bind_rows(cutoff.df, data.frame(greater.than.burden=b, num.significant.parts=num.significant, total.parts=length(p.values), percent.significant.parts=(num.significant/length(p.values))*100))
}

write.csv(cutoff.df, paste0(output.prefix, ".parts.with.at.least.burden.csv"))


not.control.data = parts.means %>% filter(burden.category != "control")

cutoff.df = data.frame()

cat("burden cutoff", "# sig> than", "% sig > than", "\n")
for (b in seq(0, 0.5, by=0.05)) {
  
  p.values = pt((b-not.control.data$burden.mean)/not.control.data$burden.sem, not.control.data$replicates-1)
  #Uncomment for multiple testing correction
  p.values = p.adjust(p.values, method="BH")
  num.significant = sum(p.values <= 0.05, na.rm=T)
  
  cutoff.df = bind_rows(cutoff.df, data.frame(greater.than.burden=b, num.significant.parts=num.significant, total.parts=length(p.values), percent.significant.parts=(num.significant/length(p.values))*100))
}

write.csv(cutoff.df, paste0(output.prefix, ".parts.with.at.least.burden.BH.adjusted.csv"))


