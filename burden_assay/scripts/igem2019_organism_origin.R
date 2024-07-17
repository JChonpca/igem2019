#!/usr/bin/env Rscript

library(tidyverse)
library(ggplot2)
library(cowplot)

input.prefix = "04-normalization/output"
output.prefix = "05-burden-final-output/output"
strain.metadata.file.string = "igem2019_strain_metadata.csv"
part.metadata.file.string = "igem2019_part_metadata.csv"

strain.metadata = read.csv(strain.metadata.file.string)
part.metadata = read.csv(part.metadata.file.string)

organism_origin_columns = c(
  "manual.curation.source.organism.all.classification",
  "manual.curation.source.organism.CDS.classification",
  "manual.curation.source.organism.all.no.FP.classification",
  "manual.curation.source.organism.CDS.no.FP.classification"
)

# Remove the virus/phage part 
for (ooc in organism_origin_columns) {
  part.metadata[[ooc]] = str_replace(part.metadata[[ooc]], coll(" (virus)"), "")
  part.metadata[[ooc]] = str_replace(part.metadata[[ooc]], coll(" (phage)"), "")
}


#final.all.measurements = read_csv(paste0(output.prefix, ".final.all.measurements.csv"))
#final.all.measurements$accession = toupper(final.all.measurements$accession)

#Analyze just the means of the replicate measurements
final.measurements = read_csv(paste0(output.prefix, ".final.csv"))
final.measurements$accession = toupper(final.measurements$accession)


# Add the organism origin columns
annotated.final.measurements = final.measurements %>% left_join(part.metadata, by="accession")

# Order things
ordered_organism_origins = c(
  "Escherichia coli",
  "Other-Enterobacterales",
  "Other-Gammaproteobacteria",
  "Other-Proteobacteria",
  "Other-Bacteria",
  "Fungi",
  "Plant",
  "Animal",
  "Eukaryotes",
  "Other",
  "Synthetic",
  "None"
)

for (ooc in organism_origin_columns) {
  annotated.final.measurements[[ooc]] = factor(annotated.final.measurements[[ooc]], levels=ordered_organism_origins)
}

#Create some columns that summarize things a bit more

classification.to.fewer.categories = data.frame(
  from = c(
    "Escherichia coli",
    "Other-Enterobacterales",
    "Other-Gammaproteobacteria",
    "Other-Proteobacteria",
    "Other-Bacteria",
    "Fungi",
    "Plant",
    "Animal",
    "Synthetic",
    "None"
  ),
  to = c(
    "Escherichia coli",
    "Other-Gammaproteobacteria",
    "Other-Gammaproteobacteria",
    "Other-Bacteria",
    "Other-Bacteria",
    "Eukaryotes",
    "Eukaryotes",
    "Eukaryotes",
    "None",
    "None"
  )
)

for (ooc in organism_origin_columns) {
  
  to.column.name = paste0(ooc, ".fewer")
  from.column.name = ooc
  translator = classification.to.fewer.categories %>% rename( !!to.column.name := to, !!from.column.name := from)
  annotated.final.measurements = annotated.final.measurements %>% left_join(translator, by=from.column.name)
}

classification.to.fewest.categories = data.frame(
  from = c(
    "Escherichia coli",
    "Other-Enterobacterales",
    "Other-Gammaproteobacteria",
    "Other-Proteobacteria",
    "Other-Bacteria",
    "Fungi",
    "Plant",
    "Animal",
    "Synthetic",
    "None"
  ),
  to = c(
    "Escherichia coli",
    "Other",
    "Other",    
    "Other",
    "Other",
    "Other",
    "Other",
    "Other",
    "None",
    "None"
  )
)

for (ooc in organism_origin_columns) {
  
  to.column.name = paste0(ooc, ".fewest")
  from.column.name = ooc
  translator = classification.to.fewest.categories %>% rename( !!to.column.name := to, !!from.column.name := from)
  annotated.final.measurements = annotated.final.measurements %>% left_join(translator, by=from.column.name)
}


write_csv(annotated.final.measurements, paste0(output.prefix, ".organism.annotated.final.csv"))


# Ditch the controls before analyzing and check that we have the expected number of BioBricks
annotated.final.measurements = annotated.final.measurements %>% filter(burden.category != "control")
cat("Number of (non-control) BiooBricks analyzed: ", nrow(annotated.final.measurements), "\n")

annotated.final.measurements = annotated.final.measurements %>% arrange(burden.category)

#########################################################################

###### Helper function

okabe_ito <- c("#555555", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
               "#0072B2", "#D55E00", "#CC79A7", "#6666FF")

plot_and_analyze <- function(df, category_column, output_name, dropping_none=T) {
  
  cat("=========================================\n")
  cat("Kruskal-Wallis test and distribution plot\n")
  cat("Analyzing by column: ", category_column, "\n")
  cat("Dropping None: ", dropping_none, "\n")
  
  X = df
  
  X$categories = factor(X[[category_column]], levels=ordered_organism_origins)
  
  if (dropping_none) {
    X = X %>% filter(categories != "None")
  }
  
  cat("Number of BioBricks: ", nrow(X), "\n")
  
  Y = model1 = lm( normalized.growth.rate.mean~categories+0, X)
  model2 = lm( normalized.growth.rate.mean~1, X)
  anova(model1, model2)
  
  kt = kruskal.test( normalized.growth.rate.mean~categories+0, X)
  print(kt)
  
  Xm = X %>% group_by(categories) %>% summarize(grand.mean.normalized.growth.rate = mean(normalized.growth.rate.mean))
  # Calc means
  
  X$color.categories = factor(X$categories, levels=c("AAAA",ordered_organism_origins))
  X$color.categories[X$burden.category!="significant"] = "AAAA"

  ggplot(X, aes(x=categories, y=normalized.growth.rate.mean)) +
    geom_jitter(width=0.25, alpha=0.5, stroke=NA, size=2.5, aes(color=color.categories, fill=color.categories)) + 
    geom_boxplot(aes(y = grand.mean.normalized.growth.rate, ymin=grand.mean.normalized.growth.rate, ymax=grand.mean.normalized.growth.rate, lower=grand.mean.normalized.growth.rate, middle=grand.mean.normalized.growth.rate, upper=grand.mean.normalized.growth.rate), stat = "identity", data=Xm, color="black", width=0.7) +
    geom_hline(yintercept=1, linetype='dashed') +
    scale_color_manual(values = okabe_ito) +
    theme_cowplot(12) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position="none", ) +
    panel_border(color = "black") +
    scale_y_continuous(breaks=(0:10)*0.1) +
    labs(x = NULL, y="Normalized growth rate")
  
  ggsave(paste0(output.prefix, ".organism.origin.", output_name, ".pdf"), width=5, height=5)
}


#########################################################################
# Kruskal-Wallis test on various subsets

plot_and_analyze(annotated.final.measurements, "manual.curation.source.organism.all.classification", "all_categories")
plot_and_analyze(annotated.final.measurements, "manual.curation.source.organism.CDS.classification", "CDS_categories")
plot_and_analyze(annotated.final.measurements, "manual.curation.source.organism.all.no.FP.classification", "all_categories.no.FP")
plot_and_analyze(annotated.final.measurements, "manual.curation.source.organism.CDS.no.FP.classification", "CDS_categories.no.FP")

plot_and_analyze(annotated.final.measurements, "manual.curation.source.organism.all.classification.fewer", "all_categories.fewer")
plot_and_analyze(annotated.final.measurements, "manual.curation.source.organism.CDS.classification.fewer", "CDS_categories.fewer")
plot_and_analyze(annotated.final.measurements, "manual.curation.source.organism.all.no.FP.classification.fewer", "all_categories.no.FP.fewer")
plot_and_analyze(annotated.final.measurements, "manual.curation.source.organism.CDS.no.FP.classification.fewer", "CDS_categories.no.FP.fewer")

plot_and_analyze(annotated.final.measurements, "manual.curation.source.organism.all.classification.fewest", "all_categories.fewest")
plot_and_analyze(annotated.final.measurements, "manual.curation.source.organism.CDS.classification.fewest", "CDS_categories.fewest")
plot_and_analyze(annotated.final.measurements, "manual.curation.source.organism.all.no.FP.classification.fewest", "all_categories.no.FP.fewest")
plot_and_analyze(annotated.final.measurements, "manual.curation.source.organism.CDS.no.FP.classification.fewest", "CDS_categories.no.FP.fewest")

#########################################################################


run_binomial_test <- function(df, category_column, output_name, dropping_none=T) {
  cat("=========================================\n")
  cat("Binomial test\n")
  cat("Analyzing by column: ", category_column, "\n")
  cat("Dropping None: ", dropping_none, "\n")
  
  X = df
  X$categories = X[[category_column]]
  
  if (dropping_none) {
    X = X %>% filter(categories != "None")
  }
  
  cat("Number of BioBricks: ", nrow(X), "\n")
  cat("Number of burdensome BioBricks: ", nrow(X %>% filter(burden.category == "significant")), "\n")
  
  X$significant.burden = 0
  X$significant.burden[X$burden.category == "significant"] = 1

  g = glm(significant.burden ~ categories+0, data=X, family = binomial)
  g1 = glm(significant.burden ~ 1, data=X, family = binomial)
  print(anova(g, g1,test="Chisq" ))
}

#########################################################################
# Binomial test of burden versus no burden on various subsets


run_binomial_test(annotated.final.measurements, "manual.curation.source.organism.all.classification", "all_categories")
run_binomial_test(annotated.final.measurements, "manual.curation.source.organism.CDS.classification", "CDS_categories")
run_binomial_test(annotated.final.measurements, "manual.curation.source.organism.all.no.FP.classification", "all_categories.no.FP")
run_binomial_test(annotated.final.measurements, "manual.curation.source.organism.CDS.no.FP.classification", "CDS_categories.no.FP")

run_binomial_test(annotated.final.measurements, "manual.curation.source.organism.all.classification.fewer", "all_categories.fewer")
run_binomial_test(annotated.final.measurements, "manual.curation.source.organism.CDS.classification.fewer", "CDS_categories.fewer")
run_binomial_test(annotated.final.measurements, "manual.curation.source.organism.all.no.FP.classification.fewer", "all_categories.no.FP.fewer")
run_binomial_test(annotated.final.measurements, "manual.curation.source.organism.CDS.no.FP.classification.fewer", "CDS_categories.no.FP.fewer")

run_binomial_test(annotated.final.measurements, "manual.curation.source.organism.all.classification.fewest", "all_categories.fewest")
run_binomial_test(annotated.final.measurements, "manual.curation.source.organism.CDS.classification.fewest", "CDS_categories.fewest")
run_binomial_test(annotated.final.measurements, "manual.curation.source.organism.all.no.FP.classification.fewest", "all_categories.no.FP.fewest")
run_binomial_test(annotated.final.measurements, "manual.curation.source.organism.CDS.no.FP.classification.fewest", "CDS_categories.no.FP.fewest")

write_csv(annotated.final.measurements, file.path("05-burden-final-output", "output.organism.origin.csv"))

