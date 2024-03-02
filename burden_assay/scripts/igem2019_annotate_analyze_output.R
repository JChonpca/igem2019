#!/usr/bin/env Rscript

### This script adds part information to the final output files
### and looks for promoter/RBS sequences associated with high burden


library(tidyverse)
library(ggplot2)
library(cowplot)

input.file.string = "05-burden-final-output/output.burden_confidence.csv"
part.metadata.file.string = "igem2019_part_metadata.csv"
output.prefix = "05-burden-final-output/output"

parts = read.csv(input.file.string)
parts.annotation = read.csv(part.metadata.file.string )

parts = parts %>% left_join(parts.annotation, by="accession")

####################################################################################
### Calculate Fisher's Exact Test on parts containing or not containing Promoter/RBS


significant.burden.parts = parts %>% filter(burden.category=="significant")
no.burden.parts = parts %>% filter(burden.category=="no burden")

#First construct the list of what we observed and in how many parts. Some have multiple promoters
parts$automated.curation.common.promoter = as.character(parts$automated.curation.common.promoter)
parts$automated.curation.common.RBS = as.character(parts$automated.curation.common.RBS)

list.promoters = c()
list.RBSs = c()

for (i in 1:nrow(parts)) {
  
  if (parts$automated.curation.common.promoter[i] == "") {
    this.promoters = c()
  } else {
    this.promoters = unlist(strsplit(parts$automated.curation.common.promoter[i], "; "))
  }
  list.promoters = c(list.promoters, this.promoters)
  
  if (parts$automated.curation.common.RBS[i] == "") {
    this.RBSs = c()
  } else {
    this.RBSs = unlist(strsplit(parts$automated.curation.common.RBS[i], "; "))
  }
  list.RBSs = c(list.RBSs, this.RBSs)
  
}

df.promoters = data.frame(accession=list.promoters)
promoter.table = count(df.promoters, accession)

promoter.strengths = data.frame (
  accession=c("J23112",
              "J23103",
              "J23113",
              "J23109",
              "J23117",
              "J23114",
              "J23115",
              "J23116",
              "J23105",
              "J23110",
              "J23107",
              "J23106",
              "J23108",
              "J23118",
              "J23111",
              "J23101",
              "J23104",
              "J23102",
              "J23100",
              "J23119"
              ),
  au=c(1,
       17,
       21,
       106,
       162,
       256,
       387,
       396,
       623,
       844,
       908,
       1185,
       1303,
       1429,
       1487,
       1791,
       1831,
       2179,
       2547,
       2000 # not really from table but it is the strong reference
  )
)

promoter.table = promoter.table %>% left_join(promoter.strengths, by="accession")

df.RBSs = data.frame(accession=list.RBSs)
RBS.table = count(df.RBSs, accession)

RBS.strengths = data.frame (
  accession=c("B0030", "B0032", "B0034"),
  au=c(0.6, 0.3, 1.0)
)

RBS.table = RBS.table %>% left_join(RBS.strengths, by="accession")

promoter.table$category="none"
promoter.table$category[promoter.table$au > 0] = "weak"
promoter.table$category[promoter.table$au >= 500] = "medium"
promoter.table$category[promoter.table$au >= 1500] = "strong"

promoter.table = promoter.table %>% arrange(au)
promoter.table$accession = factor(promoter.table$accession, levels = promoter.table$accession)


RBS.table$category="none"
RBS.table$category[RBS.table$au >= 0.3] = "weak"
RBS.table$category[RBS.table$au >= 0.6] = "medium"
RBS.table$category[RBS.table$au >= 1.0] = "strong"

RBS.table = RBS.table %>% arrange(au)
RBS.table$accession = factor(RBS.table$accession, levels = RBS.table$accession)

# Graph the relative strengths

ggplot(promoter.table, aes(x=accession, y=au, fill=category)) +
  geom_bar(stat="identity") +
  scale_y_continuous(expand=expansion(mult=c(0,0.05))) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_text(aes(label = n, y=2700)) +
  NULL
ggsave(paste0(output.prefix, ".promoter.au.pdf"))

ggplot(RBS.table, aes(x=accession, y=au, fill=category)) +
  geom_bar(stat="identity") +
  scale_y_continuous(expand=expansion(mult=c(0,0.05))) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_text(aes(label = n, y=1.02)) +
  NULL
ggsave(paste0(output.prefix, ".RBS.au.pdf"))

## Now assign a promoter.category to all rows in the parts table based
## on the strongest one that is present

parts$max.promoter.au = 0
parts$max.RBS.au = 0

for (i in 1:nrow(parts)) {
  
  if (parts$automated.curation.common.promoter[i] == "") {
    this.promoters = c()
  } else {
    this.promoters = unlist(strsplit(parts$automated.curation.common.promoter[i], "; "))
  }
  
  for (this.promoter in this.promoters) {
    this.promoter.info  = promoter.table %>% filter(accession == this.promoter)
    parts$max.promoter.au[i] = max(parts$max.promoter.au[i], this.promoter.info$au)
  }
  
  if (parts$automated.curation.common.RBS[i] == "") {
    this.RBSs = c()
  } else {
    this.RBSs = unlist(strsplit(parts$automated.curation.common.RBS[i], "; "))
  }
  
  for (this.RBS in this.RBSs) {
    this.RBS.info  = RBS.table %>% filter(accession == this.RBS)
    parts$max.RBS.au[i] = max(parts$max.RBS.au[i], this.RBS.info$au)
  }
}


## categorize based on the strongest ones present...

parts$promoter.category="none"
parts$promoter.category[parts$max.promoter.au > 0] = "weak"
parts$promoter.category[parts$max.promoter.au >= 500] = "medium"
parts$promoter.category[parts$max.promoter.au >= 1500] = "strong"

parts$promoter.any = F
parts$promoter.any[parts$promoter.category != "none"] = T

parts$RBS.category="none"
parts$RBS.category[parts$max.RBS.au >= 0.3] = "weak"
parts$RBS.category[parts$max.RBS.au >= 0.6] = "medium"
parts$RBS.category[parts$max.RBS.au >= 1.0] = "strong"

#parts$strong.promoter.and.strong.RBS = FALSE
#parts$strong.promoter.and.strong.RBS[(parts$RBS.category == "strong") & (parts$promoter.category == "strong")] = TRUE

parts.promoter.summary = parts %>% group_by(promoter.category, burden.category) %>% count()
parts.promoter.summary = parts.promoter.summary %>% filter(burden.category != "control")

parts.promoter.summary.totals = parts.promoter.summary %>% group_by(promoter.category) %>% summarize(total=sum(n))
parts.promoter.summary = parts.promoter.summary %>% left_join(parts.promoter.summary.totals, by="promoter.category")
parts.promoter.summary$fraction = parts.promoter.summary$n / parts.promoter.summary$total
parts.promoter.summary$promoter.category = factor(parts.promoter.summary$promoter.category, levels = c("none", "weak", "medium", "strong"))

#graph these
ggplot(parts.promoter.summary, aes(x=promoter.category, y=fraction, fill=burden.category)) +
  geom_bar(stat="identity") +
  scale_y_continuous(expand=expansion(mult=c(0,0.05))) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_text(aes(label = total, y=1.02)) +
  NULL

ggsave(paste0(output.prefix, ".promoter.association.pdf"))


parts.promoter.summary.any = parts %>% group_by(promoter.any, burden.category) %>% count()
parts.promoter.summary.any = parts.promoter.summary.any %>% filter(burden.category != "control")

cat("Are parts with any annotated constitutive promoter more likely to be burdensome?\n")

fisher.test(matrix(parts.promoter.summary.any$n,
                   nrow = 2,
                   dimnames = list(Burden.Category = c("no burden", "significant"),
                                   Any.Promoter = c("False", "True"))))

parts.RBS.summary = parts %>% group_by(RBS.category, burden.category) %>% count()
parts.RBS.summary = parts.RBS.summary %>% filter(burden.category != "control")

parts.RBS.summary.totals = parts.RBS.summary %>% group_by(RBS.category) %>% summarize(total=sum(n))

parts.RBS.summary = parts.RBS.summary %>% left_join(parts.RBS.summary.totals, by="RBS.category")
parts.RBS.summary$fraction = parts.RBS.summary$n / parts.RBS.summary$total
parts.RBS.summary$RBS.category = factor(parts.RBS.summary$RBS.category, levels = c("none", "weak", "medium", "strong"))

#graph these
ggplot(parts.RBS.summary, aes(x=RBS.category, y=fraction, fill=burden.category)) +
  geom_bar(stat="identity") +
  scale_y_continuous(expand=expansion(mult=c(0,0.05))) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_text(aes(label = total, y=1.02)) +
  NULL

ggsave(paste0(output.prefix, ".RBS.association.pdf"))

parts.RBS.summary = parts.RBS.summary %>% filter(RBS.category != "weak")
parts.RBS.summary = parts.RBS.summary %>% filter(RBS.category != "medium")

cat("Are parts with any annotated strong RBS more likely to be burdensome?\n")

fisher.test(matrix(parts.RBS.summary$n,
                   nrow = 2,
                   dimnames = list(Burden.Category = c("no burden", "significant"),
                                   Strong.RBS = c("False", "True"))))

# Now iterate over categories and perform FET

for (promoter in unique(parts.promoter.summary)) {
  cat(promoter, "\n")
}

for (RBS in sort(unique(list.RBSs))) {
  cat(RBS, "\n")
}

# Remove some columns before output
parts = parts %>% select(
  -sequenced,
  -measured, 
  -manual.curation.description,
  -manual.curation.total.components,
  -manual.curation.has.promoter,
  -manual.curation.has.RBS,
  -manual.curation.has.CDS,
  -manual.curation.operon.1,
  -manual.curation.operon.2,
  -manual.curation.notes,
  -manual.curation.subparts,
  -max.promoter.au,
  -max.RBS.au,
  -promoter.any
)
write_csv(parts, paste0(output.prefix, ".final.csv"))
