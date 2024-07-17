#!/usr/bin/env Rscript

library(tidyverse)
library(ggplot2)
library(cowplot)


unmutated_BFP_1 = read_csv("12-BFP-unmutated-output/exp059.rates.all.csv")
unmutated_BFP_2 = read_csv("12-BFP-unmutated-output/exp060.rates.all.csv")
unmutated_BFP = unmutated_BFP_1 %>% bind_rows(unmutated_BFP_2)

burden_assay_normalized_BFP = read_csv("04-normalization/output.no.burden.and.control.normalized.all.wells.csv")
burden_assay_normalized_BFP = burden_assay_normalized_BFP %>% filter(strain %in% c("JEB1204", "JEB1205", "JEB1206", "JEB1207", "JEB1208"))

strain_correspondence = data.frame(
  strain = c(
    "JEB1204",
    "JEB1205",
    "JEB1206",
    "JEB1207",
    "JEB1208",
    "JEB1209",
    "JEB1210",
    "JEB1211",
    "JEB1212",
    "JEB1213"
  ), 
  BFP_construct = c(
    "K3174002",
    "K3174003",
    "K3174004",
    "K3174006",
    "K3174007",
    "K3174002",
    "K3174003",
    "K3174004",
    "K3174006",
    "K3174007"
  )
)

unmutated_BFP = unmutated_BFP %>% left_join(strain_correspondence, by="strain")
burden_assay_normalized_BFP = burden_assay_normalized_BFP %>% left_join(strain_correspondence, by="strain")

#now figure out a coefficient that will adjust the unmutated BFP measurements to the same scale for JEB1206, JEB1207, and JEB1208

unmutated_BFP_unmutated_subset = unmutated_BFP %>% filter(BFP_construct %in% c("K3174004", "K3174006", "K3174007"))
burden_assay_normalized_BFP_unmutated_subset = burden_assay_normalized_BFP %>% filter(BFP_construct %in% c("K3174004", "K3174006", "K3174007"))


unmutated_BFP_unmutated_subset_summary = unmutated_BFP_unmutated_subset %>% group_by(BFP_construct)  %>% summarize(unmutated.growth.rate = mean(growth.rate))

burden_assay_normalized_BFP_unmutated_subset_summary = burden_assay_normalized_BFP_unmutated_subset %>% group_by(BFP_construct)  %>% summarize(normalized.growth.rate = mean(normalized.growth.rate))

combined = unmutated_BFP_unmutated_subset_summary %>% left_join(burden_assay_normalized_BFP_unmutated_subset_summary, by="BFP_construct")

fit_for_unmutated = lm(normalized.growth.rate ~ unmutated.growth.rate + 0, combined )

unmutated_BFP$normalized.growth.rate = unmutated_BFP$growth.rate * coef(fit_for_unmutated)[1]

## Now plot a comparison of the normalized

unmutated_BFP$type="unmutated"
burden_assay_normalized_BFP$type="burden"

combined_for_plot = burden_assay_normalized_BFP %>% select(BFP_construct, normalized.growth.rate, type)  %>% bind_rows(unmutated_BFP %>% select(BFP_construct, normalized.growth.rate, type))

combined_for_plot$group=paste0(combined_for_plot$BFP_construct,combined_for_plot$type) 

combined_for_plot_stats = combined_for_plot %>% group_by(BFP_construct, type) %>% 
  summarize(replicates=n(),
            normalized.growth.rate.mean = mean(normalized.growth.rate), 
            normalized.growth.rate.sd = sd(normalized.growth.rate),
            normalized.growth.rate.sem = normalized.growth.rate.sd/sqrt(replicates),
            normalized.growth.rate.95CI.range =  normalized.growth.rate.sem*qt(0.975, df=replicates-1),
            normalized.growth.rate.95CI.max =  normalized.growth.rate.mean+normalized.growth.rate.95CI.range,
            normalized.growth.rate.95CI.min =  normalized.growth.rate.mean-normalized.growth.rate.95CI.range,
  )

combined_for_plot_stats$group=paste0(combined_for_plot_stats$BFP_construct,combined_for_plot_stats$type) 

ggplot(combined_for_plot, aes(x=group,y=normalized.growth.rate, color=type)) + 
  geom_jitter(width=0.2, alpha=0.6, stroke=NA, size=2.5) + 
  geom_boxplot(data=combined_for_plot_stats, aes(x=group, y=normalized.growth.rate.mean, lower=normalized.growth.rate.mean, upper=normalized.growth.rate.mean, middle=normalized.growth.rate.mean, ymin=normalized.growth.rate.95CI.min, ymax=normalized.growth.rate.95CI.max ), stat = "identity", color="black") +
  theme_cowplot(12) + 
  scale_y_continuous(lim=c(0.4, 1), breaks=4:10*0.1) + 
  scale_colour_manual(name = "Strain",values = c("#0072B2", "#009E73")) +
  panel_border(color = "black")

ggsave("12-BFP-unmutated-output/normalized-burden-to-unmutated-comparison.pdf", width=8, height=5)

the_fit = lm( normalized.growth.rate~BFP_construct+type:BFP_construct+0, data = combined_for_plot)
summary(the_fit)

# Now do t-tests to compare for each strain in each condition
for (this.BFP.construct in unique(strain_correspondence$BFP_construct)) {
  
  the_fit = lm( normalized.growth.rate~type, data = combined_for_plot %>% filter(BFP_construct==this.BFP.construct))
  cat("testing: ", this.BFP.construct, "\n")
  print(summary(the_fit))
  
}

write_csv(combined_for_plot, "12-BFP-unmutated-output/combined_data_for_plot.csv")
write_csv(combined_for_plot_stats, "12-BFP-unmutated-output/combined_data_for_plot_stats.csv")

write_csv(unmutated_BFP, "12-BFP-unmutated-output/unmutated_BFP_stats.csv")

#Conclusion: Only the two suspected ones are different after correcting for multiple testing.


