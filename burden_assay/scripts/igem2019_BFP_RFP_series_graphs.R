#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(plotly))
suppressMessages(library(gridExtra))
suppressMessages(library(cowplot))
suppressMessages(library(optparse))
suppressMessages(library(xtable))
suppressMessages(library(ggrepel))

library(chemCal)
library(deming)

# RFP ##########################  

metadata = read.csv("igem2019_strain_metadata.csv")
metadata = metadata %>% select(strain,accession)

RFP.file.string = "10-RFP-series-output/exp057.rates.all.csv"
RFP.all.data = read_csv(RFP.file.string)

RFP.all.data = RFP.all.data %>% left_join(metadata, by="strain")

## Remove BB0 and controls that are not in same plasmid backbone
RFP.all.data = RFP.all.data %>% filter(! strain %in% c("BB0", "JEB1203", "JEB1204", "JEB1205", "JEB1206", "JEB1207", "JEB1208") )


#### Begin old regression analysis
#### This fitting all assumes error is only in the y-direction (GFP rate)
#### So it is not the best way to analyze - see Deming treatment below for what we use
RFP.fit.fixed.zero = lm(GFP.rate~growth.rate + 0, RFP.all.data)
RFP.slope.fixed.zero = coef(RFP.fit.fixed.zero)

RFP.fit.nonzero.intercept = lm(GFP.rate~growth.rate, RFP.all.data)
RFP.intercept.significantly.nonzero = anova(RFP.fit.fixed.zero, RFP.fit.nonzero.intercept)
RFP.y.intercept = coef(RFP.fit.nonzero.intercept)[1]
RFP.slope = coef(RFP.fit.nonzero.intercept)[2]

cat("RFP y-intercept is significantly different from zero?")
print(RFP.intercept.significantly.nonzero)

#What is the error on the x-intercept
RFP.y.intercept.model = inverse.predict(RFP.fit.nonzero.intercept, 0)
RFP.x.intercept.CL = RFP.y.intercept.model$`Confidence Limits`
RFP.x.intercept = RFP.y.intercept.model$Prediction
#### End old regression analysis

### Make graphs

RFP.file.string = "10-RFP-series-output/exp057.rates.summary.csv"
RFP.data = read_csv(RFP.file.string)

RFP.data = RFP.data %>% left_join(metadata, by="strain")

## Remove BB0 and controls that are not in same plasmid backbone
RFP.data = RFP.data %>% filter(! strain %in% c("BB0", "JEB1203", "JEB1204", "JEB1205", "JEB1206", "JEB1207", "JEB1208") )

# Correlation (uses means only)
cor(RFP.data$growth.rate, RFP.data$GFP.rate)

### Fit with maximum in a way that accounts for error in X and Y
RFP.deming.fit = deming(GFP.rate~growth.rate, RFP.data, xstd=RFP.data$growth.rate.sd, ystd=RFP.data$GFP.rate.sd,)
RFP.deming.y.intercept = coef(RFP.deming.fit)[1]
RFP.deming.slope = coef(RFP.deming.fit)[2]

RFP.promoter.strengths.file.string = "input-plate-data-RFP-series/RFP_promoter_strengths.csv"
RFP.promoter.strengths = read_csv(RFP.promoter.strengths.file.string)

RFP.data = RFP.data %>% left_join(RFP.promoter.strengths, by="strain")

# plot showing all strains burden vs growth 
burdenVsgrowthPlot = ggplot(RFP.data, aes(x=growth.rate, y=GFP.rate, color=other.rate))  +
  geom_errorbarh(aes(xmin=growth.rate.95L, xmax=growth.rate.95U), height=0) +
  geom_errorbar(aes(ymin=GFP.rate.95L, ymax=GFP.rate.95U), width=0) +
  geom_point(size=4, shape=16)  +
  scale_x_continuous(limits = c(0, 1.5), breaks = (0:7)*0.25 ) + 
  scale_y_continuous(limits = c(0, 30000), breaks = (0:6)*5000) +  
#  geom_abline(intercept=RFP.y.intercept, slope = RFP.slope, linetype=2) +
  geom_abline(intercept=RFP.deming.y.intercept, slope = RFP.deming.slope, linetype=1) +
  scale_color_gradient(low = "lightpink", high = "red3", na.value = "grey", limits=c(-200,6200), breaks=0:6*1000) +
#  scale_color_gradient(low = "lightpink", high = "red3", na.value = "gray") +
  labs(color = "RFP Rate") + 
#  geom_label_repel(aes(label = strain),size = 2.5,box.padding   = 0.35,point.padding = 0.5,segment.color = 'grey50') +
  ggtitle("RFP Series Burden Plot") + xlab("Growth Rate") + ylab("GFP Rate") + 
  theme_bw() + theme(panel.grid.minor = element_blank()) +
  coord_cartesian(expand = F) +
  NULL
burdenVsgrowthPlot
ggsave("10-RFP-series-output/growth_rate_versus_GFP_rate.pdf", width=7, height=6)


# growth BAR PLOT 
RFP.data = RFP.data %>% arrange(desc(growth.rate))
RFP.data$strain = factor(RFP.data$strain, levels=(RFP.data$strain))
growthRatePlot = ggplot(RFP.data, aes(x = strain, y=growth.rate, fill=other.rate)) +  
  geom_bar(size=3, stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=growth.rate.95L, ymax=growth.rate.95U), position=position_dodge()) + 
  scale_y_continuous(limits = c(0, max(RFP.data$growth.rate.95U)), expand=expansion(add = c(0, 0.1))) + 
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 12)) +
  scale_fill_gradient(low = "lightpink", high = "red3", na.value = "gray") +
  labs(fill = "RFP Rate") + 
  ggtitle("RFP Series Growth Rates") + xlab("Strain") + ylab("Growth Rate") +
  theme_bw() +
  NULL	
growthRatePlot 
ggsave("10-RFP-series-output/growth_rate_colored_by_RFP_rate.pdf")

#RFP Plot
RFP.data = RFP.data %>% arrange(desc(other.rate))
RFP.data$strain = factor(RFP.data$strain, levels=(RFP.data$strain))
rfpRatePlot = ggplot(RFP.data, aes(x =strain, y=other.rate, fill=promoter.strength)) +  
  geom_bar(size=3, stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=other.rate.95L, ymax=other.rate.95U), position=position_dodge()) + 
  scale_y_continuous(limits = c(-200, max(RFP.data$other.rate.95U)), expand=expansion(add = c(0, 500))) + 
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 12)) +
  scale_fill_gradient(low = "lightpink", high = "red3", na.value = "gray") +
  labs(fill = "Promoter Strength") + 
  ggtitle("RFP Rates") + xlab("Strain") + ylab("RFP Rate") +
  theme_bw() +
  NULL	
rfpRatePlot 
ggsave("10-RFP-series-output/RFP_rate_colored_by_promoter_strength.pdf")


# Do some stats on these...
RFP.file.string = "10-RFP-series-output/exp057.rates.all.csv"
no.plasmid.data = read_csv(RFP.file.string) %>% filter(strain=="BB0")
filtered.RFP.data = read_csv(RFP.file.string) %>% filter(! strain %in% c("BB0", "JEB1203", "JEB1204", "JEB1205", "JEB1206", "JEB1207", "JEB1208") )

cat("RFP strains have reduced growth rate?\n")

pvals = c()
pvals.more.than.40.percent = c()
control.growth.rates = no.plasmid.data$growth.rate
for(this.strain in unique(filtered.RFP.data$strain)) {
  cat("Testing strain:", this.strain, "\n")
  this.strain.data = filtered.RFP.data %>% filter(strain==this.strain)
  this.strain.growth.rates = this.strain.data$growth.rate
  res = t.test(control.growth.rates, this.strain.growth.rates, alternative="greater")
  print(res)
  pvals = c(pvals, res$p.value)
  
  this.strain.relative.growth.rates = this.strain.growth.rates / mean(control.growth.rates) - 0.6
  res = t.test(this.strain.relative.growth.rates, alternative="less")
  print(res)
  pvals.more.than.40.percent = c(pvals.more.than.40.percent, res$p.value)
}

#Correct for multiple testing
pvals = p.adjust(pvals, method="BH")
cat(as.character(sum(pvals<0.05)), "strains do after correction for multiple testing.\n")


# Individual value plots

# Growth rate plot
RFP.data = RFP.data %>% arrange(desc(growth.rate))
RFP.data$strain = factor(RFP.data$strain, levels=RFP.data$strain)
growthRatePlot = ggplot(RFP.data, aes(x = strain, y=growth.rate, fill="lightpink")) +  
  geom_bar(size=3, stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=growth.rate.95L, ymax=growth.rate.95U), position=position_dodge(), width=0.4) + 
  geom_jitter(data=RFP.all.data, alpha=0.5, size=1, width=0.24) +
  scale_x_discrete(labels=RFP.data$accession) +
  scale_y_continuous(limits = c(0, max(RFP.data$growth.rate.95U, RFP.all.data$growth.rate)), expand=expansion(mult = c(0, 0.05))) + 
  ggtitle("RFP Series Growth Rates") + xlab("Strain") + ylab("Growth Rate") +
  theme_bw() + 
  theme(legend.position = "None", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust=.5, hjust=1, size = 12)) +
  NULL
growthRatePlot 
ggsave("10-RFP-series-output/growth_rate_individual_points.pdf")


# Growth rate plot
growthRatePlot = ggplot(RFP.data, aes(x = strain, y=GFP.rate, fill="lightpink")) +  
  geom_bar(size=3, stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=GFP.rate.95L, ymax=GFP.rate.95U), position=position_dodge(), width=0.4) +
  geom_jitter(data=RFP.all.data, alpha=0.5, size=1, width=0.24) +
  scale_x_discrete(labels=RFP.data$accession) +
  scale_y_continuous(limits = c(0, max(RFP.data$GFP.rate.95U, RFP.all.data$GFP.rate)), expand=expansion(mult = c(0, 0.05))) + 
  ggtitle("RFP Series GFP Rates") + xlab("Strain") + ylab("GFP Rate") +
  theme_bw() + 
  theme(legend.position = "None", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust=.5, hjust=1, size = 12)) +
  NULL
growthRatePlot 
ggsave("10-RFP-series-output/GFP_individual_points.pdf")

# RFP rate plot
growthRatePlot = ggplot(RFP.data, aes(x = strain, y=other.rate, fill="lightpink")) +  
  geom_bar(size=3, stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=other.rate.95L, ymax=other.rate.95U), position=position_dodge(), width=0.4) +
  geom_jitter(data=RFP.all.data, alpha=0.5, size=1, width=0.24) +
  scale_x_discrete(labels=RFP.data$accession) +
  scale_y_continuous(limits = c(0, max(RFP.data$other.rate.95U, RFP.all.data$other.rate)), expand=expansion(mult = c(0, 0.05))) + 
  ggtitle("RFP Series RFP Rates") + xlab("Strain") + ylab("RFP Rate") +
  theme_bw() + 
  theme(legend.position = "None", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust=.5, hjust=1, size = 12)) +
  NULL
growthRatePlot 
ggsave("10-RFP-series-output/other_individual_points.pdf")

# BFP ##########################  


BFP.file.string = "11-BFP-series-output/exp061.rates.all.csv"
BFP.all.data = read_csv(BFP.file.string)

BFP.metadata = data.frame(
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
  accession = c(
    "K3174002-1",
    "K3174003-1",
    "K3174004-1",
    "K3174006-1",
    "K3174007-1",
    "K3174002-2",
    "K3174003-2",
    "K3174004-2",
    "K3174006-2",
    "K3174007-2"
  )
)

BFP.all.data = BFP.all.data %>% left_join(BFP.metadata, by="strain")


#### Begin old regression analysis
#### This fitting all assumes error is only in the y-direction (GFP rate)
#### So it is not the best way to analyze - see Demin treatment below for what we use

BFP.fit.fixed.zero = lm(GFP.rate~growth.rate + 0, BFP.all.data)
BFP.slope.fixed.zero = coef(RFP.fit.fixed.zero)

BFP.fit.nonzero.intercept = lm(GFP.rate~growth.rate, BFP.all.data)
BFP.intercept.significantly.nonzero = anova(BFP.fit.fixed.zero, BFP.fit.nonzero.intercept)
BFP.y.intercept = coef(BFP.fit.nonzero.intercept)[1]
BFP.slope = coef(BFP.fit.nonzero.intercept)[2]

cat("BFP y-intercept is significantly different from zero?")
print(BFP.intercept.significantly.nonzero)

#What is the error on the x-intercept
BFP.y.intercept.model = inverse.predict(BFP.fit.nonzero.intercept, 0)
BFP.x.intercept.CL = BFP.y.intercept.model$`Confidence Limits`
BFP.x.intercept = BFP.y.intercept.model$Prediction
#### End old regression

### Make graphs
BFP.file.string = "11-BFP-series-output/exp061.rates.summary.csv"
BFP.data = read_csv(BFP.file.string)
BFP.data$isolate=factor(BFP.data$isolate)
BFP.data = BFP.data %>% left_join(BFP.metadata, by="strain")


## Remove BB0 and controls that are not in same plasmid backbone
BFP.data = BFP.data %>% filter(! strain %in% c("JEB1203", "JEB1303") )

##Filter second set of strains, ones not used
#BFP.data = BFP.data %>% filter(! strain %in% c("JEB1304", "JEB1305", "JEB1306", "JEB1307", "JEB1308") )

# Correlation (uses means only)
cor(BFP.data$growth.rate, BFP.data$GFP.rate)

BFP.deming.fit = deming(GFP.rate~growth.rate, BFP.data, xstd=BFP.data$growth.rate.sd, ystd=BFP.data$GFP.rate.sd,)
BFP.deming.y.intercept = coef(BFP.deming.fit)[1]
BFP.deming.slope = coef(BFP.deming.fit)[2]

# plot showing all strains burden vs growth 
burdenVsgrowthPlot = ggplot(BFP.data, aes(x=growth.rate, y=GFP.rate, shape=isolate, color=other.rate))  +
  geom_errorbarh(aes(xmin=growth.rate.95L, xmax=growth.rate.95U), height=0) +
  geom_errorbar(aes(ymin=GFP.rate.95L, ymax=GFP.rate.95U), width=0) +
  geom_point(size=4)  +
  scale_x_continuous(limits = c(0, 1.5), breaks = (0:7)*0.25 ) + 
  scale_y_continuous(limits = c(0, 8000), breaks = (0:8)*1000) +  
#  geom_abline(intercept=BFP.y.intercept, slope = BFP.slope, linetype=2) +
  geom_abline(intercept=BFP.deming.y.intercept, slope = BFP.deming.slope, linetype=1) +
  scale_color_gradient(low = "lightblue", high = "blue3", na.value = "grey", limits=c(-200,5200), breaks=0:5*1000) +
  labs(color = "BFP Rate") + 
  #  geom_label_repel(aes(label = strain),size = 2.5,box.padding   = 0.35,point.padding = 0.5,segment.color = 'grey50') +
  ggtitle("BFP Series Burden Plot") + xlab("Growth Rate") + ylab("GFP Rate") + 
  theme_bw() + theme(panel.grid.minor = element_blank()) +
  coord_cartesian(expand = F) +
  NULL
burdenVsgrowthPlot
ggsave("11-BFP-series-output/growth_rate_versus_GFP_rate.pdf", width=7, height=6)


# growth BAR PLOT 
BFP.data = BFP.data %>% arrange(desc(growth.rate))
BFP.data$strain = factor(BFP.data$strain, levels=(BFP.data$strain))
growthRatePlot = ggplot(BFP.data, aes(x = strain, y=growth.rate, fill=other.rate)) +  
  geom_bar(size=3, stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=growth.rate.95L, ymax=growth.rate.95U), position=position_dodge()) + 
  scale_y_continuous(limits = c(0, max(BFP.data$growth.rate.95U)), expand=expansion(add = c(0, 0.1))) + 
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 12)) +
  scale_color_gradient(low = "lightblue", high = "blue3", na.value = "grey", limits=c(-50,850), breaks=0:4*200) +
  labs(fill = "BFP Rate") + 
  ggtitle("BFP Series Growth Rates") + xlab("Strain") + ylab("Growth Rate") +
  theme_bw() +
  NULL	
growthRatePlot 
ggsave("11-BFP-series-output/growth_rate_colored_by_BFP_rate.pdf")

# Individual value plots

# Growth rate plot
BFP.data = BFP.data %>% arrange(accession)
  
growthRatePlot = ggplot(BFP.data, aes(x = strain, y=growth.rate)) +  
  geom_bar(size=3, stat="identity", position=position_dodge(), fill="lightblue") +
  geom_errorbar(aes(ymin=growth.rate.95L, ymax=growth.rate.95U), position=position_dodge(), width=0.4) + 
  geom_jitter(data=BFP.all.data, alpha=0.5, size=1, width=0.24) +
  scale_x_discrete(labels=BFP.data$accession) +
  scale_y_continuous(limits = c(0, max(BFP.data$growth.rate.95U, BFP.all.data$growth.rate)), expand=expansion(mult = c(0, 0.05))) + 
  ggtitle("BFP Series Growth Rates") + xlab("Strain") + ylab("Growth Rate") +
  theme_bw() + 
  theme(legend.position = "None", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust=.5, hjust=1, size = 12)) +
  NULL
growthRatePlot 
ggsave("11-BFP-series-output/growth_rate_individual_points.pdf")


# Growth rate plot
growthRatePlot = ggplot(BFP.data, aes(x = strain, y=GFP.rate)) +  
  geom_bar(size=3, stat="identity", position=position_dodge(), fill="lightblue") +
  geom_errorbar(aes(ymin=GFP.rate.95L, ymax=GFP.rate.95U), position=position_dodge(), width=0.4) +
  geom_jitter(data=BFP.all.data, alpha=0.5, size=1, width=0.24) +
  scale_x_discrete(labels=BFP.data$accession) +
  scale_y_continuous(limits = c(0, max(BFP.data$GFP.rate.95U, BFP.all.data$GFP.rate)), expand=expansion(mult = c(0, 0.05))) + 
  ggtitle("BFP Series GFP Rates") + xlab("Strain") + ylab("GFP Rate") +
  theme_bw() + 
  theme(legend.position = "None", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust=.5, hjust=1, size = 12)) +
  NULL
growthRatePlot 
ggsave("11-BFP-series-output/GFP_individual_points.pdf")

# RFP rate plot
growthRatePlot = ggplot(BFP.data, aes(x = strain, y=other.rate)) +  
  geom_bar(size=3, stat="identity", position=position_dodge(), fill="lightblue") +
  geom_errorbar(aes(ymin=other.rate.95L, ymax=other.rate.95U), position=position_dodge(), width=0.4) +
  geom_jitter(data=BFP.all.data, alpha=0.5, size=1, width=0.24) +
  scale_x_discrete(labels=BFP.data$accession) +
  scale_y_continuous(limits = c(0, max(BFP.data$other.rate.95U, BFP.all.data$other.rate)), expand=expansion(mult = c(0, 0.05))) + 
  ggtitle("BFP Series BFP Rates") + xlab("Strain") + ylab("BFP Rate") +
  theme_bw() + 
  theme(legend.position = "None", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust=.5, hjust=1, size = 12)) +
  NULL
growthRatePlot 
ggsave("11-BFP-series-output/other_individual_points.pdf")

