#!/usr/bin/env Rscript

##############################################################
########### Summary of all constructs tested #################
########### *Code Tailored to iGEM 2019 data* ################
##############################################################

library(tidyverse)
library(plotly)
library(ggpubr)
library(cowplot)
library(xtable)
library(forcats)
library(plyr)
library(mgcv)


#tidy data
all.data<-read.csv("burden.sheet.csv")
control.data<-subset(all.data, category=="control")
strain.data<-subset(all.data, category!="control")

#clean data
strain1.data<-subset(strain.data, replicates>1&
                       normalized.GFP.rate.mean>0.5&normalized.GFP.rate.mean<1.5&
                       normalized.growth.rate.mean>0.5&normalized.growth.rate.mean<1.5&
                       normalized.growth.rate.sd<0.2&normalized.GFP.rate.sd<0.2)
badpoints<-subset(strain1.data, strain=="BB246")
clean.data<-anti_join(strain1.data,badpoints,by = "strain")


#Make control strain grand mean lines
grandmeans<- aggregate(control.data[, 14], list(control.data$strain), mean)
grandmeans$Controls<-grandmeans$Group.1
myvector<-as.data.frame(c("BBa_K3174002","BBa_K3174003","BBa_K3174004","BBa_K3174006","BBa_K3174007"))
myvector[2]=grandmeans$x
myvector$Part<-myvector$`c("BBa_K3174002", "BBa_K3174003", "BBa_K3174004", "BBa_K3174006", "BBa_K3174007")`
control.data$strain<- recode(control.data$strain, JEB1204="BBa_K3174002")
control.data$strain<- recode(control.data$strain, JEB1205="BBa_K3174003")
control.data$strain<- recode(control.data$strain, JEB1206="BBa_K3174004")
control.data$strain<- recode(control.data$strain, JEB1207="BBa_K3174006")
control.data$strain<- recode(control.data$strain, JEB1208="BBa_K3174007")

#make duplicates distinct for plotting
clean.data$strain <- sub('[.]', '_', make.names(clean.data$strain, unique=TRUE))
others<-subset(clean.data, clean.data$translational.burden.fraction.p.value<0.05)
clean.data[match(others$strain, clean.data$strain), ] <- others
count(clean.data$category)

#generate list of burdensome parts
sigs<-subset(clean.data, category == "significant")
write.csv(sigs, "burdensome.parts.csv")
noburden<-as.data.frame(clean.data$category=="no burden")

#confidence bands
burdenbro<-lm(normalized.GFP.rate.95CI.range~normalized.growth.rate.95CI.range, clean.data)
confint(burdenbro)

# Plot showing all strains burden vs growth 
burdenVsGrowthPlot = ggplot(clean.data, aes(x=normalized.growth.rate.mean, y=normalized.GFP.rate.mean, color =category))  +
  geom_errorbarh(aes(xmin=normalized.growth.rate.mean-normalized.growth.rate.sd, xmax=normalized.growth.rate.mean+normalized.growth.rate.sd), height=0, alpha = 0.25) +
  geom_errorbar(aes(ymin=normalized.GFP.rate.mean-normalized.GFP.rate.sd, ymax=normalized.GFP.rate.mean+normalized.GFP.rate.sd), width=0, alpha = 0.25) + 
  geom_point(size=1.5)  +
  scale_color_manual(values = c("#68c377","#F43A19", "#618AF9"))+
  scale_x_continuous(limits = c(0, max(clean.data$normalized.growth.rate.mean+clean.data$normalized.growth.rate.sd))) + 
  scale_y_continuous(limits = c(0, max(clean.data$normalized.GFP.rate.mean+clean.data$normalized.GFP.rate.sd))) + 
  geom_abline(intercept=0, slope = 1, alpha = 0.75) +
  geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.5) +
  geom_abline(intercept = -0.05840967, slope = 1, alpha = 0.25, size =1)+
  geom_abline(intercept = 0.08375753, slope = 1, alpha = 0.25, size =1)+
  ylab("Normalized GFP Expression Rate (RFU/hr)")+
  xlab("Normalized Growth Rate (hr^-1)")+
  labs(color = "Category")+
  labs(title = "Burden Measurements for Each Construct")+
  NULL
burdenVsGrowthPlot
ggsave("GFP.v.growth.rates.pdf")
burdenVsGrowthPlotly <- ggplotly(burdenVsGrowthPlot)
htmlwidgets::saveWidget(as_widget(burdenVsGrowthPlotly), "GFP.v.growth.html")

# Plot growth rate for every strain
growthRatePlot = ggplot(clean.data, aes(x= reorder(strain, -normalized.growth.rate.mean), y=normalized.growth.rate.mean, fill = normalized.growth.rate.mean)) + 
  geom_bar(size=3, stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=normalized.growth.rate.mean-normalized.growth.rate.sd, ymax=normalized.growth.rate.mean+normalized.growth.rate.sd), position=position_dodge()) + 
  scale_y_continuous(limits = c(0, max(clean.data$normalized.growth.rate.mean+clean.data$normalized.growth.rate.sd))) +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  xlab("Constructs")+
  ylab("Normalized Growth Rate (hr^-1)")+
  theme(legend.position = "none")+
  labs(title = "Growth Rate Reduction")+
  NULL
growthRatePlot
ggsave("growth.rates.pdf")
growthRatePlotly <- ggplotly(growthRatePlot)
htmlwidgets::saveWidget(as_widget(growthRatePlotly), "growth.rates.html")

#convert burden means to percentages
clean.data$burden.mean<-clean.data$burden.mean*100
myvector$V2<-myvector$V2*100

#plot burden distribution
burdenDensity<- ggdensity(data = clean.data$burden.mean, rug = TRUE, fill = "light gray")+
  geom_vline(xintercept=0, color = "black", linetype = "dashed") +
  geom_vline(data=myvector, aes(xintercept = V2, color = Part))+
  scale_color_manual(values=rainbow(5))+
  scale_x_continuous(limits = c(-20, 45))+
  scale_y_continuous(limits = c(0,0.085))+
  xlab("Burden (% Growth Rate Reduction)")+
  ylab("Relative Density")+
  labs(title = "Burden Value Distribution")+
  theme(legend.position = "right")+
  NULL
burdenDensity
ggsave("burden.density.pdf")
burdenDensityPlotly <- ggplotly(burdenDensity)
htmlwidgets::saveWidget(as_widget(burdenDensityPlotly),"burden.density.html")

