#Burden statistics

library(tidyverse)
library(plotly)
library(plyr)
require(ggiraph)
require(ggiraphExtra)
library(ggplot2)
library(data.table)

#read in file
sumdat<- read.csv("rates_summary_0.1_merged.csv")

#clean data
cleansumdat<- sumdat[complete.cases(sumdat[ ,6:7]),]%>%mutate(sumdat, CoefVar = growth.rate.sd/mean(sumdat$growth.rate, na.rm = TRUE))
stlim2<-(sd(cleansumdat$CoefVar, na.rm=TRUE))*2
cleansumdat<- subset(cleansumdat, growth.rate.sd<0.2)
cleansumdat$experiment<-sub("exp", "", cleansumdat[,1])


#check noise for each experiment based on variation coefficient
CoefVar0.1plot = ggplot(cleansumdat, aes_(x= reorder(cleansumdat$experiment, -cleansumdat$CoefVar), y=as.name("CoefVar"), fill=as.name("strain")))  +  
  geom_bar(size=3, stat="identity", position = position_dodge()) +
  scale_y_continuous(limits = c(0, max(cleansumdat$CoefVar))) +
  xlab("Strains by Experiment") +
  ylab("Variation Coefficient") +
  theme(legend.position = "none")+
  geom_hline(yintercept = stlim2, linetype = "dashed", alpha = 0.5)+
  geom_hline(yintercept = mean(cleansumdat$CoefVar), linetype = "dashed")+
  NULL
CoefVar0.1plot

ggsave("CoefVar.0.1.plot.pdf")
CoefVar0.1plotly<-ggplotly(CoefVar0.1plot)
htmlwidgets::saveWidget(as_widget(CoefVar0.1plotly), ("CoefVar0.1.plotly.html"))

#read file and subset data
alldat<- read.csv("rates_all_0.1_merged.csv")
controls <- subset(alldat, strain == "JEB1204"| strain == "JEB1205"|
                         strain =="JEB1206"| strain == "JEB1207"| strain =="JEB1208")

#create model
datdat<-data.table(x=controls$growth.rate, y=controls$GFP.rate, grp=controls$experiment)
thisdat<-datdat[,list(intercept=coef(lm(y~x))[1], coef=coef(lm(y~x))[2]),by=grp]

#check regression variation for each experiment's controls
control.lines<- ggplot() +
  geom_point(
    mapping = aes(x = growth.rate, y = GFP.rate, col = experiment, shape = strain),
    data = controls
  ) +
  geom_abline(
    mapping = aes(intercept = intercept, slope = coef, col = grp, alpha = 0.5),
    data = thisdat
  ) +
  scale_x_continuous(limits = c(0, max(controls$growth.rate))) + 
  scale_y_continuous(limits = c(0, max(controls$GFP.rate))) + 
  NULL
control.lines
ggsave("control.lines.plot.pdf")
displotly<-ggplotly(control.lines)
htmlwidgets::saveWidget(as_widget(displotly), ("control.lines.html"))


