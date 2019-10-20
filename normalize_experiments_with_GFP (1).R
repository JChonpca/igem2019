#Burden statistics

#Find parts that are burdensome from a statistical standpoint.

#Take u-test p-value of all wells for growth rate, apply bonferroni correction, keep p<0.05 wells
#Keep strains with 2/3 significant wells. These are the burdensome strains.

#Calculate growth rate coordinate difference between these strains and their corresponding point
#on the regression line.

library(tidyverse)
library(ggplot2)
library(cowplot)

##############  Read in the input files

all.data = read.csv("rates_all_minOD_0.1_maxmethod_2_merged.csv")
all.data$experiment.strain = paste0(all.data$experiment, "_", all.data$strain)

############## Subset to the experiments that we want to include

metadata = read.csv("igem_2019_experiment_tracker.csv")
include.experiments = (metadata %>% filter(Include==T))$Experiment

all.data = all.data %>% filter(experiment %in% include.experiments)
all.data$experiment = droplevels(all.data$experiment)

##############  Make a general QC graph that highlights the controls
control.data=all.data %>% filter(strain %in% c("JEB1204", "JEB1205", "JEB1206", "JEB1207", "JEB1208"))
not.control.data = all.data %>% filter(!(strain %in% c("JEB1204", "JEB1205", "JEB1206", "JEB1207", "JEB1208")))


ggplot(not.control.data , aes(x=experiment, y=growth.rate, color=experiment)) + 
   geom_jitter( width=0.2, size=0.3, color="gray") + 
   geom_jitter(data=control.data, width=0.2, size=1, aes(x=experiment, y=growth.rate, color=strain)) + 
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
   ylim(0.4,1.6)

ggplot(not.control.data , aes(x=experiment, y=GFP.rate, color=experiment)) + 
   geom_jitter( width=0.2, size=0.3, color="gray") + 
   geom_jitter(data=control.data, width=0.2, size=1, aes(x=experiment, y=GFP.rate, color=strain)) + 
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
   ylim(200,800)

#Show control strains only
#ggplot(control.data, aes(x=experiment, y=growth.rate, color=strain)) + geom_point()


no.burden.quantile = not.control.data %>% group_by(experiment, strain) %>% summarize(growth.rate.mean=mean(growth.rate), GFP.rate.mean=mean(GFP.rate)) %>% group_by(experiment) %>% summarize(growth.rate.quantile=quantile(growth.rate.mean, probs=0.8), GFP.rate.quantile=quantile(GFP.rate.mean, probs=0.8, na.rm=T))

##############  Make a  QC graph that connects the controls

control.data.means = control.data %>% group_by(experiment, strain) %>% summarize(growth.rate.mean=mean(growth.rate), growth.rate.sd = sd(growth.rate), GFP.rate.mean=mean(GFP.rate), GFP.rate.sd = sd(GFP.rate))

ggplot(control.data.means, aes(x=experiment, y=growth.rate.mean, color=strain)) + 
   geom_jitter(data=control.data, width=0.2, size=1, aes(x=experiment, y=growth.rate, color=strain)) +
   geom_line(data=control.data.means, aes(x=experiment, y=growth.rate.mean, color=strain, group=strain)) +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
   geom_errorbar(aes(ymin=growth.rate.mean-growth.rate.sd, ymax=growth.rate.mean+growth.rate.sd), width=0.2) +
   geom_point(data=no.burden.quantile, aes(x=experiment, y=growth.rate.quantile, group="new"), color="black")


ggplot(control.data.means, aes(x=experiment, y=GFP.rate.mean, color=strain)) + 
   geom_jitter(data=control.data, width=0.2, size=1, aes(x=experiment, y=GFP.rate, color=strain)) +
   geom_line(data=control.data.means, aes(x=experiment, y=GFP.rate.mean, color=strain, group=strain)) +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
   geom_errorbar(aes(ymin=GFP.rate.mean-GFP.rate.sd, ymax=GFP.rate.mean+GFP.rate.sd), width=0.2) +
   geom_point(data=no.burden.quantile, aes(x=experiment, y=GFP.rate.quantile, group="new"), color="black")

##############  Determine the no-burden growth rate for each run 
##############  by looking for the maximum density of points

metadata = metadata %>% filter(Experiment %in% include.experiments)
droplevels(metadata$Experiment)
metadata$no.burden.growth.rate = c()
metadata$no.burden.GFP.rate = c()

for (i in 1:nrow(metadata)) {
   this.experiment = as.character(metadata$Experiment[i])
   experiment.no.control.data = not.control.data %>% filter(experiment==this.experiment)
   
   this.density = density(experiment.no.control.data$growth.rate)
   max.y.index = which.max(this.density$y)
   no.burden.growth.rate = this.density$x[max.y.index]
   
   p = ggplot(experiment.no.control.data, aes(x=growth.rate)) + 
      geom_density() + scale_x_continuous(breaks = seq(0, 2, by = 0.1), limits=c(0,2)) + geom_vline(xintercept=no.burden.growth.rate) 

   metadata$no.burden.growth.rate[i] = no.burden.growth.rate
   
   #Remove some that have no GFP measurements
   experiment.no.control.data.has.GFP = experiment.no.control.data %>% filter(!is.na(GFP.rate))
   
   this.density = density(experiment.no.control.data.has.GFP$GFP.rate)
   max.y.index = which.max(this.density$y)
   no.burden.GFP.rate = this.density$x[max.y.index]
   
   ggplot(experiment.no.control.data.has.GFP, aes(x=GFP.rate)) + 
      geom_density() + scale_x_continuous(breaks = seq(200, 800, by = 50), limits=c(200,800)) + geom_vline(xintercept=no.burden.GFP.rate) 
   
   metadata$no.burden.GFP.rate[i] = no.burden.GFP.rate
}


## Normalize to the no.burden.growth.rate
no.burden.growth.rate.data = data.frame(experiment=as.character(metadata$Experiment), no.burden.growth.rate = metadata$no.burden.growth.rate, no.burden.GFP.rate = metadata$no.burden.GFP.rate)

no.burden.normalized.data = all.data %>% left_join(no.burden.growth.rate.data, by="experiment")
no.burden.normalized.data$normalized.growth.rate = no.burden.normalized.data$growth.rate / no.burden.normalized.data$no.burden.growth.rate
no.burden.normalized.data$normalized.GFP.rate = no.burden.normalized.data$GFP.rate / no.burden.normalized.data$no.burden.GFP.rate


no.burden.normalized.control.data = no.burden.normalized.data %>% filter(strain %in% c("JEB1204", "JEB1205", "JEB1206", "JEB1207", "JEB1208"))
no.burden.normalized.not.control.data = no.burden.normalized.data %>% filter(!(strain %in% c("JEB1204", "JEB1205", "JEB1206", "JEB1207", "JEB1208")))

ggplot(no.burden.normalized.not.control.data , aes(x=experiment, y=normalized.growth.rate, color=experiment)) + 
   geom_jitter( width=0.2, size=0.3, color="gray") + 
   geom_jitter(data=no.burden.normalized.control.data, width=0.2, size=1, aes(x=experiment, y=normalized.growth.rate, color=strain)) + 
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
   ylim(0.4,1.2)

ggplot(no.burden.normalized.not.control.data , aes(x=experiment, y=normalized.GFP.rate, color=experiment)) + 
   geom_jitter( width=0.2, size=0.3, color="gray") + 
   geom_jitter(data=no.burden.normalized.control.data, width=0.2, size=1, aes(x=experiment, y=normalized.GFP.rate, color=strain)) + 
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
   ylim(0.2,1.4)



##############  Further normalize values between runs based on controls. 


## determine the grand means for the controls across all experiments

control.means = no.burden.normalized.control.data %>% group_by(strain) %>% summarize(normalized.growth.rate.grand.mean=mean(normalized.growth.rate), normalized.GFP.rate.grand.mean=mean(normalized.GFP.rate))

#now fit the best line going through 1,1 for each run individually

adjustments = data.frame()
for (this.experiment in unique(all.data$experiment)) {
   
   this.normalized.control.data.means = no.burden.normalized.data %>% filter(experiment==this.experiment)
   this.normalized.control.data.means = this.normalized.control.data.means %>% left_join(control.means, by="strain")
   
   to.fit = this.normalized.control.data.means %>% select(normalized.growth.rate, normalized.growth.rate.grand.mean, normalized.GFP.rate, normalized.GFP.rate.grand.mean)
   to.fit$normalized.growth.rate = 1-to.fit$normalized.growth.rate
   to.fit$normalized.growth.rate.grand.mean = 1-to.fit$normalized.growth.rate.grand.mean
   
   to.fit$normalized.GFP.rate = 1-to.fit$normalized.GFP.rate
   to.fit$normalized.GFP.rate.grand.mean = 1-to.fit$normalized.GFP.rate.grand.mean
   
   the.growth.rate.fit = lm(normalized.growth.rate.grand.mean~normalized.growth.rate+0, data=to.fit)
   the.GFP.rate.fit = lm(normalized.GFP.rate.grand.mean~normalized.GFP.rate+0, data=to.fit)
   
   adjustments = bind_rows(adjustments, data.frame(experiment=this.experiment, growth.rate.normalization.factor = coef(the.growth.rate.fit)[1], GFP.rate.normalization.factor = coef(the.GFP.rate.fit)[1]))
}

final.normalized.data = no.burden.normalized.data

final.normalized.data = final.normalized.data %>% left_join(adjustments, by="experiment")

final.normalized.data$normalized.growth.rate = 1-(1-final.normalized.data$normalized.growth.rate)*final.normalized.data$growth.rate.normalization.factor
final.normalized.data$normalized.GFP.rate = 1-(1-final.normalized.data$normalized.GFP.rate)*final.normalized.data$GFP.rate.normalization.factor


final.normalized.control.data=final.normalized.data %>% filter(strain %in% c("JEB1204", "JEB1205", "JEB1206", "JEB1207", "JEB1208"))
final.normalized.not.control.data = final.normalized.data %>% filter(!(strain %in% c("JEB1204", "JEB1205", "JEB1206", "JEB1207", "JEB1208")))

ggplot(final.normalized.not.control.data , aes(x=experiment, y=normalized.growth.rate, color=experiment)) + 
   geom_jitter( width=0.2, size=0.3, color="gray") + 
   geom_jitter(data=final.normalized.control.data, width=0.2, size=1, aes(x=experiment, y=normalized.growth.rate, color=strain)) + 
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
   ylim(0.4,1.2)

ggplot(final.normalized.not.control.data , aes(x=experiment, y=normalized.GFP.rate, color=experiment)) + 
   geom_jitter( width=0.2, size=0.3, color="gray") + 
   geom_jitter(data=final.normalized.control.data, width=0.2, size=1, aes(x=experiment, y=normalized.GFP.rate, color=strain)) + 
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
   ylim(0.4,1.2)

##################################

#Now judge which are significantly different from one

parts = final.normalized.not.control.data
p.values = data.frame()
for (this.experiment.strain in unique(parts$experiment.strain)) {

   this.strain.experiment.data = parts %>% filter(experiment.strain == this.experiment.strain)
   
   cat(this.experiment.strain, "  ", length(this.strain.experiment.data$normalized.growth.rate), "\n")
   
   if (length(this.strain.experiment.data$normalized.growth.rate) > 1) {
      
      test.result = t.test(this.strain.experiment.data$normalized.growth.rate-1, alternative = "less")
      cat(this.experiment.strain, "  ", test.result$p.value, "\n")
   } else {
      test.result = data.frame(p.value=NA)
   }
   p.values = bind_rows(p.values, data.frame(experiment.strain=this.experiment.strain, p.value = test.result$p.value ))
}
p.values$p.value = p.adjust(p.values$p.value, method="BH")

## GFP

parts = final.normalized.not.control.data
translational.burden.fraction.p.values = data.frame()
for (this.experiment.strain in unique(parts$experiment.strain)) {
   
   this.strain.experiment.data = parts %>% filter(experiment.strain == this.experiment.strain) %>% filter(!is.na(normalized.GFP.rate))
   
   cat(this.experiment.strain, "  ", length(this.strain.experiment.data$normalized.GFP.rate), "\n")
   
   if ( (length(this.strain.experiment.data$normalized.GFP.rate) > 1)) {
      
      test.result = t.test(this.strain.experiment.data$normalized.GFP.rate, this.strain.experiment.data$normalized.growth.rate, alternative = "greater")
      cat(this.experiment.strain, "  ", test.result$p.value, "\n")
   } else {
      test.result = data.frame(p.value=NA)
   }
   translational.burden.fraction.p.values = bind_rows(translational.burden.fraction.p.values, data.frame(experiment.strain=this.experiment.strain, translational.burden.fraction.p.value = test.result$p.value ))
}
translational.burden.fraction.p.values$translational.burden.fraction.p.value = p.adjust(translational.burden.fraction.p.values$translational.burden.fraction.p.value, method="BH")


parts = final.normalized.data


parts.means = parts %>% group_by(experiment, strain, experiment.strain) %>% summarize(replicates=n(), normalized.growth.rate.mean=mean(normalized.growth.rate), normalized.growth.rate.sd=sd(normalized.growth.rate), normalized.growth.rate.sem = normalized.growth.rate.sd/sqrt(replicates), normalized.growth.rate.95CI.range =  normalized.growth.rate.sem*qt(0.975, df=replicates-1), normalized.GFP.rate.mean=mean(normalized.GFP.rate), normalized.GFP.rate.sd=sd(normalized.GFP.rate), normalized.GFP.rate.sem = normalized.GFP.rate.sd/sqrt(replicates), normalized.GFP.rate.95CI.range =  normalized.GFP.rate.sem*qt(0.975, df=replicates-1))

parts.means = parts.means %>%left_join(p.values, by="experiment.strain")
parts.means = parts.means %>%left_join(translational.burden.fraction.p.values, by="experiment.strain")

parts.means$category = "no burden";
parts.means$category[parts.means$p.value < 0.05] = "significant"
parts.means$category[parts.means$strain %in% c("JEB1204", "JEB1205", "JEB1206", "JEB1207", "JEB1208")] = "control"

parts.means = parts.means %>% group_by() %>% select(-experiment.strain)

parts.means$burden.mean = 1-parts.means$normalized.growth.rate.mean
parts.means$burden.sd = parts.means$normalized.growth.rate.sd
parts.means$burden.sem = parts.means$normalized.growth.rate.sem
parts.means$burden.95CI.range = parts.means$normalized.growth.rate.95CI.range
parts.means$burden.cv = parts.means$burden.sd  / parts.means$burden.mean


parts.means$translational.burden.fraction.mean =  (1 - parts.means$normalized.GFP.rate.mean) / (1 - parts.means$normalized.growth.rate.mean)
parts.means$translational.burden.fraction.sd = parts.means$normalized.GFP.rate.sd / (1 - parts.means$normalized.growth.rate.mean)
parts.means$translational.burden.fraction.sem = parts.means$normalized.growth.rate.sem / (1 - parts.means$normalized.growth.rate.mean)
parts.means$translational.burden.fraction.95CI.range = parts.means$normalized.growth.rate.95CI.range / (1 - parts.means$normalized.growth.rate.mean)

#parts.means = parts.means %>% select(-growth.rate.mean,-growth.rate.sd)

write.csv(parts.means, "part_burden.csv", row.names=F)

ggplot(parts.means, aes(x=experiment, y=normalized.growth.rate.mean, color=category)) + geom_jitter() +  geom_errorbar(aes(ymin=normalized.growth.rate.mean-normalized.growth.rate.sd, ymax=normalized.growth.rate.mean+normalized.growth.rate.sd), width=.2,) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


### Makes the overall distribution graph

control.burdens = control.means
control.burdens$burden.mean = 1 - control.burdens$normalized.growth.rate.grand.mean
control.burdens = control.burdens %>% select(-normalized.growth.rate.grand.mean, -normalized.GFP.rate.grand.mean)
control.burdens = rbind(data.frame(strain="no burden", burden.mean=0), control.burdens)

ggplot(parts.means %>% filter(category!="control"), aes(x=burden.mean)) + geom_density() + xlim(-0.2, 0.6) + geom_vline(data = control.burdens, aes(xintercept=burden.mean, color=strain)) 

