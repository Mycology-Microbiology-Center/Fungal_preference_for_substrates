library(RColorBrewer)
library(colorspace)
library(emmeans)
library(MuMIn)
library(marginaleffects)
library(lmerTest)
library(rcompanion)
library(dplyr)
metadata<-read.csv("metadata.final2.csv",row.names = 1)
fungi<-read.csv("table.nolichenandhost.csv",row.names = 1)
fungi<-as.data.frame(t(fungi))
r <- read.csv("all.r2.tax.csv", row.names = 1)
r2 <- r[r$p.value < 0.05, ]
###
fungi$sample_names<-row.names(fungi)
fungi2<-merge(fungi,metadata,by="sample_names",all=T)
row.names(fungi2)<-fungi2$sample_names
fungi2<-fungi2[,-1]
fungi2<-fungi2[!is.na(fungi2$ff04207f066a2b2769b796c7e8da74227ea12fe4),]
all2<-aggregate(fungi2[,1:(ncol(fungi2)-2)],by=list(substrate=fungi2$substrate,site=fungi2$site2),sum)
library(dplyr)
all <- all2 %>%
  group_by(site,substrate) %>%
  mutate(site_total = sum(across(where(is.numeric)))) %>%  # total per site
  mutate(across(where(is.numeric), ~ . / site_total)) %>%
  select(-site_total)
##calculate the frequency
fre <- all2 %>%
  group_by(site) %>%
  summarise(across(where(is.numeric), ~ sum(. > 0)))
##
library(tidyverse)
all2_long <- all %>%
  pivot_longer(cols = where(is.numeric), names_to = "OTU", values_to = "Ratio")

all2_fre <- fre %>%
  pivot_longer(cols = where(is.numeric), names_to = "OTU", values_to = "frequency")

all.try<-merge(all2_long,all2_fre,by=c("site","OTU"),all=T)
all.try<-all.try[all.try$Ratio>0,]
all.try<-all.try[all.try$frequency>0,]
###merging the lifestyle
trait<-read.csv("FungalTraits.csv")
tax<-read.csv("taxonomy.final.csv")
tax$GENUS<-gsub(".*__","",tax$genus)
table<-merge(tax[,c(1,10)],trait,by="GENUS",all=T)
table<-table[!is.na(table$qseqid),]
table<-table[!is.na(table$jrk_template),]
table<-table[!is.na(table$qseqid),]
table<-table[table$GENUS!="",]
table$primary_lifestyle[is.na(table$primary_lifestyle)]<-"unknown"
table$primary_lifestyle[table$primary_lifestyle==""]<-"unknown"
row.names(table)<-table$qseqid
###
table<-table[,c(1,2,9)]
###
all<-merge(table,all.try,by.x="qseqid",by.y="OTU",all=T)
all<-all[!is.na(all$frequency),]
all$primary_lifestyle[is.na(all$primary_lifestyle)]<-"unknown"
all$primary_lifestyle[all$primary_lifestyle==""]<-"unknown"
##
r[r$BestHost %in% "fruit.body",]$BestHost <- "fruit body"
r[r$BestHost %in% "leaf.litter",]$BestHost <- "leaf litter"
try<-all[all$qseqid %in% unique(r$OTU), ]
names(try)[c(1,5)]<-c("OTU","BestHost")
try<-merge(try,r[,c(1:2,4,19)],by=c("OTU","BestHost"),all=T)
###
aaa<-try
try[is.na(try$sign),]$sign<-FALSE
all2 <- try %>%
  arrange(primary_lifestyle, OTU,site) %>%
  mutate(qseqid_ordered = factor(OTU, levels = unique(OTU)))
all2 <- all2 %>%
  mutate(color_group = ifelse(sign, primary_lifestyle, "non-indicator"))

colors <- c(
  colorRampPalette(brewer.pal(9, "Set1"))(9),
  colorRampPalette(brewer.pal(12, "Set3"))(12),
  colorRampPalette(brewer.pal(8, "Set2"))(8))
names(colors)<-c( unique(all2$site.x))
all3<-all2[!is.na(all2$stat),]
####
all2[!all2$color_group %in% "non-indicator",]$color_group <- 'indicator'
all2[all2$color_group %in% "non-indicator",]$color_group <- 'non_indicator'
model <- lmer(frequency ~ BestHost + (1|site), data = all2) ##compare the substrate range of fungal species at particular substrate
b<-avg_comparisons(model, variables = list(BestHost = "pairwise")) 
performance::performance(model)
fre.substrate<-parameters::model_parameters(model)
write.csv(fre.substrate,"fre.substrate.parameters.csv")
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group[letter$Group == "fruitbody"] <- "fruit body"
letter$Group[letter$Group == "leaflitter"] <- "leaf litter"
difference <- all2 %>%
  group_by(BestHost) %>%
  summarise(max=mean(frequency))
y.site<-merge(letter,difference,by.x="Group",by.y="BestHost")
y.site<-data.frame(BestHost=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)
plot_predictions(model,condition = c("BestHost"))+
  labs(x="substrate",y="frequency")+
  theme_light() +
scale_x_discrete(limits = c("leaves","moss","snow","leaf litter",
                            "subsoil","topsoil","bark","fruit body",
                            "wood","feces","lichen"))+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom")+
  geom_text(data = y.site, aes(x = BestHost , y = ymax, label = letter,hjust=-0.5))
######only indicator
a<-all2[!all2$color_group %in% 'non_indicator', ]
###substrate range of guild 
substrate_range<- a %>%
  group_by(BestHost, site) %>%
  summarise(min.s=min(frequency), max.s=max(frequency))
write.csv(substrate_range,"substrate_range.csv")
###
model <- lm(frequency ~ BestHost  +  site, data = a)##compare the substrate range of indicator species at particular substrate
performance::performance(model)
fre.substrate<-parameters::model_parameters(model)
write.csv(fre.substrate,"fre.indi.substrate.parameters.csv")
b<-avg_comparisons(model, variables = list(BestHost  = "pairwise")) 
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group[letter$Group == "fruitbody"] <- "fruit body"
letter$Group[letter$Group == "leaflitter"] <- "leaf litter"
difference <- a %>%
  group_by(BestHost ) %>%
  summarise(max=mean(frequency))

y.site<-merge(letter,difference,by.x="Group",by.y="BestHost")
y.site<-data.frame(BestHost =factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)

plot_predictions(model,condition = c("BestHost"))+
  labs(x="BestHost",y="frequency")+
  scale_y_continuous(breaks = c(1,2, 2.5, 3))+
  scale_x_discrete(limits = c("leaves","moss","snow","leaf litter",
                              "subsoil","topsoil","bark","fruit body",
                              "wood","feces","lichen"))+
  scale_y_continuous(breaks = c(1,2,2.5,3), limits = c(1, 3))+
theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom")+
  geom_text(data = y.site, aes(x = BestHost  , y = ymax, label = letter,hjust=-0.5))
###indicator compare
data<-all3
data$sign<-ifelse(data$sign,'indicator','non-indicator')
data2<-as.data.frame(data)
model <- lmer(frequency ~ sign + BestHost + (1|site), data = all2)
b<-avg_comparisons(model, variables = list(sign = "pairwise"))##substrate range comparison between indicator and non-indicators
##substrate range comparison
result<-data.frame()
for (i in unique(all2$BestHost)) {
  indicator.min<-min(all2[all2$sign & all2$BestHost %in% i,]$frequency)
  indicator.max<-max(all2[all2$sign & all2$BestHost %in% i,]$frequency)
  indicator.mean<-mean(all2[all2$sign & all2$BestHost %in% i,]$frequency)
  
  nonindicator.min<-min(all2[!all2$sign  & all2$BestHost %in% i,]$frequency)
  nonindicator.max<-max(all2[!all2$sign & all2$BestHost %in% i,]$frequency)
  nonindicator.mean<-mean(all2[!all2$sign & all2$BestHost %in% i,]$frequency)
  
  c<-data.frame(indicator.min=indicator.min,indicator.max=indicator.max,
                nonindicator.min=nonindicator.min,nonindicator.max=nonindicator.max,
                nonindicator.meann=nonindicator.mean,indicator.mean=indicator.mean)
  
  result<-rbind(result,c)
}
###
indicator.rs<-aaa[!is.na(aaa$sign),]
indicator.rs<-indicator.rs[(indicator.rs$sign),]
model <- lm(stat ~ Ratio + BestHost + site, data = indicator.rs)
fre.ratio<-parameters::model_parameters(model)
performance::performance(model)
write.csv(fre.ratio,"stat.abundance.parameter.csv")
####compare the phi
aaaa.p <- r %>%
  group_by(BestHost) %>%
  summarise(fre=mean(stat),sd=sd(stat))

###compare the association strength among substrates
indicator.all<-r[r$p.value<0.05,]
model <- lm(stat ~ BestHost, data = indicator.all)
b<-avg_comparisons(model, variables = list(BestHost = "pairwise")) 
fre.ratio<-parameters::model_parameters(model)
performance::performance(model)
write.csv(fre.ratio,"stat.substrate.parameter.csv")
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group[letter$Group == "fruitbody"] <- "fruit body"
letter$Group[letter$Group == "leaflitter"] <- "leaf litter"
difference <- indicator.all %>%
  group_by(BestHost) %>%
  summarise(max=mean(stat))
y.site<-merge(letter,difference,by.x="Group",by.y="BestHost")
y.site<-data.frame(substrate=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)
plot_predictions(model,condition = c("BestHost"))+
  labs(x="substrate",y="preference strength")+
  scale_x_discrete(limits = c("leaves","moss","snow","leaf litter",
                              "subsoil","topsoil","bark","fruit body",
                              "wood","feces","lichen"))+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom")+
  geom_text(data = y.site, aes(x = substrate , y = ymax, label = letter,hjust=-0.5))



