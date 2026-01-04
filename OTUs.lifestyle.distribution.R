setwd("C:/Users/meirong/Desktop/PhD project/second preference/final.table/organised code/1.phi")
r <- read.csv("fungi2.r2.tax.csv", row.names = 1)
r2 <- r[r$p.value < 0.05, ]
setwd("C:/Users/meirong/Desktop/PhD project/second preference/final.table/third/1.calculate the result/")
metadata<-read.csv("metadata.final2.csv",row.names = 1)
setwd("C:/Users/meirong/Desktop/PhD project/second preference/final.table/organised code/table/")
fungi<-read.csv("table.nolichenandhost.csv",row.names = 1)
fungi<-as.data.frame(t(fungi))
setwd("C:/Users/meirong/Desktop/PhD project/second preference/final.table/organised code/1.phi")
###
fungi$sample_names<-row.names(fungi)
fungi2<-merge(fungi,metadata,by="sample_names",all=T)
row.names(fungi2)<-fungi2$sample_names
fungi2<-fungi2[,-1]
fungi2<-fungi2[!is.na(fungi2$ff04207f066a2b2769b796c7e8da74227ea12fe4),]
all2<-aggregate(fungi2[,1:(ncol(fungi2)-2)],by=list(substrate=fungi2$substrate,site=fungi2$site2),sum)
library(dplyr)
all <- all2 %>%
  group_by(site) %>%
  mutate(site_total = sum(across(where(is.numeric)))) %>%  # total per site
  mutate(across(where(is.numeric), ~ . / site_total)) %>%
  select(-site_total)  # remove temporary column

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

all.try<-all2_fre
all.try<-all.try[all.try$frequency>0,]
###merging the lifestyle
setwd("C:/Users/meirong/Desktop/PhD project/second preference/final.table/forth/7.living type")
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
table2<-table[,c(1,2,9,10)]
table2$Secondary_lifestyle[is.na(table2$Secondary_lifestyle)]<-"unknown"
table2$Secondary_lifestyle[table2$Secondary_lifestyle==""]<-"unknown"
table2$qseqid<-row.names(table2)
all<-merge(table2,all.try,by.x="qseqid",by.y="OTU",all=T)
all<-all[!is.na(all$frequency),]
all$primary_lifestyle[is.na(all$primary_lifestyle)]<-"unknown"
all$primary_lifestyle[all$primary_lifestyle==""]<-"unknown"
all$Secondary_lifestyle[is.na(all$Secondary_lifestyle)]<-"unknown"
all$Secondary_lifestyle[all$Secondary_lifestyle==""]<-"unknown"
###
all$Secondary_lifestyle<-ifelse(all$Secondary_lifestyle=="unknown",0,1)
all$primary_lifestyle<-ifelse(all$primary_lifestyle=="unknown",0,1)
all$type<-all$primary_lifestyle+all$Secondary_lifestyle
all<-all[all$type>0,]
all$type<-ifelse(all$type==1,"one","second")
library(lmerTest)
library(marginaleffects)
model <- lm(frequency ~ type+site, data = all)
###
b<-avg_comparisons(model, variables = list(type  = "pairwise")) 
performance::performance(model)
