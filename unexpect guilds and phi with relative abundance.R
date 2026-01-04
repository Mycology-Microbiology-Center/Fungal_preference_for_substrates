library(dplyr)
library(tidyverse)
setwd("C:/Users/meirong/Desktop/PhD project/second preference/final.table/third/1.calculate the result/")
metadata<-read.csv("metadata.final2.csv",row.names = 1)
setwd("C:/Users/meirong/Desktop/PhD project/second preference/final.table/organised code/table")
fungi<-read.csv("table.nolichenandhost.csv",row.names = 1)
table<-as.data.frame(t(fungi))
table2<-table
table2$sample_names<-row.names(table2)
table2<-merge(table2,metadata,by="sample_names")
table3 <- table2 %>%
  group_by(site2,sample_names) %>%
  mutate(group_sum = sum(across(where(is.numeric)), na.rm = TRUE)) %>%
  mutate(across(where(is.numeric), ~ .x / group_sum)) %>%
  select(-group_sum)
table4 <- aggregate(table3[,2:(ncol(table3)-2)],by=list(substrate=table3$substrate,site=table3$site2),mean)
###mean relative abundance
table_long <- table4 %>%
  pivot_longer(
    cols = -c(substrate, site),    
    names_to = "OTU",             
    values_to = "reads"            
  ) %>%
  group_by(substrate,site,OTU)  %>%
  summarise(relative.sum=sum(reads))
###reads
reads<-aggregate(table2[,2:(ncol(table2)-2)],by=list(substrate=table2$substrate,site=table2$site2),mean)
read_long <- reads %>%
  pivot_longer(
    cols = -c(substrate, site),    
    names_to = "OTU",             
    values_to = "reads"            
  ) %>%
  group_by(substrate,site,OTU)  %>%
  summarise(reads.sum=sum(reads))

read_long <- read_long %>%
  group_by(substrate,OTU)  %>%
  summarise(reads = mean(reads.sum))
###
table_long <- table_long %>%
  group_by(substrate,OTU)  %>%
  summarise(relative = mean(relative.sum))
###
setwd("C:/Users/meirong/Desktop/PhD project/second preference/final.table/organised code/1.phi")
r <- read.csv("all.r2.tax.csv", row.names = 1)
r2 <- r[r$p.value < 0.05, ]
names(r2)[2]<-"substrate"
r2$substrate<-gsub("\\."," ",r2$substrate)
read_long2<-merge(r2[,c(1,2,4)],read_long,by=c("substrate","OTU"))
table_long2<-merge(r2[,c(1,2,4)],table_long,by=c("substrate","OTU"))
###
model <- lm(stat ~ relative + substrate, data = table_long2)
performance::performance(model)
parameters::parameters(model)
####indicators belongs to unexpected guilds
setwd("C:/Users/meirong/Desktop/PhD project/second preference/final.table/forth/7.living type")
trait<-read.csv("FungalTraits.csv")
tax<-read.csv("taxonomy.final.csv")
tax$GENUS<-gsub(".*__","",tax$genus)
r<-r[,c(1,2,5,25)]
r2 <- r[r$p.value < 0.05, ]
names(r2)[4]<-"GENUS"
table<-merge(r2,trait[,c(6,8,9)],by="GENUS",all=T)
table<-table[!is.na(table$BestHost),]
table[is.na(table$primary_lifestyle),]$primary_lifestyle<-"unknown"
table[is.na(table$Secondary_lifestyle ),]$Secondary_lifestyle<-"unknown"
setwd("C:/Users/meirong/Desktop/PhD project/second preference/final.table/organised code/6.lifestyle distribution/")
write.csv(table,"unexpected.guilds2.csv")
