library(adegenet)
library(tidyverse)
library(gridExtra)
library(dplyr)
setwd("C:/Users/meirong/Desktop/PhD project/second preference/final.table/third/1.calculate the result/")
metadata<-read.csv("metadata.final2.csv",row.names = 1)
setwd("C:/Users/meirong/Desktop/PhD project/second preference/final.table/organised code/table/")
fungi<-read.csv("table.nolichenandhost.csv",row.names = 1)
###
table<-as.data.frame(t(fungi))
table2<-table
table2$sample_names<-row.names(table2)
table2<-merge(table2,metadata,by="sample_names")
##
###
fungi.distri <- aggregate(table2[,2:(ncol(table2)-2)],by=list(substrate=table2$substrate,site=table2$site2),sum)
fungi.distri2 <-  fungi.distri %>%
  pivot_longer(cols = where(is.numeric), names_to = "OTU", values_to = "Ratio")%>%
  filter(Ratio>0)
fungi.distri3 <-  fungi.distri2 %>%
  group_by(site,OTU) %>%
  summarise(number=length(substrate))
fungi.distri4 <-  fungi.distri3 %>%
  group_by(site) %>%
  summarise(number2=sum(number<3),all=sum(number>0))
##
table3 <- table2 %>%
  group_by(site2,sample_names) %>%
  mutate(group_sum = sum(across(where(is.numeric)), na.rm = TRUE)) %>%
  mutate(across(where(is.numeric), ~ .x / group_sum)) %>%
  select(-group_sum)
table4 <- aggregate(table3[,2:(ncol(table3)-2)],by=list(substrate=table3$substrate,site=table3$site2),mean)
table3 <- aggregate(table4[,3:(ncol(table4))],by=list(substrate=table4$substrate),mean)

table<-table[rowSums(table)>0,]
metadata <- metadata[match(row.names(table),metadata$sample_names),]
###
fungi<-table
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
table<-table[,c(1,2,9)]
###
table3<-as.data.frame(table3)
row.names(table3)<-table3$substrate
table3<-table3[,-1]
table_long <- table3 %>%
  rownames_to_column(var = "substrate") %>%
  pivot_longer(cols = where(is.numeric), names_to = "OTU", values_to = "Ratio")%>%
  filter(Ratio>0)
##
table42 <- table4 %>%
  pivot_longer(cols = where(is.numeric), names_to = "OTU", values_to = "Ratio")%>%
  filter(Ratio>0)
###use the mean relative abundance of OTUs as the baseline to show the key OTUs changing in relative abundance
row.names(table2)<-table2$sample_names
table2<-table2[,-1]
table_long2 <- table2 %>%
  rownames_to_column(var = "sample_names") %>%
  pivot_longer(cols = where(is.numeric), names_to = "OTU", values_to = "Ratio")%>%
  filter(Ratio>0)
ge<-merge(tax[,c(1,5)],table_long2,by.x="qseqid",by.y="OTU")
try<-ge
try$genus[(try$order=="")]<-"o__unclassified"
aaa<-try
try<- try %>%
  group_by(substrate,order,site2) %>%
  summarise(Ratio=sum(Ratio))
try<- try %>%
  group_by(site2) %>%
  mutate(site.ratio=sum(Ratio))
try$final<-try$Ratio/try$site.ratio # get the most abundant order
