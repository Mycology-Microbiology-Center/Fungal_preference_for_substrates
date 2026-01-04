library(adegenet)
library(tidyverse)
library(gridExtra)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
metadata<-read.csv("metadata.final2.csv",row.names = 1)
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
table42 <- table4 %>%
  pivot_longer(cols = where(is.numeric), names_to = "OTU", values_to = "Ratio")%>%
  filter(Ratio>0)
###
ge<-merge(tax[,c(1,7)],table_long,by.x="qseqid",by.y="OTU")
try<-ge
try$genus[(try$genus=="")]<-"g__unclassified"
aaa<-try
try<- try %>%
  group_by(substrate,genus) %>%
  summarise(Ratio=sum(Ratio))
try2<-try[!try$genus %in% "g__unclassified",]
try2 <- try2 %>%
  group_by(substrate) %>%
  mutate(across(where(is.numeric), ~ rank(desc(.), ties.method = "first")))
try2<-try2[try2$Ratio %in% c(1:10),]
try[!try$genus %in% try2$genus,]$genus<-"other"
try<- try %>%
  group_by(substrate,genus) %>%
  summarise(Ratio=sum(Ratio))
try$genus<-gsub("g__","",try$genus)
try<-try[!try$genus %in% "unclassified",]
heat_mat <- try %>%
  pivot_wider(names_from = substrate, values_from = Ratio, values_fill = 0) %>%
  column_to_rownames("genus") %>%
  as.matrix()
ht1<-Heatmap(heat_mat, name = "Ratio",
        cluster_rows = T,
        cluster_columns = TRUE, 
        col = colorRamp2(c(0,0.005,0.01,0.05, 0.15, 0.4), c("white","pink","lightblue","#FFEEA0","#F4BB44","#C04000")),
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 9))
#######
table_long2<-merge(table_long,table,by.x="OTU",by.y="qseqid",all=T)
table_long2<-table_long2[!is.na(table_long2$Ratio),]
table_long2$GENUS[is.na(table_long2$GENUS)]<-"unknown"
table_long2$primary_lifestyle[is.na(table_long2$primary_lifestyle)]<-"unknown"
try<-table_long2
##
table422<-merge(table42,table,by.x="OTU",by.y="qseqid",all=T)
table422<-table422[!is.na(table422$Ratio),]
table422$GENUS[is.na(table422$GENUS)]<-"unknown"
table422$primary_lifestyle[is.na(table422$primary_lifestyle)]<-"unknown"
consistence <- table422 %>%
  group_by(substrate,primary_lifestyle,site) %>%
  summarise(Ratio=sum(Ratio))
consistence <- consistence %>%
  group_by(substrate,site) %>%
  mutate(all=sum(Ratio))
consistence$percentage<- consistence$Ratio/consistence$all
consistence<-consistence[!consistence$primary_lifestyle %in% "unknown",]
consistence2 <- consistence %>%
  group_by(substrate, site) %>%
  slice_max(order_by = percentage, n = 1, with_ties = FALSE) %>%
  select(substrate, site, rank = percentage, primary_lifestyle) %>%
  ungroup()
##
try<- try %>%
  group_by(substrate,primary_lifestyle) %>%
  summarise(Ratio=sum(Ratio))
try<-try[!try$primary_lifestyle %in% "unknown",]
##
heat_mat <- try %>%
  pivot_wider(names_from = substrate, values_from = Ratio, values_fill = 0) %>%
  column_to_rownames("primary_lifestyle") %>%
  as.matrix()
ht2<-Heatmap(heat_mat, name = "Ratio",
        cluster_rows = T,
        cluster_columns = TRUE, 
        col = colorRamp2(c(0,0.01,0.05, 0.1, 0.3), c("white","lightblue","#FFEEA0","#F4BB44","#C04000")),
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 9))

p1 <- grid.grabExpr(draw(ht1, newpage = FALSE))
p2 <- grid.grabExpr(draw(ht2, newpage = FALSE))

grid.arrange(p1, p2, ncol = 2)

