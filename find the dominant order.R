library(adegenet)
library(tidyverse)
library(gridExtra)
library(dplyr)
metadata<-read.csv("metadata.final2.csv",row.names = 1)
fungi<-read.csv("table.nolichenandhost.csv",row.names = 1)
table<-as.data.frame(t(fungi))
table2<-table
table2$sample_names<-row.names(table2)
table2<-merge(table2,metadata,by="sample_names")
metadata <- metadata[match(row.names(table),metadata$sample_names),]
###
tax<-read.csv("taxonomy.final.csv")
###
table3<-as.data.frame(table2)

table3<-table3[,-1]
table3<- aggregate(table3[,1:(ncol(table3)-2) ],by=list(substrate=table3$substrate,site=table3$site2),mean)
  
table_long <- table3 %>%
  pivot_longer(cols = where(is.numeric), names_to = "OTU", values_to = "Ratio")%>%
  filter(Ratio>0)

###use the mean relative abundance of OTUs as the baseline to show the key OTUs changing in relative abundance
ge<-merge(tax[,c(1,5)],table_long,by.x="qseqid",by.y="OTU")
try<-ge
try<- try %>%
  group_by(substrate,order,site) %>%
  summarise(Ratio=sum(Ratio))

