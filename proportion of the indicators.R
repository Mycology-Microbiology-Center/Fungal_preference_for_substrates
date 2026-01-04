setwd("C:/Users/meirong/Desktop/PhD project/second preference/final.table/third/1.calculate the result/")
metadata<-read.csv("metadata.final2.csv",row.names = 1)
setwd("C:/Users/meirong/Desktop/PhD project/second preference/final.table/organised code/table/")
fungi<-read.csv("table.nolichenandhost.csv",row.names = 1)
fungi<-as.data.frame(t(fungi))
setwd("C:/Users/meirong/Desktop/PhD project/second preference/final.table/organised code/1.phi")
r <- read.csv("all.r2.tax.csv", row.names = 1)
r2 <- r[r$p.value < 0.05, ]
###
fungi$sample_names<-row.names(fungi)
fungi2<-merge(fungi,metadata,by="sample_names",all=T)
row.names(fungi2)<-fungi2$sample_names
fungi2<-fungi2[,-1]
fungi2<-fungi2[!is.na(fungi2$ff04207f066a2b2769b796c7e8da74227ea12fe4),]
all2<-aggregate(fungi2[,1:(ncol(fungi2)-2)],by=list(substrate=fungi2$substrate),sum)
row.names(all2)<-all2$substrate
all2<-all2[,-1]
richness<-data.frame(richness=rowSums(all2>0),BestHost=row.names(all2))
indicator<-r2 %>%
  group_by(BestHost) %>%
  summarise(number=length(OTU))
indicator$BestHost<-gsub("\\.", " ", indicator$BestHost)
all<-merge(richness,indicator,by="BestHost")
all$ratio<-all$number/all$richness
