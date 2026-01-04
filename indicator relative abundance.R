library(adegenet)
library(tidyverse)
setwd("C:/Users/meirong/Desktop/PhD project/second preference/final.table/third/1.calculate the result/")
metadata<-read.csv("metadata.final2.csv",row.names = 1)
setwd("C:/Users/meirong/Desktop/PhD project/second preference/final.table/organised code/table/")
fungi<-read.csv("table.nolichenandhost.csv",row.names = 1)
table<-as.data.frame(t(fungi))
table2<-table
table2$sample_names<-row.names(table2)
table2<-merge(table2,metadata,by="sample_names")
proportion <- table2 %>%
  mutate(across(2:(ncol(table2)-2), ~ as.numeric(.x > 0)))
# proportion <- aggregate(proportion[,2:(ncol(proportion)-2)],by=list(site=proportion$site2,substrate=proportion$substrate),sum)
# proportion <- proportion %>%
#   mutate(across(3:(ncol(proportion)), ~ as.numeric(.x > 0)))
# tablep<-data.frame(reads=rowSums(proportion[,3:(ncol(proportion))]),substrate=proportion$substrate, site=proportion$site)
table<-data.frame(reads=rowSums(table2[,2:(ncol(table2)-2)]),site=table2$site,substrate=table2$substrate,sample=table2$sample_names)
###
table_long <- table2 %>%
  pivot_longer(
    cols = -c(substrate,sample_names, site2),    # all OTU columns
    names_to = "OTU",              # OTU names go here
    values_to = "reads"            # abundance values
  ) %>%
  group_by(substrate,sample_names,site2,OTU)  %>%
  summarise(reads.sum=sum(reads))  %>%
  filter(reads.sum>0)
###
proportion_long <- proportion %>%
  pivot_longer(
    cols = -c(sample_names,substrate,site2),    # all OTU columns
    names_to = "OTU",              # OTU names go here
    values_to = "reads"            # abundance values
  ) %>%
  # group_by(substrate,site,OTU)  %>%
  # summarise(reads.sum=sum(reads))  %>%
  filter(reads>0)
#######indicator
setwd("C:/Users/meirong/Desktop/PhD project/second preference/final.table/organised code/1.phi")
r <- read.csv("fungi2.r2.tax.csv", row.names = 1)
r <- r[r$p.value < 0.05, ]
###calculateing the relative abundance of all indicators
need<-r
names(need)[2]<-"substrate"
need$substrate<-gsub("\\."," ",need$substrate)
get<-merge(table_long,need[,1:2],by=c("OTU","substrate"))
get2 <- get %>%
  distinct(OTU, substrate, sample_names, site2, .keep_all = TRUE)
####
try<-table_long
proportion2<-proportion_long
names(r)[2]<-"substrate"
r$substrate <-gsub("fruit.body","fruit body",r$substrate)
r$substrate <-gsub("leaf.litter","leaf litter",r$substrate)

all<-merge(r[,c(1:2,4,5)],try, by=c("substrate","OTU"))
pall<-merge(r[,c(1:2,4,5)],proportion2, by=c("substrate","OTU"))
pall2<-pall %>%
  distinct(substrate,site2,sample_names,OTU)
indicator.reads<- all  %>%
  group_by(substrate,sample_names,site2)  %>%
  summarise(reads=sum(reads.sum),stat = first(stat),
            p.value = first(p.value),
            .groups = "drop")

indicator.p<- pall  %>%
  group_by(substrate,site2,sample_names)  %>%
  summarise(reads=sum(reads),
            .groups = "drop")

names(indicator.reads)[4]<-"indicator.reads"
names(table)[4]<-"sample_names"
all<-merge(indicator.reads, table[,c(1,4)],by="sample_names")
result <- all %>%
  group_by(substrate,site2)  %>%
  summarise(indicator.reads2=sum(indicator.reads),reads2=sum(reads))
result$ratio<-result$indicator.reads2/result$reads2
indicator.r<-result
result <- result %>% 
  summarise(mean.v=mean(ratio),sd=sd(ratio))
#
# names(indicator.p)[3]<-"indicator.reads"
# allp<-merge(indicator.p, tablep,by=c("substrate","site"))
# 
# allp$ratio<-allp$indicator.reads/allp$reads
# presult <- allp %>% 
#   group_by(substrate)  %>%
#   summarise(mean.v=mean(ratio),sd=sd(ratio))
##relative abundance of the indicators
indicator.r2<-indicator.r %>%
  group_by(substrate) %>%
  summarise(sd=sd(ratio),mean=mean(ratio))

indicator.r2 <- indicator.r2 %>%
  mutate(substrate = factor(substrate, levels = substrate[order(-mean)]))
###
indicator.r3<-all
indicator.r3$ratio<-indicator.r3$indicator.reads/indicator.r3$reads
library(marginaleffects)
model <- lm(ratio ~ substrate + site2, data = indicator.r3)
b<-avg_comparisons(model, variables = list(substrate = "pairwise")) 
library(rcompanion)
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group[letter$Group == "fruitbody"] <- "fruit body"
letter$Group[letter$Group == "leaflitter"] <- "leaf litter"
difference <- indicator.r %>%
  group_by(substrate) %>%
  summarise(max=mean(ratio))

y.site<-merge(letter,difference,by.x="Group",by.y="substrate")
y.site<-data.frame(substrate=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)
library(ggplot2)
plot_predictions(model,condition = c("substrate"))+
  labs(x="substrate",y="relative abundance")+
  theme_light() +
  scale_x_discrete(limits = as.character(levels(indicator.r2$substrate)))+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom")+
  geom_text(data = y.site, aes(x = substrate , y = ymax, label = letter,hjust=-0.5))
parameter<-parameters::parameters(model)
write.csv(parameter,"indicator.relative.abundance.parameters.csv")
performance::performance(model)
###
# names(allp)[5]<-"indicator.proportion"
# names(indicator.r)[2]<-"site"
try<-merge(indicator.p[,c(3,4)],indicator.r3[,c(1:3,8)],by="sample_names")
names(try)[2]<-"indicator.reads"
setwd("C:/Users/meirong/Desktop/PhD project/second preference/final.table/organised code/table/")
fungi<-read.csv("table.nolichenandhost.csv",row.names = 1)
fungi<-as.data.frame(t(fungi))
r<-hill(fungi, q = 0, raref = 0, runs = 999)
r2<-as.data.frame(r)
r2$sample_names<-row.names(r2)
r2<-merge(r2,metadata,by="sample_names",all=T)
r2<-r2[!is.na(r2$`Hill 0`),]
names(r2)[2]<-'richness'
try<-merge(try, r2[,1:2],by="sample_names")
mod<-lm(ratio ~ richness + substrate + site2,data=try) ##relative abundance and richness
mod3<-lm(indicator.reads ~ richness + substrate + site2,data=try)##indicator number and richness

relative<-parameters::parameters(mod)
performance::performance(mod)
write.csv(relative,"indicator.relative.abundance and richness.parameters.csv")

number<-parameters::parameters(mod3)
performance::performance(mod3)
write.csv(number,"indicator.number.and.richness.parameters.csv")


