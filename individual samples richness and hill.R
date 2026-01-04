library(BAT)
setwd("C:/Users/meirong/Desktop/PhD project/second preference/final.table/third/1.calculate the result/")
metadata<-read.csv("metadata.final2.csv",row.names = 1)
setwd("C:/Users/meirong/Desktop/PhD project/second preference/final.table/organised code/table/")
fungi<-read.csv("table.nolichenandhost.csv",row.names = 1)
fungi<-as.data.frame(t(fungi))
r<-hill(fungi, q = 1, raref = 0, runs = 999)
r2<-as.data.frame(r)
r2$sample_names<-row.names(r2)
r2<-merge(r2,metadata,by="sample_names",all=T)
r2<-r2[!is.na(r2$`Hill 1`),]
names(r2)[2]<-'hill.1'
##
##
library(lmerTest)
library(marginaleffects)
model <- lm(hill.1 ~  substrate + site2, data = r2)
###
b<-avg_comparisons(model, variables = list(substrate = "pairwise")) 
performance::performance(model)
substrate.richness<-parameters::model_parameters(model)
library(rcompanion)
library(dplyr)
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group[letter$Group == "fruitbody"] <- "fruit body"
letter$Group[letter$Group == "leaflitter"] <- "leaf litter"
difference <- r2 %>%
  group_by(substrate) %>%
  summarise(max=mean(hill.1))

y.site<-merge(letter,difference,by.x="Group",by.y="substrate")
y.site<-data.frame(substrate=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)
library(ggplot2)

plot_predictions(model,condition = c("substrate"))+
  labs(x="substrate",y="hill number")+
  theme_light() +
  scale_x_discrete(limits = c("snow","leaves","topsoil","moss","bark"
                              ,"lichen", "subsoil","leaf litter",
                              "fruit body",
                              "wood","feces"))+
  scale_y_continuous(breaks = c(0,50))+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom")+
  geom_text(data = y.site, aes(x = substrate , y = ymax, label = letter,hjust=-0.5))

write.csv(substrate.richness,"q1.richnessparameter.csv")
###
r<-hill(fungi, q = 0, raref = 0, runs = 999)
r2<-as.data.frame(r)
r2$sample_names<-row.names(r2)
r2<-merge(r2,metadata,by="sample_names",all=T)
r2<-r2[!is.na(r2$`Hill 0`),]
names(r2)[2]<-'hill.0'
write.csv(r2,"richness.csv")
##
model <- lm(hill.0 ~  substrate + site2, data = r2)
###
b<-avg_comparisons(model, variables = list(substrate = "pairwise")) 
performance::performance(model)
# options(scipen = 999)
substrate.richness<-parameters::model_parameters(model)
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group[letter$Group == "fruitbody"] <- "fruit body"
letter$Group[letter$Group == "leaflitter"] <- "leaf litter"
difference <- r2 %>%
  group_by(substrate) %>%
  summarise(max=mean(hill.0))
y.site<-merge(letter,difference,by.x="Group",by.y="substrate")
y.site<-data.frame(substrate=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)
library(ggplot2)
plot_predictions(model,condition = c("substrate"))+
  labs(x="substrate",y="richness(OTU)")+
  theme_light() +
  scale_x_discrete(limits = c("snow","leaves","topsoil","moss","bark"
                              ,"lichen", "subsoil","leaf litter",
                              "fruit body",
                              "wood","feces"))+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom")+
  geom_text(data = y.site, aes(x = substrate , y = ymax, label = letter,hjust=-0.5))
write.csv(substrate.richness,"q0.richnessparameter.csv")
