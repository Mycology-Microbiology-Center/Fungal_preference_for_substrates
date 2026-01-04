##compare total richness
library(dplyr)
setwd("C:/Users/meirong/Desktop/PhD project/second preference/final.table/organised code/")
setwd("C:/Users/meirong/Desktop/PhD project/second preference/final.table/third/1.calculate the result/")
metadata<-read.csv("metadata.final2.csv",row.names = 1)
setwd("C:/Users/meirong/Desktop/PhD project/second preference/final.table/organised code/table/")
fungi<-read.csv("table.nolichenandhost.csv",row.names = 1)
fungi<-as.data.frame(t(fungi))
fungi$sample_names<-row.names(fungi)
table<-merge(fungi,metadata,by="sample_names")
table<-aggregate(table[,2:(ncol(table)-2)],by=list(substrate=table$substrate,site=table$site2),sum)
richness2 <- table %>%
  mutate(richness = rowSums(across(where(is.numeric)) > 0)) %>%
  group_by(substrate, site) %>%
  select(substrate, site, richness)
##
library(lmerTest)
library(marginaleffects)
model <- lm(richness ~  substrate + site, data = richness2)
###
b<-avg_comparisons(model, variables = list(substrate = "pairwise")) 
performance::performance(model)
# options(scipen = 999)
substrate.richness<-parameters::model_parameters(model)
library(rcompanion)
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group[letter$Group == "fruitbody"] <- "fruit body"
letter$Group[letter$Group == "leaflitter"] <- "leaf litter"
difference <- richness2 %>%
  group_by(substrate) %>%
  summarise(max=mean(richness))

y.site<-merge(letter,difference,by.x="Group",by.y="substrate")
y.site<-data.frame(substrate=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)
library(ggplot2)
plot_predictions(model,condition = c("substrate"))+
  labs(x="substrate",y="richness")+
  theme_light() +
  scale_x_discrete(limits = c("topsoil","snow","moss","bark"
                              ,"leaves","lichen","leaf litter","subsoil",
                              "fruit body",
                              "wood","feces"))+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom")+
  geom_text(data = y.site, aes(x = substrate , y = ymax, label = letter,hjust=-0.5))
write.csv(substrate.richness,"total.richnessparameter.csv")
