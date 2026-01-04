library(picante)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ranger)
library(shapviz)
library(rcompanion)
library(marginaleffects)
library(lme4)
###
metadata<-read.csv("metadata.final2.csv",row.names = 1)
fungi<-read.csv("table.nolichenandhost.csv",row.names = 1)
fungi<-as.data.frame(t(fungi))
tax<-read.csv("tax.final.csv",header = T,row.names = 1)
tax<-tax[tax$qseqid %in% names(fungi),]
write.csv(tax[,2:6],"fungi.csv")
fungi$sample_names<-row.names(fungi)
###
fungi.tree<-read.tree("fungi.tree.txt")
###
rename.tr<-function(tree){
  edge_matrix <- tree$edge
  tip_labels <- tree$tip.label
  node_labels <- tree$node.label
  num_tips <- length(tip_labels)
  all_labels <- c(tip_labels, node_labels)
  parent_child_df <- data.frame(
    Parent_Node = all_labels[edge_matrix[,1]],  # Convert parent node number to name
    Child_Node  = all_labels[edge_matrix[,2]],   # Convert child node number to name
    Parent_Node1 = edge_matrix[,1],  # Convert parent node number to name
    Child_Node1  = edge_matrix[,2]
  )
  
  family1<-parent_child_df[grepl("f__",parent_child_df$Child_Node),]
  names(family1)[4]<-"Child_Nodefamily"
  names(family1)[3]<-"Parent_Nodefamily"
  order<-parent_child_df[grepl("o__",parent_child_df$Child_Node),]
  try<-merge(order,family1,by.x="Child_Node1",by.y="Parent_Nodefamily",all=T)
  names(try)[4]<-"Parent_Nodeorder"
  class<-parent_child_df[grepl("c__",parent_child_df$Child_Node),]
  try<-merge(class,try,by.x="Child_Node1",by.y="Parent_Nodeorder",all=T)
  try<-try[,-1]
  names(try)[3]<-"Parent_Nodeclass"
  phylum<-parent_child_df[grepl("p__",parent_child_df$Child_Node),]
  try<-merge(phylum,try,by.x="Child_Node1",by.y="Parent_Nodeclass",all=T)
  ##change the label names as table names
  try<-try[,c(2,3,6,9,11)]
  tree$tip.label<-paste0(try$Parent_Node.x,try$Child_Node.x,try$Child_Node.y,try$Child_Node.x.1,
                         try$Child_Node.y.1)
  return(tree)
}
fungi.tree<-rename.tr(fungi.tree)
fungi.distance<-cophenetic(fungi.tree)
##table
tax<-read.csv("tax.final.csv",header = T,row.names = 1)
table<-fungi[,-ncol(fungi)]
table<-as.data.frame(t(table))
table$OTU<-row.names(table)
###indicator
r <- read.csv("fungi2.r2.tax.csv", row.names = 1)
r <- r[r$p.value < 0.05, ]
indicator.table<-table[row.names(table) %in% r$OTU,]
###
table<-merge(table,tax,by.x="OTU",by.y = "qseqid")
table$agg<-paste0(table$Kingdom,table$Phylum,table$Class,table$Order,table$Family)
table<-aggregate(table[,2:(ncol(table)-8)],by=list(table$agg),sum)
table<-table[match(fungi.tree$tip.label,table$Group.1),]
row.names(table)<-table$Group.1
table<-table[,-1]
##LZ
LV<-table[, names(table) %in% metadata[metadata$site2 %in% "LV",]$sample_names]
LV<-LV[,colSums(LV)>0]
LV<-LV[rowSums(LV)>0,]
LV.distance<-fungi.distance[row.names(LV),row.names(LV)]

LW<-table[, names(table) %in% metadata[metadata$site2 %in% "LW",]$sample_names]
LW<-LW[,colSums(LW)>0]
LW<-LW[rowSums(LW)>0,]
LW.distance<-fungi.distance[row.names(LW),row.names(LW)]

LZ<-table[, names(table) %in% metadata[metadata$site2 %in% "LZ",]$sample_names]
LZ<-LZ[,colSums(LZ)>0]
LZ<-LZ[rowSums(LZ)>0,]
LZ.distance<-fungi.distance[row.names(LZ),row.names(LZ)]
###indicator
indicator.table<-merge(indicator.table,tax,by.x="OTU",by.y = "qseqid")
indicator.table$agg<-paste0(indicator.table$Kingdom,indicator.table$Phylum,indicator.table$Class,indicator.table$Order,indicator.table$Family)
indicator.table<-aggregate(indicator.table[,2:(ncol(indicator.table)-8)],by=list(indicator.table$agg),sum)
row.names(indicator.table)<-indicator.table$Group.1
indicator.table<-indicator.table[,-1]
##LZ
LV.indicator<-indicator.table[, names(indicator.table) %in% metadata[metadata$site2 %in% "LV",]$sample_names]
LV.indicator<-LV.indicator[,colSums(LV.indicator)>0]
LV.indicator<-LV.indicator[rowSums(LV.indicator)>0,]
LV.indicator.distance<-fungi.distance[row.names(LV.indicator),row.names(LV.indicator)]

LW.indicator<-indicator.table[, names(indicator.table) %in% metadata[metadata$site2 %in% "LW",]$sample_names]
LW.indicator<-LW.indicator[,colSums(LW.indicator)>0]
LW.indicator<-LW.indicator[rowSums(LW.indicator)>0,]
LW.indicator.distance<-fungi.distance[row.names(LW.indicator),row.names(LW.indicator)]

LZ.indicator<-indicator.table[, names(indicator.table) %in% metadata[metadata$site2 %in% "LZ",]$sample_names]
LZ.indicator<-LZ.indicator[,colSums(LZ.indicator)>0]
LZ.indicator<-LZ.indicator[rowSums(LZ.indicator)>0,]
LZ.indicator.distance<-fungi.distance[row.names(LZ.indicator),row.names(LZ.indicator)]
###
LZ<-t(LZ)
LZ<-LZ[,colnames(LZ.distance)]
LZ.ses<-ses.mpd((LZ), LZ.distance,null.model ="richness")
LV<-t(LV)
LV<-LV[,colnames(LV.distance)]
LV.ses<-ses.mpd((LV), LV.distance,null.model ="richness")
LW<-t(LW)
LW<-LW[,colnames(LW.distance)]
LW.ses<-ses.mpd((LW), LW.distance,null.model ="richness")
##indicator
LZ.indicator<-t(LZ.indicator)
LZ.indicator<-LZ.indicator[,colnames(LZ.indicator.distance)]
LZ.indicator.ses<-ses.mpd((LZ.indicator), LZ.indicator.distance,null.model ="richness")
LV.indicator<-t(LV.indicator)
LV.indicator<-LV.indicator[,colnames(LV.indicator.distance)]
LV.indicator.ses<-ses.mpd((LV.indicator), LV.indicator.distance,null.model ="richness")
LW.indicator<-t(LW.indicator)
LW.indicator<-LW.indicator[,colnames(LW.indicator.distance)]
LW.indicator.ses<-ses.mpd((LW.indicator), LW.indicator.distance,null.model ="richness")
#regression between richness and mpd
richness<-read.csv("richness.csv" ,header = T,row.names = 1,sep = ",")
metadata<-read.csv("metadata.final2.csv",row.names = 1)
metadata2<-metadata[match(richness$sample,metadata$sample_names),]
ses<-rbind(LW.ses,LV.ses,LZ.ses)
ses$sample<-row.names(ses)
try<-merge(metadata,ses[,c(2,6,7,9)],by.x = "sample_names",by.y = "sample",all=T)
try<-try[!is.na(try$mpd.obs),]
try<-merge(try,richness[,1:2],by = "sample_names")
##indicator
ses.indicator<-rbind(LW.indicator.ses,LV.indicator.ses,LZ.indicator.ses)
ses.indicator$sample<-row.names(ses.indicator)
try.indicator<-merge(metadata,ses.indicator[,c(2,6,7,9)],by.x = "sample_names",by.y = "sample",all=T)
try.indicator<-try.indicator[!is.na(try.indicator$mpd.obs),]
##
names(try)[7]<-"richness"
try2<-try
try2$substrate<-factor(try2$substrate,levels=c("subsoil","topsoil","wood","feces","fruit body",
                                               "lichen","moss","bark","leaves","leaf litter",
                                               "snow"))
mod01<-lmerTest::lmer(mpd.obs.z ~ substrate+(1|site2),data = try2)
b<-avg_comparisons(mod01, variables = list(substrate = "pairwise")) 
performance::performance(mod01)
substrate.ses<-parameters::model_parameters(mod01)
write.csv(substrate.ses,"substrate.ses.parameter.csv")
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group[letter$Group == "fruitbody"] <- "fruit body"
letter$Group[letter$Group == "leaflitter"] <- "leaf litter"
difference <- try2 %>%
  group_by(substrate) %>%
  summarise(max=mean(mpd.obs.z ))

y.site<-merge(letter,difference,by.x="Group",by.y="substrate")
y.site<-data.frame(substrate=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)


plot_predictions(mod01,condition = c("substrate"))+
  labs(x="substrate",y="mpd.obs.z ")+
  scale_x_discrete(breaks=c("subsoil","topsoil","wood","feces","fruit body",
                            "lichen","moss","bark","leaves","leaf litter",
                            "snow"))+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom")+
  geom_text(data = y.site, aes(x = substrate , y = ymax, label = letter,hjust=-0.5))

###
try$substrate<-factor(try$substrate,levels=c("subsoil","topsoil","wood","feces","fruit body",
                                               "lichen","moss","bark","leaves","leaf litter",
                                               "snow"))
mod01<-lmerTest::lmer(mpd.obs ~ substrate+(1|site2),data = try)
performance::performance(mod01)
substrate.obs<-parameters::model_parameters(mod01)
write.csv(substrate.obs,"substrate.obs.parameter.csv")
b<-avg_comparisons(mod01, variables = list(substrate = "pairwise")) 
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group[letter$Group == "fruitbody"] <- "fruit body"
letter$Group[letter$Group == "leaflitter"] <- "leaf litter"
difference <- try %>%
  group_by(substrate) %>%
  summarise(max=mean(mpd.obs))

y.site<-merge(letter,difference,by.x="Group",by.y="substrate")
y.site<-data.frame(substrate=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)

plot_predictions(mod01,condition = c("substrate"))+
  labs(x="substrate",y="mpd.obs")+
  scale_x_discrete(breaks=c("subsoil","topsoil","wood","feces","fruit body",
                            "lichen","moss","bark","leaves","leaf litter",
                            "snow"))+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom")+
  geom_text(data = y.site, aes(x = substrate , y = ymax, label = letter,hjust=-0.5))
###indicator
mod01<-lmerTest::lmer(mpd.obs ~ substrate+(1|site2),data = try.indicator)
performance::performance(mod01)
b<-avg_comparisons(mod01, variables = list(substrate = "pairwise")) 
substrate.ses<-parameters::model_parameters(mod01)
write.csv(substrate.ses,"indi.substrate.ses.parameter.csv")

letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group[letter$Group == "fruitbody"] <- "fruit body"
letter$Group[letter$Group == "leaflitter"] <- "leaf litter"
difference <- try %>%
  group_by(substrate) %>%
  summarise(max=mean(mpd.obs))

y.site<-merge(letter,difference,by.x="Group",by.y="substrate")
y.site<-data.frame(substrate=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)

plot_predictions(mod01,condition = c("substrate"))+
  labs(x="substrate",y="mpd.obs")+
  scale_x_discrete(breaks=c("subsoil","topsoil","wood","feces","fruit body",
                            "lichen","moss","bark","leaves","leaf litter",
                            "snow"))+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom")+
  geom_text(data = y.site, aes(x = substrate , y = ymax, label = letter,hjust=-0.5))
###
try.indicator2<-try.indicator
mod01<-lmerTest::lmer(mpd.obs.z ~ substrate+(1|site2),data = try.indicator2)
performance::performance(mod01)
substrate.obs<-parameters::model_parameters(mod01)
write.csv(substrate.obs,"indi.substrate.obs.parameter.csv")
b<-avg_comparisons(mod01, variables = list(substrate = "pairwise")) 
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group[letter$Group == "fruitbody"] <- "fruit body"
letter$Group[letter$Group == "leaflitter"] <- "leaf litter"
difference <- try.indicator %>%
  group_by(substrate) %>%
  summarise(max=mean(mpd.obs.z))

y.site<-merge(letter,difference,by.x="Group",by.y="substrate")
y.site<-data.frame(substrate=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)

plot_predictions(mod01,condition = c("substrate"))+
  labs(x="substrate",y="mpd.obs.z")+
  scale_x_discrete(breaks=c("subsoil","topsoil","wood","feces","fruit body",
                            "lichen","moss","bark","leaves","leaf litter",
                            "snow"))+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom")+
  geom_text(data = y.site, aes(x = substrate , y = ymax, label = letter,hjust=-0.5))

#
try<-try[,-7]
try$type<-"whole"
try.indicator$type<-"indicator"
all.s<-rbind(try, try.indicator)
##
fungi_overlap <- bind_rows(try, try.indicator)
fungi_overlap$substrate <- factor(
  fungi_overlap$substrate,
  levels = c("subsoil","topsoil","wood","feces","fruit body",
             "lichen","moss","bark","leaves","leaf litter","snow")
)
# plot overlayed points, faceted by lifestyle
out <- ggplot(fungi_overlap, aes(x = substrate, y = mpd.obs,color=type)) +
  geom_boxplot() +
  scale_color_manual(values = c("whole" = "#008080", "indicator" = "black")) +
  scale_x_discrete(breaks=c("subsoil","topsoil","wood","feces","fruit body",
                            "lichen","moss","bark","leaves","leaf litter",
                            "snow"))+
  labs(
    color = "type"
  ) +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )

out


out <- ggplot(fungi_overlap, aes(x = substrate, y = mpd.obs.z,color=type)) +
  geom_boxplot() +
  scale_color_manual(values = c("whole" = "#008080", "indicator" = "black")) +
  labs(
    color = "type"
  ) +
  scale_x_discrete(breaks=c("subsoil","topsoil","wood","feces","fruit body",
                   "lichen","moss","bark","leaves","leaf litter",
                   "snow"))+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )

out


