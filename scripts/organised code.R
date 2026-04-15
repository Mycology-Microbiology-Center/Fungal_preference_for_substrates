#=======package library=======
library(data.table)
library(dplyr)
library(tidyr)
library(tidyverse)
library(lme4)
library(marginaleffects)
library(rcompanion)
library(ranger)

library(vegan)
library(phyloseq)
library(picante)
library(indicspecies)
library(permute)
library(pairwiseAdonis)

library(igraph)
library(ggraph)
library(tidygraph)

library(ape)
library(adegenet)

library(ggplot2)
library(ggalluvial)
library(ggforce)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)
library(gridExtra)

library(shapviz)

library(future.apply)
library(doFuture)
library(progressr)

library(reshape2)
library(metagMisc)
#####=======rarefaction=======
table2<-read.csv("table2.csv",header = T,row.names = 1)
tax<-read.csv("taxonomy.final.csv",header = T,row.names = 1)
row.names(metadata)<-metadata$sample_names
row.names(tax)<-tax$qseqid
fungi2<-table2
fungi2$qseqid<-row.names(fungi2)
tax2<-merge(tax,fungi2,by="qseqid",all=T)
##taxonomy assignment filtering
tax2<-tax2[tax2$pident>85,]
tax2 <- tax2[tax2$evalue < 1e-20, ]
tax2<-tax2[tax2$qcovs>90,]
fungi2<-tax2[grepl("k__Fungi",tax2$stitle),]
fungi2<-fungi2[!is.na(fungi2$TUE127696),]
fungi2<-fungi2[,c(1,19:ncol(fungi2))]
tax2<-tax2[,1:2]
row.names(tax2)<-tax2$qseqid
row.names(fungi2)<-fungi2$qseqid
fungi2<-fungi2[,-1]
table2<-fungi2
table3<-table2[,colSums(table2)>=960]
table3<-table3[rowSums(table3)!=0,]
table3<-table3[,!names(table3) %in% names(sort(colSums(table3[,metadata[paste0(metadata$substrate,metadata$site2)=="topsoilLV",]$sample_names])))[c(1:(length(names(sort(colSums(table3[,metadata[paste0(metadata$substrate,metadata$site2)=="topsoilLV",]$sample_names]))))-15))]]
table3<-table3[,!names(table3) %in% names(sort(colSums(table3[,metadata[paste0(metadata$substrate,metadata$site2)=="topsoilLZ",]$sample_names])))[c(1:(length(names(sort(colSums(table3[,metadata[paste0(metadata$substrate,metadata$site2)=="topsoilLZ",]$sample_names]))))-15))]]
table3<-table3[,!names(table3) %in% names(sort(colSums(table3[,metadata[paste0(metadata$substrate,metadata$site2)=="topsoilLW",]$sample_names])))[c(1:(length(names(sort(colSums(table3[,metadata[paste0(metadata$substrate,metadata$site2)=="topsoilLW",]$sample_names]))))-15))]]

metadata2<-metadata[match(names(table3),metadata$sample_names),]
row.names(metadata2)<-metadata2$sample_names
tax3<-tax2[match(row.names(table3),tax2$qseqid),]

#get the rarefy table
lichen<-table3[,names(table3) %in% metadata2[metadata2$substrate=="lichen",]$sample_names]
lichen<-lichen[(rowSums(lichen)!=0),]
lichen[] <- lapply(lichen, function(x) {
  x[x == max(x, na.rm = TRUE)] <- 0
  return(x)
})
lichen<-lichen[,colSums(lichen)>0]

fruit.body<-table3[,names(table3) %in% metadata2[metadata2$substrate=="fruit body",]$sample_names]
fruit.body<-fruit.body[(rowSums(fruit.body)!=0),]
fruit.body[] <- lapply(fruit.body, function(x) {
  x[x == max(x, na.rm = TRUE)] <- 0
  return(x)
})
fruit.body<-fruit.body[,colSums(fruit.body)>0]


table2<-table3[,names(table3) %in% metadata2[metadata2$substrate!="fruit body"&metadata2$substrate!="lichen",]$sample_names]
table2$OTU<-row.names(table2)
fruit.body$OTU<-row.names(fruit.body)
lichen$OTU<-row.names(lichen)
me <- function(x,y){
  a<-merge(x,y,by="OTU",all=T)
  return(a)
}
tableall <- list(fruit.body, lichen, table2)
table<-Reduce(me, tableall)
table[is.na(table)]<-0
write.csv(table,"table.remove.host.csv")
######
row.names(table) <- table$OTU
table <- table[,-1]
###alpha diversity of raw table
library(vegan)
ttable<-as.data.frame(t(table))
richness1<-data.frame(richness=rowSums((ttable)>0), sample_names= row.names(ttable), reads=rowSums(ttable))
hill <- data.frame(shannon=exp(diversity(ttable, index = "shannon")),sample_names= names(exp(diversity(ttable, index = "shannon"))))
r.diversity<-merge(hill, richness1, by="sample_names")
write.csv(r.diversity,"diversity.of.raw.table.csv")
####after filtering out the host reads, some lichen or fruiting body samples have quite low sequencing depths
#### so separate them from rarefaction
table2 <- table [, colSums(table) > 976]
table2 <- table2[rowSums(table2)>0,]
fruit.lichen <- table [, colSums(table) < 976] 
metadata3 <- metadata2[metadata2$sample_names %in% names(table2), ]
fungi <- phyloseq(
  otu_table = otu_table(as.matrix(table2), taxa_are_rows = TRUE),
  tax_table = tax_table(as.matrix(tax)),
  sample_data = sample_data(metadata3)
)
min<-min(colSums(table2))
d2 <- fungi

fungirarefy<-phyloseq_mult_raref(d2,trimOTUs = F, replace = T, MinSizeThreshold = min, SampSize = min, iter = 1000,multithread=T,seeds = NULL)
####divide list into four parts to speed up
fungirarefy1<-fungirarefy[c(1:500)]
fungirarefy12<-fungirarefy1[c(1:250)]
fungirarefy11<-fungirarefy1[c(251:500)]
fungirarefy2<-fungirarefy[c(501:1000)]
fungirarefy22<-fungirarefy2[c(1:250)]
fungirarefy21<-fungirarefy2[c(251:500)]

table1<-list()
table1 <- lapply(fungirarefy11, function(x) {
  otu <- x@otu_table
  
  mat <- otu@.Data
  
  if (otu@taxa_are_rows) {
    mat <- t(mat)
  }
  
  Matrix(mat, sparse = TRUE)
})

table2<-list()
table2 <- lapply(fungirarefy12, function(x) {
  otu <- x@otu_table
  
  mat <- otu@.Data
  
  if (otu@taxa_are_rows) {
    mat <- t(mat)
  }
  
  Matrix(mat, sparse = TRUE)
})

table3<-list()
table3 <- lapply(fungirarefy21, function(x) {
  otu <- x@otu_table
  
  mat <- otu@.Data
  
  if (otu@taxa_are_rows) {
    mat <- t(mat)
  }
  
  Matrix(mat, sparse = TRUE)
})


table4<-list()
table4 <- lapply(fungirarefy22, function(x) {
  otu <- x@otu_table
  
  mat <- otu@.Data
  
  if (otu@taxa_are_rows) {
    mat <- t(mat)
  }
  
  Matrix(mat, sparse = TRUE)
})

###
plan(multisession, workers = 4)
options(future.globals.maxSize = 5e9)
handlers(global = TRUE)
handlers("txtprogressbar")
calc_phi <- function(x, perm = 999){
  
  # setDF(x)
  
  res <- list()
  
  res$PHI <- try(multipatt(
    x = x[, -ncol(x)],
    cluster = x[, ncol(x)],
    func = "r.g",
    control = how(nperm = perm),
    duleg = TRUE
  ))
  
  res$III <- try(multipatt(
    x = x[, -ncol(x)],
    cluster = x[, ncol(x)],
    func = "IndVal.g",
    control = how(nperm = perm),
    duleg = TRUE
  ))
  
  return(res)
}

fruit.lichen$OTU <- row.names(fruit.lichen)

cal_phi2 <- function(x){
  p <- progressor(along = x)
  PIRES <- future_lapply(x, function(obj){
    p()  #progressing
    a <- as.data.frame(t(as.matrix(obj)))
    a$OTU <- row.names(a)
    fungi.ind <- merge(a, fruit.lichen, by = "OTU", all = TRUE)
    fungi.ind[is.na(fungi.ind)] <- 0
    row.names(fungi.ind) <- fungi.ind$OTU
    fungi.ind <- fungi.ind[, -1]
    
    fungi.ind <- (fungi.ind > 0) * 1
    fungi.ind <- as.data.frame(t(fungi.ind))
    
    fungi.ind <- fungi.ind[, colSums(fungi.ind > 0) > 2]
    
    fungi.ind$sample_names <- row.names(fungi.ind)
    fungi2 <- merge(fungi.ind, metadata, by = "sample_names", all = TRUE)
    row.names(fungi2) <- fungi2$sample_names
    fungi2 <- fungi2[, -1]
    
    fungi2 <- fungi2[
      !is.na(fungi2$e57bceb79df4a766f8fb9d049798c662ba89dc9e),
    ]
    
    fungi2 <- fungi2[, -ncol(fungi2)]
    res <- calc_phi(fungi2)
    rm(a, fungi.ind, fungi2)
    gc()
    return(res)
  })
  return(PIRES)
}
with_progress({
  PIRES <- cal_phi2(table1)
})

with_progress({
  PIRES2 <- cal_phi2(table2)
})

with_progress({
  PIRES3 <- cal_phi2(table3)
})

with_progress({
  PIRES4 <- cal_phi2(table4)
})
###extract the IndVal.g
extract<-function(fungi2.PIRES){
  IIIs <- llply(.data = fungi2.PIRES, .fun = function(x){
    
    iii <- metagMisc::dfRowName(x = x$III$sign, name = "OTU")
    setDT(iii)
    
    treesp <- colnames(iii)
    treesp <- treesp[ ! treesp %in% c("OTU", "index", "stat", "p.value") ]
    
    iii[, BestHost := treesp[ index ] ]
    setcolorder(x = iii, neworder = c("OTU", "BestHost", "index", "stat", "p.value"))
    
    return(iii)
  })
}
result1<-extract(PIRES)
result2<-extract(PIRES2)
result3<-extract(PIRES3)
result4<-extract(PIRES4)
####taxonomy assignment
tax<-read.csv("taxonomy.final.csv",row.names = 1)
tax$phylumn[is.na(tax$phylumn)]<-"p__unclassified"
tax$class[is.na(tax$class)]<-"c__unclassified"
tax$order[is.na(tax$order)]<-"o__unclassified"
tax$family[is.na(tax$family)]<-"f__unclassified"
tax$genus[is.na(tax$genus)]<-"g__unclassified"
tax$species[is.na(tax$species)]<-"s__unclassified"
tax2<-tax
tax2$genus<-gsub(".*_","",tax2$genus)
tax2$OTU<-row.names(tax2)

for(i in 1:250){
  result1[[i]]<-merge(result1[[i]],tax2,by="OTU",all=T)
  result1[[i]]<-result1[[i]][!is.na(result1[[i]]$BestHost),]
}

for(i in 1:250){
  result2[[i]]<-merge(result2[[i]],tax2,by="OTU",all=T)
  result2[[i]]<-result2[[i]][!is.na(result2[[i]]$BestHost),]
}

for(i in 1:250){
  result3[[i]]<-merge(result3[[i]],tax2,by="OTU",all=T)
  result3[[i]]<-result3[[i]][!is.na(result3[[i]]$BestHost),]
}

for(i in 1:250){
  result4[[i]]<-merge(result4[[i]],tax2,by="OTU",all=T)
  result4[[i]]<-result4[[i]][!is.na(result4[[i]]$BestHost),]
}
save.image("phi.table.RData")
##################### individual richness and Hill number (q = 1)
fruit.lichen$OTU <- row.names(fruit.lichen)
prep_data <- function(x){
  lapply(x, function(obj){
    a <- as.data.frame(t(as.matrix(obj)))
    a$OTU <- row.names(a)
    
    fungi.ind <- merge(a, fruit.lichen, by = "OTU", all = TRUE)
    fungi.ind[is.na(fungi.ind)] <- 0
    
    row.names(fungi.ind) <- fungi.ind$OTU
    fungi.ind <- fungi.ind[, -1]
    fungi.ind <- as.data.frame(t(fungi.ind))
    fungi.ind <- fungi.ind[, colSums(fungi.ind > 0) > 2]
    fungi.ind <- fungi.ind[rowSums(fungi.ind > 0) > 0,]
    fungi.ind <- Matrix::Matrix(as.matrix(fungi.ind), sparse = TRUE)
    
    return(fungi.ind)
  })
}
table12 <- prep_data(table1)
table22 <- prep_data(table2)
table32 <- prep_data(table3)
table42 <- prep_data(table4)
save.image("tableall.RData")
richness1<-list()
hill1<-list()
for (i in 1:250) {
  richness1[[i]]<-data.frame(richness=rowSums(as.data.frame(as.matrix(table12[[i]]))>0), sample_names= row.names(as.data.frame(as.matrix(table12[[i]]))))
  hill1[[i]] <- data.frame(shannon=exp(diversity(as.data.frame(as.matrix(table12[[i]])), index = "shannon")),sample_names= names(exp(diversity(as.data.frame(as.matrix(table12[[i]])), index = "shannon"))))
}

richness2<-list()
hill2<-list()
for (i in 1:250) {
  richness2[[i]]<-data.frame(richness=rowSums(as.data.frame(as.matrix(table22[[i]]))>0), sample_names= row.names(as.data.frame(as.matrix(table22[[i]]))))
  hill2[[i]] <- data.frame(shannon=exp(diversity(as.data.frame(as.matrix(table22[[i]])), index = "shannon")),sample_names= names(exp(diversity(as.data.frame(as.matrix(table22[[i]])), index = "shannon"))))
}

richness3<-list()
hill3<-list()
for (i in 1:250) {
  richness3[[i]]<-data.frame(richness=rowSums(as.data.frame(as.matrix(table32[[i]]))>0), sample_names= row.names(as.data.frame(as.matrix(table32[[i]]))))
  hill3[[i]] <- data.frame(shannon=exp(diversity(as.data.frame(as.matrix(table32[[i]])), index = "shannon")),sample_names= names(exp(diversity(as.data.frame(as.matrix(table32[[i]])), index = "shannon"))))
}

richness4<-list()
hill4<-list()
for (i in 1:250) {
  richness4[[i]]<-data.frame(richness=rowSums(as.data.frame(as.matrix(table42[[i]]))>0), sample_names= row.names(as.data.frame(as.matrix(table42[[i]]))))
  hill4[[i]] <- data.frame(shannon=exp(diversity(as.data.frame(as.matrix(table42[[i]])), index = "shannon")),sample_names= names(exp(diversity(as.data.frame(as.matrix(table42[[i]])), index = "shannon"))))
}

richness <- bind_rows(richness1, richness2, richness3, richness4)
hill <- bind_rows(hill1, hill2, hill3, hill4)
write.table(richness,"rarefy.richness.allsamples.csv",sep=",")
write.table(hill,"rarefy.hill.allsamples.csv",sep=",")

richness2 <- richness %>%
  group_by(sample_names) %>%
  summarise(richness2=mean(richness))
hill2 <- hill %>%
  group_by(sample_names) %>%
  summarise(hill2=mean(shannon))
richness2 <- merge(richness2, metadata, by="sample_names")
hill2 <- merge(hill2, metadata, by="sample_names")
####

model <- lm(hill2~  substrate + site2, data =hill2)
###
b<-avg_comparisons(model, variables = list(substrate = "pairwise")) 
performance::performance(model)
substrate.richness<-parameters::model_parameters(model)

letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group[letter$Group == "fruitbody"] <- "fruit body"
letter$Group[letter$Group == "leaflitter"] <- "leaf litter"
difference <- hill2 %>%
  group_by(substrate) %>%
  summarise(max=mean(hill2))

y.site<-merge(letter,difference,by.x="Group",by.y="substrate")
y.site<-data.frame(substrate=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)

plot_predictions(model,condition = c("substrate"))+
  labs(x="substrate",y="hill number")+
  theme_light() +
  scale_x_discrete(limits = c("snow","leaves","topsoil","moss"
                              ,"lichen","bark", "subsoil","leaf litter",
                              "fruit body",
                              "wood","feces"))+
  scale_y_continuous(breaks = c(0,50))+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom")+
  geom_text(data = y.site, aes(x = substrate , y = ymax, label = letter,hjust=-0.5))

write.csv(substrate.richness,"hillparameter.csv")
##
model <- lm(richness2 ~  substrate + site2, data = richness2)
###
b<-avg_comparisons(model, variables = list(substrate = "pairwise")) 
performance::performance(model)
# options(scipen = 999)
substrate.richness<-parameters::model_parameters(model)
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group[letter$Group == "fruitbody"] <- "fruit body"
letter$Group[letter$Group == "leaflitter"] <- "leaf litter"
difference <- richness2 %>%
  group_by(substrate) %>%
  summarise(max=mean(richness2))
y.site<-merge(letter,difference,by.x="Group",by.y="substrate")
y.site<-data.frame(substrate=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)
plot_predictions(model,condition = c("substrate"))+
  labs(x="substrate",y="richness(OTU)")+
  theme_light() +
  scale_x_discrete(limits = c("snow","leaves","topsoil","moss"
                              ,"lichen","bark", "subsoil","leaf litter","feces",
                              "fruit body",
                              "wood"))+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom")+
  geom_text(data = y.site, aes(x = substrate , y = ymax, label = letter,hjust=-0.5))
write.csv(substrate.richness,"individual.richnessparameter.csv")
######compare total substrate richness
load("tableall.RData")
total_richness_calculation<- function(fungi){
  fungi$sample_names<-row.names(fungi)
  metadata2 <- metadata[match(fungi$sample_names, metadata$sample_names), ]
  table <- cbind(fungi[,-ncol(fungi)], metadata2[,2:3])
  table2 <- rowsum(
    table[,1:(ncol(table)-2)],
    group = interaction(table$substrate, table$site2)
  ) 
  richness2 <- data.frame(richness=rowSums(table2>0), group=row.names(table2)) %>%
    separate(group, into = c("substrate", "site"), sep = "\\.")
  return(richness2)
  }

richness1<-list()
for(i in 1:250) {
  a<-as.data.frame(as.matrix(table12[[i]]))
  richness1[[i]] <- total_richness_calculation(a) 
}

richness2<-list()
for(i in 1:250) {
  a<-as.data.frame(as.matrix(table22[[i]]))
  richness2[[i]] <- total_richness_calculation(a) 
}


richness3<-list()
for(i in 1:250) {
  a<-as.data.frame(as.matrix(table32[[i]]))
  richness3[[i]] <- total_richness_calculation(a) 
}


richness4<-list()
for(i in 1:250) {
  a<-as.data.frame(as.matrix(table42[[i]]))
  richness4[[i]] <- total_richness_calculation(a) 
}

richness <- bind_rows(richness1, richness2, richness3, richness4)
richness2 <- richness %>%
  group_by(site, substrate) %>%
  dplyr::summarise(richness2=mean(richness))
#####
model <- lm(richness2 ~  substrate + site, data = richness2)
b<-avg_comparisons(model, variables = list(substrate = "pairwise")) 
performance::performance(model)
substrate.richness<-parameters::model_parameters(model)
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group[letter$Group == "fruitbody"] <- "fruit body"
letter$Group[letter$Group == "leaflitter"] <- "leaf litter"
difference <- richness2 %>%
  group_by(substrate) %>%
  summarise(max=mean(richness2))

y.site<-merge(letter,difference,by.x="Group",by.y="substrate")
y.site<-data.frame(substrate=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)
library(ggplot2)
plot_predictions(model,condition = c("substrate"))+
  labs(x="substrate",y="richness")+
  theme_light() +
  scale_x_discrete(limits = c("topsoil","snow","moss","lichen","bark"
                              ,"leaves","leaf litter",
                              "fruit body","subsoil",
                              "wood","feces"))+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom")+
  geom_text(data = y.site, aes(x = substrate , y = ymax, label = letter,hjust=-0.5))
write.csv(substrate.richness,"total.richnessparameter.csv")
#######relative abundance within substrate

tables <- c(table12, table22, table32, table42)
rel_tables <- lapply(tables, function(x) {
  rs <- Matrix::rowSums(x)
  x / rs
})

all_cols <- Reduce(union, lapply(tables, colnames))
align_matrix <- function(x, all_cols) {
  missing <- setdiff(all_cols, colnames(x))
  
  if (length(missing) > 0) {
    extra <- Matrix(0, nrow(x), length(missing), sparse = TRUE)
    colnames(extra) <- missing
    x <- cbind(x, extra)
  }
  
  x[, all_cols]
}

tables_aligned <- lapply(tables, align_matrix, all_cols = all_cols)

rel_tables <- lapply(tables_aligned, function(x) {
  rs <- Matrix::rowSums(x)
  x / rs
})
######get the mean relative abundance of 1,000 tables
mean_table <- Reduce(`+`, rel_tables) / length(rel_tables)
mean_table <- as.data.frame(as.matrix(mean_table))

table_2<- mean_table
table_2$sample_names<-row.names(table_2)
table_2<-merge(table_2,metadata,by="sample_names")
table_4 <- aggregate(table_2[,2:(ncol(table_2)-2)],by=list(substrate=table_2$substrate,site=table_2$site2),mean)
table_3 <- aggregate(table_4[,3:(ncol(table_4))],by=list(substrate=table_4$substrate),mean)

mean_table<-mean_table[rowSums(mean_table)>0,]
metadata <- metadata[match(row.names(mean_table),metadata$sample_names),]
###
fungi<-mean_table
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
row.names(table_3)<-table_3$substrate
table_3<-table_3[,-1]
table_long <- table_3 %>%
  rownames_to_column(var = "substrate") %>%
  pivot_longer(cols = where(is.numeric), names_to = "OTU", values_to = "Ratio")%>%
  filter(Ratio>0)
##
table_42 <- table_4 %>%
  pivot_longer(cols = where(is.numeric), names_to = "OTU", values_to = "Ratio")%>%
  filter(Ratio>0)
###use the mean relative abundance of OTUs as the baseline to show the key OTUs changing in relative abundance
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
##calculate the consistence of dominant lifestyle
table422<-merge(table_42,table,by.x="OTU",by.y="qseqid",all=T)
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
grid.arrange( p2,p1, ncol = 2)
###
unique(table_long2$GENUS)
a<- table_long2 %>%
  group_by(GENUS) %>%
  summarise(Ratio=length(unique(substrate))/11)
n<-table_long2[table_long2$GENUS %in% a$GENUS[order(a$Ratio)[1:10]],]

a<- try %>%
  group_by(primary_lifestyle) %>%
  summarise(Ratio=length(unique(substrate))/11)
n<-try[try$primary_lifestyle %in% a$primary_lifestyle[order(a$Ratio)[1:10]],]

####Simpson diddimilarity#####
simpson_abund <- function(x){
  # x = otu table, taxa = columns
  
  ## Function for Simpson dissimilarity for a pair of samples
  abund_simps_pair <- function(x){
    # x = data.frame or matrix with two rows (samples), relative species abundances
    
    ## Remove double zeros and subset to shared species
    xx <- x[, colSums(x > 0) == 2, drop = FALSE]
    
    if(ncol(xx) > 0){
      
      ## Estimate sum of total [relative] abundances of shared species for each community
      U <- rowSums(xx[1,, drop = FALSE])
      V <- rowSums(xx[2,, drop = FALSE])
      
      ## Estimate Simpson's / Lennon's dissimilarity coefficient (see Chao et al., 2006 Biometrics; Table 2)
      UV <- U*V
      dd <- UV / (UV + min(U - UV, V - UV))
      dd <- 1 - dd
    } else {
      # No shared species
      dd <- 1
    }
    return(dd)
  }
  
  ## Convert to relative abundances
  x <- decostand(x, method = "total", MARGIN = 1)
  
  ## Initialize distance matrix
  dd <- matrix(data = NA, nrow = nrow(x), ncol = nrow(x),
               dimnames = list(rownames(x), rownames(x)))
  
  ## Estimate pairwise dissimilarities
  cmb <- combn(x = 1:nrow(x), m = 2)
  
  for(i in 1:ncol(cmb)){
    s1 <- cmb[1, i]
    s2 <- cmb[2, i]
    ds <- abund_simps_pair(x = x[c(s1, s2), ])
    dd[s1, s2] <- ds
    dd[s2, s1] <- ds
    rm(ds)
  }
  
  diag(dd) <- 0
  dd <- as.dist(dd)
  return(dd)
}


cal_dist2 <- function(x){
  p <- progressor(along = x)
  
  dist <- future_lapply(x, function(obj){
    p()
    if(inherits(obj, "dgCMatrix")){
      obj <- as.data.frame(as.matrix(obj))
    }
    
    res <- simpson_abund_fast(obj)
    return(res)
  })
  
  return(dist)
}


plan(multisession, workers = 4)
options(future.globals.maxSize = 5e9)
handlers(global = TRUE)

dist41 <- with_progress({
  cal_dist2(table42)
})


dist31 <- with_progress({
  cal_dist2(table32)
})

dist21 <- with_progress({
  cal_dist2(table22)
})

dist11 <- with_progress({
  cal_dist2(table12)
})
#######======continue
files <- list.files("simpson dis", full.names = TRUE) ###here we exported the Simpson dissimilarity to a separated directory.
###
data_list <- c(dist11, dist21, dist31, dist41)
data_list <- lapply(data_list, as.matrix)
all(sapply(data_list, function(x)
  identical(rownames(x), rownames(data_list[[1]]))
))  ##check if the order is identical
###
arr <- simplify2array(data_list)
###get the median the Simpson dissimilarity within a same pair in 1,000 distance matrix
median_mat <- apply(arr, c(1, 2), median, na.rm = TRUE)
metadata<-read.csv("metadata.final2.csv",row.names = 1)
row.names(metadata)<-metadata$sample_names
median_mat<-as.dist(median_mat)
metadata3<-metadata[attr(median_mat,"Labels"),]
####
permanova_result <- adonis2(median_mat ~ substrate,strata = metadata3$site2, data = metadata3,permutations = 999)
pairwise_results <- pairwise.adonis(median_mat,factors = metadata3$substrate, perm = 999, p.adjust.m = "fdr")
substrate_shapes <- c("feces" = 15,"leaves" = 7,  "snow" = 8,  "fruit body" = 18,
                      "lichen" = 16,  "topsoil"=1,  "subsoil"=13,  "leaf litter"=22,  "wood"=10,
                      "bark"=3,  "moss"=4)
distance<-median_mat
pcoa1<-pcoa(distance)
data.scores <- as.data.frame(pcoa1$vectors)
data.scores$sample_names <- row.names(data.scores)
aaa<-merge(metadata3,data.scores,by="sample_names")
##
out1<-ggplot(aaa, aes(x = Axis.1, y = Axis.2, )) + 
  geom_point(aes(color = substrate,shape=substrate)) +
  scale_shape_manual(values = substrate_shapes) +
  stat_ellipse(geom = "polygon", aes(group = substrate, color = substrate, fill = substrate),
               alpha = 0.01, level = 0.95)+
  labs(x = paste("PCoA1 (", round(pcoa1$values$Relative_eig[1] * 100, 1), "%)", sep = ""),
       y = paste("PCoA2 (", round(pcoa1$values$Relative_eig[2] * 100, 1), "%)", sep = "")) +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )

colors<-c("#882233","#E691C9", "#EDBB00","#E31864",
          "#0005A4", "#99999D", "#0FF001" , "#B35EEE", "#009666",
          "#996677","#00CCFF")
colors <- setNames(colors, c("topsoil", "subsoil", "leaf litter","bark","wood",       
                             "lichen","moss","fruit body","leaves","feces",      
                             "snow"))
out1<-out1 + scale_colour_manual(values = colors) + 
  scale_fill_manual(values = colors)
###
out2<-ggplot(aaa, aes(x = Axis.3, y = Axis.4, )) + 
  geom_point(aes(color = substrate,shape=substrate)) +
  scale_shape_manual(values = substrate_shapes) +
  stat_ellipse(geom = "polygon", aes(group = substrate, color = substrate, fill = substrate),
               alpha = 0.01, level = 0.95)+
  labs(x = paste("PCoA3 (", round(pcoa1$values$Relative_eig[3] * 100, 1), "%)", sep = ""),
       y = paste("PCoA4 (", round(pcoa1$values$Relative_eig[4] * 100, 1), "%)", sep = "")) +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )
out2<-out2 + scale_colour_manual(values = colors) + 
  scale_fill_manual(values = colors)
#####similarity heat map
fungi.d2 <- dist2list(median_mat)
fungi.d2<-merge(fungi.d2,metadata3,by.x="col",by.y = "sample_names",all=T )
fungi.d2<-merge(fungi.d2,metadata3,by.x="row",by.y = "sample_names",all=T )
fungi.d2<-fungi.d2[fungi.d2$site2.x == fungi.d2$site2.y,]
fungi.d2<-fungi.d2[!is.na(fungi.d2$row),]
fungi.d2<-fungi.d2[!is.na(fungi.d2$col),]
fungi.d2$value <- 1-fungi.d2$value
fungi.d2<-aggregate(fungi.d2$value, by=list(substrate.x=fungi.d2$substrate.x, substrate.y=fungi.d2$substrate.y,site=fungi.d2$site2.x), data = fungi.d2, FUN = mean, na.rm = TRUE)
fungi.d2<-aggregate(x ~ substrate.x + substrate.y, data = fungi.d2, FUN = mean, na.rm = TRUE)
##
names(fungi.d2)<-c("row","col","value")
samples <- sort(unique(c(fungi.d2$row, fungi.d2$col)))
similarity_matrix <- matrix(NA, nrow = length(samples), ncol = length(samples),
                            dimnames = list(samples, samples))
for (i in seq_len(nrow(fungi.d2))) {
  r <- fungi.d2$row[i]
  c <- fungi.d2$col[i]
  v <- fungi.d2$value[i]
  similarity_matrix[r, c] <- v
  similarity_matrix[c, r] <- v  # Fill symmetric cell with the same value
}
similarity_matrix <- as.matrix(similarity_matrix)

od =  hclust(dist(similarity_matrix))$order
m2 = similarity_matrix[od, od]
h<-Heatmap(similarity_matrix, name = "Ratio",
        cluster_rows = T,
        cluster_columns = TRUE,
        col = colorRamp2(c(0, 0.1, 0.2, 0.3, 0.5), c("white", "lightblue", "#FFEEA0", "#F4BB44", "#C04000")),
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 9),
        clustering_method_rows = "average",
        clustering_method_columns = "average",
        cell_fun = function(j, i, x, y, width, height, fill) {
          val <- similarity_matrix[i, j]
          txt <- sub("\\.?0+$", "", format(val, trim = TRUE, digits = 2))  # remove trailing zeros
          grid.text(txt, x, y, gp = gpar(fontsize = 7))
        })
write.csv(pairwise_results,"pairwise_results.csv")
write.csv(permanova_result,"permanova_result.csv")
write.csv(fungi.d2,"mean.similarity.result.across.site.csv")
##relative abundance of indicator
table_1<-list()
table_long1<-list()
proportion_long1<-list()
for (i in 1:250) {
  fungi<-table12[[i]]
  if(inherits(fungi, "dgCMatrix")){
    table2 <- as.data.frame(as.matrix(fungi))
  }
table2$sample_names<-row.names(table2)
table2<-merge(table2,metadata,by="sample_names")
proportion <- table2 %>%
  mutate(across(2:(ncol(table2)-2), ~ as.numeric(.x > 0)))
table_1[[i]]<-data.frame(reads=rowSums(table2[,2:(ncol(table2)-2)]),site=table2$site,substrate=table2$substrate,sample=table2$sample_names)
###
table_long1[[i]] <- table2 %>%
  pivot_longer(
    cols = -c(substrate,sample_names, site2),    # all OTU columns
    names_to = "OTU",              # OTU names go here
    values_to = "reads"            # abundance values
  ) %>%
  group_by(substrate,sample_names,site2,OTU)  %>%
  summarise(reads.sum=sum(reads))  %>%
  filter(reads.sum>0)
###
proportion_long1[[i]] <- proportion %>%
  pivot_longer(
    cols = -c(sample_names,substrate,site2),    # all OTU columns
    names_to = "OTU",              # OTU names go here
    values_to = "reads"            # abundance values
  ) %>%
  filter(reads>0)
}

table_2<-list()
table_long2<-list()
proportion_long2<-list()
for (i in 1:250) {
  fungi<-table22[[i]]
  if(inherits(fungi, "dgCMatrix")){
    table2 <- as.data.frame(as.matrix(fungi))
  }
  table2$sample_names<-row.names(table2)
  table2<-merge(table2,metadata,by="sample_names")
  proportion <- table2 %>%
    mutate(across(2:(ncol(table2)-2), ~ as.numeric(.x > 0)))
  table_2[[i]]<-data.frame(reads=rowSums(table2[,2:(ncol(table2)-2)]),site=table2$site,substrate=table2$substrate,sample=table2$sample_names)
  ###
  table_long2[[i]] <- table2 %>%
    pivot_longer(
      cols = -c(substrate,sample_names, site2),    # all OTU columns
      names_to = "OTU",              # OTU names go here
      values_to = "reads"            # abundance values
    ) %>%
    group_by(substrate,sample_names,site2,OTU)  %>%
    summarise(reads.sum=sum(reads))  %>%
    filter(reads.sum>0)
  ###
  proportion_long2[[i]] <- proportion %>%
    pivot_longer(
      cols = -c(sample_names,substrate,site2),    # all OTU columns
      names_to = "OTU",              # OTU names go here
      values_to = "reads"            # abundance values
    ) %>%
    filter(reads>0)
}

table_3<-list()
table_long3<-list()
proportion_long3<-list()
for (i in 1:250) {
  fungi<-table32[[i]]
  if(inherits(fungi, "dgCMatrix")){
    table2 <- as.data.frame(as.matrix(fungi))
  }
  table2$sample_names<-row.names(table2)
  table2<-merge(table2,metadata,by="sample_names")
  proportion <- table2 %>%
    mutate(across(2:(ncol(table2)-2), ~ as.numeric(.x > 0)))
  table_3[[i]]<-data.frame(reads=rowSums(table2[,2:(ncol(table2)-2)]),site=table2$site,substrate=table2$substrate,sample=table2$sample_names)
  ###
  table_long3[[i]] <- table2 %>%
    pivot_longer(
      cols = -c(substrate,sample_names, site2),    # all OTU columns
      names_to = "OTU",              # OTU names go here
      values_to = "reads"            # abundance values
    ) %>%
    group_by(substrate,sample_names,site2,OTU)  %>%
    summarise(reads.sum=sum(reads))  %>%
    filter(reads.sum>0)
  ###
  proportion_long3[[i]] <- proportion %>%
    pivot_longer(
      cols = -c(sample_names,substrate,site2),    # all OTU columns
      names_to = "OTU",              # OTU names go here
      values_to = "reads"            # abundance values
    ) %>%
    filter(reads>0)
}

table_4<-list()
table_long4<-list()
proportion_long4<-list()
for (i in 1:250) {
  fungi<-table42[[i]]
  if(inherits(fungi, "dgCMatrix")){
    table2 <- as.data.frame(as.matrix(fungi))
  }
  table2$sample_names<-row.names(table2)
  table2<-merge(table2,metadata,by="sample_names")
  proportion <- table2 %>%
    mutate(across(2:(ncol(table2)-2), ~ as.numeric(.x > 0)))
  table_4[[i]]<-data.frame(reads=rowSums(table2[,2:(ncol(table2)-2)]),site=table2$site,substrate=table2$substrate,sample=table2$sample_names)
  ###
  table_long4[[i]] <- table2 %>%
    pivot_longer(
      cols = -c(substrate,sample_names, site2),    # all OTU columns
      names_to = "OTU",              # OTU names go here
      values_to = "reads"            # abundance values
    ) %>%
    group_by(substrate,sample_names,site2,OTU)  %>%
    summarise(reads.sum=sum(reads))  %>%
    filter(reads.sum>0)
  proportion_long4[[i]] <- proportion %>%
    pivot_longer(
      cols = -c(sample_names,substrate,site2),    # all OTU columns
      names_to = "OTU",              # OTU names go here
      values_to = "reads"            # abundance values
    ) %>%
    filter(reads>0)
}
#######indicator
for (i in 1:250) {
  result2[[i]]$BestHost <- gsub("^s\\.", "", result2[[i]]$BestHost ) 
  result2[[i]]$BestHost <- gsub("\\.", " ", result2[[i]]$BestHost )
  result2[[i]] <- result2[[i]][result2[[i]]$p.value <0.05,]
}

for (i in 1:250) {
  result3[[i]]$BestHost <- gsub("^s\\.", "", result3[[i]]$BestHost ) 
  result3[[i]]$BestHost <- gsub("\\.", " ", result3[[i]]$BestHost ) 
  result3[[i]] <- result3[[i]][result3[[i]]$p.value <0.05,]
}

for (i in 1:250) {
  result4[[i]]$BestHost <- gsub("^s\\.", "", result4[[i]]$BestHost ) 
  result4[[i]]$BestHost <- gsub("\\.", " ", result4[[i]]$BestHost ) 
  result4[[i]] <- result4[[i]][result4[[i]]$p.value <0.05,]
}

for (i in 1:250) {
  result1[[i]]$BestHost <- gsub("^s\\.", "", result1[[i]]$BestHost ) 
  result1[[i]]$BestHost <- gsub("\\.", " ", result1[[i]]$BestHost ) 
  result1[[i]] <- result1[[i]][result1[[i]]$p.value <0.05,]
}

###calculateing the relative abundance of all indicators
indicator.p1<-list()
indicator.r2_1 <- list()
all1<-list()
for(i in 1:250){
  r<-result1[[i]]
  try<-table_long1[[i]]
  proportion2<-proportion_long1[[i]]
  names(r)[2]<-"substrate"
  
  all<-merge(r[,c(1:2,4,5)],try, by=c("substrate","OTU"))

  pall<-merge(r[,c(1:2,4,5)],proportion2, by=c("substrate","OTU"))
  pall2<-pall %>%
    distinct(substrate,site2,sample_names,OTU)
  indicator.reads<- all  %>%
    group_by(substrate,sample_names,site2)  %>%
    summarise(reads=sum(reads.sum),stat = first(stat),
              p.value = first(p.value),
              .groups = "drop")
  
  indicator.p1[[i]]<- pall  %>%
    group_by(substrate,site2,sample_names)  %>%
    summarise(reads=sum(reads),
              .groups = "drop")
  
  names(indicator.reads)[4]<-"indicator.reads"
  names(table_1[[i]])[4]<-"sample_names"
  all<-merge(indicator.reads, table_1[[i]][,c(1,4)],by="sample_names")
  all1[[i]]<-all
  result <- all %>%
    group_by(substrate,site2)  %>%
    summarise(indicator.reads2=sum(indicator.reads),reads2=sum(reads))
  result$ratio<-result$indicator.reads2/result$reads2
  indicator.r<-result
  indicator.r2_1[[i]]<-indicator.r %>%
    group_by(substrate) %>%
    summarise(sd=sd(ratio),mean=mean(ratio))
  
}


indicator.p2<-list()
indicator.r2_2 <- list()
all2<-list()
for(i in 1:250){
  r<-result2[[i]]
  ####
  try<-table_long2[[i]]
  proportion2<-proportion_long2[[i]]
  names(r)[2]<-"substrate"
  
  all<-merge(r[,c(1:2,4,5)],try, by=c("substrate","OTU"))
  pall<-merge(r[,c(1:2,4,5)],proportion2, by=c("substrate","OTU"))
  pall2<-pall %>%
    distinct(substrate,site2,sample_names,OTU)
  indicator.reads<- all  %>%
    group_by(substrate,sample_names,site2)  %>%
    summarise(reads=sum(reads.sum),stat = first(stat),
              p.value = first(p.value),
              .groups = "drop")
  
  indicator.p2[[i]]<- pall  %>%
    group_by(substrate,site2,sample_names)  %>%
    summarise(reads=sum(reads),
              .groups = "drop")
  
  names(indicator.reads)[4]<-"indicator.reads"
  names(table_2[[i]])[4]<-"sample_names"
  all<-merge(indicator.reads, table_2[[i]][,c(1,4)],by="sample_names")
  all2[[i]]<-all
  result <- all %>%
    group_by(substrate,site2)  %>%
    summarise(indicator.reads2=sum(indicator.reads),reads2=sum(reads))
  result$ratio<-result$indicator.reads2/result$reads2
  indicator.r<-result
  indicator.r2_2[[i]]<-indicator.r %>%
    group_by(substrate) %>%
    summarise(sd=sd(ratio),mean=mean(ratio))
  
}

indicator.p3<-list()
indicator.r2_3 <- list()
all3<-list()
for(i in 1:250){
  r<-result3[[i]]
  try<-table_long3[[i]]
  proportion2<-proportion_long3[[i]]
  names(r)[2]<-"substrate"
  
  all<-merge(r[,c(1:2,4,5)],try, by=c("substrate","OTU"))

  pall<-merge(r[,c(1:2,4,5)],proportion2, by=c("substrate","OTU"))
  pall2<-pall %>%
    distinct(substrate,site2,sample_names,OTU)
  indicator.reads<- all  %>%
    group_by(substrate,sample_names,site2)  %>%
    summarise(reads=sum(reads.sum),stat = first(stat),
              p.value = first(p.value),
              .groups = "drop")
  
  indicator.p3[[i]]<- pall  %>%
    group_by(substrate,site2,sample_names)  %>%
    summarise(reads=sum(reads),
              .groups = "drop")
  
  names(indicator.reads)[4]<-"indicator.reads"
  names(table_3[[i]])[4]<-"sample_names"
  all<-merge(indicator.reads, table_3[[i]][,c(1,4)],by="sample_names")
  all3[[i]]<-all
  result <- all %>%
    group_by(substrate,site2)  %>%
    summarise(indicator.reads2=sum(indicator.reads),reads2=sum(reads))
  result$ratio<-result$indicator.reads2/result$reads2
  indicator.r<-result
  indicator.r2_3[[i]]<-indicator.r %>%
    group_by(substrate) %>%
    summarise(sd=sd(ratio),mean=mean(ratio))
  
}

indicator.p4<-list()
indicator.r2_4<- list()
all4<-list()
for(i in 1:250){
  r<-result3[[i]]
  try<-table_long4[[i]]
  proportion2<-proportion_long4[[i]]
  names(r)[2]<-"substrate"
  all<-merge(r[,c(1:2,4,5)],try, by=c("substrate","OTU"))

  pall<-merge(r[,c(1:2,4,5)],proportion2, by=c("substrate","OTU"))
  pall2<-pall %>%
    distinct(substrate,site2,sample_names,OTU)
  indicator.reads<- all  %>%
    group_by(substrate,sample_names,site2)  %>%
    summarise(reads=sum(reads.sum),stat = first(stat),
              p.value = first(p.value),
              .groups = "drop")
  
  indicator.p4[[i]]<- pall  %>%
    group_by(substrate,site2,sample_names)  %>%
    summarise(reads=sum(reads),
              .groups = "drop")
  
  names(indicator.reads)[4]<-"indicator.reads"
  names(table_4[[i]])[4]<-"sample_names"
  all<-merge(indicator.reads, table_4[[i]][,c(1,4)],by="sample_names")
  all4[[i]]<-all
  result <- all %>%
    group_by(substrate,site2)  %>%
    summarise(indicator.reads2=sum(indicator.reads),reads2=sum(reads))
  result$ratio<-result$indicator.reads2/result$reads2
  indicator.r<-result
  indicator.r2_4[[i]]<-indicator.r %>%
    group_by(substrate) %>%
    summarise(sd=sd(ratio),mean=mean(ratio))
  
}

###
indicator.r2 <- c(indicator.r2_1, indicator.r2_2, indicator.r2_3, indicator.r2_4)
indicator.r2<- bind_rows(indicator.r2)
indicator.r2 <- indicator.r2 %>%
  group_by(substrate) %>%
  summarise(mean2=mean(mean))
indicator.r2 <- indicator.r2 %>%
  mutate(substrate = factor(substrate, levels = substrate[order(-mean2)]))

indicator.r3<-c(all1,all2,all3,all4)
indicator.r3<- bind_rows(indicator.r3)
indicator.r3$ratio <- with(indicator.r3, indicator.reads / reads)
indicator.r3<-indicator.r3 %>%
  group_by(sample_names, substrate, site2) %>%
  summarise(ratio2=mean(ratio))

model <- lm(ratio2 ~ substrate + site2, data = indicator.r3)
b<-avg_comparisons(model, variables = list(substrate = "pairwise")) 
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group[letter$Group == "fruitbody"] <- "fruit body"
letter$Group[letter$Group == "leaflitter"] <- "leaf litter"
difference <- indicator.r3 %>%
  group_by(substrate) %>%
  summarise(max=mean(ratio2))

y.site<-merge(letter,difference,by.x="Group",by.y="substrate")
y.site<-data.frame(substrate=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)

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
indicator.p<-c(indicator.p1, indicator.p2, indicator.p3, indicator.p4)
indicator.p<- bind_rows(indicator.p)
indicator.p <- indicator.p %>%
  group_by(sample_names, substrate, site2) %>%
  summarise(reads2=mean(reads))

try<-merge(indicator.p[,c(1,4)],indicator.r3,by="sample_names")
names(try)[2]<-"indicator.reads"

r2<-richness %>%
  group_by (sample_names)%>%
  dplyr::summarise(richness2=mean(richness))
try<-merge(try, r2,by="sample_names")
mod<-lm(ratio2 ~ richness2 + substrate + site2,data=try) ##relative abundance and richness
relative<-parameters::parameters(mod)
performance::performance(mod)
write.csv(relative,"indicator.relative.abundance and richness.parameters.csv")
mod3<-lm(indicator.reads ~ richness2 + substrate + site2,data=try)##indicator number and richness
number<-parameters::parameters(mod3)
performance::performance(mod3)
write.csv(number,"indicator.number.and.richness.parameters.csv")
######sankey plot##
try1<-list()
try12<-list()
for (i in 1:250) {
  fungi<-table12[[i]]
  if(inherits(fungi, "dgCMatrix")){
    table2 <- as.data.frame(as.matrix(fungi))
  }
  
  table2$sample_names<-row.names(table2)
  table2<-merge(table2,metadata,by="sample_names")
  table4 <- aggregate(table2[,2:(ncol(table2)-2)],by=list(substrate=table2$substrate,site=table2$site2),mean)
  ###
  table_long <- table4 %>%
    pivot_longer(
      cols = -c(substrate, site),    
      names_to = "OTU",             
      values_to = "reads"            
    ) %>%
    group_by(substrate,OTU)  %>%
    summarise(relative.sum=mean(reads))  %>%
    filter(relative.sum>0)
  indicator<-result1[[i]][,c(1,2)]
  names(indicator)<-c("OTU","substrate")
  indicator$substrate<-gsub("\\.", " ", indicator$substrate)
  try<-merge(indicator,table_long,by=c("OTU","substrate"))
  try2 <- aggregate(try$relative.sum,by=list(substrate=try$substrate,OTU=try$OTU),sum)
  try3 <- table_long[table_long$OTU %in% indicator$OTU,] ##the indicators go into other substrate that are not indicators
  indicator$group<-"indicator"
  try<-merge(indicator,try3,by=c("OTU","substrate"),all=T)
  try[is.na(try$group),]$group<-"nonindicator"
  try1[[i]]<-try
  try12[[i]]<-try2
}

try_2<-list()
try22<-list()
for (i in 1:250) {
  fungi<-table22[[i]]
  if(inherits(fungi, "dgCMatrix")){
    table2 <- as.data.frame(as.matrix(fungi))
  }
  
  table2$sample_names<-row.names(table2)
  table2<-merge(table2,metadata,by="sample_names")
  table4 <- aggregate(table2[,2:(ncol(table2)-2)],by=list(substrate=table2$substrate,site=table2$site2),mean)
  ###
  table_long <- table4 %>%
    pivot_longer(
      cols = -c(substrate, site),    
      names_to = "OTU",             
      values_to = "reads"            
    ) %>%
    group_by(substrate,OTU)  %>%
    summarise(relative.sum=mean(reads))  %>%
    filter(relative.sum>0)
  indicator<-result2[[i]][,c(1,2)]
  names(indicator)<-c("OTU","substrate")
  indicator$substrate<-gsub("\\.", " ", indicator$substrate)
  try<-merge(indicator,table_long,by=c("OTU","substrate"))
  try2 <- aggregate(try$relative.sum,by=list(substrate=try$substrate,OTU=try$OTU),sum)
  try3 <- table_long[table_long$OTU %in% indicator$OTU,] ##the indicators go into other substrate that are not indicators
  indicator$group<-"indicator"
  try<-merge(indicator,try3,by=c("OTU","substrate"),all=T)
  try[is.na(try$group),]$group<-"nonindicator"
  try_2[[i]]<-try
  try22[[i]]<-try2
}

try_3<-list()
try32<-list()
for (i in 1:250) {
  fungi<-table32[[i]]
  if(inherits(fungi, "dgCMatrix")){
    table2 <- as.data.frame(as.matrix(fungi))
  }
  
  table2$sample_names<-row.names(table2)
  table2<-merge(table2,metadata,by="sample_names")
  table4 <- aggregate(table2[,2:(ncol(table2)-2)],by=list(substrate=table2$substrate,site=table2$site2),mean)
  ###
  table_long <- table4 %>%
    pivot_longer(
      cols = -c(substrate, site),    
      names_to = "OTU",             
      values_to = "reads"            
    ) %>%
    group_by(substrate,OTU)  %>%
    summarise(relative.sum=mean(reads))  %>%
    filter(relative.sum>0)
  indicator<-result3[[i]][,c(1,2)]
  names(indicator)<-c("OTU","substrate")
  indicator$substrate<-gsub("\\.", " ", indicator$substrate)
  try<-merge(indicator,table_long,by=c("OTU","substrate"))
  try2 <- aggregate(try$relative.sum,by=list(substrate=try$substrate,OTU=try$OTU),sum)
  try3 <- table_long[table_long$OTU %in% indicator$OTU,] ##the indicators go into other substrate that are not indicators
  indicator$group<-"indicator"
  try<-merge(indicator,try3,by=c("OTU","substrate"),all=T)
  try[is.na(try$group),]$group<-"nonindicator"
  try_3[[i]]<-try
  try32[[i]]<-try2
}

try4<-list()
try42<-list()
for (i in 1:250) {
  fungi<-table42[[i]]
  if(inherits(fungi, "dgCMatrix")){
    table2 <- as.data.frame(as.matrix(fungi))
  }
  
  table2$sample_names<-row.names(table2)
  table2<-merge(table2,metadata,by="sample_names")
  table4 <- aggregate(table2[,2:(ncol(table2)-2)],by=list(substrate=table2$substrate,site=table2$site2),mean)
  ###
  table_long <- table4 %>%
    pivot_longer(
      cols = -c(substrate, site),    
      names_to = "OTU",             
      values_to = "reads"            
    ) %>%
    group_by(substrate,OTU)  %>%
    summarise(relative.sum=mean(reads))  %>%
    filter(relative.sum>0)
  indicator<-result4[[i]][,c(1,2)]
  names(indicator)<-c("OTU","substrate")
  indicator$substrate<-gsub("\\.", " ", indicator$substrate)
  try<-merge(indicator,table_long,by=c("OTU","substrate"))
  try2 <- aggregate(try$relative.sum,by=list(substrate=try$substrate,OTU=try$OTU),sum)
  try3 <- table_long[table_long$OTU %in% indicator$OTU,] ##the indicators go into other substrate that are not indicators
  indicator$group<-"indicator"
  try<-merge(indicator,try3,by=c("OTU","substrate"),all=T)
  try[is.na(try$group),]$group<-"nonindicator"
  try4[[i]]<-try
  try42[[i]]<-try2
}
########
df_sankey1<-list()
flow_edges1<-list()
for (i in 1:250) {
  ### ----------------------------------------------------
  ### 1. Sample data
  ### ----------------------------------------------------
  df <- try1[[i]]
  ###
  sankey<-try
  sankey$relative.sum<-ifelse(sankey$relative.sum>0,1,0)
  ### ----------------------------------------------------
  ### 2. Define indicator OTUs for each substrate
  ### Replace with your own indicator-OTU list if needed
  ### ----------------------------------------------------
  indicator_list <- df %>% 
    group_by(substrate) %>%
    summarise(indicators = list(unique(OTU)))
  indicator_df <- try2[,1:2]
  ### ----------------------------------------------------
  ### 3. Build indicator-flow edges
  ### ----------------------------------------------------
  flow_edges <- data.frame()
  
  substrates <- sort(unique(df$substrate))
  
  for (A in substrates) {
    
    # indicator OTUs OF A
    inds_A <- indicator_df$OTU[indicator_df$substrate == A]
    
    # records for indicator OTUs
    df_i <- df %>% filter(OTU %in% inds_A)
    
    # find other substrates
    other_subs <- setdiff(unique(df_i$substrate), A)
    
    for (B in other_subs) {
      
      df_AB <- df_i %>% filter(substrate == B)
      otus_AB <- unique(df_AB$OTU)
      
      # indicator OTUs of B
      inds_B <- indicator_df$OTU[indicator_df$substrate == B]
      
      # keep only OTUs that are indicator of A
      # AND not indicator of B
      exchanged <- setdiff(otus_AB, inds_B)
      
      if (length(exchanged) > 0) {
        flow_edges <- rbind(
          flow_edges,
          data.frame(
            from = A,
            to = B,
            weight = length(exchanged)
          )
        )
      }
    }
  }
  ###consistent indicators
  df2<-df[duplicated(df$OTU),]
  df3<-df[!df$OTU %in% unique(df2$OTU), ]
  consistency<-df3
  consistency$relative.sum<-ifelse(consistency$relative.sum>0,1,0)
  consistency <- consistency %>%
    group_by(substrate) %>%
    summarise(weight=sum(relative.sum))
  
  consistency2 <- data.frame(from= consistency$substrate, to= consistency$substrate, weight= consistency $weight)
  flow_edges <- rbind(
    flow_edges,consistency2)
  flow_edges1[[i]] <-flow_edges
  df_sankey <- flow_edges %>%
    rename(source = from,
           target = to,
           value = weight)
  df_sankey1[[i]] <- df_sankey
  
}
flow_edges_2 <-list()
df_sankey2<-list()
for (i in 1:250) {
  ### ----------------------------------------------------
  ### 1. Sample data
  ### ----------------------------------------------------
  df <- try_2[[i]]
  ###
  sankey<-try
  sankey$relative.sum<-ifelse(sankey$relative.sum>0,1,0)
  ### ----------------------------------------------------
  ### 2. Define indicator OTUs for each substrate
  ### Replace with your own indicator-OTU list if needed
  ### ----------------------------------------------------
  indicator_list <- df %>% 
    group_by(substrate) %>%
    summarise(indicators = list(unique(OTU)))
  indicator_df <- try2[,1:2]
  ### ----------------------------------------------------
  ### 3. Build indicator-flow edges
  ### ----------------------------------------------------
  flow_edges <- data.frame()
  
  substrates <- sort(unique(df$substrate))
  
  for (A in substrates) {
    
    # indicator OTUs OF A
    inds_A <- indicator_df$OTU[indicator_df$substrate == A]
    
    # records for indicator OTUs
    df_i <- df %>% filter(OTU %in% inds_A)
    
    # find other substrates
    other_subs <- setdiff(unique(df_i$substrate), A)
    
    for (B in other_subs) {
      
      df_AB <- df_i %>% filter(substrate == B)
      otus_AB <- unique(df_AB$OTU)
      
      # indicator OTUs of B
      inds_B <- indicator_df$OTU[indicator_df$substrate == B]
      
      # keep only OTUs that are indicator of A
      # AND not indicator of B
      exchanged <- setdiff(otus_AB, inds_B)
      
      if (length(exchanged) > 0) {
        flow_edges <- rbind(
          flow_edges,
          data.frame(
            from = A,
            to = B,
            weight = length(exchanged)
          )
        )
      }
    }
  }
  ###consistent indicators
  df2<-df[duplicated(df$OTU),]
  df3<-df[!df$OTU %in% unique(df2$OTU), ]
  consistency<-df3
  consistency$relative.sum<-ifelse(consistency$relative.sum>0,1,0)
  consistency <- consistency %>%
    group_by(substrate) %>%
    summarise(weight=sum(relative.sum))
  
  consistency2 <- data.frame(from= consistency$substrate, to= consistency$substrate, weight= consistency $weight)
  flow_edges <- rbind(
    flow_edges,consistency2)
  flow_edges_2[[i]] <-flow_edges
  df_sankey <- flow_edges %>%
    rename(source = from,
           target = to,
           value = weight)
  df_sankey2[[i]] <- df_sankey
  
}
flow_edges_3 <-list()
df_sankey3<-list()
for (i in 1:250) {
  ### ----------------------------------------------------
  ### 1. Sample data
  ### ----------------------------------------------------
  df <- try_3[[i]]
  ###
  sankey<-try
  sankey$relative.sum<-ifelse(sankey$relative.sum>0,1,0)
  ### ----------------------------------------------------
  ### 2. Define indicator OTUs for each substrate
  ### Replace with your own indicator-OTU list if needed
  ### ----------------------------------------------------
  indicator_list <- df %>% 
    group_by(substrate) %>%
    summarise(indicators = list(unique(OTU)))
  indicator_df <- try2[,1:2]
  ### ----------------------------------------------------
  ### 3. Build indicator-flow edges
  ### ----------------------------------------------------
  flow_edges <- data.frame()
  
  substrates <- sort(unique(df$substrate))
  
  for (A in substrates) {
    
    # indicator OTUs OF A
    inds_A <- indicator_df$OTU[indicator_df$substrate == A]
    
    # records for indicator OTUs
    df_i <- df %>% filter(OTU %in% inds_A)
    
    # find other substrates
    other_subs <- setdiff(unique(df_i$substrate), A)
    
    for (B in other_subs) {
      
      df_AB <- df_i %>% filter(substrate == B)
      otus_AB <- unique(df_AB$OTU)
      
      # indicator OTUs of B
      inds_B <- indicator_df$OTU[indicator_df$substrate == B]
      
      # keep only OTUs that are indicator of A
      # AND not indicator of B
      exchanged <- setdiff(otus_AB, inds_B)
      
      if (length(exchanged) > 0) {
        flow_edges <- rbind(
          flow_edges,
          data.frame(
            from = A,
            to = B,
            weight = length(exchanged)
          )
        )
      }
    }
  }
  ###consistent indicators
  df2<-df[duplicated(df$OTU),]
  df3<-df[!df$OTU %in% unique(df2$OTU), ]
  consistency<-df3
  consistency$relative.sum<-ifelse(consistency$relative.sum>0,1,0)
  consistency <- consistency %>%
    group_by(substrate) %>%
    summarise(weight=sum(relative.sum))
  
  consistency2 <- data.frame(from= consistency$substrate, to= consistency$substrate, weight= consistency $weight)
  flow_edges <- rbind(
    flow_edges,consistency2)
  flow_edges_3[[i]] <-flow_edges
  df_sankey <- flow_edges %>%
    rename(source = from,
           target = to,
           value = weight)
  df_sankey3[[i]] <- df_sankey
  
}
flow_edges4 <-list()
df_sankey4<-list()
for (i in 1:250) {
  ### ----------------------------------------------------
  ### 1. Sample data
  ### ----------------------------------------------------
  df <- try4[[i]]
  ###
  sankey<-try
  sankey$relative.sum<-ifelse(sankey$relative.sum>0,1,0)
  ### ----------------------------------------------------
  ### 2. Define indicator OTUs for each substrate
  ### Replace with your own indicator-OTU list if needed
  ### ----------------------------------------------------
  indicator_list <- df %>% 
    group_by(substrate) %>%
    summarise(indicators = list(unique(OTU)))
  indicator_df <- try2[,1:2]
  ### ----------------------------------------------------
  ### 3. Build indicator-flow edges
  ### ----------------------------------------------------
  flow_edges <- data.frame()
  
  substrates <- sort(unique(df$substrate))
  
  for (A in substrates) {
    
    # indicator OTUs OF A
    inds_A <- indicator_df$OTU[indicator_df$substrate == A]
    
    # records for indicator OTUs
    df_i <- df %>% filter(OTU %in% inds_A)
    
    # find other substrates
    other_subs <- setdiff(unique(df_i$substrate), A)
    
    for (B in other_subs) {
      
      df_AB <- df_i %>% filter(substrate == B)
      otus_AB <- unique(df_AB$OTU)
      
      # indicator OTUs of B
      inds_B <- indicator_df$OTU[indicator_df$substrate == B]
      
      # keep only OTUs that are indicator of A
      # AND not indicator of B
      exchanged <- setdiff(otus_AB, inds_B)
      
      if (length(exchanged) > 0) {
        flow_edges <- rbind(
          flow_edges,
          data.frame(
            from = A,
            to = B,
            weight = length(exchanged)
          )
        )
      }
    }
  }
  ###consistent indicators
  df2<-df[duplicated(df$OTU),]
  df3<-df[!df$OTU %in% unique(df2$OTU), ]
  consistency<-df3
  consistency$relative.sum<-ifelse(consistency$relative.sum>0,1,0)
  consistency <- consistency %>%
    group_by(substrate) %>%
    summarise(weight=sum(relative.sum))
  
  consistency2 <- data.frame(from= consistency$substrate, to= consistency$substrate, weight= consistency $weight)
  flow_edges <- rbind(
    flow_edges,consistency2)
  flow_edges4[[i]] <-flow_edges
  df_sankey <- flow_edges %>%
    rename(source = from,
           target = to,
           value = weight)
  df_sankey4[[i]] <- df_sankey
  
}
######
df_sankey<-c(df_sankey1,df_sankey2,df_sankey3,df_sankey4)
df_sankey<-bind_rows(df_sankey)
df_sankey2 <- df_sankey %>%
  group_by(source, target) %>%
  summarise(value = mean(value))
df_sankey2 <- df_sankey2  %>%
  gather_set_data(1:2) %>%        # <- ggforce helper function
  arrange(x,source,desc(target))
df_sankey2 <- df_sankey2 %>%
  mutate_at(vars(source, target), 
            funs(factor(., levels = c("topsoil", "subsoil",
                                      "leaf litter","bark","wood", 
                                      "lichen","moss","fruit body","leaves","feces", "snow"))))
df_sankey2$self <- df_sankey2$source == df_sankey2$target
df_sankey2$order_rank <- ifelse(df_sankey2$self, 2, 1)
df_sankey2$flow_color <- df_sankey2$source
colors <- c("#882233","#E691C9", "#EDBB00","#E31864",
            "#0005A4", "#99999D", "#0FF001" , "#B35EEE",
            "#009666", "#996677","#00CCFF")
colors <- setNames(colors, c("topsoil", "subsoil",
                             "leaf litter","bark","wood", 
                             "lichen","moss","fruit body","leaves","feces", "snow"))

# Ensure self flows are TRUE/FALSE
df_sankey2$self <- as.logical(df_sankey2$self)

# For ordering on axes: move self flows to last
df_sankey2$flow_color <- factor(
  ifelse(df_sankey2$source == df_sankey2$target, "self", "indicator"),
  levels = c("indicator", "self")  # indicator on top, self at bottom
)
df_sankey2$flow_color <- ifelse(df_sankey2$flow_color=="indicator",as.character(df_sankey2$source),"self")
df_sankey2$flow_color <- factor (df_sankey2$flow_color, levels = c(unique(as.character(df_sankey2$source)),"self"))

df_sankey2$y <- factor (df_sankey2$y, levels = c("leaves","lichen","bark","moss", "snow",
                                                 "fruit body","wood",
                                                 "feces",
                                                 "leaf litter", 
                                                 "topsoil", "subsoil"))

ggplot(df_sankey2, aes(x = x, id = id, split = y, value = value))+
  geom_parallel_sets(aes(fill = flow_color, alpha = self), axis.width = 0.2,
                     n=110, strength = 0.5) +
  geom_parallel_sets_axes(axis.width = 0.01, fill = "gray95",
                          color = "gray80", size = 0.15) +
  geom_parallel_sets_labels(colour = 'gray35', size = 4.5, angle = 0, fontface="bold") +
  scale_alpha_manual(values = c(`TRUE` = 0, `FALSE` = 1),
                     guide = "none") +
  scale_fill_manual(values = colors, guide = "none") +
  scale_x_discrete(limits = c("source", "target"),
                   expand = c(0, 1)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

###
OTU.number <-r2 %>%
  group_by(BestHost) %>%
  summarise(number=length(OTU))

flow_edges<-c(flow_edges1, flow_edges_2, flow_edges_3, flow_edges4)
flow_edges2<-bind_rows(flow_edges)
flow_edges2<-flow_edges2 %>%
  group_by(from, to)%>%
  summarise(weight=mean(weight))

##inflow calculating
flow_edges3<-flow_edges2[flow_edges2$from!=flow_edges2$to,]
flow_edges3<-flow_edges3 %>%
  group_by(to) %>%
  summarise(weight2=sum(weight))
###get the biggest sources
flow_edges4<-flow_edges2[flow_edges2$from!=flow_edges2$to,]
flow_edges4<-flow_edges4 %>%
  group_by(from) %>%
  summarise(weight2=sum(weight))
flow_edges5<-flow_edges2 %>%
  group_by(from) %>%
  summarise(weight2=sum(weight))
all<-merge(flow_edges4,flow_edges5,by = "from")
all$ratio<-all$weight2.x/all$weight2.y
###
flow_edges4$ratio<-flow_edges4$weight2/sum(flow_edges4$weight2)
###get the sources
flow_edges6<-flow_edges2 %>%
  group_by(from) %>%
  summarise(weight2=sum(weight))
##outflow/inflow calculating
flow_edges6<-flow_edges2[flow_edges2$from!=flow_edges2$to,]
outflow<-flow_edges6 %>%
  group_by(from) %>%
  summarise(outflow=sum(weight))
inflow<-flow_edges6 %>%
  group_by(to) %>%
  summarise(inflow=sum(weight))
all<-merge(inflow,outflow,by.x="to",by.y = "from")
all$ratio<-all$outflow/all$inflow

##====mntd and ses.mntd calculation
fruit.lichen$OTU <- row.names(fruit.lichen)
prep_data <- function(x){
  
  LV_list <- vector("list", length(x))
  LW_list <- vector("list", length(x))
  LZ_list <- vector("list", length(x))
  
  for(i in seq_along(x)){
    
    a <- as.data.frame(t(as.matrix(x[[i]])))
    a$OTU <- row.names(a)
    
    fungi.ind <- merge(a, fruit.lichen, by = "OTU", all = TRUE)
    fungi.ind[is.na(fungi.ind)] <- 0
    
    row.names(fungi.ind) <- fungi.ind$OTU
    fungi.ind <- fungi.ind[, -1]
    
    fungi.ind <- as.data.frame(t(fungi.ind))
    fungi.ind$sample_names <- row.names(fungi.ind)
    
    fungi2 <- merge(fungi.ind, metadata, by = "sample_names", all = TRUE)
    row.names(fungi2) <- fungi2$sample_names
    fungi2 <- fungi2[,-1]
    
    fungi2 <- fungi2[!is.na(fungi2$e57bceb79df4a766f8fb9d049798c662ba89dc9e),]
    
    LV <- fungi2[fungi2$site %in% "LV",]
    LV <- LV[, c(colSums(LV[, -((ncol(LV)-1):ncol(LV))]>0) > 2, TRUE, TRUE)]
    LV <- LV[,-ncol(LV)]
    
    LW <- fungi2[fungi2$site %in% "LW",]
    LW <- LW[, c(colSums(LW[, -((ncol(LW)-1):ncol(LW))]>0) > 2, TRUE, TRUE)]
    LW <- LW[,-ncol(LW)]
    
    LZ <- fungi2[fungi2$site %in% "LZ",]
    LZ <- LZ[, c(colSums(LZ[, -((ncol(LZ)-1):ncol(LZ))]>0) > 2, TRUE, TRUE)]
    LZ <- LZ[,-ncol(LZ)]
    
    LV_list[[i]] <- LV
    LW_list[[i]] <- LW
    LZ_list[[i]] <- LZ
  }
  
  return(list(LV = LV_list, LW = LW_list, LZ = LZ_list))
}
res <- prep_data(table1)
res2 <- prep_data(table2)
res3 <- prep_data(table3)
res4 <- prep_data(table4)

######
fungi.tree<-read.tree("fungi.tree2.txt")
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
taxonomy<-read.csv("tax2.csv",header = T,row.names = 1)
###indicator
LV1.ses<- list()
LV.indicator.ses1<-list()
LZ1.ses<- list()
LZ.indicator.ses1<-list()
LW1.ses<- list()
LW.indicator.ses1<-list()

for (i in 1:250) {
  LV1<-res[["LV"]][[i]]
  tax.LV1<-taxonomy[(taxonomy$qseqid) %in% names(LV1),]
  LV1$sample_names<-row.names(LV1)
  table.LV1<-LV1[,-c((ncol(LV1)-1):ncol(LV1))]
  table.LV1<-as.data.frame(t(table.LV1))
  table.LV1$OTU<-row.names(table.LV1)
  table.LV1<-merge(table.LV1,tax.LV1,by="OTU", by.y="qseqid")
  table.LV1$agg<-paste0(table.LV1$kindom,table.LV1$phylum,table.LV1$class,table.LV1$order,table.LV1$family)
  table.LV1<-aggregate(table.LV1[,2:(ncol(table.LV1)-8)],by=list(table.LV1$agg),sum)
  table.LV1<-table.LV1[table.LV1$Group.1 %in% fungi.tree$tip.label,]
  row.names(table.LV1)<-table.LV1$Group.1
  table.LV1<-table.LV1[,-1]
  table.LV1<-table.LV1[,colSums(table.LV1)>0]
  table.LV1<-table.LV1[rowSums(table.LV1)>0,]
  LV.distance<-fungi.distance[row.names(table.LV1),row.names(table.LV1)]
  table.LV1<-t(table.LV1)
  table.LV1<-table.LV1[,colnames(LV.distance)]
  LV1.ses[[i]]<-ses.mpd((table.LV1), LV.distance,null.model ="richness")
  LW1<-res[["LW"]][[i]]
  tax.LW1<-taxonomy[(taxonomy$qseqid) %in% names(LW1),]
  LW1$sample_names<-row.names(LW1)
  table.LW1<-LW1[,-c((ncol(LW1)-1):ncol(LW1))]
  table.LW1<-as.data.frame(t(table.LW1))
  table.LW1$OTU<-row.names(table.LW1)
  table.LW1<-merge(table.LW1,tax.LW1,by="OTU", by.y="qseqid")
  table.LW1$agg<-paste0(table.LW1$kindom,table.LW1$phylum,table.LW1$class,table.LW1$order,table.LW1$family)
  table.LW1<-aggregate(table.LW1[,2:(ncol(table.LW1)-8)],by=list(table.LW1$agg),sum)
  table.LW1<-table.LW1[table.LW1$Group.1 %in% fungi.tree$tip.label,]
  row.names(table.LW1)<-table.LW1$Group.1
  table.LW1<-table.LW1[,-1]
  table.LW1<-table.LW1[,colSums(table.LW1)>0]
  table.LW1<-table.LW1[rowSums(table.LW1)>0,]
  LW.distance<-fungi.distance[row.names(table.LW1),row.names(table.LW1)]
  table.LW1<-t(table.LW1)
  table.LW1<-table.LW1[,colnames(LW.distance)]
  LW1.ses[[i]]<-ses.mpd((table.LW1), LW.distance,null.model ="richness")
  LZ1<-res[["LZ"]][[i]]
  tax.LZ1<-taxonomy[(taxonomy$qseqid) %in% names(LZ1),]
  LZ1$sample_names<-row.names(LZ1)
  table.LZ1<-LZ1[,-c((ncol(LZ1)-1):ncol(LZ1))]
  table.LZ1<-as.data.frame(t(table.LZ1))
  table.LZ1$OTU<-row.names(table.LZ1)
  table.LZ1<-merge(table.LZ1,tax.LZ1,by="OTU", by.y="qseqid")
  table.LZ1$agg<-paste0(table.LZ1$kindom,table.LZ1$phylum,table.LZ1$class,table.LZ1$order,table.LZ1$family)
  table.LZ1<-aggregate(table.LZ1[,2:(ncol(table.LZ1)-8)],by=list(table.LZ1$agg),sum)
  table.LZ1<-table.LZ1[table.LZ1$Group.1 %in% fungi.tree$tip.label,]
  row.names(table.LZ1)<-table.LZ1$Group.1
  table.LZ1<-table.LZ1[,-1]
  table.LZ1<-table.LZ1[,colSums(table.LZ1)>0]
  table.LZ1<-table.LZ1[rowSums(table.LZ1)>0,]
  LZ.distance<-fungi.distance[row.names(table.LZ1),row.names(table.LZ1)]
  table.LZ1<-t(table.LZ1)
  table.LZ1<-table.LZ1[,colnames(LZ.distance)]
  LZ1.ses[[i]]<-ses.mpd((table.LZ1), LZ.distance,null.model ="richness")
  #####indicator
  r <- result1[[i]][result1[[i]]$p.value < 0.05, ]
  LZ1<-LZ1[,-c((ncol(LZ1)-1):ncol(LZ1))]
  LZ1<-as.data.frame(t((LZ1)))
  indicator.table.LZ1<-LZ1[row.names(LZ1) %in% r$OTU,]
  indicator.table.LZ1$OTU<-row.names(indicator.table.LZ1)
  indicator.table.LZ1<-merge(indicator.table.LZ1,tax.LZ1,by.x="OTU",by.y = "qseqid")
  indicator.table.LZ1$agg<-paste0(indicator.table.LZ1$kindom,indicator.table.LZ1$phylum,indicator.table.LZ1$class,indicator.table.LZ1$order,indicator.table.LZ1$family)
  indicator.table.LZ1<-aggregate(indicator.table.LZ1[,2:(ncol(indicator.table.LZ1)-8)],by=list(indicator.table.LZ1$agg),sum)
  row.names(indicator.table.LZ1)<-indicator.table.LZ1$Group.1
  indicator.table.LZ1<-indicator.table.LZ1[,-1]
  indicator.table.LZ1<-  indicator.table.LZ1[,colSums(indicator.table.LZ1)>0]
  indicator.table.LZ1<-  indicator.table.LZ1[rowSums(indicator.table.LZ1)>0,]
  indicator.table.LZ1<-indicator.table.LZ1[row.names(indicator.table.LZ1) %in% fungi.tree$tip.label,]
  LZ.indicator.distance<-fungi.distance[row.names(indicator.table.LZ1),row.names(indicator.table.LZ1)]
  indicator.table.LZ1<-t(indicator.table.LZ1)
  indicator.table.LZ1<-indicator.table.LZ1[,colnames(LZ.indicator.distance)]
  LZ.indicator.ses1[[i]]<-ses.mpd((indicator.table.LZ1), LZ.indicator.distance,null.model ="richness")
  
  #####indicator
  r <- result1[[i]][result1[[i]]$p.value < 0.05, ]
  LW1<-LW1[,-c((ncol(LW1)-1):ncol(LW1))]
  LW1<-as.data.frame(t((LW1)))
  indicator.table.LW1<-LW1[row.names(LW1) %in% r$OTU,]
  indicator.table.LW1$OTU<-row.names(indicator.table.LW1)
  indicator.table.LW1<-merge(indicator.table.LW1,tax.LW1,by.x="OTU",by.y = "qseqid")
  indicator.table.LW1$agg<-paste0(indicator.table.LW1$kindom,indicator.table.LW1$phylum,indicator.table.LW1$class,indicator.table.LW1$order,indicator.table.LW1$family)
  indicator.table.LW1<-aggregate(indicator.table.LW1[,2:(ncol(indicator.table.LW1)-8)],by=list(indicator.table.LW1$agg),sum)
  row.names(indicator.table.LW1)<-indicator.table.LW1$Group.1
  indicator.table.LW1<-indicator.table.LW1[,-1]
  indicator.table.LW1<-  indicator.table.LW1[,colSums(indicator.table.LW1)>0]
  indicator.table.LW1<-  indicator.table.LW1[rowSums(indicator.table.LW1)>0,]
  indicator.table.LW1<-indicator.table.LW1[row.names(indicator.table.LW1) %in% fungi.tree$tip.label,]
  LW.indicator.distance<-fungi.distance[row.names(indicator.table.LW1),row.names(indicator.table.LW1)]
  indicator.table.LW1<-t(indicator.table.LW1)
  indicator.table.LW1<-indicator.table.LW1[,colnames(LW.indicator.distance)]
  LW.indicator.ses1[[i]]<-ses.mpd((indicator.table.LW1), LW.indicator.distance,null.model ="richness")
  #####indicator
  r <- result1[[i]][result1[[i]]$p.value < 0.05, ]
  LV1<-LV1[,-c((ncol(LV1)-1):ncol(LV1))]
  LV1<-as.data.frame(t((LV1)))
  indicator.table.LV1<-LV1[row.names(LV1) %in% r$OTU,]
  indicator.table.LV1$OTU<-row.names(indicator.table.LV1)
  indicator.table.LV1<-merge(indicator.table.LV1,tax.LV1,by.x="OTU",by.y = "qseqid")
  indicator.table.LV1$agg<-paste0(indicator.table.LV1$kindom,indicator.table.LV1$phylum,indicator.table.LV1$class,indicator.table.LV1$order,indicator.table.LV1$family)
  indicator.table.LV1<-aggregate(indicator.table.LV1[,2:(ncol(indicator.table.LV1)-8)],by=list(indicator.table.LV1$agg),sum)
  row.names(indicator.table.LV1)<-indicator.table.LV1$Group.1
  indicator.table.LV1<-indicator.table.LV1[,-1]
  indicator.table.LV1<-  indicator.table.LV1[,colSums(indicator.table.LV1)>0]
  indicator.table.LV1<-  indicator.table.LV1[rowSums(indicator.table.LV1)>0,]
  indicator.table.LV1<-indicator.table.LV1[row.names(indicator.table.LV1) %in% fungi.tree$tip.label,]
  LV.indicator.distance<-fungi.distance[row.names(indicator.table.LV1),row.names(indicator.table.LV1)]
  indicator.table.LV1<-t(indicator.table.LV1)
  indicator.table.LV1<-indicator.table.LV1[,colnames(LV.indicator.distance)]
  LV.indicator.ses1[[i]]<-ses.mpd((indicator.table.LV1), LV.indicator.distance,null.model ="richness")
}

LV2.ses<- list()
LV.indicator.ses2<-list()
LZ2.ses<- list()
LZ.indicator.ses2<-list()
LW2.ses<- list()
LW.indicator.ses2<-list()

for (i in 1:250) {
  LV1<-res2[["LV"]][[i]]
  tax.LV1<-taxonomy[(taxonomy$qseqid) %in% names(LV1),]
  LV1$sample_names<-row.names(LV1)
  table.LV1<-LV1[,-c((ncol(LV1)-1):ncol(LV1))]
  table.LV1<-as.data.frame(t(table.LV1))
  table.LV1$OTU<-row.names(table.LV1)
  table.LV1<-merge(table.LV1,tax.LV1,by="OTU", by.y="qseqid")
  table.LV1$agg<-paste0(table.LV1$kindom,table.LV1$phylum,table.LV1$class,table.LV1$order,table.LV1$family)
  table.LV1<-aggregate(table.LV1[,2:(ncol(table.LV1)-8)],by=list(table.LV1$agg),sum)
  table.LV1<-table.LV1[table.LV1$Group.1 %in% fungi.tree$tip.label,]
  row.names(table.LV1)<-table.LV1$Group.1
  table.LV1<-table.LV1[,-1]
  table.LV1<-table.LV1[,colSums(table.LV1)>0]
  table.LV1<-table.LV1[rowSums(table.LV1)>0,]
  LV.distance<-fungi.distance[row.names(table.LV1),row.names(table.LV1)]
  table.LV1<-t(table.LV1)
  table.LV1<-table.LV1[,colnames(LV.distance)]
  LV2.ses[[i]]<-ses.mpd((table.LV1), LV.distance,null.model ="richness")
  LW1<-res2[["LW"]][[i]]
  tax.LW1<-taxonomy[(taxonomy$qseqid) %in% names(LW1),]
  LW1$sample_names<-row.names(LW1)
  table.LW1<-LW1[,-c((ncol(LW1)-1):ncol(LW1))]
  table.LW1<-as.data.frame(t(table.LW1))
  table.LW1$OTU<-row.names(table.LW1)
  table.LW1<-merge(table.LW1,tax.LW1,by="OTU", by.y="qseqid")
  table.LW1$agg<-paste0(table.LW1$kindom,table.LW1$phylum,table.LW1$class,table.LW1$order,table.LW1$family)
  table.LW1<-aggregate(table.LW1[,2:(ncol(table.LW1)-8)],by=list(table.LW1$agg),sum)
  table.LW1<-table.LW1[table.LW1$Group.1 %in% fungi.tree$tip.label,]
  row.names(table.LW1)<-table.LW1$Group.1
  table.LW1<-table.LW1[,-1]
  table.LW1<-table.LW1[,colSums(table.LW1)>0]
  table.LW1<-table.LW1[rowSums(table.LW1)>0,]
  LW.distance<-fungi.distance[row.names(table.LW1),row.names(table.LW1)]
  table.LW1<-t(table.LW1)
  table.LW1<-table.LW1[,colnames(LW.distance)]
  LW2.ses[[i]]<-ses.mpd((table.LW1), LW.distance,null.model ="richness")
  LZ1<-res2[["LZ"]][[i]]
  tax.LZ1<-taxonomy[(taxonomy$qseqid) %in% names(LZ1),]
  LZ1$sample_names<-row.names(LZ1)
  table.LZ1<-LZ1[,-c((ncol(LZ1)-1):ncol(LZ1))]
  table.LZ1<-as.data.frame(t(table.LZ1))
  table.LZ1$OTU<-row.names(table.LZ1)
  table.LZ1<-merge(table.LZ1,tax.LZ1,by="OTU", by.y="qseqid")
  table.LZ1$agg<-paste0(table.LZ1$kindom,table.LZ1$phylum,table.LZ1$class,table.LZ1$order,table.LZ1$family)
  table.LZ1<-aggregate(table.LZ1[,2:(ncol(table.LZ1)-8)],by=list(table.LZ1$agg),sum)
  table.LZ1<-table.LZ1[table.LZ1$Group.1 %in% fungi.tree$tip.label,]
  row.names(table.LZ1)<-table.LZ1$Group.1
  table.LZ1<-table.LZ1[,-1]
  table.LZ1<-table.LZ1[,colSums(table.LZ1)>0]
  table.LZ1<-table.LZ1[rowSums(table.LZ1)>0,]
  LZ.distance<-fungi.distance[row.names(table.LZ1),row.names(table.LZ1)]
  table.LZ1<-t(table.LZ1)
  table.LZ1<-table.LZ1[,colnames(LZ.distance)]
  LZ2.ses[[i]]<-ses.mpd((table.LZ1), LZ.distance,null.model ="richness")
  #####indicator
  r <- result2[[i]][result2[[i]]$p.value < 0.05, ]
  LZ1<-LZ1[,-c((ncol(LZ1)-1):ncol(LZ1))]
  LZ1<-as.data.frame(t((LZ1)))
  indicator.table.LZ1<-LZ1[row.names(LZ1) %in% r$OTU,]
  indicator.table.LZ1$OTU<-row.names(indicator.table.LZ1)
  indicator.table.LZ1<-merge(indicator.table.LZ1,tax.LZ1,by.x="OTU",by.y = "qseqid")
  indicator.table.LZ1$agg<-paste0(indicator.table.LZ1$kindom,indicator.table.LZ1$phylum,indicator.table.LZ1$class,indicator.table.LZ1$order,indicator.table.LZ1$family)
  indicator.table.LZ1<-aggregate(indicator.table.LZ1[,2:(ncol(indicator.table.LZ1)-8)],by=list(indicator.table.LZ1$agg),sum)
  row.names(indicator.table.LZ1)<-indicator.table.LZ1$Group.1
  indicator.table.LZ1<-indicator.table.LZ1[,-1]
  indicator.table.LZ1<-  indicator.table.LZ1[,colSums(indicator.table.LZ1)>0]
  indicator.table.LZ1<-  indicator.table.LZ1[rowSums(indicator.table.LZ1)>0,]
  indicator.table.LZ1<-indicator.table.LZ1[row.names(indicator.table.LZ1) %in% fungi.tree$tip.label,]
  LZ.indicator.distance<-fungi.distance[row.names(indicator.table.LZ1),row.names(indicator.table.LZ1)]
  indicator.table.LZ1<-t(indicator.table.LZ1)
  indicator.table.LZ1<-indicator.table.LZ1[,colnames(LZ.indicator.distance)]
  LZ.indicator.ses2[[i]]<-ses.mpd((indicator.table.LZ1), LZ.indicator.distance,null.model ="richness")
  
  #####indicator
  r <- result2[[i]][result2[[i]]$p.value < 0.05, ]
  LW1<-LW1[,-c((ncol(LW1)-1):ncol(LW1))]
  LW1<-as.data.frame(t((LW1)))
  indicator.table.LW1<-LW1[row.names(LW1) %in% r$OTU,]
  indicator.table.LW1$OTU<-row.names(indicator.table.LW1)
  indicator.table.LW1<-merge(indicator.table.LW1,tax.LW1,by.x="OTU",by.y = "qseqid")
  indicator.table.LW1$agg<-paste0(indicator.table.LW1$kindom,indicator.table.LW1$phylum,indicator.table.LW1$class,indicator.table.LW1$order,indicator.table.LW1$family)
  indicator.table.LW1<-aggregate(indicator.table.LW1[,2:(ncol(indicator.table.LW1)-8)],by=list(indicator.table.LW1$agg),sum)
  row.names(indicator.table.LW1)<-indicator.table.LW1$Group.1
  indicator.table.LW1<-indicator.table.LW1[,-1]
  indicator.table.LW1<-  indicator.table.LW1[,colSums(indicator.table.LW1)>0]
  indicator.table.LW1<-  indicator.table.LW1[rowSums(indicator.table.LW1)>0,]
  indicator.table.LW1<-indicator.table.LW1[row.names(indicator.table.LW1) %in% fungi.tree$tip.label,]
  LW.indicator.distance<-fungi.distance[row.names(indicator.table.LW1),row.names(indicator.table.LW1)]
  indicator.table.LW1<-t(indicator.table.LW1)
  indicator.table.LW1<-indicator.table.LW1[,colnames(LW.indicator.distance)]
  LW.indicator.ses2[[i]]<-ses.mpd((indicator.table.LW1), LW.indicator.distance,null.model ="richness")
  #####indicator
  r <- result2[[i]][result2[[i]]$p.value < 0.05, ]
  LV1<-LV1[,-c((ncol(LV1)-1):ncol(LV1))]
  LV1<-as.data.frame(t((LV1)))
  indicator.table.LV1<-LV1[row.names(LV1) %in% r$OTU,]
  indicator.table.LV1$OTU<-row.names(indicator.table.LV1)
  indicator.table.LV1<-merge(indicator.table.LV1,tax.LV1,by.x="OTU",by.y = "qseqid")
  indicator.table.LV1$agg<-paste0(indicator.table.LV1$kindom,indicator.table.LV1$phylum,indicator.table.LV1$class,indicator.table.LV1$order,indicator.table.LV1$family)
  indicator.table.LV1<-aggregate(indicator.table.LV1[,2:(ncol(indicator.table.LV1)-8)],by=list(indicator.table.LV1$agg),sum)
  row.names(indicator.table.LV1)<-indicator.table.LV1$Group.1
  indicator.table.LV1<-indicator.table.LV1[,-1]
  indicator.table.LV1<-  indicator.table.LV1[,colSums(indicator.table.LV1)>0]
  indicator.table.LV1<-  indicator.table.LV1[rowSums(indicator.table.LV1)>0,]
 indicator.table.LV1<-indicator.table.LV1[row.names(indicator.table.LV1) %in% fungi.tree$tip.label,]
  LV.indicator.distance<-fungi.distance[row.names(indicator.table.LV1),row.names(indicator.table.LV1)]
  indicator.table.LV1<-t(indicator.table.LV1)
  indicator.table.LV1<-indicator.table.LV1[,colnames(LV.indicator.distance)]
  LV.indicator.ses2[[i]]<-ses.mpd((indicator.table.LV1), LV.indicator.distance,null.model ="richness")
}

LV3.ses<- list()
LV.indicator.ses3<-list()
LZ3.ses<- list()
LZ.indicator.ses3<-list()
LW3.ses<- list()
LW.indicator.ses3<-list()

for (i in 1:250) {
  LV1<-res3[["LV"]][[i]]
  tax.LV1<-taxonomy[(taxonomy$qseqid) %in% names(LV1),]
  LV1$sample_names<-row.names(LV1)
  table.LV1<-LV1[,-c((ncol(LV1)-1):ncol(LV1))]
  table.LV1<-as.data.frame(t(table.LV1))
  table.LV1$OTU<-row.names(table.LV1)
  table.LV1<-merge(table.LV1,tax.LV1,by="OTU", by.y="qseqid")
  table.LV1$agg<-paste0(table.LV1$kindom,table.LV1$phylum,table.LV1$class,table.LV1$order,table.LV1$family)
  table.LV1<-aggregate(table.LV1[,2:(ncol(table.LV1)-8)],by=list(table.LV1$agg),sum)
  table.LV1<-table.LV1[table.LV1$Group.1 %in% fungi.tree$tip.label,]
  row.names(table.LV1)<-table.LV1$Group.1
  table.LV1<-table.LV1[,-1]
  table.LV1<-table.LV1[,colSums(table.LV1)>0]
  table.LV1<-table.LV1[rowSums(table.LV1)>0,]
  LV.distance<-fungi.distance[row.names(table.LV1),row.names(table.LV1)]
  table.LV1<-t(table.LV1)
  table.LV1<-table.LV1[,colnames(LV.distance)]
  LV3.ses[[i]]<-ses.mpd((table.LV1), LV.distance,null.model ="richness")
  LW1<-res3[["LW"]][[i]]
  tax.LW1<-taxonomy[(taxonomy$qseqid) %in% names(LW1),]
  LW1$sample_names<-row.names(LW1)
  table.LW1<-LW1[,-c((ncol(LW1)-1):ncol(LW1))]
  table.LW1<-as.data.frame(t(table.LW1))
  table.LW1$OTU<-row.names(table.LW1)
  table.LW1<-merge(table.LW1,tax.LW1,by="OTU", by.y="qseqid")
  table.LW1$agg<-paste0(table.LW1$kindom,table.LW1$phylum,table.LW1$class,table.LW1$order,table.LW1$family)
  table.LW1<-aggregate(table.LW1[,2:(ncol(table.LW1)-8)],by=list(table.LW1$agg),sum)
  table.LW1<-table.LW1[table.LW1$Group.1 %in% fungi.tree$tip.label,]
  row.names(table.LW1)<-table.LW1$Group.1
  table.LW1<-table.LW1[,-1]
  table.LW1<-table.LW1[,colSums(table.LW1)>0]
  table.LW1<-table.LW1[rowSums(table.LW1)>0,]
  LW.distance<-fungi.distance[row.names(table.LW1),row.names(table.LW1)]
  table.LW1<-t(table.LW1)
  table.LW1<-table.LW1[,colnames(LW.distance)]
  LW3.ses[[i]]<-ses.mpd((table.LW1), LW.distance,null.model ="richness")
  LZ1<-res3[["LZ"]][[i]]
  tax.LZ1<-taxonomy[(taxonomy$qseqid) %in% names(LZ1),]
  LZ1$sample_names<-row.names(LZ1)
  table.LZ1<-LZ1[,-c((ncol(LZ1)-1):ncol(LZ1))]
  table.LZ1<-as.data.frame(t(table.LZ1))
  table.LZ1$OTU<-row.names(table.LZ1)
  table.LZ1<-merge(table.LZ1,tax.LZ1,by="OTU", by.y="qseqid")
  table.LZ1$agg<-paste0(table.LZ1$kindom,table.LZ1$phylum,table.LZ1$class,table.LZ1$order,table.LZ1$family)
  table.LZ1<-aggregate(table.LZ1[,2:(ncol(table.LZ1)-8)],by=list(table.LZ1$agg),sum)
  table.LZ1<-table.LZ1[table.LZ1$Group.1 %in% fungi.tree$tip.label,]
  row.names(table.LZ1)<-table.LZ1$Group.1
  table.LZ1<-table.LZ1[,-1]
  table.LZ1<-table.LZ1[,colSums(table.LZ1)>0]
  table.LZ1<-table.LZ1[rowSums(table.LZ1)>0,]
  LZ.distance<-fungi.distance[row.names(table.LZ1),row.names(table.LZ1)]
  table.LZ1<-t(table.LZ1)
  table.LZ1<-table.LZ1[,colnames(LZ.distance)]
  LZ3.ses[[i]]<-ses.mpd((table.LZ1), LZ.distance,null.model ="richness")
  #####indicator
  r <- result3[[i]][result3[[i]]$p.value < 0.05, ]
  LZ1<-LZ1[,-c((ncol(LZ1)-1):ncol(LZ1))]
  LZ1<-as.data.frame(t((LZ1)))
  indicator.table.LZ1<-LZ1[row.names(LZ1) %in% r$OTU,]
  indicator.table.LZ1$OTU<- row.names(indicator.table.LZ1)
  indicator.table.LZ1<-merge(indicator.table.LZ1,tax.LZ1,by.x="OTU",by.y = "qseqid")
  indicator.table.LZ1$agg<-paste0(indicator.table.LZ1$kindom,indicator.table.LZ1$phylum,indicator.table.LZ1$class,indicator.table.LZ1$order,indicator.table.LZ1$family)
  indicator.table.LZ1<-aggregate(indicator.table.LZ1[,2:(ncol(indicator.table.LZ1)-8)],by=list(indicator.table.LZ1$agg),sum)
  row.names(indicator.table.LZ1)<-indicator.table.LZ1$Group.1
  indicator.table.LZ1<-indicator.table.LZ1[,-1]
  indicator.table.LZ1<-  indicator.table.LZ1[,colSums(indicator.table.LZ1)>0]
  indicator.table.LZ1<-  indicator.table.LZ1[rowSums(indicator.table.LZ1)>0,]
  indicator.table.LZ1<-indicator.table.LZ1[row.names(indicator.table.LZ1) %in% fungi.tree$tip.label,]
  LZ.indicator.distance<-fungi.distance[row.names(indicator.table.LZ1),row.names(indicator.table.LZ1)]
  indicator.table.LZ1<-t(indicator.table.LZ1)
  indicator.table.LZ1<-indicator.table.LZ1[,colnames(LZ.indicator.distance)]
  LZ.indicator.ses3[[i]]<-ses.mpd((indicator.table.LZ1), LZ.indicator.distance,null.model ="richness")
  
  #####indicator
  r <- result3[[i]][result3[[i]]$p.value < 0.05, ]
  LW1<-LW1[,-c((ncol(LW1)-1):ncol(LW1))]
  LW1<-as.data.frame(t((LW1)))
  indicator.table.LW1<-LW1[row.names(LW1) %in% r$OTU,]
  indicator.table.LW1$OTU<- row.names(indicator.table.LW1)
  indicator.table.LW1<-merge(indicator.table.LW1,tax.LW1,by.x="OTU",by.y = "qseqid")
  indicator.table.LW1$agg<-paste0(indicator.table.LW1$kindom,indicator.table.LW1$phylum,indicator.table.LW1$class,indicator.table.LW1$order,indicator.table.LW1$family)
  indicator.table.LW1<-aggregate(indicator.table.LW1[,2:(ncol(indicator.table.LW1)-8)],by=list(indicator.table.LW1$agg),sum)
  row.names(indicator.table.LW1)<-indicator.table.LW1$Group.1
  indicator.table.LW1<-indicator.table.LW1[,-1]
  indicator.table.LW1<-  indicator.table.LW1[,colSums(indicator.table.LW1)>0]
  indicator.table.LW1<-  indicator.table.LW1[rowSums(indicator.table.LW1)>0,]
  indicator.table.LW1<-indicator.table.LW1[row.names(indicator.table.LW1) %in% fungi.tree$tip.label,]
  LW.indicator.distance<-fungi.distance[row.names(indicator.table.LW1),row.names(indicator.table.LW1)]
  indicator.table.LW1<-t(indicator.table.LW1)
  indicator.table.LW1<-indicator.table.LW1[,colnames(LW.indicator.distance)]
  LW.indicator.ses3[[i]]<-ses.mpd((indicator.table.LW1), LW.indicator.distance,null.model ="richness")
  #####indicator
  r <- result3[[i]][result3[[i]]$p.value < 0.05, ]
  LV1<-LV1[,-c((ncol(LV1)-1):ncol(LV1))]
  LV1<-as.data.frame(t((LV1)))
  indicator.table.LV1<-LV1[row.names(LV1) %in% r$OTU,]
  indicator.table.LV1$OTU<- row.names(indicator.table.LV1)
  indicator.table.LV1<-merge(indicator.table.LV1,tax.LV1,by.x="OTU",by.y = "qseqid")
  indicator.table.LV1$agg<-paste0(indicator.table.LV1$kindom,indicator.table.LV1$phylum,indicator.table.LV1$class,indicator.table.LV1$order,indicator.table.LV1$family)
  indicator.table.LV1<-aggregate(indicator.table.LV1[,2:(ncol(indicator.table.LV1)-8)],by=list(indicator.table.LV1$agg),sum)
  row.names(indicator.table.LV1)<-indicator.table.LV1$Group.1
  indicator.table.LV1<-indicator.table.LV1[,-1]
  indicator.table.LV1<-  indicator.table.LV1[,colSums(indicator.table.LV1)>0]
  indicator.table.LV1<-  indicator.table.LV1[rowSums(indicator.table.LV1)>0,]
  indicator.table.LV1<-indicator.table.LV1[row.names(indicator.table.LV1) %in% fungi.tree$tip.label,]
  LV.indicator.distance<-fungi.distance[row.names(indicator.table.LV1),row.names(indicator.table.LV1)]
  indicator.table.LV1<-t(indicator.table.LV1)
  indicator.table.LV1<-indicator.table.LV1[,colnames(LV.indicator.distance)]
  LV.indicator.ses3[[i]]<-ses.mpd((indicator.table.LV1), LV.indicator.distance,null.model ="richness")
}

LV4.ses<- list()
LV.indicator.ses4<-list()
LZ4.ses<- list()
LZ.indicator.ses4<-list()
LW4.ses<- list()
LW.indicator.ses4<-list()

for (i in 1:250) {
  LV1<-res4[["LV"]][[i]]
  tax.LV1<-taxonomy[(taxonomy$qseqid) %in% names(LV1),]
  LV1$sample_names<-row.names(LV1)
  table.LV1<-LV1[,-c((ncol(LV1)-1):ncol(LV1))]
  table.LV1<-as.data.frame(t(table.LV1))
  table.LV1$OTU<-row.names(table.LV1)
  table.LV1<-merge(table.LV1,tax.LV1,by="OTU", by.y="qseqid")
  table.LV1$agg<-paste0(table.LV1$kindom,table.LV1$phylum,table.LV1$class,table.LV1$order,table.LV1$family)
  table.LV1<-aggregate(table.LV1[,2:(ncol(table.LV1)-8)],by=list(table.LV1$agg),sum)
  table.LV1<-table.LV1[table.LV1$Group.1 %in% fungi.tree$tip.label,]
  row.names(table.LV1)<-table.LV1$Group.1
  table.LV1<-table.LV1[,-1]
  table.LV1<-table.LV1[,colSums(table.LV1)>0]
  table.LV1<-table.LV1[rowSums(table.LV1)>0,]
  LV.distance<-fungi.distance[row.names(table.LV1),row.names(table.LV1)]
  table.LV1<-t(table.LV1)
  table.LV1<-table.LV1[,colnames(LV.distance)]
  LV4.ses[[i]]<-ses.mpd((table.LV1), LV.distance,null.model ="richness")
  LW1<-res4[["LW"]][[i]]
  tax.LW1<-taxonomy[(taxonomy$qseqid) %in% names(LW1),]
  LW1$sample_names<-row.names(LW1)
  table.LW1<-LW1[,-c((ncol(LW1)-1):ncol(LW1))]
  table.LW1<-as.data.frame(t(table.LW1))
  table.LW1$OTU<-row.names(table.LW1)
  table.LW1<-merge(table.LW1,tax.LW1,by="OTU", by.y="qseqid")
  table.LW1$agg<-paste0(table.LW1$kindom,table.LW1$phylum,table.LW1$class,table.LW1$order,table.LW1$family)
  table.LW1<-aggregate(table.LW1[,2:(ncol(table.LW1)-8)],by=list(table.LW1$agg),sum)
  table.LW1<-table.LW1[table.LW1$Group.1 %in% fungi.tree$tip.label,]
  row.names(table.LW1)<-table.LW1$Group.1
  table.LW1<-table.LW1[,-1]
  table.LW1<-table.LW1[,colSums(table.LW1)>0]
  table.LW1<-table.LW1[rowSums(table.LW1)>0,]
  LW.distance<-fungi.distance[row.names(table.LW1),row.names(table.LW1)]
  table.LW1<-t(table.LW1)
  table.LW1<-table.LW1[,colnames(LW.distance)]
  LW4.ses[[i]]<-ses.mpd((table.LW1), LW.distance,null.model ="richness")
  LZ1<-res4[["LZ"]][[i]]
  tax.LZ1<-taxonomy[(taxonomy$qseqid) %in% names(LZ1),]
  LZ1$sample_names<-row.names(LZ1)
  table.LZ1<-LZ1[,-c((ncol(LZ1)-1):ncol(LZ1))]
  table.LZ1<-as.data.frame(t(table.LZ1))
  table.LZ1$OTU<-row.names(table.LZ1)
  table.LZ1<-merge(table.LZ1,tax.LZ1,by="OTU", by.y="qseqid")
  table.LZ1$agg<-paste0(table.LZ1$kindom,table.LZ1$phylum,table.LZ1$class,table.LZ1$order,table.LZ1$family)
  table.LZ1<-aggregate(table.LZ1[,2:(ncol(table.LZ1)-8)],by=list(table.LZ1$agg),sum)
  table.LZ1<-table.LZ1[table.LZ1$Group.1 %in% fungi.tree$tip.label,]
  row.names(table.LZ1)<-table.LZ1$Group.1
  table.LZ1<-table.LZ1[,-1]
  table.LZ1<-table.LZ1[,colSums(table.LZ1)>0]
  table.LZ1<-table.LZ1[rowSums(table.LZ1)>0,]
  LZ.distance<-fungi.distance[row.names(table.LZ1),row.names(table.LZ1)]
  table.LZ1<-t(table.LZ1)
  table.LZ1<-table.LZ1[,colnames(LZ.distance)]
  LZ4.ses[[i]]<-ses.mpd((table.LZ1), LZ.distance,null.model ="richness")
  #####indicator
  r <- result4[[i]][result4[[i]]$p.value < 0.05, ]
  LZ1<-LZ1[,-c((ncol(LZ1)-1):ncol(LZ1))]
  LZ1<-as.data.frame(t((LZ1)))
  indicator.table.LZ1<-LZ1[row.names(LZ1) %in% r$OTU,]
  indicator.table.LZ1$OTU<-row.names(indicator.table.LZ1)
  indicator.table.LZ1<-merge(indicator.table.LZ1,tax.LZ1,by.x="OTU",by.y = "qseqid")
  indicator.table.LZ1$agg<-paste0(indicator.table.LZ1$kindom,indicator.table.LZ1$phylum,indicator.table.LZ1$class,indicator.table.LZ1$order,indicator.table.LZ1$family)
  indicator.table.LZ1<-aggregate(indicator.table.LZ1[,2:(ncol(indicator.table.LZ1)-8)],by=list(indicator.table.LZ1$agg),sum)
  row.names(indicator.table.LZ1)<-indicator.table.LZ1$Group.1
  indicator.table.LZ1<-indicator.table.LZ1[,-1]
  indicator.table.LZ1<-  indicator.table.LZ1[,colSums(indicator.table.LZ1)>0]
  indicator.table.LZ1<-  indicator.table.LZ1[rowSums(indicator.table.LZ1)>0,]
  indicator.table.LZ1<-indicator.table.LZ1[row.names(indicator.table.LZ1) %in% fungi.tree$tip.label,]
  LZ.indicator.distance<-fungi.distance[row.names(indicator.table.LZ1),row.names(indicator.table.LZ1)]
  indicator.table.LZ1<-t(indicator.table.LZ1)
  indicator.table.LZ1<-indicator.table.LZ1[,colnames(LZ.indicator.distance)]
  LZ.indicator.ses4[[i]]<-ses.mpd((indicator.table.LZ1), LZ.indicator.distance,null.model ="richness")
  
  #####indicator
  r <- result4[[i]][result4[[i]]$p.value < 0.05, ]
  LW1<-LW1[,-c((ncol(LW1)-1):ncol(LW1))]
  LW1<-as.data.frame(t((LW1)))
  indicator.table.LW1<-LW1[row.names(LW1) %in% r$OTU,]
  indicator.table.LW1$OTU<-row.names(indicator.table.LW1)
  indicator.table.LW1<-merge(indicator.table.LW1,tax.LW1,by.x="OTU",by.y = "qseqid")
  indicator.table.LW1$agg<-paste0(indicator.table.LW1$kindom,indicator.table.LW1$phylum,indicator.table.LW1$class,indicator.table.LW1$order,indicator.table.LW1$family)
  indicator.table.LW1<-aggregate(indicator.table.LW1[,2:(ncol(indicator.table.LW1)-8)],by=list(indicator.table.LW1$agg),sum)
  row.names(indicator.table.LW1)<-indicator.table.LW1$Group.1
  indicator.table.LW1<-indicator.table.LW1[,-1]
  indicator.table.LW1<-  indicator.table.LW1[,colSums(indicator.table.LW1)>0]
  indicator.table.LW1<-  indicator.table.LW1[rowSums(indicator.table.LW1)>0,]
  indicator.table.LW1<-indicator.table.LW1[row.names(indicator.table.LW1) %in% fungi.tree$tip.label,]
  LW.indicator.distance<-fungi.distance[row.names(indicator.table.LW1),row.names(indicator.table.LW1)]
  indicator.table.LW1<-t(indicator.table.LW1)
  indicator.table.LW1<-indicator.table.LW1[,colnames(LW.indicator.distance)]
  LW.indicator.ses4[[i]]<-ses.mpd((indicator.table.LW1), LW.indicator.distance,null.model ="richness")
  #####indicator
  r <- result4[[i]][result4[[i]]$p.value < 0.05, ]
  LV1<-LV1[,-c((ncol(LV1)-1):ncol(LV1))]
  LV1<-as.data.frame(t((LV1)))
  indicator.table.LV1<-LV1[row.names(LV1) %in% r$OTU,]
  indicator.table.LV1$OTU<-row.names(indicator.table.LV1)
  indicator.table.LV1<-merge(indicator.table.LV1,tax.LV1,by.x="OTU",by.y = "qseqid")
  indicator.table.LV1$agg<-paste0(indicator.table.LV1$kindom,indicator.table.LV1$phylum,indicator.table.LV1$class,indicator.table.LV1$order,indicator.table.LV1$family)
  indicator.table.LV1<-aggregate(indicator.table.LV1[,2:(ncol(indicator.table.LV1)-8)],by=list(indicator.table.LV1$agg),sum)
  row.names(indicator.table.LV1)<-indicator.table.LV1$Group.1
  indicator.table.LV1<-indicator.table.LV1[,-1]
  indicator.table.LV1<-  indicator.table.LV1[,colSums(indicator.table.LV1)>0]
  indicator.table.LV1<-  indicator.table.LV1[rowSums(indicator.table.LV1)>0,]
  indicator.table.LV1<-indicator.table.LV1[row.names(indicator.table.LV1) %in% fungi.tree$tip.label,]
  LV.indicator.distance<-fungi.distance[row.names(indicator.table.LV1),row.names(indicator.table.LV1)]
  indicator.table.LV1<-t(indicator.table.LV1)
  indicator.table.LV1<-indicator.table.LV1[,colnames(LV.indicator.distance)]
  LV.indicator.ses4[[i]]<-ses.mpd((indicator.table.LV1), LV.indicator.distance,null.model ="richness")
}

#regression between richness and mpd
LZ<- c(LZ1.ses,LZ2.ses,LZ3.ses,LZ4.ses)
for (i in 1:1000) {
 LZ[[i]]$sample_names <- row.names(LZ[[i]])
}
LZ<-bind_rows(LZ)
LZ <- LZ %>%
  group_by(sample_names) %>%
  summarise(mpd.obs2 = mean(mpd.obs, na.rm = TRUE), mpd.obs.z2=mean(mpd.obs.z, na.rm = TRUE))

LW<- c(LW1.ses,LW2.ses,LW3.ses,LW4.ses)
for (i in 1:1000) {
  LW[[i]]$sample_names <- row.names(LW[[i]])
}
LW<-bind_rows(LW)
LW <- LW %>%
  group_by(sample_names) %>%
  summarise(mpd.obs2 = mean(mpd.obs, na.rm = TRUE), mpd.obs.z2=mean(mpd.obs.z, na.rm = TRUE))


LV<- c(LV1.ses,LV2.ses,LV3.ses,LV4.ses)
for (i in 1:1000) {
  LV[[i]]$sample_names <- row.names(LV[[i]])
}
LV<-bind_rows(LV)

LV <- LV %>%
  group_by(sample_names) %>%
  summarise(mpd.obs2 = mean(mpd.obs, na.rm = TRUE), mpd.obs.z2=mean(mpd.obs.z, na.rm = TRUE))

ses<-rbind(LW,LV,LZ)
try<-merge(metadata,ses,by= "sample_names",all=T)
try<-try[!is.na(try$mpd.obs2),]

try2<-try
try2$substrate<-factor(try2$substrate,levels=c("subsoil","topsoil","wood","fruit body","feces",
                                               "lichen","moss","bark","leaves","leaf litter",
                                               "snow"))
mod01<-lmerTest::lmer(mpd.obs.z2 ~ substrate+(1|site2),data = try2)
b<-avg_comparisons(mod01, variables = list(substrate = "pairwise")) 
performance::performance(mod01)
substrate.ses<-parameters::model_parameters(mod01)

write.csv(substrate.ses,"substrate.ses.parameter.csv")
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group[letter$Group == "fruitbody"] <- "fruit body"
letter$Group[letter$Group == "leaflitter"] <- "leaf litter"
difference <- try2 %>%
  group_by(substrate) %>%
  summarise(max=mean(mpd.obs.z2 ))

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
mod01<-lmerTest::lmer(mpd.obs2 ~ substrate+(1|site2),data = try2)
performance::performance(mod01)
substrate.obs<-parameters::model_parameters(mod01)
write.csv(substrate.obs,"substrate.obs.parameter.csv")
b<-avg_comparisons(mod01, variables = list(substrate = "pairwise")) 
library(rcompanion)
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group[letter$Group == "fruitbody"] <- "fruit body"
letter$Group[letter$Group == "leaflitter"] <- "leaf litter"
difference <- try %>%
  group_by(substrate) %>%
  summarise(max=mean(mpd.obs2))

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
