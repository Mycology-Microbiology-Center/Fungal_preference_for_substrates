load("tableall.RData")
load("phi.table.RData")
library(tidyverse)
library(marginaleffects)
library(lme4)
library(rcompanion)
library(dplyr)
library(marginaleffects)
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
######
all2_1<-list()
all3_1<-list()
for (i in 1:250) {
  fungi<-table12[[i]]
  if(inherits(fungi, "dgCMatrix")){
    fungi <- as.data.frame(as.matrix(fungi))
  }
  fungi$sample_names<-row.names(fungi)
  fungi2<-merge(fungi,metadata,by="sample_names",all=T)
  row.names(fungi2)<-fungi2$sample_names
  fungi2<-fungi2[,-1]
  fungi2<-fungi2[!is.na(fungi2$ff04207f066a2b2769b796c7e8da74227ea12fe4),]
  all2<-aggregate(fungi2[,1:(ncol(fungi2)-2)],by=list(substrate=fungi2$substrate,site=fungi2$site2),sum)
  all <- all2 %>%
    group_by(site,substrate) %>%
    mutate(site_total = sum(across(where(is.numeric)))) %>%  # total per site
    mutate(across(where(is.numeric), ~ . / site_total)) %>%
    select(-site_total)
  ##calculate the frequency
  fre <- all2 %>%
    group_by(site) %>%
    summarise(across(where(is.numeric), ~ sum(. > 0)))
  ##
  all2_long <- all %>%
    pivot_longer(cols = where(is.numeric), names_to = "OTU", values_to = "Ratio")
  
  all2_fre <- fre %>%
    pivot_longer(cols = where(is.numeric), names_to = "OTU", values_to = "frequency")
  
  all.try<-merge(all2_long,all2_fre,by=c("site","OTU"),all=T)
  all.try<-all.try[all.try$Ratio>0,]
  all.try<-all.try[all.try$frequency>0,]
  ###merging the lifestyle
  all<-merge(table,all.try,by.x="qseqid",by.y="OTU",all=T)
  all<-all[!is.na(all$frequency),]
  all$primary_lifestyle[is.na(all$primary_lifestyle)]<-"unknown"
  all$primary_lifestyle[all$primary_lifestyle==""]<-"unknown"
  r<-result1[[i]]
  r$BestHost <- gsub("^s\\.","", r$BestHost)
  r[r$BestHost %in% "fruit.body",]$BestHost <- "fruit body"
  r[r$BestHost %in% "leaf.litter",]$BestHost <- "leaf litter"
  try<-all[all$qseqid %in% unique(r$OTU), ]
  names(try)[c(1,5)]<-c("OTU","BestHost")
  try<-merge(try,r[,c(1:2,4,5,22)],by=c("OTU","BestHost"),all=T)
  try$sign <-ifelse(try$p.value < 0.05, TRUE, FALSE)
  ###
  aaa<-try
  try[is.na(try$sign),]$sign<-FALSE
  all2 <- try %>%
    arrange(primary_lifestyle, OTU,site) %>%
    mutate(qseqid_ordered = factor(OTU, levels = unique(OTU)))
  
  all2 <- all2 %>%
    mutate(color_group = ifelse(sign, primary_lifestyle, "non-indicator"))
  all3<-all2[!is.na(all2$stat),]
  all2_1[[i]]<-all2
  all3_1[[i]]<-all3
}

all2_2<-list()
all3_2<-list()
for (i in 1:250) {
  fungi<-table22[[i]]
  if(inherits(fungi, "dgCMatrix")){
    fungi <- as.data.frame(as.matrix(fungi))
  }
  fungi$sample_names<-row.names(fungi)
  fungi2<-merge(fungi,metadata,by="sample_names",all=T)
  row.names(fungi2)<-fungi2$sample_names
  fungi2<-fungi2[,-1]
  fungi2<-fungi2[!is.na(fungi2$ff04207f066a2b2769b796c7e8da74227ea12fe4),]
  all2<-aggregate(fungi2[,1:(ncol(fungi2)-2)],by=list(substrate=fungi2$substrate,site=fungi2$site2),sum)
  all <- all2 %>%
    group_by(site,substrate) %>%
    mutate(site_total = sum(across(where(is.numeric)))) %>%  # total per site
    mutate(across(where(is.numeric), ~ . / site_total)) %>%
    select(-site_total)
  ##calculate the frequency
  fre <- all2 %>%
    group_by(site) %>%
    summarise(across(where(is.numeric), ~ sum(. > 0)))
  ##
  all2_long <- all %>%
    pivot_longer(cols = where(is.numeric), names_to = "OTU", values_to = "Ratio")
  
  all2_fre <- fre %>%
    pivot_longer(cols = where(is.numeric), names_to = "OTU", values_to = "frequency")
  
  all.try<-merge(all2_long,all2_fre,by=c("site","OTU"),all=T)
  all.try<-all.try[all.try$Ratio>0,]
  all.try<-all.try[all.try$frequency>0,]
  ###merging the lifestyle
  all<-merge(table,all.try,by.x="qseqid",by.y="OTU",all=T)
  all<-all[!is.na(all$frequency),]
  all$primary_lifestyle[is.na(all$primary_lifestyle)]<-"unknown"
  all$primary_lifestyle[all$primary_lifestyle==""]<-"unknown"
  r<-result2[[i]]
  r$BestHost <- gsub("^s\\.","", r$BestHost)
  r[r$BestHost %in% "fruit.body",]$BestHost <- "fruit body"
  r[r$BestHost %in% "leaf.litter",]$BestHost <- "leaf litter"
  try<-all[all$qseqid %in% unique(r$OTU), ]
  names(try)[c(1,5)]<-c("OTU","BestHost")
  try<-merge(try,r[,c(1:2,4,5,22)],by=c("OTU","BestHost"),all=T)
  try$sign <-ifelse(try$p.value < 0.05, TRUE, FALSE)
  ###
  aaa<-try
  try[is.na(try$sign),]$sign<-FALSE
  all2 <- try %>%
    arrange(primary_lifestyle, OTU,site) %>%
    mutate(qseqid_ordered = factor(OTU, levels = unique(OTU)))
  
  all2 <- all2 %>%
    mutate(color_group = ifelse(sign, primary_lifestyle, "non-indicator"))
  all3<-all2[!is.na(all2$stat),]
  all2_2[[i]]<-all2
  all3_2[[i]]<-all3
}

all2_3<-list()
all3_3<-list()
for (i in 1:250) {
  fungi<-table32[[i]]
  if(inherits(fungi, "dgCMatrix")){
    fungi <- as.data.frame(as.matrix(fungi))
  }
  fungi$sample_names<-row.names(fungi)
  fungi2<-merge(fungi,metadata,by="sample_names",all=T)
  row.names(fungi2)<-fungi2$sample_names
  fungi2<-fungi2[,-1]
  fungi2<-fungi2[!is.na(fungi2$ff04207f066a2b2769b796c7e8da74227ea12fe4),]
  all2<-aggregate(fungi2[,1:(ncol(fungi2)-2)],by=list(substrate=fungi2$substrate,site=fungi2$site2),sum)
  library(dplyr)
  all <- all2 %>%
    group_by(site,substrate) %>%
    mutate(site_total = sum(across(where(is.numeric)))) %>%  # total per site
    mutate(across(where(is.numeric), ~ . / site_total)) %>%
    select(-site_total)
  ##calculate the frequency
  fre <- all2 %>%
    group_by(site) %>%
    summarise(across(where(is.numeric), ~ sum(. > 0)))
  ##
  all2_long <- all %>%
    pivot_longer(cols = where(is.numeric), names_to = "OTU", values_to = "Ratio")
  
  all2_fre <- fre %>%
    pivot_longer(cols = where(is.numeric), names_to = "OTU", values_to = "frequency")
  
  all.try<-merge(all2_long,all2_fre,by=c("site","OTU"),all=T)
  all.try<-all.try[all.try$Ratio>0,]
  all.try<-all.try[all.try$frequency>0,]
  ###merging the lifestyle
  all<-merge(table,all.try,by.x="qseqid",by.y="OTU",all=T)
  all<-all[!is.na(all$frequency),]
  all$primary_lifestyle[is.na(all$primary_lifestyle)]<-"unknown"
  all$primary_lifestyle[all$primary_lifestyle==""]<-"unknown"
  r<-result3[[i]]
  r$BestHost <- gsub("^s\\.","", r$BestHost)
  r[r$BestHost %in% "fruit.body",]$BestHost <- "fruit body"
  r[r$BestHost %in% "leaf.litter",]$BestHost <- "leaf litter"
  try<-all[all$qseqid %in% unique(r$OTU), ]
  names(try)[c(1,5)]<-c("OTU","BestHost")
  try<-merge(try,r[,c(1:2,4,5,22)],by=c("OTU","BestHost"),all=T)
  try$sign <-ifelse(try$p.value < 0.05, TRUE, FALSE)
  ###
  aaa<-try
  try[is.na(try$sign),]$sign<-FALSE
  all2 <- try %>%
    arrange(primary_lifestyle, OTU,site) %>%
    mutate(qseqid_ordered = factor(OTU, levels = unique(OTU)))
  
  all2 <- all2 %>%
    mutate(color_group = ifelse(sign, primary_lifestyle, "non-indicator"))
  all3<-all2[!is.na(all2$stat),]
  all2_3[[i]]<-all2
  all3_3[[i]]<-all3
}


all2_4<-list()
all3_4<-list()
for (i in 1:250) {
  fungi<-table32[[i]]
  if(inherits(fungi, "dgCMatrix")){
    fungi <- as.data.frame(as.matrix(fungi))
  }
  fungi$sample_names<-row.names(fungi)
  fungi2<-merge(fungi,metadata,by="sample_names",all=T)
  row.names(fungi2)<-fungi2$sample_names
  fungi2<-fungi2[,-1]
  fungi2<-fungi2[!is.na(fungi2$ff04207f066a2b2769b796c7e8da74227ea12fe4),]
  all2<-aggregate(fungi2[,1:(ncol(fungi2)-2)],by=list(substrate=fungi2$substrate,site=fungi2$site2),sum)
  all <- all2 %>%
    group_by(site,substrate) %>%
    mutate(site_total = sum(across(where(is.numeric)))) %>%  # total per site
    mutate(across(where(is.numeric), ~ . / site_total)) %>%
    select(-site_total)
  ##calculate the substrate range
  fre <- all2 %>%
    group_by(site) %>%
    summarise(across(where(is.numeric), ~ sum(. > 0)))
  ##
  all2_long <- all %>%
    pivot_longer(cols = where(is.numeric), names_to = "OTU", values_to = "Ratio")
  
  all2_fre <- fre %>%
    pivot_longer(cols = where(is.numeric), names_to = "OTU", values_to = "frequency")
  
  all.try<-merge(all2_long,all2_fre,by=c("site","OTU"),all=T)
  all.try<-all.try[all.try$Ratio>0,]
  all.try<-all.try[all.try$frequency>0,]
  ###merging the lifestyle
  all<-merge(table,all.try,by.x="qseqid",by.y="OTU",all=T)
  all<-all[!is.na(all$frequency),]
  all$primary_lifestyle[is.na(all$primary_lifestyle)]<-"unknown"
  all$primary_lifestyle[all$primary_lifestyle==""]<-"unknown"
  r<-result4[[i]]
  r$BestHost <- gsub("^s\\.","", r$BestHost)
  r[r$BestHost %in% "fruit.body",]$BestHost <- "fruit body"
  r[r$BestHost %in% "leaf.litter",]$BestHost <- "leaf litter"
  try<-all[all$qseqid %in% unique(r$OTU), ]
  names(try)[c(1,5)]<-c("OTU","BestHost")
  try<-merge(try,r[,c(1:2,4,5,22)],by=c("OTU","BestHost"),all=T)
  try$sign <-ifelse(try$p.value < 0.05, TRUE, FALSE)
  ###
  aaa<-try
  try[is.na(try$sign),]$sign<-FALSE
  all2 <- try %>%
    arrange(primary_lifestyle, OTU,site) %>%
    mutate(qseqid_ordered = factor(OTU, levels = unique(OTU)))
  
  all2 <- all2 %>%
    mutate(color_group = ifelse(sign, primary_lifestyle, "non-indicator"))
  all3<-all2[!is.na(all2$stat),]
  all2_4[[i]]<-all2
  all3_4[[i]]<-all3
}

all2<-c(all2_1,all2_2,all2_3,all2_4)
all2<-bind_rows(all2)
all3<-all2[!is.na(all2$frequency),]
all3<-all3[!is.na(all3$site),]
all2_final <- all3 %>%
  group_by(OTU, BestHost, site) %>%
  summarise(frequency = mean(frequency), .groups = "drop")
####substrate range comparison across substrates
model <- lmer(frequency ~ BestHost + (1|site), data = all2_final)
###calculate the OTU +- SD
aaaa <- all2_final %>%
 group_by(BestHost,site) %>%
  summarise(fre=mean(frequency))%>%
  group_by(BestHost) %>%
  summarise(frequency=mean(fre),sd=sd(fre))
###
b<-avg_comparisons(model, variables = list(BestHost = "pairwise")) 
performance::performance(model)
fre.substrate<-parameters::model_parameters(model)
write.csv(fre.substrate,"fre.substrate.parameters.csv")
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group[letter$Group == "fruitbody"] <- "fruit body"
letter$Group[letter$Group == "leaflitter"] <- "leaf litter"
difference <- all2_final %>%
  group_by(BestHost) %>%
  summarise(max=mean(frequency))

y.site<-merge(letter,difference,by.x="Group",by.y="BestHost")
y.site<-data.frame(BestHost=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)

plot_predictions(model,condition = c("BestHost"))+
  labs(x="substrate",y="frequency")+
  theme_light() +
scale_x_discrete(limits = c("leaves","moss","snow","leaf litter",
                            "subsoil","topsoil","bark","fruit body",
                            "wood","feces","lichen"))+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom")+
  geom_text(data = y.site, aes(x = BestHost , y = ymax, label = letter,hjust=-0.5))
######only indicators
all3<-all2[!is.na(all2$frequency),]
all3<-all3[!is.na(all3$site),]
a<-all3[all3$sign, ]
###substrate range of guild 
substrate_range <- a %>%
  group_by(OTU, BestHost, site) %>%
  summarise(frequency = mean(frequency), .groups = "drop")
substrate_range2 <- substrate_range %>%
  group_by(BestHost, site) %>%
  summarise(min.s=min(frequency), max.s=max(frequency))

write.csv(substrate_range2,"substrate_range.csv")
########substrate range comparison across substrates
model <- lm(frequency ~ BestHost  +  site, data = substrate_range)
performance::performance(model)
fre.substrate<-parameters::model_parameters(model)
write.csv(fre.substrate,"fre.indi.substrate.parameters.csv")

###
b<-avg_comparisons(model, variables = list(BestHost  = "pairwise")) 
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group[letter$Group == "fruitbody"] <- "fruit body"
letter$Group[letter$Group == "leaflitter"] <- "leaf litter"
difference <- substrate_range %>%
  group_by(BestHost ) %>%
  summarise(max=mean(frequency))

y.site<-merge(letter,difference,by.x="Group",by.y="BestHost")
y.site<-data.frame(BestHost =factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)

plot_predictions(model,condition = c("BestHost"))+
  labs(x="BestHost",y="frequency")+
  scale_y_continuous(breaks = c(1,2, 2.5, 3))+
  scale_x_discrete(limits = c("leaves","moss","snow","leaf litter",
                              "subsoil","topsoil","bark","fruit body",
                              "wood","feces","lichen"))+
  scale_y_continuous(breaks = c(1,2,2.5,3), limits = c(1, 3))+
theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom")+
  geom_text(data = y.site, aes(x = BestHost  , y = ymax, label = letter,hjust=-0.5))

###indicator comparison
#======comparison based on 1000 time regression
for(i in 1:250){
  data<-all3_1[[i]]
  data$sign<-ifelse(data$sign,'indicator','non-indicator')
  data$iteration<-i
  all3_1[[i]]<-as.data.frame(data) 
}

for(i in 1:250){
  data<-all3_2[[i]]
  data$sign<-ifelse(data$sign,'indicator','non-indicator')
  data$iteration<-i+250
  all3_2[[i]]<-as.data.frame(data) 
}

for(i in 1:250){
  data<-all3_3[[i]]
  data$sign<-ifelse(data$sign,'indicator','non-indicator')
  data$iteration<-i+500
  all3_3[[i]]<-as.data.frame(data) 
}

for(i in 1:250){
  data<-all3_4[[i]]
  data$sign<-ifelse(data$sign,'indicator','non-indicator')
  data$iteration<-i+750
  all3_4[[i]]<-as.data.frame(data) 
}
#####substrate range comparison between indicators and non-indicators
try<-c(all3_1,all3_2, all3_3, all3_4 )
try<-bind_rows(try)

try<-try[,-9]
models <- lapply(1:1000, function(i) {
  df <- try[try$iteration == i, ]
  lm(frequency ~ sign + BestHost + site, data = df)
})


out<-list()
for (i in 1:1000) {
  out[[i]] <- avg_comparisons(models[[i]], variables =  list(sign  = "pairwise"))#substrate range comparison between indicator and non-indicators
}

df <- do.call(rbind, out)
df_sig <- df[df$p.value < 0.05, ]
p_robust <- nrow(df_sig) / nrow(df)
p_robust
quantile(df$estimate)
out<-list()
for (i in 1:1000) {
  out[[i]] <- summary(models[[i]])$r.squared
}
r <- do.call(rbind, out)
quantile(r, c(0.025, 0.5, 0.975))

##relationship between association stength and relative abundance
result<-data.frame()
for (j in 1:1000) {
  try2<-try[try$iteration ==j,]
  try2$sign<-ifelse(try2$sign == "non-indicator", FALSE, TRUE)
  for (i in unique(try2$BestHost)) {
    indicator.min<-min(try2[try2$sign & try2$BestHost %in% i,]$frequency)
    indicator.max<-max(try2[try2$sign & try2$BestHost %in% i,]$frequency)
    indicator.mean<-mean(try2[try2$sign & try2$BestHost %in% i,]$frequency)
    
    nonindicator.min<-min(try2[!try2$sign  & try2$BestHost %in% i,]$frequency)
    nonindicator.max<-max(try2[!try2$sign & try2$BestHost %in% i,]$frequency)
    nonindicator.mean<-mean(try2[!try2$sign & try2$BestHost %in% i,]$frequency)
    
    c<-data.frame(indicator.min=indicator.min,indicator.max=indicator.max,
                  nonindicator.min=nonindicator.min,nonindicator.max=nonindicator.max,
                  nonindicator.meann=nonindicator.mean,indicator.mean=indicator.mean)
    c$substrate <-i
    
    result<-rbind(result,c)
    }
}


result2 <- result %>%
  group_by(substrate) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))
min(result2$nonindicator.meann/result2$indicator.mean)
max(result2$nonindicator.meann/result2$indicator.mean)
##all2
try3<- try[!is.na(try$frequency),]
try3 <- try3 %>%
  group_by(OTU, BestHost, site) %>%
  dplyr::slice_min(abs(frequency - mean(frequency, na.rm = TRUE)), n = 1) %>%
  ungroup()
try3<- try3 %>%
  distinct(OTU, BestHost, site, .keep_all = T)###because last step keep more than one OTUs having the same frequency

indicator.rs<-try3[!is.na(try3$sign),]
indicator.rs<-indicator.rs[(indicator.rs$sign=="indicator" ),]
model <- lm(stat ~ Ratio + BestHost + site, data = indicator.rs)
fre.ratio<-parameters::model_parameters(model)
performance::performance(model)
write.csv(fre.ratio,"stat.abundance.parameter.csv")
####compare the association strength across substrates###
substrate_stat <- a %>%
  group_by(OTU, BestHost) %>%
  summarise(stat = median(stat), .groups = "drop")
aaaa.p <- substrate_stat %>%
  group_by(BestHost) %>%
  summarise(fre=mean(stat),sd=sd(stat))
model <- lm(stat ~ BestHost, data = substrate_stat)
b<-avg_comparisons(model, variables = list(BestHost = "pairwise")) 
fre.ratio<-parameters::model_parameters(model)
performance::performance(model)
write.csv(fre.ratio,"stat.substrate.parameter.csv")

letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group[letter$Group == "fruitbody"] <- "fruit body"
letter$Group[letter$Group == "leaflitter"] <- "leaf litter"
difference <- substrate_stat %>%
  group_by(BestHost) %>%
  summarise(max=mean(stat))

y.site<-merge(letter,difference,by.x="Group",by.y="BestHost")
y.site<-data.frame(substrate=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)

plot_predictions(model,condition = c("BestHost"))+
  labs(x="substrate",y="preference strength")+
  scale_x_discrete(limits = c("leaves","moss","snow","leaf litter",
                              "subsoil","topsoil","bark","fruit body",
                              "wood","feces","lichen"))+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom")+
  geom_text(data = y.site, aes(x = substrate , y = ymax, label = letter,hjust=-0.5))

#######get the indicators within site
load("rarefy.RData")
##
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
    fungi.ind <- ifelse(fungi.ind > 0, 1, 0)
    
    fungi.ind <- as.data.frame(t(fungi.ind))
    fungi.ind$sample_names <- row.names(fungi.ind)
    
    fungi2 <- merge(fungi.ind, metadata, by = "sample_names", all = TRUE)
    row.names(fungi2) <- fungi2$sample_names
    fungi2 <- fungi2[,-1]
    
    fungi2 <- fungi2[!is.na(fungi2$e57bceb79df4a766f8fb9d049798c662ba89dc9e),]
    
    LV <- fungi2[fungi2$site %in% "LV",]
    LV <- LV[, c(colSums(LV[, -((ncol(LV)-1):ncol(LV))]) > 2, TRUE, TRUE)]
    LV <- LV[,-ncol(LV)]
    
    LW <- fungi2[fungi2$site %in% "LW",]
    LW <- LW[, c(colSums(LW[, -((ncol(LW)-1):ncol(LW))]) > 2, TRUE, TRUE)]
    LW <- LW[,-ncol(LW)]
    
    LZ <- fungi2[fungi2$site %in% "LZ",]
    LZ <- LZ[, c(colSums(LZ[, -((ncol(LZ)-1):ncol(LZ))]) > 2, TRUE, TRUE)]
    LZ <- LZ[,-ncol(LZ)]
    
    LV_list[[i]] <- LV
    LW_list[[i]] <- LW
    LZ_list[[i]] <- LZ
  }
  
  return(list(LV = LV_list, LW = LW_list, LZ = LZ_list))
}
######
calc_phi <- function(x, perm = 999){
  setDF(x)
  res <- list()
  res$PHI <- try( multipatt(
    x = x[, -ncol(x)],
    cluster = x[,ncol(x)],
    func = "r.g",
    control = how(nperm = perm),
    duleg = TRUE) )
  ## Calculate complementary fidelity values for the tree species (UbinB metric)
  res$III <- try( multipatt(
    x = x[, -ncol(x)],
    cluster = x[,ncol(x)],
    func = "IndVal.g",
    control = how(nperm = perm),
    duleg = TRUE) )
  
  return(res)
}
res <- prep_data(table1)
res2 <- prep_data(table2)
res3 <- prep_data(table3)
res4 <- prep_data(table4)

LV.ind1<-list(res$LV)
LV.ind2<-list(res2$LV)
LV.ind3<-list(res3$LV)
LV.indk4<-list(res4$LV)

LW.ind1<-list(res$LW)
LW.ind2<-list(res2$LW)
LW.ind3<-list(res3$LW)
LW.indk4<-list(res4$LW)

LZ.ind1<-list(res$LZ)
LZ.ind2<-list(res2$LZ)
LZ.ind3<-list(res3$LZ)
LZ.indk4<-list(res4$LZ)
rm(res, res2, res3, res4)
library(plyr)
library(doFuture)
library(data.table)
library(indicspecies)
registerDoFuture()
plan(multisession, workers = 4)      # for RStudio
options(future.globals.maxSize = 5e9)  # 5GB; default = 500 * 1024 ^ 2 = 500 MiB
####
LW.PIRES1 <- llply(.data = LW.ind1[[1]], .fun = function(x){
  calc_phi(x)
}, .parallel = TRUE)
saveRDS(LW.PIRES1, "LW_PIRES1.rds")

LW.PIRES2 <- llply(.data = LW.ind2[[1]], .fun = function(x){
  calc_phi(x)
}, .parallel = TRUE)
saveRDS(LW.PIRES2, "LW_PIRES2.rds")


LW.PIRES3 <- llply(.data = LW.ind3[[1]], .fun = function(x){
  calc_phi(x)
}, .parallel = TRUE)
saveRDS(LW.PIRES3, "LW_PIRES3.rds")


LW.PIRES4 <- llply(.data = LW.indk4[[1]], .fun = function(x){
  calc_phi(x)
}, .parallel = TRUE)
saveRDS(LW.PIRES4, "LW_PIRES4.rds")

####
LV.PIRES1 <- llply(.data = LV.ind1[[1]], .fun = function(x){
  calc_phi(x)
}, .parallel = TRUE)
saveRDS(LV.PIRES1, "LV_PIRES1.rds")


LV.PIRES2 <- llply(.data = LV.ind2[[1]], .fun = function(x){
  calc_phi(x)
}, .parallel = TRUE)
saveRDS(LV.PIRES2, "LV_PIRES2.rds")


LV.PIRES3 <- llply(.data = LV.ind3[[1]], .fun = function(x){
  calc_phi(x)
}, .parallel = TRUE)
saveRDS(LV.PIRES3, "LV_PIRES3.rds")


LV.PIRES4 <- llply(.data = LV.indk4[[1]], .fun = function(x){
  calc_phi(x)
}, .parallel = TRUE)
saveRDS(LV.PIRES4, "LV_PIRES4.rds")

####
LZ.PIRES1 <- llply(.data = LZ.ind1[[1]], .fun = function(x){
  calc_phi(x)
}, .parallel = TRUE)
saveRDS(LZ.PIRES1, "LZ_PIRES1.rds")


LZ.PIRES2 <- llply(.data = LZ.ind2[[1]], .fun = function(x){
  calc_phi(x)
}, .parallel = TRUE)
saveRDS(LZ.PIRES2, "LZ_PIRES2.rds")

LZ.PIRES3 <- llply(.data = LZ.ind3[[1]], .fun = function(x){
  calc_phi(x)
}, .parallel = TRUE)
saveRDS(LZ.PIRES3, "LZ_PIRES3.rds")

LZ.PIRES4 <- llply(.data = LZ.indk4[[1]], .fun = function(x){
  calc_phi(x)
}, .parallel = TRUE)
saveRDS(LZ.PIRES4, "LZ_PIRES4.rds")
##indicators
##using phi to analysis
LZ_PIRES1<-readRDS("LZ_PIRES1.rds")
LZ_PIRES2<-readRDS("LZ_PIRES2.rds")
LZ_PIRES3<-readRDS("LZ_PIRES3.rds")
LZ_PIRES4<-readRDS("LZ_PIRES4.rds")
LW_PIRES1<-readRDS("LW_PIRES1.rds")
LW_PIRES2<-readRDS("LW_PIRES2.rds")
LW_PIRES3<-readRDS("LW_PIRES3.rds")
LW_PIRES4<-readRDS("LW_PIRES4.rds")
LV_PIRES1<-readRDS("LV_PIRES1.rds")
LV_PIRES2<-readRDS("LV_PIRES2.rds")
LV_PIRES3<-readRDS("LV_PIRES3.rds")
LV_PIRES4<-readRDS("LV_PIRES4.rds")
###extract phi
library(metagMisc)
library(plyr)
library(data.table)
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
LW.result1<-extract(LW_PIRES1)
LW.result2<-extract(LW_PIRES2)
LW.result3<-extract(LW_PIRES3)
LW.result4<-extract(LW_PIRES4)
LZ.result1<-extract(LZ_PIRES1)
LZ.result2<-extract(LZ_PIRES2)
LZ.result3<-extract(LZ_PIRES3)
LZ.result4<-extract(LZ_PIRES4)
LV.result1<-extract(LV_PIRES1)
LV.result2<-extract(LV_PIRES2)
LV.result3<-extract(LV_PIRES3)
LV.result4<-extract(LV_PIRES4)
##growth forms
tax<-read.csv("taxonomy.final.csv",row.names = 1)
tax$phylumn[is.na(tax$phylumn)]<-"p__unclassified"
tax$class[is.na(tax$class)]<-"c__unclassified"
tax$order[is.na(tax$order)]<-"o__unclassified"
tax$family[is.na(tax$family)]<-"f__unclassified"
tax$genus[is.na(tax$genus)]<-"g__unclassified"
tax$species[is.na(tax$species)]<-"s__unclassified"
tax2<-tax
tax2$genus<-gsub(".*_","",tax2$genus)
##
tax2$OTU<-row.names(tax2)
for(i in 1:250){
  LZ.result1[[i]]<-merge(LZ.result1[[i]],tax2,by="OTU",all=T)
  LZ.result1[[i]]$BestHost<-gsub("^s.","",LZ.result1[[i]]$BestHost)
  LZ.result1[[i]]<-LZ.result1[[i]][!is.na(LZ.result1[[i]]$BestHost),]
}

for(i in 1:250){
  LZ.result2[[i]]<-merge(LZ.result2[[i]],tax2,by="OTU",all=T)
  LZ.result2[[i]]<-LZ.result2[[i]][!is.na(LZ.result2[[i]]$BestHost),]
}

for(i in 1:250){
  LZ.result3[[i]]<-merge(LZ.result3[[i]],tax2,by="OTU",all=T)
  LZ.result3[[i]]<-LZ.result3[[i]][!is.na(LZ.result3[[i]]$BestHost),]
}

for(i in 1:250){
  LZ.result4[[i]]<-merge(LZ.result4[[i]],tax2,by="OTU",all=T)
  LZ.result4[[i]]<-LZ.result4[[i]][!is.na(LZ.result4[[i]]$BestHost),]
}

for(i in 1:250){
  LV.result1[[i]]<-merge(LV.result1[[i]],tax2,by="OTU",all=T)
  LV.result1[[i]]$BestHost<-gsub("^s.","",LV.result1[[i]]$BestHost)
  LV.result1[[i]]<-LV.result1[[i]][!is.na(LV.result1[[i]]$BestHost),]
}

for(i in 1:250){
  LV.result2[[i]]<-merge(LV.result2[[i]],tax2,by="OTU",all=T)
  LV.result2[[i]]<-LV.result2[[i]][!is.na(LV.result2[[i]]$BestHost),]
}

for(i in 1:250){
  LV.result3[[i]]<-merge(LV.result3[[i]],tax2,by="OTU",all=T)
  LV.result3[[i]]<-LV.result3[[i]][!is.na(LV.result3[[i]]$BestHost),]
}

for(i in 1:250){
  LV.result4[[i]]<-merge(LV.result4[[i]],tax2,by="OTU",all=T)
  LV.result4[[i]]<-LV.result4[[i]][!is.na(LV.result4[[i]]$BestHost),]
}

for(i in 1:250){
  LW.result1[[i]]<-merge(LW.result1[[i]],tax2,by="OTU",all=T)
  LW.result1[[i]]$BestHost<-gsub("^s.","",LW.result1[[i]]$BestHost)
  LW.result1[[i]]<-LW.result1[[i]][!is.na(LW.result1[[i]]$BestHost),]
}

for(i in 1:250){
  LW.result2[[i]]<-merge(LW.result2[[i]],tax2,by="OTU",all=T)
  LW.result2[[i]]<-LW.result2[[i]][!is.na(LW.result2[[i]]$BestHost),]
}

for(i in 1:250){
  LW.result3[[i]]<-merge(LW.result3[[i]],tax2,by="OTU",all=T)
  LW.result3[[i]]<-LW.result3[[i]][!is.na(LW.result3[[i]]$BestHost),]
}

for(i in 1:250){
  LW.result4[[i]]<-merge(LW.result4[[i]],tax2,by="OTU",all=T)
  LW.result4[[i]]<-LW.result4[[i]][!is.na(LW.result4[[i]]$BestHost),]
}
save.image("separate.phi.table.RData")
######calculate the propogule flow among substrates
library(tidyr)
library(tibble)
library(purrr)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(performance)
load("separate.phi.table.RData")
for (i in 1:250) {
  LZ.result1[[i]]$BestHost <- gsub("^s\\.", "", LZ.result1[[i]]$BestHost ) 
}

for (i in 1:250) {
  LZ.result2[[i]]$BestHost <- gsub("^s\\.", "", LZ.result2[[i]]$BestHost ) 
}

for (i in 1:250) {
  LZ.result3[[i]]$BestHost <- gsub("^s\\.", "", LZ.result3[[i]]$BestHost ) 
}

for (i in 1:250) {
  LZ.result4[[i]]$BestHost <- gsub("^s\\.", "", LZ.result4[[i]]$BestHost ) 
}

for (i in 1:250) {
  LV.result1[[i]]$BestHost <- gsub("^s\\.", "", LV.result1[[i]]$BestHost ) 
}

for (i in 1:250) {
  LV.result2[[i]]$BestHost <- gsub("^s\\.", "", LV.result2[[i]]$BestHost ) 
}

for (i in 1:250) {
  LV.result3[[i]]$BestHost <- gsub("^s\\.", "", LV.result3[[i]]$BestHost ) 
}

for (i in 1:250) {
  LV.result4[[i]]$BestHost <- gsub("^s\\.", "", LV.result4[[i]]$BestHost ) 
}
for (i in 1:250) {
  LW.result1[[i]]$BestHost <- gsub("^s\\.", "", LW.result1[[i]]$BestHost ) 
}

for (i in 1:250) {
  LW.result2[[i]]$BestHost <- gsub("^s\\.", "", LW.result2[[i]]$BestHost ) 
}

for (i in 1:250) {
  LW.result3[[i]]$BestHost <- gsub("^s\\.", "", LW.result3[[i]]$BestHost ) 
}

for (i in 1:250) {
  LW.result4[[i]]$BestHost <- gsub("^s\\.", "", LW.result4[[i]]$BestHost ) 
}

switch_new1<-list()
consistence31 <-list()
switch31<-list()
try1<-list()
for (i in 1:250) {
  LV.r2<-LV.result1[[i]]
  LW.r2<-LW.result1[[i]]
  LZ.r2<-LZ.result1[[i]]
  LZ.r2$site <- "LZ"
  LV.r2$site <- "LV"
  LW.r2$site <- "LW"
  all<-rbind(LZ.r2,LW.r2,LV.r2)
  non.indicator<-all[!all$p.value<0.05,]
  all<- all[all$p.value<0.05,]
  switch<- all[all$OTU %in% unique(all$OTU[duplicated(all$OTU)]),]
  switch <- switch %>%
    add_count(OTU, BestHost) %>%   
    filter(n == 1) %>%             
    select(OTU, BestHost, site)    
  
  switch2<-all[all$OTU %in% unique(switch$OTU),]
  consistence<-all[!all$OTU %in% unique(switch$OTU),] 
  consistence3<-select(consistence,c("OTU", "BestHost", 'site','stat'))
  consistence31[[i]]<-consistence3
  ########calculate the unstable bias in consistence
  stable <- consistence %>%
    add_count(OTU, BestHost) %>%   
    filter(n == 3)
  unstable <- consistence %>%
    add_count(OTU, BestHost) %>%   
    filter(n != 3)
  unstable <-unstable[,-24]
  unstable <-unstable[,-25]
  non.indicator <-non.indicator[,-24]
  non.indicator2<-rbind(unstable,non.indicator)
  non.indicator2<-non.indicator2[non.indicator2$OTU %in% unique(unstable$OTU),]
  ###site difference
  site.d <- non.indicator2 %>%
    add_count(OTU) %>%   
    filter(n == 1)
  ####non.indicator
  non.i.d <- non.indicator2 %>%
    add_count(OTU) %>%   
    filter(n != 1)
  non.i.d2<-non.i.d[non.i.d$p.value<0.05,]
  non.ind2 <- non.i.d2 %>%
    add_count(OTU,BestHost) %>%   
    filter(nn >=2)
  non.ind3 <- non.i.d2[!non.i.d2$OTU %in% unique(non.ind2$OTU),]#get the indicator that changed into non-indicator more than once
  non.ind4 <- non.i.d2[non.i.d2$OTU %in% unique(non.ind2$OTU),]#get the indicator that changed into non-indicator once
  #####select the unstable consistent
  switch3<-select(switch2,c("OTU", "BestHost", 'site','stat'))
  switch31[[i]]<-switch3
  switch2<-select(switch2,c("OTU", "BestHost", 'site'))
  switch2$BestHost<-paste0(switch2$BestHost,switch2$site)
  switch_new <- switch2 %>%
    group_by(OTU) %>%
    summarise(
      combos = list(combn(unique(BestHost), 2, simplify = FALSE)),
      .groups = "drop"
    ) %>%
    unnest(combos) %>%
    mutate(
      combos = map(combos, sort)  # sort each pair alphabetically
    ) %>%
    mutate(
      from = map_chr(combos, 1),
      to   = map_chr(combos, 2)
    ) %>%
    select(OTU, from, to)
  
  switch_new$from_site<-substr(switch_new$from, (nchar(switch_new $from)-1),nchar(switch_new $from))
  switch_new$to_site<-substr(switch_new$to, (nchar(switch_new$to)-1),nchar(switch_new $to))
  switch_new$from<-gsub("LZ|LW|LV","",switch_new$from)
  switch_new$to<-gsub("LZ|LW|LV","",switch_new$to)
  ###
  switch_new2 <- switch_new %>%
    filter(from!=to)
  switch_new2$number<-1
  switch_new2 <- switch_new2 %>%
    group_by(from, to, from_site,to_site) %>%
    summarise(number.c=sum(number))
  ##
  consistence2<-select(stable,c("OTU", "BestHost", 'site'))
  consistence2 <- consistence2[consistence2$OTU %in% unique(consistence2$OTU[duplicated(consistence2$OTU)]),]
  consistence2$BestHost<-paste0(consistence2$BestHost,consistence2$site)
  consistence_new <- consistence2 %>%
    group_by(OTU) %>%
    summarise(
      combos = list(combn(unique(BestHost), 2, simplify = FALSE)),
      .groups = "drop"
    ) %>%
    unnest(combos) %>%
    mutate(
      combos = map(combos, sort)  # sort each pair alphabetically
    ) %>%
    mutate(
      from = map_chr(combos, 1),
      to   = map_chr(combos, 2)
    ) %>%
    select(OTU, from, to)
  
  consistence_new$from_site<-substr(consistence_new$from, (nchar(consistence_new $from)-1),nchar(consistence_new $from))
  consistence_new$to_site<-substr(consistence_new$to, (nchar(consistence_new$to)-1),nchar(consistence_new $to))
  consistence_new$from<-gsub("LZ|LW|LV","",consistence_new$from)
  consistence_new$to<-gsub("LZ|LW|LV","",consistence_new$to)
  ###
  consistence_new2 <- consistence_new %>%
    filter(from!=to)
  consistence_new2 <- consistence_new
  consistence_new2$number<-1
  consistence_new2 <- consistence_new2 %>%
    group_by(from, to, from_site,to_site) %>%
    summarise(number.c=sum(number))
  ####changed into the non-indicator once
  tnon.ind<-select(non.ind4,c("OTU", "BestHost", 'site'))
  tnon.ind$BestHost<-paste0(tnon.ind$BestHost,tnon.ind$site)
  tnon.ind_new <- tnon.ind %>%
    group_by(OTU) %>%
    summarise(
      combos = list(combn(unique(BestHost), 2, simplify = FALSE)),
      .groups = "drop"
    ) %>%
    unnest(combos) %>%
    mutate(
      combos = map(combos, sort)  # sort each pair alphabetically
    ) %>%
    mutate(
      from = map_chr(combos, 1),
      to   = map_chr(combos, 2)
    ) %>%
    select(OTU, from, to)
  
  tnon.ind_new$from_site<-substr(tnon.ind_new$from, (nchar(tnon.ind_new $from)-1),nchar(tnon.ind_new $from))
  tnon.ind_new$to_site<-substr(tnon.ind_new$to, (nchar(tnon.ind_new$to)-1),nchar(tnon.ind_new $to))
  tnon.ind_new$from<-gsub("LZ|LW|LV","",tnon.ind_new$from)
  tnon.ind_new$to<-gsub("LZ|LW|LV","",tnon.ind_new$to)
  ###
  tnon.ind_new2 <- tnon.ind_new %>%
    filter(from!=to)
  tnon.ind_new2 <- tnon.ind_new 
  tnon.ind_new2$number<-1
  tnon.ind_new2 <- tnon.ind_new2 %>%
    group_by(from, to, from_site,to_site) %>%
    summarise(number.c=sum(number))
  ###
  switch_new<-rbind(switch_new2,consistence_new2,tnon.ind_new2)
  ##
  switch_new <- switch_new %>%
    group_by(from, to) %>%
    summarise(number.c=mean(number.c))
  switch_new1[[i]]<-switch_new
  ####
  switch2$number<-1
  switch2$BestHost<-gsub("LZ|LV|LW","",switch2$BestHost)
  switch.p<-switch2 %>%
    group_by(BestHost,site) %>%
    summarise(number.s=sum(number))
  stable<-stable[,-24]
  consistence<-rbind(stable,non.ind4)
  consistence$number<-1
  consistence.p<-consistence %>%
    group_by(BestHost,site) %>%
    summarise(number.c=sum(number))
  all2<-select(all,c("OTU", "BestHost", 'site'))
  all2<-all2 %>%
    distinct(OTU, BestHost,site)
  all2$number<-1
  all.p<-all2 %>%
    group_by(BestHost,site) %>%
    summarise(number.a=sum(number))
  ###
  non.i.d2<-non.ind3
  non.i.d2$number<-1
  non.i.d.p<-non.i.d2 %>%
    group_by(BestHost,site) %>%
    summarise(number.is=sum(number))
  site.d2<-site.d
  site.d2<-site.d2[site.d2$p.value<0.05,]
  site.d2 $number<-1
  site.d.p<-site.d2 %>%
    group_by(BestHost,site) %>%
    summarise(number.sd=sum(number))
  ##
  me<-function(x,y){
    a<-merge(x,y,by=c("BestHost",'site'),all=T)
    return(a)
  }
  try<-list(switch.p,consistence.p,non.i.d.p,site.d.p,all.p)
  try<-Reduce(me,try)
  try[, 3:7][is.na(try[, 3:7])] <- 0
  try$s.r<-try$number.s/try$number.a
  try$c.r<-try$number.c/try$number.a
  #try$ic.r<-try$number.ic/try$number.a
  try$is.r<-try$number.is/try$number.a##non-indicator
  try$sd.r<-try$number.sd/try$number.a###site unique
  try<-try[!is.na(try$s.r),]##filter out data without significant indicator
  try1[[i]]<-try
}

switch_new22<-list()
consistence32 <-list()
switch32<-list()
try2<-list()
for (i in 1:250) {
  LV.r2<-LV.result2[[i]]
  LW.r2<-LW.result2[[i]]
  LZ.r2<-LZ.result2[[i]]
  LZ.r2$site <- "LZ"
  LV.r2$site <- "LV"
  LW.r2$site <- "LW"
  all<-rbind(LZ.r2,LW.r2,LV.r2)
  non.indicator<-all[!all$p.value<0.05,]
  all<- all[all$p.value<0.05,]
  switch<- all[all$OTU %in% unique(all$OTU[duplicated(all$OTU)]),]
  switch <- switch %>%
    add_count(OTU, BestHost) %>%   
    filter(n == 1) %>%             
    select(OTU, BestHost, site)    
  
  switch2<-all[all$OTU %in% unique(switch$OTU),]
  consistence<-all[!all$OTU %in% unique(switch$OTU),] 
  consistence3<-select(consistence,c("OTU", "BestHost", 'site','stat'))
  consistence32[[i]]<-consistence3
  ########calculate the unstable bias in consistence
  stable <- consistence %>%
    add_count(OTU, BestHost) %>%   
    filter(n == 3)
  unstable <- consistence %>%
    add_count(OTU, BestHost) %>%   
    filter(n != 3)
  unstable <-unstable[,-24]
  unstable <-unstable[,-25]
  non.indicator <-non.indicator[,-24]
  non.indicator2<-rbind(unstable,non.indicator)
  non.indicator2<-non.indicator2[non.indicator2$OTU %in% unique(unstable$OTU),]
  ###site difference
  site.d <- non.indicator2 %>%
    add_count(OTU) %>%   
    filter(n == 1)
  ####non.indicator
  non.i.d <- non.indicator2 %>%
    add_count(OTU) %>%   
    filter(n != 1)
  non.i.d2<-non.i.d[non.i.d$p.value<0.05,]
  non.ind2 <- non.i.d2 %>%
    add_count(OTU,BestHost) %>%   
    filter(nn >=2)
  non.ind3 <- non.i.d2[!non.i.d2$OTU %in% unique(non.ind2$OTU),]#get the indicator that changed into non-indicator more than once
  non.ind4 <- non.i.d2[non.i.d2$OTU %in% unique(non.ind2$OTU),]#get the indicator that changed into non-indicator once
  #####select the unstable consistent
  switch3<-select(switch2,c("OTU", "BestHost", 'site','stat'))
  switch32[[i]]<-switch3
  switch2<-select(switch2,c("OTU", "BestHost", 'site'))
  switch2$BestHost<-paste0(switch2$BestHost,switch2$site)
  switch_new <- switch2 %>%
    group_by(OTU) %>%
    summarise(
      combos = list(combn(unique(BestHost), 2, simplify = FALSE)),
      .groups = "drop"
    ) %>%
    unnest(combos) %>%
    mutate(
      combos = map(combos, sort)  # sort each pair alphabetically
    ) %>%
    mutate(
      from = map_chr(combos, 1),
      to   = map_chr(combos, 2)
    ) %>%
    select(OTU, from, to)
  
  switch_new$from_site<-substr(switch_new$from, (nchar(switch_new $from)-1),nchar(switch_new $from))
  switch_new$to_site<-substr(switch_new$to, (nchar(switch_new$to)-1),nchar(switch_new $to))
  switch_new$from<-gsub("LZ|LW|LV","",switch_new$from)
  switch_new$to<-gsub("LZ|LW|LV","",switch_new$to)
  ###
  switch_new2 <- switch_new %>%
    filter(from!=to)
  switch_new2$number<-1
  switch_new2 <- switch_new2 %>%
    group_by(from, to, from_site,to_site) %>%
    summarise(number.c=sum(number))
  ##
  consistence2<-select(stable,c("OTU", "BestHost", 'site'))
  consistence2 <- consistence2[consistence2$OTU %in% unique(consistence2$OTU[duplicated(consistence2$OTU)]),]
  consistence2$BestHost<-paste0(consistence2$BestHost,consistence2$site)
  consistence_new <- consistence2 %>%
    group_by(OTU) %>%
    summarise(
      combos = list(combn(unique(BestHost), 2, simplify = FALSE)),
      .groups = "drop"
    ) %>%
    unnest(combos) %>%
    mutate(
      combos = map(combos, sort)  # sort each pair alphabetically
    ) %>%
    mutate(
      from = map_chr(combos, 1),
      to   = map_chr(combos, 2)
    ) %>%
    select(OTU, from, to)
  
  consistence_new$from_site<-substr(consistence_new$from, (nchar(consistence_new $from)-1),nchar(consistence_new $from))
  consistence_new$to_site<-substr(consistence_new$to, (nchar(consistence_new$to)-1),nchar(consistence_new $to))
  consistence_new$from<-gsub("LZ|LW|LV","",consistence_new$from)
  consistence_new$to<-gsub("LZ|LW|LV","",consistence_new$to)
  ###
  consistence_new2 <- consistence_new %>%
    filter(from!=to)
  consistence_new2 <- consistence_new
  consistence_new2$number<-1
  consistence_new2 <- consistence_new2 %>%
    group_by(from, to, from_site,to_site) %>%
    summarise(number.c=sum(number))
  ####changed into the non-indicator once
  tnon.ind<-select(non.ind4,c("OTU", "BestHost", 'site'))
  tnon.ind$BestHost<-paste0(tnon.ind$BestHost,tnon.ind$site)
  tnon.ind_new <- tnon.ind %>%
    group_by(OTU) %>%
    summarise(
      combos = list(combn(unique(BestHost), 2, simplify = FALSE)),
      .groups = "drop"
    ) %>%
    unnest(combos) %>%
    mutate(
      combos = map(combos, sort)  # sort each pair alphabetically
    ) %>%
    mutate(
      from = map_chr(combos, 1),
      to   = map_chr(combos, 2)
    ) %>%
    select(OTU, from, to)
  
  tnon.ind_new$from_site<-substr(tnon.ind_new$from, (nchar(tnon.ind_new $from)-1),nchar(tnon.ind_new $from))
  tnon.ind_new$to_site<-substr(tnon.ind_new$to, (nchar(tnon.ind_new$to)-1),nchar(tnon.ind_new $to))
  tnon.ind_new$from<-gsub("LZ|LW|LV","",tnon.ind_new$from)
  tnon.ind_new$to<-gsub("LZ|LW|LV","",tnon.ind_new$to)
  ###
  tnon.ind_new2 <- tnon.ind_new %>%
    filter(from!=to)
  tnon.ind_new2 <- tnon.ind_new 
  tnon.ind_new2$number<-1
  tnon.ind_new2 <- tnon.ind_new2 %>%
    group_by(from, to, from_site,to_site) %>%
    summarise(number.c=sum(number))
  ###
  switch_new<-rbind(switch_new2,consistence_new2,tnon.ind_new2)
  ##
  switch_new <- switch_new %>%
    group_by(from, to) %>%
    summarise(number.c=mean(number.c))
  switch_new22[[i]]<-switch_new
  ####
  switch2$number<-1
  switch2$BestHost<-gsub("LZ|LV|LW","",switch2$BestHost)
  switch.p<-switch2 %>%
    group_by(BestHost,site) %>%
    summarise(number.s=sum(number))
  stable<-stable[,-24]
  consistence<-rbind(stable,non.ind4)
  consistence$number<-1
  consistence.p<-consistence %>%
    group_by(BestHost,site) %>%
    summarise(number.c=sum(number))
  all2<-select(all,c("OTU", "BestHost", 'site'))
  all2<-all2 %>%
    distinct(OTU, BestHost,site)
  all2$number<-1
  all.p<-all2 %>%
    group_by(BestHost,site) %>%
    summarise(number.a=sum(number))
  ###
  non.i.d2<-non.ind3
  non.i.d2$number<-1
  non.i.d.p<-non.i.d2 %>%
    group_by(BestHost,site) %>%
    summarise(number.is=sum(number))
  site.d2<-site.d
  site.d2<-site.d2[site.d2$p.value<0.05,]
  site.d2 $number<-1
  site.d.p<-site.d2 %>%
    group_by(BestHost,site) %>%
    summarise(number.sd=sum(number))
  ##
  me<-function(x,y){
    a<-merge(x,y,by=c("BestHost",'site'),all=T)
    return(a)
  }
  try<-list(switch.p,consistence.p,non.i.d.p,site.d.p,all.p)
  try<-Reduce(me,try)
  try[, 3:7][is.na(try[, 3:7])] <- 0
  try$s.r<-try$number.s/try$number.a
  try$c.r<-try$number.c/try$number.a
  #try$ic.r<-try$number.ic/try$number.a
  try$is.r<-try$number.is/try$number.a##non-indicator
  try$sd.r<-try$number.sd/try$number.a###site unique
  try<-try[!is.na(try$s.r),]##filter out data without significant indicator
  try2[[i]]<-try
}

switch_new23<-list()
consistence33 <-list()
switch33<-list()
try3<-list()
for (i in 1:250) {
  LV.r2<-LV.result3[[i]]
  LW.r2<-LW.result3[[i]]
  LZ.r2<-LZ.result3[[i]]
  LZ.r2$site <- "LZ"
  LV.r2$site <- "LV"
  LW.r2$site <- "LW"
  all<-rbind(LZ.r2,LW.r2,LV.r2)
  non.indicator<-all[!all$p.value<0.05,]
  all<- all[all$p.value<0.05,]
  switch<- all[all$OTU %in% unique(all$OTU[duplicated(all$OTU)]),]
  library(dplyr)
  switch <- switch %>%
    add_count(OTU, BestHost) %>%   
    filter(n == 1) %>%             
    select(OTU, BestHost, site)    
  
  switch2<-all[all$OTU %in% unique(switch$OTU),]
  consistence<-all[!all$OTU %in% unique(switch$OTU),] 
  consistence3<-select(consistence,c("OTU", "BestHost", 'site','stat'))
  consistence33[[i]]<-consistence3
  ########calculate the unstable bias in consistence
  stable <- consistence %>%
    add_count(OTU, BestHost) %>%   
    filter(n == 3)
  unstable <- consistence %>%
    add_count(OTU, BestHost) %>%   
    filter(n != 3)
  unstable <-unstable[,-24]
  unstable <-unstable[,-25]
  non.indicator <-non.indicator[,-24]
  non.indicator2<-rbind(unstable,non.indicator)
  non.indicator2<-non.indicator2[non.indicator2$OTU %in% unique(unstable$OTU),]
  ###site difference
  site.d <- non.indicator2 %>%
    add_count(OTU) %>%   
    filter(n == 1)
  ####non.indicator
  non.i.d <- non.indicator2 %>%
    add_count(OTU) %>%   
    filter(n != 1)
  non.i.d2<-non.i.d[non.i.d$p.value<0.05,]
  non.ind2 <- non.i.d2 %>%
    add_count(OTU,BestHost) %>%   
    filter(nn >=2)
  non.ind3 <- non.i.d2[!non.i.d2$OTU %in% unique(non.ind2$OTU),]#get the indicator that changed into non-indicator more than once
  non.ind4 <- non.i.d2[non.i.d2$OTU %in% unique(non.ind2$OTU),]#get the indicator that changed into non-indicator once
  #####select the unstable consistent
  switch3<-select(switch2,c("OTU", "BestHost", 'site','stat'))
  switch33[[i]]<-switch3
  switch2<-select(switch2,c("OTU", "BestHost", 'site'))
  switch2$BestHost<-paste0(switch2$BestHost,switch2$site)
  switch_new <- switch2 %>%
    group_by(OTU) %>%
    summarise(
      combos = list(combn(unique(BestHost), 2, simplify = FALSE)),
      .groups = "drop"
    ) %>%
    unnest(combos) %>%
    mutate(
      combos = map(combos, sort)  # sort each pair alphabetically
    ) %>%
    mutate(
      from = map_chr(combos, 1),
      to   = map_chr(combos, 2)
    ) %>%
    select(OTU, from, to)
  
  switch_new$from_site<-substr(switch_new$from, (nchar(switch_new $from)-1),nchar(switch_new $from))
  switch_new$to_site<-substr(switch_new$to, (nchar(switch_new$to)-1),nchar(switch_new $to))
  switch_new$from<-gsub("LZ|LW|LV","",switch_new$from)
  switch_new$to<-gsub("LZ|LW|LV","",switch_new$to)
  ###
  switch_new2 <- switch_new %>%
    filter(from!=to)
  switch_new2$number<-1
  switch_new2 <- switch_new2 %>%
    group_by(from, to, from_site,to_site) %>%
    summarise(number.c=sum(number))
  ##
  consistence2<-select(stable,c("OTU", "BestHost", 'site'))
  consistence2 <- consistence2[consistence2$OTU %in% unique(consistence2$OTU[duplicated(consistence2$OTU)]),]
  consistence2$BestHost<-paste0(consistence2$BestHost,consistence2$site)
  consistence_new <- consistence2 %>%
    group_by(OTU) %>%
    summarise(
      combos = list(combn(unique(BestHost), 2, simplify = FALSE)),
      .groups = "drop"
    ) %>%
    unnest(combos) %>%
    mutate(
      combos = map(combos, sort)  # sort each pair alphabetically
    ) %>%
    mutate(
      from = map_chr(combos, 1),
      to   = map_chr(combos, 2)
    ) %>%
    select(OTU, from, to)
  
  consistence_new$from_site<-substr(consistence_new$from, (nchar(consistence_new $from)-1),nchar(consistence_new $from))
  consistence_new$to_site<-substr(consistence_new$to, (nchar(consistence_new$to)-1),nchar(consistence_new $to))
  consistence_new$from<-gsub("LZ|LW|LV","",consistence_new$from)
  consistence_new$to<-gsub("LZ|LW|LV","",consistence_new$to)
  ###
  consistence_new2 <- consistence_new %>%
    filter(from!=to)
  consistence_new2 <- consistence_new
  consistence_new2$number<-1
  consistence_new2 <- consistence_new2 %>%
    group_by(from, to, from_site,to_site) %>%
    summarise(number.c=sum(number))
  ####changed into the non-indicator once
  tnon.ind<-select(non.ind4,c("OTU", "BestHost", 'site'))
  tnon.ind$BestHost<-paste0(tnon.ind$BestHost,tnon.ind$site)
  tnon.ind_new <- tnon.ind %>%
    group_by(OTU) %>%
    summarise(
      combos = list(combn(unique(BestHost), 2, simplify = FALSE)),
      .groups = "drop"
    ) %>%
    unnest(combos) %>%
    mutate(
      combos = map(combos, sort)  # sort each pair alphabetically
    ) %>%
    mutate(
      from = map_chr(combos, 1),
      to   = map_chr(combos, 2)
    ) %>%
    select(OTU, from, to)
  
  tnon.ind_new$from_site<-substr(tnon.ind_new$from, (nchar(tnon.ind_new $from)-1),nchar(tnon.ind_new $from))
  tnon.ind_new$to_site<-substr(tnon.ind_new$to, (nchar(tnon.ind_new$to)-1),nchar(tnon.ind_new $to))
  tnon.ind_new$from<-gsub("LZ|LW|LV","",tnon.ind_new$from)
  tnon.ind_new$to<-gsub("LZ|LW|LV","",tnon.ind_new$to)
  ###
  tnon.ind_new2 <- tnon.ind_new %>%
    filter(from!=to)
  tnon.ind_new2 <- tnon.ind_new 
  tnon.ind_new2$number<-1
  tnon.ind_new2 <- tnon.ind_new2 %>%
    group_by(from, to, from_site,to_site) %>%
    summarise(number.c=sum(number))
  ###
  switch_new<-rbind(switch_new2,consistence_new2,tnon.ind_new2)
  ##
  switch_new <- switch_new %>%
    group_by(from, to) %>%
    summarise(number.c=mean(number.c))
  switch_new23[[i]]<-switch_new
  ####
  switch2$number<-1
  switch2$BestHost<-gsub("LZ|LV|LW","",switch2$BestHost)
  switch.p<-switch2 %>%
    group_by(BestHost,site) %>%
    summarise(number.s=sum(number))
  stable<-stable[,-24]
  consistence<-rbind(stable,non.ind4)
  consistence$number<-1
  consistence.p<-consistence %>%
    group_by(BestHost,site) %>%
    summarise(number.c=sum(number))
  all2<-select(all,c("OTU", "BestHost", 'site'))
  all2<-all2 %>%
    distinct(OTU, BestHost,site)
  all2$number<-1
  all.p<-all2 %>%
    group_by(BestHost,site) %>%
    summarise(number.a=sum(number))
  ###
  non.i.d2<-non.ind3
  non.i.d2$number<-1
  non.i.d.p<-non.i.d2 %>%
    group_by(BestHost,site) %>%
    summarise(number.is=sum(number))
  site.d2<-site.d
  site.d2<-site.d2[site.d2$p.value<0.05,]
  site.d2 $number<-1
  site.d.p<-site.d2 %>%
    group_by(BestHost,site) %>%
    summarise(number.sd=sum(number))
  ##
  me<-function(x,y){
    a<-merge(x,y,by=c("BestHost",'site'),all=T)
    return(a)
  }
  try<-list(switch.p,consistence.p,non.i.d.p,site.d.p,all.p)
  try<-Reduce(me,try)
  try[, 3:7][is.na(try[, 3:7])] <- 0
  try$s.r<-try$number.s/try$number.a
  try$c.r<-try$number.c/try$number.a
  #try$ic.r<-try$number.ic/try$number.a
  try$is.r<-try$number.is/try$number.a##non-indicator
  try$sd.r<-try$number.sd/try$number.a###site unique
  try<-try[!is.na(try$s.r),]##filter out data without significant indicator
  try3[[i]]<-try
}

switch_new24<-list()
consistence34 <-list()
switch34<-list()
try4<-list()
for (i in 1:250) {
  LV.r2<-LV.result4[[i]]
  LW.r2<-LW.result4[[i]]
  LZ.r2<-LZ.result4[[i]]
  LZ.r2$site <- "LZ"
  LV.r2$site <- "LV"
  LW.r2$site <- "LW"
  all<-rbind(LZ.r2,LW.r2,LV.r2)
  non.indicator<-all[!all$p.value<0.05,]
  all<- all[all$p.value<0.05,]
  switch<- all[all$OTU %in% unique(all$OTU[duplicated(all$OTU)]),]
  switch <- switch %>%
    add_count(OTU, BestHost) %>%   
    filter(n == 1) %>%             
    select(OTU, BestHost, site)    
  
  switch2<-all[all$OTU %in% unique(switch$OTU),]
  consistence<-all[!all$OTU %in% unique(switch$OTU),] 
  consistence3<-select(consistence,c("OTU", "BestHost", 'site','stat'))
  consistence34[[i]]<-consistence3
  ########calculate the unstable bias in consistence
  stable <- consistence %>%
    add_count(OTU, BestHost) %>%   
    filter(n == 3)
  unstable <- consistence %>%
    add_count(OTU, BestHost) %>%   
    filter(n != 3)
  unstable <-unstable[,-24]
  unstable <-unstable[,-25]
  non.indicator <-non.indicator[,-24]
  non.indicator2<-rbind(unstable,non.indicator)
  non.indicator2<-non.indicator2[non.indicator2$OTU %in% unique(unstable$OTU),]
  ###site difference
  site.d <- non.indicator2 %>%
    add_count(OTU) %>%   
    filter(n == 1)
  ####non.indicator
  non.i.d <- non.indicator2 %>%
    add_count(OTU) %>%   
    filter(n != 1)
  non.i.d2<-non.i.d[non.i.d$p.value<0.05,]
  non.ind2 <- non.i.d2 %>%
    add_count(OTU,BestHost) %>%   
    filter(nn >=2)
  non.ind3 <- non.i.d2[!non.i.d2$OTU %in% unique(non.ind2$OTU),]#get the indicator that changed into non-indicator more than once
  non.ind4 <- non.i.d2[non.i.d2$OTU %in% unique(non.ind2$OTU),]#get the indicator that changed into non-indicator once
  #####select the unstable consistent
  switch3<-select(switch2,c("OTU", "BestHost", 'site','stat'))
  switch34[[i]]<-switch3
  switch2<-select(switch2,c("OTU", "BestHost", 'site'))

  switch2$BestHost<-paste0(switch2$BestHost,switch2$site)
  switch_new <- switch2 %>%
    group_by(OTU) %>%
    summarise(
      combos = list(combn(unique(BestHost), 2, simplify = FALSE)),
      .groups = "drop"
    ) %>%
    unnest(combos) %>%
    mutate(
      combos = map(combos, sort)  # sort each pair alphabetically
    ) %>%
    mutate(
      from = map_chr(combos, 1),
      to   = map_chr(combos, 2)
    ) %>%
    select(OTU, from, to)
  
  switch_new$from_site<-substr(switch_new$from, (nchar(switch_new $from)-1),nchar(switch_new $from))
  switch_new$to_site<-substr(switch_new$to, (nchar(switch_new$to)-1),nchar(switch_new $to))
  switch_new$from<-gsub("LZ|LW|LV","",switch_new$from)
  switch_new$to<-gsub("LZ|LW|LV","",switch_new$to)
  ###
  switch_new2 <- switch_new %>%
    filter(from!=to)
  switch_new2$number<-1
  switch_new2 <- switch_new2 %>%
    group_by(from, to, from_site,to_site) %>%
    summarise(number.c=sum(number))
  ##
  consistence2<-select(stable,c("OTU", "BestHost", 'site'))
  consistence2 <- consistence2[consistence2$OTU %in% unique(consistence2$OTU[duplicated(consistence2$OTU)]),]
  consistence2$BestHost<-paste0(consistence2$BestHost,consistence2$site)
  consistence_new <- consistence2 %>%
    group_by(OTU) %>%
    summarise(
      combos = list(combn(unique(BestHost), 2, simplify = FALSE)),
      .groups = "drop"
    ) %>%
    unnest(combos) %>%
    mutate(
      combos = map(combos, sort)  # sort each pair alphabetically
    ) %>%
    mutate(
      from = map_chr(combos, 1),
      to   = map_chr(combos, 2)
    ) %>%
    select(OTU, from, to)
  
  consistence_new$from_site<-substr(consistence_new$from, (nchar(consistence_new $from)-1),nchar(consistence_new $from))
  consistence_new$to_site<-substr(consistence_new$to, (nchar(consistence_new$to)-1),nchar(consistence_new $to))
  consistence_new$from<-gsub("LZ|LW|LV","",consistence_new$from)
  consistence_new$to<-gsub("LZ|LW|LV","",consistence_new$to)
  ###
  consistence_new2 <- consistence_new %>%
    filter(from!=to)
  consistence_new2 <- consistence_new
  consistence_new2$number<-1
  consistence_new2 <- consistence_new2 %>%
    group_by(from, to, from_site,to_site) %>%
    summarise(number.c=sum(number))
  ####changed into the non-indicator once
  tnon.ind<-select(non.ind4,c("OTU", "BestHost", 'site'))
  tnon.ind$BestHost<-paste0(tnon.ind$BestHost,tnon.ind$site)
  tnon.ind_new <- tnon.ind %>%
    group_by(OTU) %>%
    summarise(
      combos = list(combn(unique(BestHost), 2, simplify = FALSE)),
      .groups = "drop"
    ) %>%
    unnest(combos) %>%
    mutate(
      combos = map(combos, sort)  # sort each pair alphabetically
    ) %>%
    mutate(
      from = map_chr(combos, 1),
      to   = map_chr(combos, 2)
    ) %>%
    select(OTU, from, to)
  
  tnon.ind_new$from_site<-substr(tnon.ind_new$from, (nchar(tnon.ind_new $from)-1),nchar(tnon.ind_new $from))
  tnon.ind_new$to_site<-substr(tnon.ind_new$to, (nchar(tnon.ind_new$to)-1),nchar(tnon.ind_new $to))
  tnon.ind_new$from<-gsub("LZ|LW|LV","",tnon.ind_new$from)
  tnon.ind_new$to<-gsub("LZ|LW|LV","",tnon.ind_new$to)
  ###
  tnon.ind_new2 <- tnon.ind_new %>%
    filter(from!=to)
  tnon.ind_new2 <- tnon.ind_new 
  tnon.ind_new2$number<-1
  tnon.ind_new2 <- tnon.ind_new2 %>%
    group_by(from, to, from_site,to_site) %>%
    summarise(number.c=sum(number))
  ###
  switch_new<-rbind(switch_new2,consistence_new2,tnon.ind_new2)
  ##
  switch_new <- switch_new %>%
    group_by(from, to) %>%
    summarise(number.c=mean(number.c))
  switch_new24[[i]]<-switch_new
  ####
  switch2$number<-1
  switch2$BestHost<-gsub("LZ|LV|LW","",switch2$BestHost)
  switch.p<-switch2 %>%
    group_by(BestHost,site) %>%
    summarise(number.s=sum(number))
  stable<-stable[,-24]
  consistence<-rbind(stable,non.ind4)
  consistence$number<-1
  consistence.p<-consistence %>%
    group_by(BestHost,site) %>%
    summarise(number.c=sum(number))
  all2<-select(all,c("OTU", "BestHost", 'site'))
  all2<-all2 %>%
    distinct(OTU, BestHost,site)
  all2$number<-1
  all.p<-all2 %>%
    group_by(BestHost,site) %>%
    summarise(number.a=sum(number))
  ###
  non.i.d2<-non.ind3
  non.i.d2$number<-1
  non.i.d.p<-non.i.d2 %>%
    group_by(BestHost,site) %>%
    summarise(number.is=sum(number))
  site.d2<-site.d
  site.d2<-site.d2[site.d2$p.value<0.05,]
  site.d2 $number<-1
  site.d.p<-site.d2 %>%
    group_by(BestHost,site) %>%
    summarise(number.sd=sum(number))
  ##
  me<-function(x,y){
    a<-merge(x,y,by=c("BestHost",'site'),all=T)
    return(a)
  }
  try<-list(switch.p,consistence.p,non.i.d.p,site.d.p,all.p)
  try<-Reduce(me,try)
  try[, 3:7][is.na(try[, 3:7])] <- 0
  try$s.r<-try$number.s/try$number.a
  try$c.r<-try$number.c/try$number.a
  #try$ic.r<-try$number.ic/try$number.a
  try$is.r<-try$number.is/try$number.a##non-indicator
  try$sd.r<-try$number.sd/try$number.a###site unique
  try<-try[!is.na(try$s.r),]##filter out data without significant indicator
  try4[[i]]<-try
}
try<-c(try1, try2, try3, try4)
consistence3<-c(consistence31, consistence32,consistence33, consistence34)
switch3 <- c(switch31, switch32, switch33, switch34)
switch_new2<- c(switch_new1, switch_new22, switch_new23, switch_new24)

for (i in seq_along(consistence3)) {
  consistence3[[i]]$iteration <- i
  consistence3[[i]]$BestHost <- gsub("^s\\.", "", consistence3[[i]]$BestHost)
}
consistence3 <-bind_rows(consistence3)


for (i in seq_along(switch3)) {
  switch3[[i]]$iteration <- i
  switch3[[i]]$BestHost <- gsub("^s\\.", "", switch3[[i]]$BestHost)
}
switch3<-bind_rows(switch3)


switch_new2<-bind_rows(switch_new2)
try<-bind_rows(try)
setDT(try)
try <- try[, lapply(.SD, mean), 
           by = .(BestHost, site), 
           .SDcols = 3:ncol(try)]
write.csv(try, "switching.csv")

setDT(switch_new2)
switch_new222 <- switch_new2[, 
                             .(number.c = mean(number.c)), 
                             by = .(from, to)
]
####plot
switch_new<-switch_new222
samples <- sort(unique(c(switch_new$from, switch_new$to)))
similarity_matrix <- matrix(NA, nrow = length(samples), ncol = length(samples),
                            dimnames = list(samples, samples))
for (i in seq_len(nrow(switch_new))) {
  r <- switch_new$from[i]
  c <- switch_new$to[i]
  v <- switch_new$number.c[i]
  similarity_matrix[r, c] <- v
  similarity_matrix[c, r] <- v  # Fill symmetric cell with the same value
}
similarity_matrix <- as.matrix(similarity_matrix)
similarity_matrix[is.na(similarity_matrix)]<-0
p<-Heatmap(similarity_matrix, name = "Ratio",
           cluster_rows = T,
           cluster_columns = TRUE, 
           col = colorRamp2(c(0,1,5, 10, 20), c("white","lightblue","#FFEEA0","#F4BB44","#C04000")),
           row_names_gp = gpar(fontsize = 8),
           column_names_gp = gpar(fontsize = 9))
########frequency
load("tableall.RData")
all.try1<-list()
for (i in 1:250) {
  obj<-table12[[i]]
  if(inherits(obj, "dgCMatrix")){
    fungi <- as.data.frame(as.matrix(obj))
  }
  fungi$sample_names<-row.names(fungi)
  fungi2<-merge(fungi,metadata,by="sample_names",all=T)
  row.names(fungi2)<-fungi2$sample_names
  fungi2<-fungi2[,-1]
  fungi2<-fungi2[!is.na(fungi2$ff04207f066a2b2769b796c7e8da74227ea12fe4),]
  all2<-aggregate(fungi2[,1:(ncol(fungi2)-2)],by=list(substrate=fungi2$substrate,site=fungi2$site2),sum)
  all <- all2 %>%
    group_by(site) %>%
    mutate(site_total = sum(across(where(is.numeric), sum))) %>%  # total per site
    mutate(across(where(is.numeric), ~ . / site_total)) %>%
    select(-site_total)  # remove temporary column
  ##calculate the frequency
  fre <- all %>%
    group_by(site) %>%
    summarise(across(where(is.numeric), ~ sum(. > 0)))
  ##
  library(tidyverse)
  all2_long <- all %>%
    pivot_longer(cols = where(is.numeric), names_to = "OTU", values_to = "Ratio")
  
  all2_fre <- fre %>%
    pivot_longer(cols = where(is.numeric), names_to = "OTU", values_to = "frequency")
  
  all2_long$me<-paste0(all2_long$site,all2_long$OTU)
  all2_fre$me<-paste0(all2_fre$site,all2_fre$OTU)
  all.try<-merge(all2_long,all2_fre,by="me",all=T)
  all.try<-all.try[all.try$Ratio>0,]
  all.try<-all.try[all.try$frequency>0,]
  all.try<-all.try[,c(2:5,8)]
  names(all.try)[2:3]<-c("site", "OTU") 
  all.try1[[i]]<-all.try
}

all.try2<-list()
for (i in 1:250) {
  obj<-table22[[i]]
  if(inherits(obj, "dgCMatrix")){
    fungi <- as.data.frame(as.matrix(obj))
  }
  fungi$sample_names<-row.names(fungi)
  fungi2<-merge(fungi,metadata,by="sample_names",all=T)
  row.names(fungi2)<-fungi2$sample_names
  fungi2<-fungi2[,-1]
  fungi2<-fungi2[!is.na(fungi2$ff04207f066a2b2769b796c7e8da74227ea12fe4),]
  all2<-aggregate(fungi2[,1:(ncol(fungi2)-2)],by=list(substrate=fungi2$substrate,site=fungi2$site2),sum)
  all <- all2 %>%
    group_by(site) %>%
    mutate(site_total = sum(across(where(is.numeric), sum))) %>%  # total per site
    mutate(across(where(is.numeric), ~ . / site_total)) %>%
    select(-site_total)  # remove temporary column
  ##calculate the frequency
  fre <- all %>%
    group_by(site) %>%
    summarise(across(where(is.numeric), ~ sum(. > 0)))
  ##
  all2_long <- all %>%
    pivot_longer(cols = where(is.numeric), names_to = "OTU", values_to = "Ratio")
  
  all2_fre <- fre %>%
    pivot_longer(cols = where(is.numeric), names_to = "OTU", values_to = "frequency")
  
  all2_long$me<-paste0(all2_long$site,all2_long$OTU)
  all2_fre$me<-paste0(all2_fre$site,all2_fre$OTU)
  all.try<-merge(all2_long,all2_fre,by="me",all=T)
  all.try<-all.try[all.try$Ratio>0,]
  all.try<-all.try[all.try$frequency>0,]
  all.try<-all.try[,c(2:5,8)]
  names(all.try)[2:3]<-c("site", "OTU") 
  all.try2[[i]]<-all.try
}

all.try3<-list()
for (i in 1:250) {
  obj<-table32[[i]]
  if(inherits(obj, "dgCMatrix")){
    fungi <- as.data.frame(as.matrix(obj))
  }
  fungi$sample_names<-row.names(fungi)
  fungi2<-merge(fungi,metadata,by="sample_names",all=T)
  row.names(fungi2)<-fungi2$sample_names
  fungi2<-fungi2[,-1]
  fungi2<-fungi2[!is.na(fungi2$ff04207f066a2b2769b796c7e8da74227ea12fe4),]
  all2<-aggregate(fungi2[,1:(ncol(fungi2)-2)],by=list(substrate=fungi2$substrate,site=fungi2$site2),sum)
  all <- all2 %>%
    group_by(site) %>%
    mutate(site_total = sum(across(where(is.numeric), sum))) %>%  # total per site
    mutate(across(where(is.numeric), ~ . / site_total)) %>%
    select(-site_total)  # remove temporary column
  ##calculate the frequency
  fre <- all %>%
    group_by(site) %>%
    summarise(across(where(is.numeric), ~ sum(. > 0)))
  ##
  all2_long <- all %>%
    pivot_longer(cols = where(is.numeric), names_to = "OTU", values_to = "Ratio")
  
  all2_fre <- fre %>%
    pivot_longer(cols = where(is.numeric), names_to = "OTU", values_to = "frequency")
  
  all2_long$me<-paste0(all2_long$site,all2_long$OTU)
  all2_fre$me<-paste0(all2_fre$site,all2_fre$OTU)
  all.try<-merge(all2_long,all2_fre,by="me",all=T)
  all.try<-all.try[all.try$Ratio>0,]
  all.try<-all.try[all.try$frequency>0,]
  all.try<-all.try[,c(2:5,8)]
  names(all.try)[2:3]<-c("site", "OTU") 
  all.try3[[i]]<-all.try
}

all.try4<-list()
for (i in 1:250) {
  obj<-table42[[i]]
  if(inherits(obj, "dgCMatrix")){
    fungi <- as.data.frame(as.matrix(obj))
  }
  fungi$sample_names<-row.names(fungi)
  fungi2<-merge(fungi,metadata,by="sample_names",all=T)
  row.names(fungi2)<-fungi2$sample_names
  fungi2<-fungi2[,-1]
  fungi2<-fungi2[!is.na(fungi2$ff04207f066a2b2769b796c7e8da74227ea12fe4),]
  all2<-aggregate(fungi2[,1:(ncol(fungi2)-2)],by=list(substrate=fungi2$substrate,site=fungi2$site2),sum)
  all <- all2 %>%
    group_by(site) %>%
    mutate(site_total = sum(across(where(is.numeric), sum))) %>%  # total per site
    mutate(across(where(is.numeric), ~ . / site_total)) %>%
    select(-site_total)  # remove temporary column
  ##calculate the frequency
  fre <- all %>%
    group_by(site) %>%
    summarise(across(where(is.numeric), ~ sum(. > 0)))
  ##
  all2_long <- all %>%
    pivot_longer(cols = where(is.numeric), names_to = "OTU", values_to = "Ratio")
  
  all2_fre <- fre %>%
    pivot_longer(cols = where(is.numeric), names_to = "OTU", values_to = "frequency")
  
  all2_long$me<-paste0(all2_long$site,all2_long$OTU)
  all2_fre$me<-paste0(all2_fre$site,all2_fre$OTU)
  all.try<-merge(all2_long,all2_fre,by="me",all=T)
  all.try<-all.try[all.try$Ratio>0,]
  all.try<-all.try[all.try$frequency>0,]
  all.try<-all.try[,c(2:5,8)]
  names(all.try)[2:3]<-c("site", "OTU") 
  all.try4[[i]]<-all.try
}
######
consistence2<-consistence3
names(consistence2)[2]<-"substrate"
consistence2$type<-"consistant"

switch2<-switch3
names(switch2)[2]<-"substrate"
switch2$type<-"switching"

type<-rbind(consistence2,switch2)

all.try <- c(all.try1, all.try2, all.try3, all.try4)

for (i in seq_along(all.try)) {
  all.try[[i]]$iteration <- i
}
all.try <- bind_rows(all.try)
setDT(type)
setDT(all.try)
setkey(type, substrate, site, OTU, iteration)
setkey(all.try, substrate, site, OTU, iteration)

try <- merge(type, all.try, all = TRUE, sort = FALSE)
try <- try[!is.na(type)]
try<-try[!is.na(try$frequency),]
library(marginaleffects)
models <- lapply(1:1000, function(i) {
  df <- try[try$iteration == i, ]
 lm(frequency ~ substrate + type + site, data = df)
})

out<-list()
for (i in 1:1000) {
  out[[i]] <- avg_comparisons(models[[i]], variables =  list(type  = "pairwise"))
}

df <- do.call(rbind, out)
df_sig <- df[df$p.value < 0.05, ]
p_robust <- nrow(df_sig) / nrow(df)
p_robust
quantile(df$estimate)

out<-list()
for (i in 1:1000) {
  out[[i]] <- summary(models[[i]])$r.squared
}
r <- do.call(rbind, out)
quantile(r, c(0.025, 0.5, 0.975))

#####
models <- lapply(1:1000, function(i) {
  df <- try[try$iteration == i, ]
  lm(stat ~ substrate + type + site, data = df)
})

out<-list()
for (i in 1:1000) {
  out[[i]] <- avg_comparisons(models[[i]], variables =  list(type  = "pairwise"))
}

df <- do.call(rbind, out)
df_sig <- df[df$p.value < 0.05, ]
p_robust <- nrow(df_sig) / nrow(df)
p_robust
quantile(df$estimate)
out<-list()
for (i in 1:1000) {
  out[[i]] <- summary(models[[i]])$r.squared
}
r <- do.call(rbind, out)
quantile(r, c(0.025, 0.5, 0.975))

######difference between single lifestyle and twe lifestyles in substrate range
load("tableall.RData")
for (i in 1:250) {
  result1[[i]]<-result1[[i]][result1[[i]]$p.value < 0.05, ]
 
}
for (i in 1:250) {
  result2[[i]]<-result2[[i]][result2[[i]]$p.value < 0.05, ]
  
}

for (i in 1:250) {
  result3[[i]]<-result3[[i]][result3[[i]]$p.value < 0.05, ]
  
}

for (i in 1:250) {
  result4[[i]]<-result4[[i]][result4[[i]]$p.value < 0.05, ]
  
}

###
fre1<-list()
for (i in 1:250) {
  fungi<-table12[[i]]
  if(inherits(fungi, "dgCMatrix")){
    fungi <- as.data.frame(as.matrix(fungi))
  }
  fungi$sample_names<-row.names(fungi)
  fungi2<-merge(fungi,metadata,by="sample_names",all=T)
  row.names(fungi2)<-fungi2$sample_names
  fungi2<-fungi2[,-1]
  fungi2<-fungi2[!is.na(fungi2$ff04207f066a2b2769b796c7e8da74227ea12fe4),]
  all2<-aggregate(fungi2[,1:(ncol(fungi2)-2)],by=list(substrate=fungi2$substrate,site=fungi2$site2),sum)
  all <- all2 %>%
    group_by(site) %>%
    mutate(site_total = sum(across(where(is.numeric)))) %>%  # total per site
    mutate(across(where(is.numeric), ~ . / site_total)) %>%
    select(-site_total)  # remove temporary column
  fre1[[i]] <- all2 %>%
    group_by(site) %>%
    summarise(across(where(is.numeric), ~ sum(. > 0)))
}

all2_long1<-list()
all2_fre1<-list()
for (i in 1:250) {
  fungi<-table12[[i]]
  if(inherits(fungi, "dgCMatrix")){
    fungi <- as.data.frame(as.matrix(fungi))
  }
  fungi$sample_names<-row.names(fungi)
  fungi2<-merge(fungi,metadata,by="sample_names",all=T)
  row.names(fungi2)<-fungi2$sample_names
  fungi2<-fungi2[,-1]
  fungi2<-fungi2[!is.na(fungi2$ff04207f066a2b2769b796c7e8da74227ea12fe4),]
  all2<-aggregate(fungi2[,1:(ncol(fungi2)-2)],by=list(substrate=fungi2$substrate,site=fungi2$site2),sum)
  all <- all2 %>%
    group_by(site) %>%
    mutate(site_total = sum(across(where(is.numeric)))) %>%  # total per site
    mutate(across(where(is.numeric), ~ . / site_total)) %>%
    select(-site_total)  # remove temporary column
  fre <- all2 %>%
    group_by(site) %>%
    summarise(across(where(is.numeric), ~ sum(. > 0)))
  all2_long1[[i]] <- all %>%
    pivot_longer(cols = where(is.numeric), names_to = "OTU", values_to = "Ratio")
  
  all2_fre1[[i]] <- fre %>%
    pivot_longer(cols = where(is.numeric), names_to = "OTU", values_to = "frequency")
}

all2_long2<-list()
all2_fre2<-list()
for (i in 1:250) {
  fungi<-table12[[i]]
  if(inherits(fungi, "dgCMatrix")){
    fungi <- as.data.frame(as.matrix(fungi))
  }
  fungi$sample_names<-row.names(fungi)
  fungi2<-merge(fungi,metadata,by="sample_names",all=T)
  row.names(fungi2)<-fungi2$sample_names
  fungi2<-fungi2[,-1]
  fungi2<-fungi2[!is.na(fungi2$ff04207f066a2b2769b796c7e8da74227ea12fe4),]
  all2<-aggregate(fungi2[,1:(ncol(fungi2)-2)],by=list(substrate=fungi2$substrate,site=fungi2$site2),sum)
  all <- all2 %>%
    group_by(site) %>%
    mutate(site_total = sum(across(where(is.numeric)))) %>%  # total per site
    mutate(across(where(is.numeric), ~ . / site_total)) %>%
    select(-site_total)  # remove temporary column
  fre <- all2 %>%
    group_by(site) %>%
    summarise(across(where(is.numeric), ~ sum(. > 0)))
  library(tidyverse)
  all2_long2[[i]] <- all %>%
    pivot_longer(cols = where(is.numeric), names_to = "OTU", values_to = "Ratio")
  
  all2_fre2[[i]] <- fre %>%
    pivot_longer(cols = where(is.numeric), names_to = "OTU", values_to = "frequency")
}

all2_long3<-list()
all2_fre3<-list()
for (i in 1:250) {
  fungi<-table12[[i]]
  if(inherits(fungi, "dgCMatrix")){
    fungi <- as.data.frame(as.matrix(fungi))
  }
  fungi$sample_names<-row.names(fungi)
  fungi2<-merge(fungi,metadata,by="sample_names",all=T)
  row.names(fungi2)<-fungi2$sample_names
  fungi2<-fungi2[,-1]
  fungi2<-fungi2[!is.na(fungi2$ff04207f066a2b2769b796c7e8da74227ea12fe4),]
  all2<-aggregate(fungi2[,1:(ncol(fungi2)-2)],by=list(substrate=fungi2$substrate,site=fungi2$site2),sum)
  all <- all2 %>%
    group_by(site) %>%
    mutate(site_total = sum(across(where(is.numeric)))) %>%  # total per site
    mutate(across(where(is.numeric), ~ . / site_total)) %>%
    select(-site_total)  # remove temporary column
  fre <- all2 %>%
    group_by(site) %>%
    summarise(across(where(is.numeric), ~ sum(. > 0)))
  all2_long3[[i]] <- all %>%
    pivot_longer(cols = where(is.numeric), names_to = "OTU", values_to = "Ratio")
  
  all2_fre3[[i]] <- fre %>%
    pivot_longer(cols = where(is.numeric), names_to = "OTU", values_to = "frequency")
}


all2_long4<-list()
all2_fre4<-list()
for (i in 1:250) {
  fungi<-table12[[i]]
  if(inherits(fungi, "dgCMatrix")){
    fungi <- as.data.frame(as.matrix(fungi))
  }
  fungi$sample_names<-row.names(fungi)
  fungi2<-merge(fungi,metadata,by="sample_names",all=T)
  row.names(fungi2)<-fungi2$sample_names
  fungi2<-fungi2[,-1]
  fungi2<-fungi2[!is.na(fungi2$ff04207f066a2b2769b796c7e8da74227ea12fe4),]
  all2<-aggregate(fungi2[,1:(ncol(fungi2)-2)],by=list(substrate=fungi2$substrate,site=fungi2$site2),sum)
  all <- all2 %>%
    group_by(site) %>%
    mutate(site_total = sum(across(where(is.numeric)))) %>%  # total per site
    mutate(across(where(is.numeric), ~ . / site_total)) %>%
    select(-site_total)  # remove temporary column
  fre <- all2 %>%
    group_by(site) %>%
    summarise(across(where(is.numeric), ~ sum(. > 0)))
  all2_long4[[i]] <- all %>%
    pivot_longer(cols = where(is.numeric), names_to = "OTU", values_to = "Ratio")
  
  all2_fre4[[i]] <- fre %>%
    pivot_longer(cols = where(is.numeric), names_to = "OTU", values_to = "frequency")
}
####calculate the frequency
###merging the lifestyle
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
###
table2<-table[,c(1,2,9,10)]
table2$Secondary_lifestyle[is.na(table2$Secondary_lifestyle)]<-"unknown"
table2$Secondary_lifestyle[table2$Secondary_lifestyle==""]<-"unknown"
table2$qseqid<-row.names(table2)
all1<-list()
for (i in 1:250) {
  all2_fre1[[i]] <- all2_fre1[[i]] [all2_fre1[[i]] $frequency>0,]
  all<-merge(table2,all2_fre1[[i]],by.x="qseqid",by.y="OTU",all=T)
  all<-all[!is.na(all$frequency),]
  all$primary_lifestyle[is.na(all$primary_lifestyle)]<-"unknown"
  all$primary_lifestyle[all$primary_lifestyle==""]<-"unknown"
  all$Secondary_lifestyle[is.na(all$Secondary_lifestyle)]<-"unknown"
  all$Secondary_lifestyle[all$Secondary_lifestyle==""]<-"unknown"
  all$Secondary_lifestyle<-ifelse(all$Secondary_lifestyle=="unknown",0,1)
  all$primary_lifestyle<-ifelse(all$primary_lifestyle=="unknown",0,1)
  all$type<-all$primary_lifestyle+all$Secondary_lifestyle
  all<-all[all$type>0,]
  all$type<-ifelse(all$type==1,"one","second")
  all1[[i]]<-all
}

all1<-list()
for (i in 1:250) {
  all2_fre1[[i]] <- all2_fre1[[i]] [all2_fre1[[i]] $frequency>0,]
  all<-merge(table2,all2_fre1[[i]],by.x="qseqid",by.y="OTU",all=T)
  all<-all[!is.na(all$frequency),]
  all$primary_lifestyle[is.na(all$primary_lifestyle)]<-"unknown"
  all$primary_lifestyle[all$primary_lifestyle==""]<-"unknown"
  all$Secondary_lifestyle[is.na(all$Secondary_lifestyle)]<-"unknown"
  all$Secondary_lifestyle[all$Secondary_lifestyle==""]<-"unknown"
  all$Secondary_lifestyle<-ifelse(all$Secondary_lifestyle=="unknown",0,1)
  all$primary_lifestyle<-ifelse(all$primary_lifestyle=="unknown",0,1)
  all$type<-all$primary_lifestyle+all$Secondary_lifestyle
  all<-all[all$type>0,]
  all$type<-ifelse(all$type==1,"one","second")
  all$iteration<-i
  all1[[i]]<-all
}

all2<-list()
for (i in 1:250) {
  all2_fre2[[i]] <- all2_fre2[[i]] [all2_fre2[[i]] $frequency>0,]
  all<-merge(table2,all2_fre2[[i]],by.x="qseqid",by.y="OTU",all=T)
  all<-all[!is.na(all$frequency),]
  all$primary_lifestyle[is.na(all$primary_lifestyle)]<-"unknown"
  all$primary_lifestyle[all$primary_lifestyle==""]<-"unknown"
  all$Secondary_lifestyle[is.na(all$Secondary_lifestyle)]<-"unknown"
  all$Secondary_lifestyle[all$Secondary_lifestyle==""]<-"unknown"
  all$Secondary_lifestyle<-ifelse(all$Secondary_lifestyle=="unknown",0,1)
  all$primary_lifestyle<-ifelse(all$primary_lifestyle=="unknown",0,1)
  all$type<-all$primary_lifestyle+all$Secondary_lifestyle
  all<-all[all$type>0,]
  all$type<-ifelse(all$type==1,"one","second")
  all$iteration<-i+250
  all2[[i]]<-all
}

all3<-list()
for (i in 1:250) {
  all2_fre3[[i]] <- all2_fre3[[i]] [all2_fre3[[i]] $frequency>0,]
  all<-merge(table2,all2_fre3[[i]],by.x="qseqid",by.y="OTU",all=T)
  all<-all[!is.na(all$frequency),]
  all$primary_lifestyle[is.na(all$primary_lifestyle)]<-"unknown"
  all$primary_lifestyle[all$primary_lifestyle==""]<-"unknown"
  all$Secondary_lifestyle[is.na(all$Secondary_lifestyle)]<-"unknown"
  all$Secondary_lifestyle[all$Secondary_lifestyle==""]<-"unknown"
  all$Secondary_lifestyle<-ifelse(all$Secondary_lifestyle=="unknown",0,1)
  all$primary_lifestyle<-ifelse(all$primary_lifestyle=="unknown",0,1)
  all$type<-all$primary_lifestyle+all$Secondary_lifestyle
  all<-all[all$type>0,]
  all$type<-ifelse(all$type==1,"one","second")
  all$iteration<-i+500
  all3[[i]]<-all
}

all4<-list()
for (i in 1:250) {
  all2_fre4[[i]] <- all2_fre4[[i]] [all2_fre4[[i]] $frequency>0,]
  all<-merge(table2,all2_fre4[[i]],by.x="qseqid",by.y="OTU",all=T)
  all<-all[!is.na(all$frequency),]
  all$primary_lifestyle[is.na(all$primary_lifestyle)]<-"unknown"
  all$primary_lifestyle[all$primary_lifestyle==""]<-"unknown"
  all$Secondary_lifestyle[is.na(all$Secondary_lifestyle)]<-"unknown"
  all$Secondary_lifestyle[all$Secondary_lifestyle==""]<-"unknown"
  all$Secondary_lifestyle<-ifelse(all$Secondary_lifestyle=="unknown",0,1)
  all$primary_lifestyle<-ifelse(all$primary_lifestyle=="unknown",0,1)
  all$type<-all$primary_lifestyle+all$Secondary_lifestyle
  all<-all[all$type>0,]
  all$type<-ifelse(all$type==1,"one","second")
  all$iteration<-i+750
  all4[[i]]<-all
}
all<-bind_rows(all1,all2,all3,all4)
models <- lapply(1:1000, function(i) {
  df <- all[all$iteration == i, ]
  lm(frequency ~ type + site, data = df)
})

out<-list()
for (i in 1:1000) {
  out[[i]] <- avg_comparisons(models[[i]], variables =  list(type  = "pairwise"))
}

df <- do.call(rbind, out)
df_sig <- df[df$p.value < 0.05, ]
p_robust <- nrow(df_sig) / nrow(df)
p_robust
quantile(df$estimate)
out<-list()
for (i in 1:1000) {
  out[[i]] <- summary(models[[i]])$r.squared
}
r <- do.call(rbind, out)
quantile(r, c(0.025, 0.5, 0.975))
######calculate the unexpected fungi proportion
load("tableall.RData")
library(dplyr)
library(tidyverse)
library(data.table)
table_long1 <- vector("list", 250)

for(i in seq_along(table12)) {
  
  fungi <- table12[[i]]
  
  if (inherits(fungi, "dgCMatrix")) {
    dt <- as.data.table(as.matrix(fungi))
  } else {
    dt <- as.data.table(fungi)
  }

  dt[, sample_names := rownames(fungi)]
  dt <- merge(dt, metadata, by="sample_names")
  dt2 <- dt %>%
    group_by(site2,sample_names) %>%
    mutate(group_sum = sum(across(where(is.numeric)), na.rm = TRUE)) %>%
    mutate(across(where(is.numeric), ~ .x / group_sum)) %>%
    select(-group_sum)
  otu_cols <- setdiff(names(dt2),
                      c("sample_names", "site2", "substrate"))
  dt2<-as.data.table(dt2)
  dt2 <- dt2[, lapply(.SD, mean),
             by = .(substrate, site = site2),
             .SDcols = otu_cols]
  
  dt_long <- melt(dt2, id.vars=c("substrate","site"),
                  variable.name="OTU", value.name="reads")
  
  table_long1[[i]] <- dt_long[, .(relative = mean(reads)),
                              by=.(substrate, OTU)]
}

table_long2 <- vector("list", 250)

for(i in seq_along(table22)) {
  
  fungi <- table22[[i]]
  
  if (inherits(fungi, "dgCMatrix")) {
    dt <- as.data.table(as.matrix(fungi))
  } else {
    dt <- as.data.table(fungi)
  }

  dt[, sample_names := rownames(fungi)]
  dt <- merge(dt, metadata, by="sample_names")
  dt2 <- dt %>%
    group_by(site2,sample_names) %>%
    mutate(group_sum = sum(across(where(is.numeric)), na.rm = TRUE)) %>%
    mutate(across(where(is.numeric), ~ .x / group_sum)) %>%
    select(-group_sum)
  otu_cols <- setdiff(names(dt2),
                      c("sample_names", "site2", "substrate"))
  dt2<-as.data.table(dt2)
  dt2 <- dt2[, lapply(.SD, mean),
             by = .(substrate, site = site2),
             .SDcols = otu_cols]
  
  dt_long <- melt(dt2, id.vars=c("substrate","site"),
                  variable.name="OTU", value.name="reads")
  
  
  table_long2[[i]] <- dt_long[, .(relative = mean(reads)),
                              by=.(substrate, OTU)]
}

table_long3 <- vector("list", 250)

for(i in seq_along(table32)) {
  
  fungi <- table32[[i]]
  
  if (inherits(fungi, "dgCMatrix")) {
    dt <- as.data.table(as.matrix(fungi))
  } else {
    dt <- as.data.table(fungi)
  }
  
  dt[, sample_names := rownames(fungi)]
  dt <- merge(dt, metadata, by="sample_names")
  dt2 <- dt %>%
    group_by(site2,sample_names) %>%
    mutate(group_sum = sum(across(where(is.numeric)), na.rm = TRUE)) %>%
    mutate(across(where(is.numeric), ~ .x / group_sum)) %>%
    select(-group_sum)
  otu_cols <- setdiff(names(dt2),
                      c("sample_names", "site2", "substrate"))
  dt2<-as.data.table(dt2)
  dt2 <- dt2[, lapply(.SD, mean),
               by = .(substrate, site = site2),
               .SDcols = otu_cols]
 
  dt_long <- melt(dt2, id.vars=c("substrate","site"),
                  variable.name="OTU", value.name="reads")
  
  table_long3[[i]] <- dt_long[, .(relative = mean(reads)),
                              by=.(substrate, OTU)]
}

table_long4 <- vector("list", 250)

for(i in seq_along(table42)) {
  
  fungi <- table42[[i]]
  
  if (inherits(fungi, "dgCMatrix")) {
    dt <- as.data.table(as.matrix(fungi))
  } else {
    dt <- as.data.table(fungi)
  }
  
  
  dt[, sample_names := rownames(fungi)]
  dt <- merge(dt, metadata, by="sample_names")
  dt2 <- dt %>%
    group_by(site2,sample_names) %>%
    mutate(group_sum = sum(across(where(is.numeric)), na.rm = TRUE)) %>%
    mutate(across(where(is.numeric), ~ .x / group_sum)) %>%
    select(-group_sum)
  otu_cols <- setdiff(names(dt2),
                      c("sample_names", "site2", "substrate"))
  dt2<-as.data.table(dt2)
  dt2 <- dt2[, lapply(.SD, mean),
             by = .(substrate, site = site2),
             .SDcols = otu_cols]
  
  dt_long <- melt(dt2, id.vars=c("substrate","site"),
                  variable.name="OTU", value.name="reads")
  
  table_long4[[i]] <- dt_long[, .(relative = mean(reads)),
                              by=.(substrate, OTU)]
}
###
table_long12<-list()
for (i in 1:250) {
  result1[[i]]<-result1[[i]][result1[[i]]$p.value < 0.05, ]
  names(result1[[i]])[2]<-"substrate"
  result1[[i]]$substrate<-gsub("^s\\.","",result1[[i]]$substrate)
  result1[[i]]$substrate<-gsub("\\."," ",result1[[i]]$substrate)
  table_long12[[i]]<-merge(result1[[i]][,c(1,2,4)],table_long1[[i]],by=c("substrate","OTU"))
}
table_long22<-list()
for (i in 1:250) {
  result2[[i]]<-result2[[i]][result2[[i]]$p.value < 0.05, ]
  names(result2[[i]])[2]<-"substrate"
  result2[[i]]$substrate<-gsub("^s\\.","",result2[[i]]$substrate)
  result2[[i]]$substrate<-gsub("\\."," ",result2[[i]]$substrate)
  table_long22[[i]]<-merge(result2[[i]][,c(1,2,4)],table_long2[[i]],by=c("substrate","OTU"))
}

table_long32<-list()
for (i in 1:250) {
  result3[[i]]<-result3[[i]][result3[[i]]$p.value < 0.05, ]
  names(result3[[i]])[2]<-"substrate"
  result3[[i]]$substrate<-gsub("^s\\.","",result3[[i]]$substrate)
  result3[[i]]$substrate<-gsub("\\."," ",result3[[i]]$substrate)
  table_long32[[i]]<-merge(result3[[i]][,c(1,2,4)],table_long3[[i]],by=c("substrate","OTU"))
}

table_long42<-list()
for (i in 1:250) {
  result4[[i]]<-result4[[i]][result4[[i]]$p.value < 0.05, ]
  names(result4[[i]])[2]<-"substrate"
  result4[[i]]$substrate<-gsub("^s\\.","",result4[[i]]$substrate)
  result4[[i]]$substrate<-gsub("\\."," ",result4[[i]]$substrate)
  table_long42[[i]]<-merge(result4[[i]][,c(1,2,4)],table_long4[[i]],by=c("substrate","OTU"))
}

save.image("proportion.relative.abundance.RData")
###
table_long<- c(table_long12, table_long22, table_long32,table_long42)
table_long<- bind_rows(table_long)

table_long <- table_long %>%
  group_by(substrate, OTU) %>%
  filter(abs(relative - mean(relative)) == min(abs(relative - mean(relative))))

model <- lm(stat ~ relative + substrate, data = table_long)
performance::performance(model)
parameters::parameters(model)
####indicators belongs to unexpected guilds
trait<-read.csv("FungalTraits.csv")
tax<-read.csv("taxonomy.final.csv")
tax$GENUS<-gsub(".*__","",tax$genus)
un.table1 <-list()
for (i in 1:250) {
  result1[[i]]<-result1[[i]][,c(1,2,5,22)]
  names(result1[[i]])[4]<-"GENUS"
  result1[[i]]$GENUS[(result1[[i]]$GENUS=="")]<-"unclassified"
  table<-merge(result1[[i]],trait[,c(6,8,9)],by="GENUS",all=T)
  table<-table[!is.na(table$substrate),]
  table[is.na(table$primary_lifestyle),]$primary_lifestyle<-"unknown"
  table[is.na(table$Secondary_lifestyle ),]$Secondary_lifestyle<-"unknown"
  un.table1[[i]]<-table
}

un.table2 <-list()
for (i in 1:250) {
  result2[[i]]<-result2[[i]][,c(1,2,5,22)]
  names(result2[[i]])[4]<-"GENUS"
  result2[[i]]$GENUS[(result2[[i]]$GENUS=="")]<-"unclassified"
  table<-merge(result2[[i]],trait[,c(6,8,9)],by="GENUS",all=T)
  table<-table[!is.na(table$substrate),]
  table[is.na(table$primary_lifestyle),]$primary_lifestyle<-"unknown"
  table[is.na(table$Secondary_lifestyle ),]$Secondary_lifestyle<-"unknown"
  un.table2[[i]]<-table
}

un.table3 <-list()
for (i in 1:250) {
  result3[[i]]<-result3[[i]][,c(1,2,5,22)]
  names(result3[[i]])[4]<-"GENUS"
  result3[[i]]$GENUS[(result3[[i]]$GENUS=="")]<-"unclassified"
  table<-merge(result3[[i]],trait[,c(6,8,9)],by="GENUS",all=T)
  table<-table[!is.na(table$substrate),]
  table[is.na(table$primary_lifestyle),]$primary_lifestyle<-"unknown"
  table[is.na(table$Secondary_lifestyle ),]$Secondary_lifestyle<-"unknown"
  un.table3[[i]]<-table
}

un.table4 <-list()
for (i in 1:250) {
  result4[[i]]<-result4[[i]][,c(1,2,5,22)]
  names(result4[[i]])[4]<-"GENUS"
  result4[[i]]$GENUS[(result4[[i]]$GENUS=="")]<-"unclassified"
  table<-merge(result4[[i]],trait[,c(6,8,9)],by="GENUS",all=T)
  table<-table[!is.na(table$substrate),]
  table[is.na(table$primary_lifestyle),]$primary_lifestyle<-"unknown"
  table[is.na(table$Secondary_lifestyle ),]$Secondary_lifestyle<-"unknown"
  un.table4[[i]]<-table
}

unexpected.index <- list(topsoil=c("plant_pathogen", "litter_saprotroph", "wood_saprotroph", "mycoparasite",
             "dung_saprotroph", "lichenized", "sooty_mold", "animal_parasite",
             "unspecified_saprotroph", "pollen_saprotroph", "epiphyte",
             "animal_endosymbiont", "foliar_endophyte", "algal_parasite",       
             "nectar/tap_saprotroph", "rock-inhabiting",            
             "fungal_decomposer","animal_decomposer",         
             "bryophilous",
             "unspecified_symbiotroph",
             "algal_symbiont","unsepcified_saprotroph",     
             "arthropod_parasite","resin_saprotroph"),
subsoil= c("plant_pathogen", "litter_saprotroph", "wood_saprotroph", "mycoparasite",
               "dung_saprotroph", "lichenized", "sooty_mold", "animal_parasite",
               "unspecified_saprotroph", "pollen_saprotroph", "epiphyte",
               "animal_endosymbiont", "foliar_endophyte", "algal_parasite",        
               "nectar/tap_saprotroph", "rock-inhabiting",            
           "fungal_decomposer","animal_decomposer",         
           "bryophilous",
           "unspecified_symbiotroph",
           "algal_symbiont","unsepcified_saprotroph",     
           "arthropod_parasite","resin_saprotroph" ),

lichen = c("soil_saprotroph" , "plant_pathogen", "litter_saprotroph",     
            "ectomycorrhizal", "wood_saprotroph", "mycoparasite",          
            "dung_saprotroph", "sooty_mold", "animal_parasite",       
            "unspecified_saprotroph", "pollen_saprotroph",               
            "arbuscular_mycorrhizal", "animal_endosymbiont",             
            "root_endophyte", "foliar_endophyte",        
            "nectar/tap_saprotroph",
           "root-associated","rock-inhabiting",            
           "fungal_decomposer","animal_decomposer",         
           "bryophilous","root_endophyte_dark_septate",
           "unspecified_symbiotroph","ericoid_mycorrhizal",        
           "unsepcified_saprotroph",     
           "arthropod_parasite","resin_saprotroph"),

leaf.litter= c("soil_saprotroph", "plant_pathogen",      
                   "ectomycorrhizal", "wood_saprotroph","mycoparasite",          
                   "dung_saprotroph","lichenized",            
                   "lichen_parasite", "sooty_mold","animal_parasite",      
                   "unspecified_saprotroph", "pollen_saprotroph", "epiphyte",              
                   "arbuscular_mycorrhizal", "animal_endosymbiont",                    
                   "root_endophyte", "foliar_endophyte","algal_parasite",        
                   "nectar/tap_saprotroph","root-associated","rock-inhabiting",            
               "fungal_decomposer","animal_decomposer",         
               "bryophilous","root_endophyte_dark_septate",
               "unspecified_symbiotroph","ericoid_mycorrhizal",        
               "algal_symbiont","unsepcified_saprotroph",     
               "arthropod_parasite","resin_saprotroph"),

leaves= c("soil_saprotroph", "litter_saprotroph",    
              "ectomycorrhizal","wood_saprotroph", "mycoparasite",         
              "dung_saprotroph", "lichenized",            
              "lichen_parasite","sooty_mold","animal_parasite",      
              "unspecified_saprotroph", "pollen_saprotroph",              
              "arbuscular_mycorrhizal", "animal_endosymbiont",                     
              "root_endophyte","foliar_endophyte","algal_parasite",        
              "nectar/tap_saprotroph","root-associated","rock-inhabiting",            
          "fungal_decomposer","animal_decomposer",         
          "bryophilous","root_endophyte_dark_septate",
          "unspecified_symbiotroph","ericoid_mycorrhizal",        
          "algal_symbiont","unsepcified_saprotroph",     
          "arthropod_parasite","resin_saprotroph"),
wood = c("soil_saprotroph","plant_pathogen","litter_saprotroph",     
"ectomycorrhizal", "mycoparasite",          
"dung_saprotroph","lichenized",            
"lichen_parasite", "sooty_mold", "animal_parasite",       
"unspecified_saprotroph", "pollen_saprotroph", "epiphyte",              
"arbuscular_mycorrhizal", "animal_endosymbiont",                     
"root_endophyte", "foliar_endophyte", "algal_parasite",       
"nectar/tap_saprotroph","root-associated","rock-inhabiting",            
"fungal_decomposer","animal_decomposer",         
"bryophilous","root_endophyte_dark_septate",
"unspecified_symbiotroph","ericoid_mycorrhizal",        
"algal_symbiont","unsepcified_saprotroph",     
"arthropod_parasite","resin_saprotroph"),
moss= c("soil_saprotroph", "litter_saprotroph",     
         "ectomycorrhizal", "wood_saprotroph","mycoparasite",          
         "dung_saprotroph","lichenized",            
         "lichen_parasite","sooty_mold","animal_parasite",      
         "unspecified_saprotroph", "pollen_saprotroph",              
         "arbuscular_mycorrhizal", "animal_endosymbiont",                     
         "root_endophyte", "foliar_endophyte", "algal_parasite",        
         "nectar/tap_saprotroph","root-associated","rock-inhabiting",            
        "fungal_decomposer","animal_decomposer",         
        "root_endophyte_dark_septate",
        "unspecified_symbiotroph","ericoid_mycorrhizal",        
        "algal_symbiont","unsepcified_saprotroph",     
        "arthropod_parasite","resin_saprotroph"),
bark= c("soil_saprotroph", "litter_saprotroph",     
         "ectomycorrhizal", "mycoparasite",          
         "dung_saprotroph", "lichenized",            
         "lichen_parasite", "sooty_mold", "animal_parasite",       
         "unspecified_saprotroph", "pollen_saprotroph", "epiphyte",              
         "arbuscular_mycorrhizal", "animal_endosymbiont",                     
         "root_endophyte", "foliar_endophyte","algal_parasite",        
         "nectar/tap_saprotroph","root-associated","rock-inhabiting",            
        "fungal_decomposer","animal_decomposer",         
        "bryophilous","root_endophyte_dark_septate",
        "unspecified_symbiotroph","ericoid_mycorrhizal",        
        "algal_symbiont","unsepcified_saprotroph",     
        "arthropod_parasite","resin_saprotroph"),
feces =c("soil_saprotroph", "plant_pathogen","litter_saprotroph",     
          "ectomycorrhizal","wood_saprotroph","mycoparasite",          
          "lichenized","lichen_parasite" ,"sooty_mold",      
          "unspecified_saprotroph", "pollen_saprotroph", "epiphyte",              
          "arbuscular_mycorrhizal",              
          "root_endophyte", "foliar_endophyte", "algal_parasite",        
          "nectar/tap_saprotroph","root-associated","rock-inhabiting",            
         "fungal_decomposer","animal_decomposer",         
         "bryophilous","root_endophyte_dark_septate",
         "unspecified_symbiotroph","ericoid_mycorrhizal",        
         "algal_symbiont","unsepcified_saprotroph",     
         "arthropod_parasite","resin_saprotroph"),
fruit.body = c("soil_saprotroph", "plant_pathogen", "litter_saprotroph",     
                "ectomycorrhizal","wood_saprotroph",          
                "dung_saprotroph", "lichenized",            
                "lichen_parasite","sooty_mold","animal_parasite",       
                "unspecified_saprotroph", "pollen_saprotroph","epiphyte",              
                "arbuscular_mycorrhizal", "animal_endosymbiont",            
                "root_endophyte","foliar_endophyte","algal_parasite",        
                "nectar/tap_saprotroph","root-associated","rock-inhabiting",            
               "animal_decomposer",         
               "bryophilous","root_endophyte_dark_septate",
               "unspecified_symbiotroph","ericoid_mycorrhizal",        
               "algal_symbiont","unsepcified_saprotroph",     
               "arthropod_parasite","resin_saprotroph"
),

snow = c("soil_saprotroph","plant_pathogen","litter_saprotroph",     
         "ectomycorrhizal","wood_saprotroph","mycoparasite",          
         "dung_saprotroph","lichenized",            
         "lichen_parasite","sooty_mold","animal_parasite",       
         "unspecified_saprotroph", "pollen_saprotroph","epiphyte",              
         "arbuscular_mycorrhizal", "animal_endosymbiont",                      
         "root_endophyte","foliar_endophyte", "algal_parasite",        
         "nectar/tap_saprotroph","root-associated","rock-inhabiting",            
         "fungal_decomposer","animal_decomposer",         
         "bryophilous","root_endophyte_dark_septate",
         "unspecified_symbiotroph","ericoid_mycorrhizal",        
         "algal_symbiont","unsepcified_saprotroph",     
         "arthropod_parasite","resin_saprotroph")
)

names(unexpected.index) <- gsub("\\.", " ", names(unexpected.index))

un.table<- c(un.table1, un.table2, un.table3,un.table4)
for (i in 1:1000) {
  un.table[[i]]$interation <- i
  un.table[[i]]$primary_lifestyle[un.table[[i]]$primary_lifestyle == ""]<-"unknown"
  un.table[[i]]$Secondary_lifestyle[un.table[[i]]$Secondary_lifestyle == ""]<-"unknown"
}

un.table2 <- lapply(un.table, function(df) {
  df$type  <- NA
  df$type2 <- NA
  df$final <- NA
  idx1 <- mapply(function(sub, life) {
    if (!sub %in% names(unexpected.index)) return(FALSE)
    life %in% unexpected.index[[sub]]
  }, df$substrate, df$primary_lifestyle)
  
  idx2 <- mapply(function(sub, life) {
    if (!sub %in% names(unexpected.index)) return(FALSE)
    life %in% unexpected.index[[sub]]
  }, df$substrate, df$Secondary_lifestyle)
  
  df$type[idx1]  <- "unexpected"
  df$type2[idx2] <- "unexpected"
  df$type[!idx1]  <- "expected"
  df$type2[!idx2] <- "expected"
  df
})

result<-list()
for (i in 1:1000) {
  un.table2[[i]]$final[un.table2[[i]]$type == "unexpected" & un.table2[[i]]$type2 == "unexpected"]  <- "unexpected"
  un.table2[[i]]$final[un.table2[[i]]$type == "expected" | un.table2[[i]]$type2 == "expected"]  <- "expected"
  un.table2[[i]]$final[un.table2[[i]]$primary_lifestyle == "unknown" | un.table2[[i]]$primary_lifestyle == "" ]  <- "unknown"
  un.table2[[i]]$final[un.table2[[i]]$Secondary_lifestyle == "unknown" & un.table2[[i]]$type == "unexpected" ]  <- "unexpected"
  
  un.table2[[i]]$number <-1
  result[[i]]<- un.table2[[i]] %>%
    group_by(substrate, final) %>%
    summarise(number2 = sum(number)) %>%
    group_by(substrate) %>%
    mutate(ratio = number2/sum(number2))
}

resultall<-bind_rows(result)
resultall<-resultall %>%
  group_by(substrate, final) %>%
  summarise(ratio2=mean(ratio))
write.csv(resultall, "unexpected.proportion.csv")

