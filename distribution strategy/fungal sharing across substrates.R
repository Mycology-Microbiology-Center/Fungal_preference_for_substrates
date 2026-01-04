library(lme4)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(tidyr)
library(dplyr)
library(tibble)
library(purrr)
load("phi.separate.sites.RData")
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
########calculate the unstable bias in consistence
stable <- consistence %>%
  add_count(OTU, BestHost) %>%   
  filter(n == 3)
unstable <- consistence %>%
  add_count(OTU, BestHost) %>%   
  filter(n != 3)
unstable <-unstable[,-28]
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
Heatmap(similarity_matrix, name = "Ratio",
        cluster_rows = T,
        cluster_columns = TRUE, 
        col = colorRamp2(c(0,1,5, 10, 20), c("white","lightblue","#FFEEA0","#F4BB44","#C04000")),
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 9))
####
switch2$number<-1
switch2$BestHost<-gsub("LZ|LV|LW","",switch2$BestHost)
switch.p<-switch2 %>%
  group_by(BestHost,site) %>%
  summarise(number.s=sum(number))
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
write.csv(try,"switching.csv")
result<-try %>%
  group_by(BestHost) %>%
  summarise(c.m=mean(c.r),s.m=mean(s.r),c=mean(number.c),
            s=mean(number.s),all=mean(number.a),
            is=mean(is.r),sd=mean(sd.r),
            sd.m=mean(number.sd),is.m=mean(number.is))

