library(ggalluvial)
library(ggplot2)
library(dplyr)
library(ggforce)
setwd("C:/Users/meirong/Desktop/PhD project/second preference/final.table/organised code/1.phi")
r <- read.csv("all.r2.tax.csv", row.names = 1)
r2 <- r[r$p.value < 0.05, ]
library(dplyr)
library(tidyverse)
##
setwd("C:/Users/meirong/Desktop/PhD project/second preference/final.table/third/1.calculate the result/")
metadata<-read.csv("metadata.final2.csv",row.names = 1)
setwd("C:/Users/meirong/Desktop/PhD project/second preference/final.table/organised code/table")
fungi<-read.csv("table.nolichenandhost.csv",row.names = 1)
table<-as.data.frame(t(fungi))
table2<-table
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
indicator<-r2[,c(1,2)]
names(indicator)<-c("OTU","substrate")
indicator$substrate<-gsub("\\.", " ", indicator$substrate)
try<-merge(indicator,table_long,by=c("OTU","substrate"))
try2 <- aggregate(try$relative.sum,by=list(substrate=try$substrate,OTU=try$OTU),sum)
###
try3 <- table_long[table_long$OTU %in% indicator$OTU,] ##the indicators go into other substrate that are not indicators
indicator$group<-"indicator"
try<-merge(indicator,try3,by=c("OTU","substrate"),all=T)
try[is.na(try$group),]$group<-"nonindicator"
#try3 <- aggregate(try3$relative.sum,by=list(substrate=try3$substrate,OTU=try3$OTU),sum)
####
setwd("C:/Users/meirong/Desktop/PhD project/second preference/final.table/organised code/10.indicator exchange")
write.csv(try,"indicator.among. substrates.csv")
write.csv(try2,"indicator.in.substrates.csv")
####
library(dplyr)
library(tidyr)
library(igraph)
library(ggraph)
library(tidygraph)

### ----------------------------------------------------
### 1. Sample data (REPLACE with your df)
### ----------------------------------------------------
df <- try
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
# Use your flow_edges
df_sankey <- flow_edges %>%
  rename(source = from,
         target = to,
         value = weight)

df_sankey2 <- df_sankey  %>%
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
  scale_x_discrete(limits = c("Source", "Target"),
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
##inflow calculating
flow_edges2<-flow_edges[flow_edges$from!=flow_edges$to,]
flow_edges2<-flow_edges2 %>%
  group_by(to) %>%
  summarise(weight2=sum(weight))
###get the biggest sources
flow_edges2<-flow_edges[flow_edges$from!=flow_edges$to,]
flow_edges2<-flow_edges2 %>%
  group_by(from) %>%
  summarise(weight2=sum(weight))
flow_edges3<-flow_edges %>%
  group_by(from) %>%
  summarise(weight2=sum(weight))
all<-merge(flow_edges2,flow_edges3,by = "from")
all$ratio<-all$weight2.x/all$weight2.y
###
flow_edges2$ratio<-flow_edges2$weight2/sum(flow_edges2$weight2)
###get the sources
flow_edges2<-flow_edges %>%
  group_by(from) %>%
  summarise(weight2=sum(weight))
##outflow/inflow calculating
flow_edges2<-flow_edges[flow_edges$from!=flow_edges$to,]
outflow<-flow_edges2 %>%
  group_by(from) %>%
  summarise(outflow=sum(weight))
inflow<-flow_edges2 %>%
  group_by(to) %>%
  summarise(inflow=sum(weight))
all<-merge(inflow,outflow,by.x="to",by.y = "from")
all$ratio<-all$outflow/all$inflow
