library(pairwiseAdonis)
library(ape)
library(metagMisc)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
metadata<-read.csv("metadata.final2.csv",row.names = 1)
fungi<-read.csv("table.nolichenandhost.csv",row.names = 1)
fungi<-as.data.frame(t(fungi))
fungi$sample_names<-row.names(fungi)
fungi2<-merge(fungi,metadata,by="sample_names",all=T)
row.names(fungi2)<-fungi2$sample_names
fungi2<-fungi2[,-1]
fungi2<-fungi2[!is.na(fungi2$e57bceb79df4a766f8fb9d049798c662ba89dc9e),]
fungi<-fungi2[,-c((ncol(fungi2)-1):ncol(fungi2))]
fungi<-fungi[,colSums(fungi)>0]
fungi<-fungi[rowSums(fungi)>0,]
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
fungi.d<- simpson_abund(fungi)
row.names(metadata)<-metadata$sample_names
metadata3<-metadata[attr(fungi.d,"Labels"),]
####
permanova_result <- adonis2(fungi.d ~ substrate,strata = metadata3$site2, data = metadata3,permutations = 999)
pairwise_results <- pairwise.adonis(fungi.d,factors = metadata3$substrate, perm = 999, p.adjust.m = "fdr")
substrate_shapes <- c("feces" = 15,"leaves" = 7,  "snow" = 8,  "fruit body" = 18,
                      "lichen" = 16,  "topsoil"=1,  "subsoil"=13,  "leaf litter"=22,  "wood"=10,
                      "bark"=3,  "moss"=4)
distance<-fungi.d
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
fungi.d2 <- dist2list(fungi.d)
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
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(grid)  # for grid.text
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

