rm(list=ls())

#library(kalis)
#library(ComplexHeatmap)
library(fastcluster)
library(dendextend)
library(Matrix)
library(data.table)


run_group <- "WashU_CCDG"
run_name <- "whole_genome_yale/wg_unrelated_METSIMonlyData_101phenos_Bjitter_local_10cM_clean"

grand.data.dir <- "/gpfs/gibbs/pi/ycgh/xw445/projects/"


res.storage.path <- paste0(grand.data.dir,run_group,"/screening/data/",run_name,"/")

distinguishable_colors <- c("#000000", "#D55E00","#56B4E9","#F0E442","#009E73",
                            "#CC79A7","#E69F00","#0072B2","#000000", "#D55E00","#56B4E9","#F0E442","#009E73",
                            "#CC79A7","#E69F00","#0072B2")



labels <- c(rep("*",8),rep("▲",8))

keep.percentage <- 5
# we would like to only keep a certain percentage of unimportant sprigs for viz purposes

significant.sprig.info.file <- fread(paste0(res.storage.path,"significant_sprig_info_additional.txt"))

viz.dir <- paste0(grand.data.dir,run_group,"/screening/data/whole_genome_yale/locater_screening_local_plots/",run_name,"/misc/")
if(!dir.exists(viz.dir)){ dir.create(viz.dir)}

for (row in 1:nrow(significant.sprig.info.file)){

  chr <- as.integer(significant.sprig.info.file[row,]$chr)
  local.index <- as.integer(significant.sprig.info.file[row,]$local.index)

  
  interested.pheno <- significant.sprig.info.file[row,]$interested.pheno
  significant.sprigs <- significant.sprig.info.file[row,]$significant.sprigs
  
  if (grepl(",", significant.sprigs)) {
    significant.sprigs <- as.integer(unlist(strsplit(significant.sprigs, ",")))
  } else {
    significant.sprigs <- as.integer(significant.sprigs)
  }
  # this is the index of significant sprigs 

  y <- readRDS(paste0(res.storage.path,"y_",local.index,"_chr",chr,".rds"))
  index <- which(colnames(y) == interested.pheno)
  
  sprigs <- readRDS(paste0(res.storage.path,"sprigs_",local.index,"_chr",chr,".rds"))
 
  interesting.haps <- vector(mode="list",length=length(significant.sprigs))
  for (i in 1:length(significant.sprigs)){
    interesting.haps[[i]] <- which(sprigs$assignments == significant.sprigs[i])
  }
  
  color_list <- lapply(seq_along(interesting.haps), function(i) {
    rep(distinguishable_colors[i], length(interesting.haps[[i]]))
  })
  
  label_list <- lapply(seq_along(interesting.haps), function(i) {
    rep(labels[i], length(interesting.haps[[i]]))
  })
  colors <- unlist(color_list)
  unlisted.lables <- unlist(label_list)
  interesting.haps <- unlist(interesting.haps)
  
  
  hapM <- readRDS(paste0(res.storage.path,"hapM_",local.index,"_chr",chr,".rds"))
  
  
  distance_matrix <- hapM + t(hapM)
  # # 
  # # # Perform hierarchical clustering
  method <- "average"
  hc <- fastcluster::hclust(as.dist(distance_matrix), method = method)
  # # 
  inds.in.hc <- match(interesting.haps,hc$order)
  # 
  dend <- as.dendrogram(hc)
  # 
  # # Plot the dendrogram
  labels.to.prune <- setdiff(labels(dend), interesting.haps)
  # # 
  no.other.individuals <- floor((nrow(distance_matrix)-length(interesting.haps)) * keep.percentage /100)
  # # 
  labels.to.prune <- labels.to.prune[sample(size=(nrow(distance_matrix)-length(interesting.haps)-no.other.individuals),
                                            x=1:length(labels.to.prune),replace = FALSE)]
  # #
  pruned.dend <- dendextend::prune(dend, labels.to.prune)
  # this will take a long time
  
  saveRDS(pruned.dend, paste0(viz.dir,"pruned_dend_keepextra",
                              no.other.individuals,"_",method,"_chr",chr,"_",local.index,".rds"))
 # save the object for future use 
  
  
  pruned.dend <- readRDS(paste0(viz.dir,"pruned_dend_keepextra",
                                no.other.individuals,"_",method,"_chr",chr,"_",local.index,".rds"))

  
  current.order <- match(labels(pruned.dend),interesting.haps)
  current.color <- rep("black",length(current.order))
  
  current.color[which(!is.na(current.order))] <- colors[current.order[which(!is.na(current.order))]]
  
  # Get existing labels
  existing_labels <- labels(pruned.dend)
  
  # Create a new labels vector
  new_labels <- ifelse(current.color == "black", "", existing_labels)
  new_labels[which(!is.na(current.order))] <- unlisted.lables[current.order[which(!is.na(current.order))]]
  # Apply the new labels to the dendrogram
  pruned.dend <- pruned.dend %>% dendextend::set("labels", value = new_labels)
  
  #par(mar = c(5, 4, 4, 2) + 0.1)
  
  cairo_pdf(paste0(viz.dir,"dendrogram_pruned_keepextra",no.other.individuals,"_",
                   method,"_chr",chr,"_",local.index,"_replot.pdf"), width = 16.5, height = 12)
  
  
  colored_dend <- pruned.dend %>% 
    dendextend::set("labels_cex", 0.3) %>% 
    dendextend::set("by_labels_branches_col", value = unlisted.lables[current.order[which(!is.na(current.order))]], type = "any",TF_values = c("#D55E00","#AAAAAA")) %>% 
    dendextend::set("by_labels_branches_lwd", value = unlisted.lables[current.order[which(!is.na(current.order))]],type = "any", TF_values = c(1,0.6)) %>% 
    dendextend::set("labels_col", current.color)
  
  plot(colored_dend,cex.axis = 1)
  
  
  dev.off()
}


