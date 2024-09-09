# this script is to plot manhattan plot for qform, rd and smt seperately
# input: interested location list;
# local data
# output: split manhattan

rm(list=ls())

require(locater)
require(scattermore)
library(ggplot2)
require(dplyr)
require(data.table)

experiment.name <- "WashU_CCDG"
test.name <- "putative-investigation"

all.hits <- 
  fread(paste0("/home/xw445/gibbs/WashU_CCDG/screening/to_share/locater_screening_local_data/",
               test.name,"/all_hits_new_thresh.txt"))


colpal <- c("#000000", "#D55E00","#56B4E9","#F0E442","#009E73","#CC79A7","#E69F00","#0072B2","#AAAAAA")


for ( i in 1:nrow(all.hits)){
  pheno <- all.hits$phenotype[i]
  start <- all.hits$start[i]
  end <- all.hits$end[i]
  chr <- all.hits$chr[i]
  
  local.data.dir <- paste0("/gpfs/gibbs/pi/ycgh/xw445/projects/",experiment.name,"/screening/to_share/locater_screening_local_data/",test.name,"/gc_controlled_alt/chr",chr,"/")
  # this is the directory that contains the local data for a region and a specific phenotype
  
  file.name <- paste0(local.data.dir,"local_data_",pheno,"_",start,"-",end,"_adjusted.rds")
  # this is the file for the local data
  
  local.res <- readRDS(file.name)
  # the rds file should contain an object that has the following components:
  # smt.in.plot: a table with columns 
  # locater.in.plot: a table with columns 
  
  smt.in.plot <- local.res$smt.in.plot
  locater.in.plot <- local.res$locater.in.plot
  
  rm(local.res)
  
  
  viz.dir <- paste0("/gpfs/gibbs/pi/ycgh/xw445/projects/",experiment.name,"/screening/to_share/locater_screening_local_plots/",test.name,"/gc_controlled_alt/chr",chr,"/")
  if (!dir.exists(viz.dir)) {
    dir.create(viz.dir)
  }
  ####################################################
  ##### section for Manhattan plot visualizing sub-tests of LOCATER
  ####################################################
  cairo_pdf(filename=paste0(viz.dir, "split_Manhattan_", pheno, "_", start, "-", end, "_locaters.pdf"),
            width = 10,
            height = 5,
            fallback_resolution = 600)
  
  plot.max <- max(c(locater.in.plot$tot.controlled,
                    smt.in.plot$unjitter$smt.p), na.rm = TRUE) + 2
  
  # Base ggplot with SMT points
  p <- ggplot(locater.in.plot, aes(x = locater.pos/1e6, y = smt, colour = "SMT")) +
    geom_point(size=0.6) +
    xlim(start/1e6, end/1e6) + ylim(0, plot.max) +
    ggtitle(paste0("Manhattan of ",pheno," in chr ",chr,":",start,"-", end)) + 
    theme(plot.title = element_text(hjust = 0.5))+
    xlab("Position (Mb)") + ylab("-log10(p-value)") + theme_classic()
  
  # Add RD and QForm points with respective colors
  p <- p + geom_point(data = locater.in.plot, aes(x = locater.pos/1e6, y = rd.controlled, colour = "SD"), size=0.6) +
    geom_point(data = locater.in.plot, aes(x = locater.pos/1e6, y = qform.controlled, colour = "QForm"), size=0.6)
  
  # Customize colors and add horizontal line
  p <- p + scale_color_manual(values = c("SMT" = "#56B4E9", "SD" = "#F0E442", "QForm" = "#009E73")) +
    # geom_hline(yintercept=locater.threshold, linetype="dashed", color = "black")   +
    theme(legend.key.size = unit(5, 'mm')) + 
    guides(colour = guide_legend(
      override.aes = list(shape = 19, size = 3)  # Shape 22 is a filled square
    )) + 
    labs(colour = NULL) +
    # Customizing theme
    theme(axis.text = element_text(size = 12),  # Increase axis number size
          legend.text = element_text(size = 12), axis.title.x = element_text(size = 14),  # Increase x-axis label size
          axis.title.y = element_text(size = 14))  # Increase legend text size
  
  # Print plot
  print(p)
  dev.off()
  
  ####################################################
  ####################################################
  
  
  
  ####################################################
  ##### section for bar plot visualizing sub-tests of LOCATER at lead marker
  ####################################################
  locater.peak <- locater.in.plot[which.max(locater.in.plot$tot.controlled),]
  
  effective.n.tests.locater <- 0.05/10 ^ (-7.851456)
  effective.n.tests.smt <- 0.05/10 ^ (-8.144726)
  
  
  secondary.adjustment.diff <- log10(effective.n.tests.smt/effective.n.tests.locater)
  # Extract values for plotting
  

  values <- c(locater.peak$smt, locater.peak$rd, locater.peak$rd.controlled, locater.peak$qform, 
              locater.peak$qform.controlled,locater.peak$tot,locater.peak$tot.controlled-secondary.adjustment.diff,locater.peak$tot.controlled)
  names(values) <- c("SMT", "SD", "SD\nAdjusted", "QForm", "QForm\nAdjusted","Combined","Combined\nAdjusted","Combined\nThreshold\nStandardized")
  
  # Create a vector of colors
  colors <- c("black", "grey","black","grey","black","grey","black","black")
  
  cairo_pdf(paste0(viz.dir,"bar_peak_",pheno, "_", start, "-", end,,".pdf"),width=10,height=7)
  par(mar=c(7,5,5,5))
  
  # Create the bar plot
  midpoints <- barplot(values, col = colors, main = paste("Phenotype:", locater.peak$phenotype,
                                                          " LOCATER Pos:", locater.peak$locater.pos, 
                                                          # "\nProp.Var:", round(locater.peak$prop.var,5),
                                                          "\nCombined controlled",round(locater.peak$tot.controlled-secondary.adjustment.diff,5)),
                       ylab = "-log10(p-values)", xlab = "",xaxt = "n",las=1)
  text(x = midpoints+0.25, y = -0.5,  # Adjust y for label placement below bars
       labels = names(values), srt = 45, adj = 1,
       xpd = TRUE, cex = 1)  # 'xpd = TRUE' allows text outside plot margins
  
  dev.off()
  
  ####################################################
  ####################################################
  
}
