rm(list=ls())

require(locater)
require(scattermore)
library(ggplot2)

experiment.name <- "WashU_CCDG"
run_name <- "2-wg-screening"

#data.dir <- paste0("/Users/xinxin.wang/Dropbox/finn/",experiment.name,"/screening_results/data/",run_name,"/")
locater.dir <- paste0("/gpfs/gibbs/pi/ycgh/xw445/projects/",experiment.name,
                      "/screening/to_share/locater_only/",run_name,"/")

smt.dir <- paste0("/gpfs/gibbs/pi/ycgh/xw445/projects/",experiment.name,
                  "/screening/to_share/smt_only/",run_name,"/")


chr.segments <- read.table("/gpfs/gibbs/pi/ycgh/xw445/projects/WashU_CCDG/chr_segments_4000segs_new_pass.txt",header=TRUE)


for (chr in 1:22){
  
  iters <- which(chr.segments$chr==chr)
  
  all.final.locater <- vector(mode = "list",length = length(iters))
  all.final.smt <- vector(mode = "list",length = length(iters))
  # create big list to gather data from all the iters in this chr
  
  for (iter in iters){
    
    locater.file <- paste0(locater.dir,"locater_res_",as.character(iter),".txt")
    smt.file <- paste0(smt.dir,"smt_res_",as.character(iter),".txt")
    
    if (!file.exists(locater.file) |!file.exists(smt.file) ){
      print(paste0("files for iter ",iter," doesn't exist."))
      next
    } else{ 
      locater.res <- data.table::fread(locater.file)
      smt.res <- data.table::fread(smt.file)
      # read in the locater res and smt result
      
      phenos <- locater.res$phenotype[locater.res$locus.idx==unique(locater.res$locus.idx)[1]]
      # get the phenotypes in the locater. 
      
      all.final.locater[[match(iter,iters)]]  <-  locater.res
      all.final.smt[[match(iter,iters)]]  <-  smt.res
      # assign the res to each element of big list

      # record if there's any NA
      rm(locater.res)
      rm(smt.res)
      gc()
    }
    
    
  }# this is for multiple iter loop
  
  # make the big list as a big table
  all.final.locater <- data.table::rbindlist(all.final.locater)
  all.final.smt <- data.table::rbindlist(all.final.smt)
  
  
  # Create a list with elements for each phenotype
  all.final.locater <- split(all.final.locater, all.final.locater$phenotype)
  all.final.smt <- split(all.final.smt, all.final.smt$phenotype)
  
  # Order the list elements based on the order of phenotypes in the vector
  all.final.locater <- all.final.locater[phenos]
  all.final.smt <- all.final.smt[phenos]
  
  
  whole.chr.data.dir <- paste0("/gpfs/gibbs/pi/ycgh/xw445/projects/",experiment.name,"/screening/to_share/locater_screening_whole_chr_data/",run_name,"/chr",chr,"/")
  if (!dir.exists(whole.chr.data.dir)){dir.create(whole.chr.data.dir,recursive = TRUE)}
  
  
  print(paste("chr",chr,"finished."))

  
  saveRDS(all.final.smt,paste0(whole.chr.data.dir,"whole_chr_smt_data_chr",chr,".rds"))
  saveRDS(all.final.locater,paste0(whole.chr.data.dir,"whole_chr_locater_data_chr",chr,".rds"))
  
  
  
}
# for real data: need to save the data, then get the region

##########################################################
##########################################################
##########################################################
rm(list=ls())


require(locater)
require(scattermore)
library(ggplot2)
library(dplyr)

colpal <- c("#000000", "#D55E00","#56B4E9","#F0E442","#009E73","#CC79A7","#E69F00","#0072B2","#AAAAAA")

FindInterestingRegion <- function(all_positions,buffer.region){
  # all_positions : vector with all the positions that exceed interested threshold for all the tests 
  # combined together
  # buffer.region: length (in bp ) of buffer region on both sides of the interested positions
  
  regions <- data.frame(start = all_positions - buffer.region, end = all_positions + buffer.region)
  
  if (nrow(regions)!=0) {
    regions <- regions %>% 
      arrange(start) %>%
      group_by(group = cumsum(start - lag(end, default = first(start)) > 0)) %>%
      summarise(start = min(start), end = max(end)) %>%
      select(-group)
    
    # 1. arrange the start position to make sure all the small regions are sorted
    # 2. group the lines so that overlapping regions are into the same group
    # 3. summarize the grouped lines, so that overlapping regons merge to biger regions
    # 4, only keep start and end column.
  }
  return(regions)
  
}

PvalAdjust <- function(x,intercept,slope){
  # x: -log10(p-vals)
  (x - intercept)/slope
}


experiment.name <- "WashU_CCDG"
run_name <- "2-wg-screening"

locater.threshold <- -log10(7.17e-9)-1
smt.threshold <- -log10(7.17e-9)


effective.n.tests.locater <- 0.05/10 ^ (-7.851456)
effective.n.tests.smt <- 0.05/10 ^ (-8.144726)

buffer.region <- 6e5

grand.data.dir <- "/gpfs/gibbs/pi/ycgh/xw445/projects/"
burnin.file <- paste0(grand.data.dir,experiment.name,"/phasing-vcfs.compute1.MANIFEST.burnin.table.txt")
burnin.table <- read.table(burnin.file,header=TRUE)

slope_intercept_table <- 
  read.table(paste0("/home/xw445/gibbs/",experiment.name,
                    "/screening/data/whole_genome_yale/locater_only/wg_unrelated_METSIMonlyData_101phenos_QQ_Bjitter/lambda_table_real_screened_phenos.txt"))

# replace this to the slopd and intercept information for traits that are actually screened in this experiment.

for (chr in 1:22){
  
  local.burnin.table <- burnin.table[burnin.table$chr==chr,]
  whole.chr.data.dir <- paste0("/gpfs/gibbs/pi/ycgh/xw445/projects/",experiment.name,"/screening/to_share/locater_screening_whole_chr_data/",run_name,"/chr",chr,"/")
  if (!dir.exists(whole.chr.data.dir)){
    next
  }

  all.final.smt <- readRDS(paste0(whole.chr.data.dir,"whole_chr_smt_data_chr",chr,".rds"))
  all.final.locater <- readRDS(paste0(whole.chr.data.dir,"whole_chr_locater_data_chr",chr,".rds"))

  phenos <- names(all.final.locater)
  
  local.data.dir <- paste0("/gpfs/gibbs/pi/ycgh/xw445/projects/",experiment.name,"/screening/to_share/locater_screening_local_data/",run_name,"/gc_controlled_alt/chr",chr,"/")
  viz.dir <- paste0("/gpfs/gibbs/pi/ycgh/xw445/projects/",experiment.name,"/screening/to_share/locater_screening_local_plots/",run_name,"/gc_controlled_alt/chr",chr,"/")
  
  if (!dir.exists(viz.dir)){dir.create(viz.dir,recursive = TRUE)}
  if (!dir.exists(local.data.dir)){dir.create(local.data.dir,recursive = TRUE)}
  
  for (pheno.idx in 1:length(phenos)){
    gc()
    
    pheno <- phenos[pheno.idx]
  
    slope.rd <- slope_intercept_table$rd_slope[pheno.idx]
    intercept.rd <- slope_intercept_table$rd_intercept[pheno.idx]
    
    slope.qform <- slope_intercept_table$qform_slope[pheno.idx]
    intercept.qform <- slope_intercept_table$qform_intercept[pheno.idx]
    
    final.smt <- all.final.smt[pheno.idx][[1]]
    
    final.locater <- all.final.locater[pheno.idx][[1]]
    
    ### sort the data by postion
    final.smt <- final.smt[order(final.smt$smt.pos),]
    
    final.locater <- final.locater[order(final.locater$locater.pos),]
    
    final.smt$in.locater <- 
      ifelse(final.smt$smt.pos %in% final.locater$locater.pos,1,0)
    
    
    # just control rd and Qform 
    final.locater$rd.controlled <- PvalAdjust(final.locater$rd,intercept = intercept.rd,slope = slope.rd)
    final.locater$rd.controlled[final.locater$rd.controlled < 0] <- 0
    final.locater$qform.controlled <- PvalAdjust(final.locater$qform,intercept = intercept.qform,slope = slope.qform)
    final.locater$qform.controlled[final.locater$qform.controlled < 0] <- 0
    # make sure -log10(p-values) are positive
    
    # then use the exact msse function in TestLoci
    final.locater$tot.controlled <- locater::msse.test(final.locater$smt,
                                                                final.locater$rd.controlled,
                                                                final.locater$qform.controlled,test.1.solo = TRUE)
    
  
    gc()
    
 
      final.locater$tot.controlled <- final.locater$tot.controlled + log10(effective.n.tests.smt/effective.n.tests.locater)
      

    
    filtered.smt <- subset(final.smt, smt.p > smt.threshold)
    filtered.locater <- subset(final.locater,tot.controlled > locater.threshold)

    all_positions <- unique(c(filtered.smt$smt.pos,
                              filtered.locater$locater.pos))
    
    
    regions <- FindInterestingRegion(all_positions,buffer.region)
    
    if (nrow(regions)==0){
      next
    }
    start.positions <- regions$start
    end.positions <- regions$end
    
    for (i in 1:length(start.positions)){
      smt.in.plot <-subset(final.smt,smt.pos >= start.positions[i] & smt.pos <= end.positions[i])
      locater.in.plot <- subset(final.locater,locater.pos >= start.positions[i] & locater.pos <= end.positions[i])


        saveRDS(list("smt.in.plot" = smt.in.plot,
                     "locater.in.plot" = locater.in.plot
                    
        ), paste0(local.data.dir,"local_data_",pheno,"_",start.positions[i],"-",end.positions[i],"_adjusted.rds"))
        cairo_pdf(filename=paste0(viz.dir,"Manhattan_",pheno,"_",start.positions[i],"-",end.positions[i],"_adjusted.pdf"),
                  width = 10,
                  height=5,
                  fallback_resolution=600)

      
      
      ############ start visualization ###########
      plot.max <- max(c(locater.in.plot$tot.controlled,
                        smt.in.plot$smt.p),na.rm=TRUE) +2
      
      sub.chr.segments <- subset(local.burnin.table,core.end >= start.positions[i] & core.start <= end.positions[i])
      if (nrow(sub.chr.segments)!=0){
        edges <- sub.chr.segments$burnin.start
        edges <- sort(edges)[-1]
        #small segment edges
      }
      
      smt.in.plot$in.locater <- factor(smt.in.plot$in.locater, levels = c("0", "1"))
      
      smt.in.locater <- subset(smt.in.plot,in.locater == 1)
      smt.not.in.locater <- subset(smt.in.plot,in.locater == 0)
      
    
      p <- ggplot(locater.in.plot, aes(x=locater.pos/1e6, y=tot.controlled,color="LOCATER")) +
        geom_point(size=0.6) + xlim(start/1e6, end/1e6) + ylim(0,plot.max) +
        ggtitle(paste0("Manhattan of ",pheno," in chr ",chr,":",start.positions[i],"-", end.positions[i])) + theme(plot.title = element_text(hjust = 0.5))+
        xlab("Position (Mb)") + ylab("-log10(p-value)") +
        geom_point(data=smt.in.locater,aes(x=smt.pos/1e6, y=smt.p,color="smt_in"),size=0.3) + 
        geom_point(data=smt.not.in.locater,aes(x=smt.pos/1e6, y=smt.p,color="smt_not"),size=0.3) +
        theme_classic()
      
      
   
        p <- p + geom_hline(yintercept=smt.threshold, linetype="dashed", color = "black")  +
          geom_vline(xintercept = edges,linetype="dashed", color = "grey")

      
      p <- p + scale_colour_manual(values = c("LOCATER" = "#D55E00","smt_in" = "#56B4E9", "smt_not" = "#000000"),  name = "", # Adjust the colors accordingly
                                   labels = c("LOCATER","SMT (Tested in\nLOCATER)", "SMT (Not tested\nin LOCATER)")) +
        theme(legend.key.size = unit(5, 'mm')) +
        guides(colour = guide_legend(
          override.aes = list(shape = 19, size = 3)  # Shape 22 is a filled square
        )) +
        labs(colour = NULL) + 
        theme(axis.text = element_text(size = 12),  # Increase axis number size
              legend.text = element_text(size = 12), axis.title.x = element_text(size = 14),  # Increase x-axis label size
              axis.title.y = element_text(size = 14))  # Increase legend text size
      
      

      print(p)
      dev.off()
      

    }
    rm(final.smt)
    rm(final.locater) 
    ############ end visualization ###########
    
    print(paste0("chr ",chr,", pheno ",pheno.idx," finished."))
    
  }
  print(paste0("chr ",chr," finished."))
  
  
}



