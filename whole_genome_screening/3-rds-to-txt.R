##################
# convert rds file to txt
##################

rm(list=ls())

require(data.table)

experiment.name <- "WashU_CCDG"
run_name <- "2-wg-screening"
grand.data.dir <- "/gpfs/gibbs/pi/ycgh/xw445/projects/"


data.dir <- paste0(grand.data.dir,experiment.name,"/screening/to_share/",run_name,"/")
# specify input directories


chr.segments <- read.table("/gpfs/gibbs/pi/ycgh/xw445/projects/WashU_CCDG/chr_segments_4000segs_new_pass.txt",header=TRUE)
# specify the information for small segments

iters <- 1:4000
options(scipen = 999)
# suppress scientific notation

locater.dir <- paste0("/gpfs/gibbs/pi/ycgh/xw445/projects/",experiment.name,"/screening/to_share/locater_only/",run_name,"/")
if(!dir.exists(locater.dir)){dir.create(locater.dir,recursive = TRUE)}
# specify output directories

smt.dir <- paste0("/gpfs/gibbs/pi/ycgh/xw445/projects/",experiment.name,"/screening/to_share/smt_only/",run_name,"/")
if(!dir.exists(smt.dir)){dir.create(smt.dir,recursive = TRUE)}

for (iter in iters){
  file <- paste0(data.dir,"res_",as.character(iter),".rds")
  if (!file.exists(file)){
    print(paste0("file for iter ",iter," doesn't exist."))
    next
  } else{ 
    res <- readRDS(file) 
    region <- res$sim
    
    locater.res <- res$out.list[[1]]
    
    locater.pos <- region$pos[locater.res$locus.idx]
    # add in position information for LOCATER output
    
    locater.res <- cbind(locater.pos,locater.res)
    
    data.table::fwrite(locater.res,file=paste0(locater.dir,"locater_res_",as.character(iter),".txt"),
                       quote=FALSE,sep="\t",na = "NA")
    
    
    ####### SMT ####### 
    
    phenos <- locater.res$phenotype[locater.res$locus.idx==unique(locater.res$locus.idx)[1]]
    # query phenotype information
    
    all.smt <- vector(mode = "list",length = length(phenos))
    # create a list to store SMT results
    
    for (i in 1:length(res$smt.res)){
      smt.p <- res$smt.res[[i]]$smt.results[core.idx]
      # extrat SMT results for variant in the core region
      
      smt.regional <- data.frame(smt.pos,map,MAC,smt.p)
      # combine p-value with map, minor allele count (SMT) and position
      
      smt.regional$phenotype <- rep(phenos[i],nrow(smt.regional))
      # add in the phenotype information
      
      all.smt[[i]] <- smt.regional
    }
    
    all.smt <- data.table::rbindlist(all.smt)
    
    #
    data.table::fwrite(all.smt,file=paste0(smt.dir,"smt_res_",chr,"_",start,"-",end,".txt"),
                       quote=FALSE,sep="\t",na = "NA")
    
    
  }
  
}



