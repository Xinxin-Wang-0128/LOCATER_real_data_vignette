##################
# convert rds file to txt
##################

rm(list=ls())

require(data.table)

experiment.name <- "WashU_CCDG"
test.name <- "2-pval-rank-matched"
grand.data.dir <- "/gpfs/gibbs/pi/ycgh/xw445/projects/"


data.dir <- paste0(grand.data.dir,experiment.name,"/screening/to_share/",test.name,"/")
# specify input directories


chr.segments <- read.table("/gpfs/gibbs/pi/ycgh/xw445/projects/WashU_CCDG/chr_segments_4000segs_new_pass.txt",header=TRUE)
# specify the information for small segments

iters <- 1:4000
options(scipen = 999)
# suppress scientific notation

locater.dir <- paste0("/gpfs/gibbs/pi/ycgh/xw445/projects/",experiment.name,"/screening/to_share/locater_only/",test.name,"/")
if(!dir.exists(locater.dir)){dir.create(locater.dir,recursive = TRUE)}
# specify output directories

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
                         quote=FALSE,sep="\t")
      
    }
    
  }
  


