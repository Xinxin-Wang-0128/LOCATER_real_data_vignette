rm(list=ls())

require(locater)
require(scattermore)
library(ggplot2)



getMaxRegion3 <- function(file_path,locater.threshold,smt.threshold,locater.column){
  loaded_object <- readRDS(file_path)
  smt.in.plot <- loaded_object$smt.in.plot
  locater.in.plot <- loaded_object$locater.in.plot
  rm(loaded_object)

  max.locater <- max(c(max(locater.in.plot[[locater.column]])))
  max.smt <- max(c(max(smt.in.plot$smt.p)))
  max.locater.pos <- locater.in.plot$locater.pos[which.max(locater.in.plot[[locater.column]])]
  max.smt.pos <- smt.in.plot$smt.pos[which.max(smt.in.plot$smt.p)]
  
  
  if(max.locater >= locater.threshold |  max.smt >= smt.threshold){
    
    chr_match <- regmatches(file_path, regexpr("chr(\\d+)", file_path))
    chr <- gsub("chr", "", chr_match)
    
    phenotype <- sub(".*/local_data_(.*)_\\d+-\\d+(_nojitter)?\\_adjusted.rds", "\\1", file_path)
    
    # Extract start and end using regular expression
    start_match <- regmatches(file_path, regexpr("\\d+-", file_path))
    end_match <- regmatches(file_path, regexpr("-\\d+", file_path))
    start <- as.numeric(gsub("-", "", start_match))
    end<- as.numeric(gsub("-", "", end_match))
    
    max.locater.pos <- mean(max.locater.pos)
    max.smt.pos <- mean(max.smt.pos)
    # when there are mutiple max position, use the mean
    
    return(c(chr,phenotype,start,end, max.locater,max.smt,max.locater.pos,max.smt.pos))
    
  }
  
}

experiment.name <- "WashU_CCDG"
run_name <- "2-wg-screening"

locater.threshold <- -log10(7.17e-9)-1
smt.threshold <- -log10(7.17e-9)

all.maxes <- vector(mode = "list",length = 22)


for (chr in 1:22){
  
  local.data.dir <- paste0("/gpfs/gibbs/pi/ycgh/xw445/projects/",experiment.name,"/screening/to_share/locater_screening_local_data/",run_name,"/gc_controlled_alt/chr",chr,"/")
  if (!dir.exists(local.data.dir)){next}
  rds_files <- list.files(local.data.dir, pattern = "\\_adjusted.rds$", full.names = TRUE)
  # Use lapply to apply the condition function to each file in the list
  maxes <- lapply(rds_files, getMaxRegion3,locater.threshold=locater.threshold,
                  smt.threshold=smt.threshold,locater.column = "tot.controlled")
  maxes <- maxes[!sapply(maxes, is.null)]

  if (length(maxes)==0){next}
  maxes <- do.call(rbind, maxes)
  colnames(maxes) <- c("chr","phenotype","start","end","locater_tot",
                       "smt","max.locater.pos","max.smt.pos")
  all.maxes[[chr]] <- maxes
}
all.maxes <- all.maxes[!sapply(all.maxes, is.null)]
all.maxes <- do.call(rbind,all.maxes)
all.maxes <- as.data.frame(all.maxes) 


all.maxes$locater_tot <- as.numeric(all.maxes$locater_tot)

all.maxes$smt <- as.numeric(all.maxes$smt)

all.maxes$sig.locater <- ifelse(all.maxes$locater_tot > all.maxes$smt, TRUE, FALSE)

write.table(all.maxes,file =  paste0("/gpfs/gibbs/pi/ycgh/xw445/projects/",experiment.name,"/screening/to_share/locater_screening_local_data/",run_name,"/all_hits_new_thresh.txt"),
            quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE
)
