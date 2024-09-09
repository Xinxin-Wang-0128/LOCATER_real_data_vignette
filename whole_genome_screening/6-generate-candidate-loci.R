# this is to generate regions for more precise estimate of the null
# usually k=5. This will take a longer time so we did not use it in whole genome screening

rm(list=ls())

require(locater)
require(dplyr)
library(tidyr)
library(zoo)


experiment.name <- "WashU_CCDG"
run_name <- "2-wg-screening"
grand.data.dir <- "/gpfs/gibbs/pi/ycgh/xw445/projects/"

# read in all hits table
all.hits <- read.table(paste0("/gpfs/gibbs/pi/ycgh/xw445/projects/",experiment.name,"/screening/to_share/locater_screening_local_data/",run_name,"/all_hits_new_thresh.txt"),
                       header=TRUE)



burnin.file <- paste0(grand.data.dir,experiment.name,"/phasing-vcfs.compute1.MANIFEST.burnin.table.txt")
burnin.table <- read.table(burnin.file,header=TRUE)
edges <- burnin.table[,match(c("chr","core.start"),names(burnin.table))]
edges  <- edges[edges$core.start!=0,]

# Create a function to split rows based on edges
split_row <- function(row, edges) {
  chr <- row$chr
  start <- row$start
  end <- row$end
  
  # Find edges that are between start and end
  valid_edges <- edges[which(edges > start & edges < end)]
  
  # Create new data frame based on valid edges
  new_df <- data.frame(chr = integer(0), start = integer(0), end = integer(0))
  for(edge in valid_edges) {
    new_df <- rbind(new_df, data.frame(chr = chr, start = start, end = edge))
    new_df <- rbind(new_df, data.frame(chr = chr, start = edge, end = end))
  }
  
  # If no valid edges, return the original row
  if(nrow(new_df) == 0) {
    return(row)
  } else {
    return(new_df)
  }
}

merge_intervals <- function(df) {
  df %>%
    arrange(start) %>% 
    # Sorts the data by the start column.
    mutate(group = cumsum(c(lag(end, default = first(end)) < start))) %>%
    #  checks if the current start is greater than the previous end (i.e., no overlap).
    group_by(group) %>%
    summarize(
      start = first(start),
      end = max(end)
    ) %>%
    select(-group)
}


candidate.regions <- vector(mode="list",length=22)

for (chr in 1:22){
  local.hits <- all.hits[all.hits$chr==chr,]
  
  
  if (nrow(local.hits)!=0) {
    regions <- data.frame(start=local.hits$start,end=local.hits$end)
    
    regions <- merge_intervals(regions)
    
    regions <- data.frame(chr=rep(chr,nrow(regions)),regions)
    
    local.edges <- edges$core.start[edges$chr==chr]
    
    # Split each row based on edges
    result <- do.call(rbind, lapply(1:nrow(regions), function(i) split_row(regions[i,], local.edges)))
    
    
    candidate.regions[[chr]] <- result
    
  }
  
  
}
candidate.regions <- data.table::rbindlist(candidate.regions)


write.table(candidate.regions,file =  paste0(grand.data.dir,experiment.name,"/screening/to_share/locater_screening_local_data/",run_name,"/candidate_loci.txt"),
            quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)

