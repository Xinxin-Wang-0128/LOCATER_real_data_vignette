#################################
# this script is an example of how you could collect the tuning results.
##################################

rm(list=ls())

run_group <- "WashU_CCDG"
run_name <- "METSIM1"

grand.data.dir <- "/gpfs/gibbs/pi/ycgh/xw445/projects/"
res.storage.path <- paste0(grand.data.dir,run_group,"/tuning/to_share/",run_name,"/")
# you can change this path to where you store the .rds file output from tuning.

n.reps <- 1
n.segs <- 155
# in my study. we evaluated tuning objective in 155 files.

############## start of collecting results ################

all.res <- vector(mode="list",length=n.segs)
# initialize a list for storing results

for (i in 1:n.segs){
  file <- paste0(res.storage.path,"tuning_res_",i,".rds")
  if (file.exists(file)){
    res <- readRDS(file)
    result <- res$res
  
    HMM_pars <- res$HMM_pars
    smt.res <- res$smt.res
    error <- res$error.list
    warning <- res$warning.list
    
    result$normalized.signal <- NA
    result$smt.causal <- NA
    result$smt.target <- NA
    # initialize new columns
    
    for (j in 1:nrow(result)){
      if (length(smt.res[[j]]$V1[result$causal.idx[j]])!=0){
        # if the SMT result for causal variant is not empty
        result$smt.causal[j] <- smt.res[[j]]$V1[result$causal.idx[j]]
        result$smt.target[j] <- smt.res[[j]]$V1[result$target.idx[j]]
        # store the SMT result for causal and target variants
      }

      result$normalized.signal[j] <- result$signal[j] / result$smt.causal[j]
      # normalize by causal SMT
      
    }
   
    par.count <- nrow(result) / n.reps
    # calculate how many parameters we tested
    # this will be useful when n.reps > 1
    
    # Use split to split the vector into different reps (each rep is a group, and every group contain all parameters)
    grouped_vector <- split(result$smt.causal, rep(1:n.reps, each = par.count))
    
    # Calculate the maximum value for each group using sapply
    max.of.causal.smt <- sapply(grouped_vector, max)
    
    significant.smt <- max.of.causal.smt >= 20
    # the cutoff for significant SMT -log10(P) is 20

    significant.smt <- rep(significant.smt,each=par.count)
    
    result <- result[which(significant.smt),]
    # this is to delete iterations when SMT is not significant.
    
    all.res[[i]] <- result
    
  }
}

all.res <- data.table::rbindlist(all.res)
# combine all results into a single data.table

write.table(all.res,file=paste0(res.storage.path,"all_res_remove_nosignal.txt"),
            quote=FALSE,
            col.names=TRUE,sep = "\t")

############## end of collecting results ################


