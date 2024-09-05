#### this script is a example to do tuning for METSIM
### we set the distance between causal and target loci to be 0.05 cM
#### user could change that by changing the desired.cM.dist parameter in the select.tuning.loci function

rm(list=ls())
set.seed(28)
#######################
# customized function to select tuning loci
#######################
select.tuning.loci <- function (idx, map, core.idx, Q, desired.cM.dist = 0.01, tol = 0.1){   
  
  # idx: candidate index for target variants
  # map: genetic map
  # core.idx: index for core region, if we want to set boundary for causal and target loci
  # Q: QR decomposition of the background covariate matrix
  # desired.cM.dist: desired distance between causal and target loci
  # tol: tolerance for the desired distance
  
  default.return <- list(causal.idx = NA_integer_, target.idx = NA_integer_)
  if (!length(idx)) {
    warning("candidate.target.idx passed to select.tuning.loci is empty, retuning NAs for causal and target idx")
    return(default.return)
  }
  min.dist.from.causal.var <- desired.cM.dist * (1 - tol)
  max.dist.from.causal.var <- desired.cM.dist * (1 + tol)
  perm.candidate.target.idx <- sample(idx)
  for (i in 1:length(idx)) {
    target.idx <- perm.candidate.target.idx[i]
    # choose the target index for this iteration
    candidate.causal.idx <- which((map[target.idx] + min.dist.from.causal.var <= map & map[target.idx] + max.dist.from.causal.var >= map) |
                                    (map[target.idx] - min.dist.from.causal.var >=map & map[target.idx] - max.dist.from.causal.var <=map))
    # choose candidate causal index based on the distance between causal and target loci
    if (length(candidate.causal.idx)) {
      # if there exists some options for candidate causal loci:
      
      ############### start of ensuring signal is not co-linear with background covariates ############
      
      g <- locater:::Haps2Genotypes(QueryCache(loci.idx = candidate.causal.idx),ploidy = ploidy)
      # get the genotype vector of candidate causal loci
      DACs <- rowSums(g)
      g_2 <- g^2
      g_2 <- apply(g_2,1,sum)
      g_2 <- sqrt(g_2)
      for (i in 1:nrow(g)){
        g[i,]<- g[i,] / g_2[i]
      }
      # normalize the genotype matrix
      ### then calculate sum(crossprod(g,Q)^2)
      crosprod <- apply(g,1,crossprod,Q)
      crosprod <- crosprod ^ 2
      sum_square <- apply(crosprod,2,sum)
      candidate.causal.idx <- candidate.causal.idx[sum_square<=0.02]
      # make sure the causal variant is not co-linear with background covariates, ensuring enough signal
      
      ############### end of ensuring signal is not co-linear with background covariates ############
      
      if (length(candidate.causal.idx)){
        break
        # if we find some loci that satisfy the constraints, break the loop
      }
    }
  }
  if (!length(candidate.causal.idx)) {
    # if we cannot find any loci that satisfy the constraints at the end, return NA
    warning("no causal loci could be found satisfying the desired.cM.dist constraint, adjust tol? retuning NAs for causal and target idx")
    return(default.return)
  }
  list(target.idx = target.idx,  causal.idx = sample(candidate.causal.idx,
                                                    1),sum_square=sum_square,DACs = DACs)
}


#######################
# pipe in parameter
#######################

nthreads <- as.integer(commandArgs(TRUE)[1])
job.index <- commandArgs(TRUE)[2]


##########################
# Specify Directories    #
##########################

run_group <- "WashU_CCDG"
run_name <- "METSIM1"

grand.data.dir <- "/gpfs/gibbs/pi/ycgh/xw445/projects/"

res.storage.path <- paste0(grand.data.dir,run_group,"/tuning/to_share/",run_name,"/")
if(!dir.exists(res.storage.path)){dir.create(res.storage.path,recursive = TRUE)}

input.genome.dir <- paste0(grand.data.dir,run_group,"/haplotypes/segmented_phased_merged_reannotated_ancestral_split_new_pass_singleton_filtered_6806metsim_hap/")

pc.dir <- paste0(grand.data.dir,run_group,"/pc_pheno/pc_finn_only/")
pc.prefix <-"finnish_SD8_iteration4_to_keep_table"

#####################
# Load Dependencies #
#####################

#### Load Libraries ####
.libPaths(c("/rchrist/lib/R.4",.libPaths())) # CRITICAL TO PUT IN THIS ORDER because it will prioritize the custom library
library(dplyr)
library(RSpectra)
library(kalis)
library(locater)

library(parallel)

require(glmnet)
library(tidyr)
library(data.table)
library(Matrix)


####################################################################
# Specify # run reps / job  and whether we're using forking #
####################################################################
ploidy <- 2L

use.forking <- FALSE
use.speidel <- TRUE
# if we plan to use Speidel et al. version of local ancestry inference
thresh <- 0.2
# threshold for matrix regularization
max1var <- TRUE
# the maximum number of variants on each branch is set to 1

burnin.var.count <- 5000
# leave N variatns at each end of the sampled chromosome, make sure local ancestry inference is not affected by the boundary

#################################
# generate sim object, an object storing information for tuning
#################################

burnin.file <- paste0(grand.data.dir,run_group,"/phasing-vcfs.compute1.MANIFEST.burnin.table.txt")
burnin.table <- read.table(burnin.file,header=TRUE)
file_prefix <- paste(burnin.table[job.index,1:3],collapse = "-")


start0 <- proc.time()[3]

propagation.window.var.count <- burnin.var.count*2+15000
# the number of variants for the sampled window. Here we set it to be 2*burnin.var.count + 15000. 
# 15000 is the length for core region, where we expect local ancestry inference to be stable and LOCATER p-value to be reliable

sim <- list()

# create haps
haps.path <- paste0(input.genome.dir,file_prefix,".phased_multi_merged_unique_ID_ancestral_split_sample_AC_var_metsim6806_nosingleton.hap.gz")
sim$haps <- haps.path

# create propagation window idx and map
map.dir <- paste0(grand.data.dir,run_group,"/recomb_map_new_pass_6806metsim/")
map.file <- read.table(paste0(map.dir,file_prefix,"_AF_weighted_cM.txt"),header=TRUE)
map <- map.file$cM

pos <- map.file$POS
pass <- map.file$REF != "N" & nchar(map.file$REF) + nchar(map.file$ALT) <= 2
# only consider SNPs with ancestral allele annotation

if (sum(pass)<propagation.window.var.count){stop("not enough variants in this segment for tuning.")}

propagation.window.idx <-(floor(length(map)/2)-floor(propagation.window.var.count/2)+1):(floor(length(map)/2)+floor(propagation.window.var.count/2))
# we sample the window from the middle of this segment

propagation.window.idx <- intersect(propagation.window.idx,which(pass==TRUE))

sim$propagation.window.idx <- propagation.window.idx
sim$map <- map
sim$pos <- pos

print("sim object created!")
#################################################################
# Specify list of parameters and read in A  #
#################################################################
# subset haplotypes, so that the tuning results are applicable to METSIM
# since we have decided on samples to run, I load the names of these samples and subset genome and A matrix based on that

start1 <- proc.time()[3]
samples.in.experiment <- readRDS("/gpfs/gibbs/pi/ycgh/xw445/projects/WashU_CCDG/screening/data/whole_genome_yale/wg_unrelated_METSIMonlyData_101phenos_local_10cM/samples_in_experiment.RDS")
sample.file <- paste0(input.genome.dir,file_prefix,".phased_multi_merged_unique_ID_ancestral_split_sample_AC_var_metsim6806_nosingleton.samples")
sample.table <- read.table(sample.file,header=TRUE)

sample.idx <- match(samples.in.experiment,sample.table$sample)
sample.idx <- na.omit(sample.idx)
sample.idx <- as.vector(sample.idx)


if (all.equal(sort(sample.idx),sample.idx)){
  hap.idx <- sort(c(2*sample.idx,2*sample.idx-1))
} else {
  stop("need to sort phenotype file to enable proper hap subsetting.")
}
hap.idx <- as.integer(hap.idx)

# read in background matrix 
pc.file <- paste0(pc.dir,pc.prefix,".tsv")
A <- read.table(pc.file,stringsAsFactors = FALSE,header = TRUE)

# sort A based on sample with phenos
sample.in.pc <- match(samples.in.experiment,A$ID)
sample.in.pc <- na.omit(sample.in.pc)
A <- A[sample.in.pc,]

A <- cbind(rep(1,nrow(A)),A[,2:11])
# as required by LOCATER, add a column of "1" before the background covariate matrices
A <- as.matrix(A)
Q <- qr.Q(qr(A))

# Load haplotypes and recombination map
CacheHaplotypes(sim$haps,
                loci.idx = sim$propagation.window.idx,
                hap.idx = hap.idx,
                warn.singletons = TRUE)

# since we subsetted the haplotypes. 
# need to see if we have singletons or monomorphics.
MAC <- rep(NA,L())
for (i in 1:L()){
  ac <- sum(QueryCache(loci.idx = i)) 
  MAC[i] <- min(ac,N()-ac)
}

if (sum(MAC==1 | MAC==0)>0){
  # if subsetting individuals created more singletons and monomorphics,
  ClearHaplotypeCache()
  new.loci.idx <- sim$propagation.window.idx[which(MAC>=2)]
  # get new loci index 
  CacheHaplotypes(sim$haps,
                  loci.idx = new.loci.idx,
                  hap.idx = hap.idx,
                  warn.singletons = TRUE)
  sim$map <- sim$map[new.loci.idx]
  sim$pos <- sim$pos[new.loci.idx]
} 

if (L()<=10000){stop("not enough variants to search for a target loci.")}
print(paste("region created and cached:",proc.time()[3] - start1,"seconds."))


n <- N() / ploidy

# set tuning parameters grid we wish to optimize over
HMM_pars <- as.data.table(expand.grid("neglog10Ne" = c(-1,2,6,8),"neglog10mu" = c(2,5,8,10,20),
                                      # pars$pars$mu need to be larger than 1
                                      "eff.num.clades" = NA_real_,
                                      "num.sprigs" = NA_integer_,
                                      "sd" = NA_real_,"qform"=NA_real_,"signal" = NA_real_,
                                      "DAC" = NA_real_,"causal.idx"=NA_real_,"target.idx" = NA_real_))

del.index <- c(which(HMM_pars$neglog10Ne==-1 &HMM_pars$neglog10mu >= 5),
               which(HMM_pars$neglog10Ne==2 &HMM_pars$neglog10mu >= 10))
HMM_pars <- HMM_pars[-del.index,]
# delete indices that usually cost very long time to run. 
# we don't want them in real data analysis for practicality anyways

HMM_pars$par_id <- 1:nrow(HMM_pars)


# set number of samples / tuning parameter
n.reps <- 1
res <- rbindlist(replicate(n.reps,HMM_pars,simplify = FALSE),idcol = "rep")
smt.res <- vector(mode = "list", length = nrow(res))


# specify options for selecting causal variant and target variant
candidate.target.idx <- seq(floor(L()/2)-20,len=40)
# select target variants from the middle of the region

core.idx <- (burnin.var.count+1):(L()-burnin.var.count)
print(paste0("There are ",L(), " variants in the cache."))

for(i in 1:n.reps){
  
  # make sure candidate causal variants are not in burnin region
  tuning.loci <- select.tuning.loci(candidate.target.idx,
                                    sim$map, 
                                    core.idx,
                                    Q,
                                    desired.cM.dist = 0.05,tol = 0.1) 

  target.idx <- tuning.loci$target.idx
  causal.idx <- tuning.loci$causal.idx

  row.idx <- which(res$rep==i)
  # locate the row for the current replication
  res$DAC[row.idx] <- sum(QueryCache(loci.idx = causal.idx))
  res$causal.idx[row.idx] <- causal.idx
  res$target.idx[row.idx] <- target.idx
  
  for(j in 1:nrow(HMM_pars)){
    
    neglog10Ne <- res$neglog10Ne[j]
    neglog10mu <- res$neglog10mu[j]

    pars <- Parameters(CalcRho(diff(sim$map),s = 10^(-neglog10Ne)),mu = 10^(-neglog10mu),
                       use.speidel = use.speidel) 
    # specify HMM parameters
    
    error.list <- list()
    warning.list <- list()
    tryCatch({
      if (exists("fwd")){
        rm(fwd)
      }
      if (exists("bck")){
        rm(bck)
      }
      # clean up existing fwd and bck matrices from ancestry inference machine
      gc()
      fwd <- MakeForwardTable(pars)
      bck <- MakeBackwardTable(pars)
      start1 <- proc.time()[3]
      Forward(fwd, pars, target.idx, nthreads = nthreads)
      Backward(bck, pars, target.idx, nthreads = nthreads)
      print(paste("Propagating HMM to target locus took",
                  signif(proc.time()[3] - start1,
                         digits = 3), "seconds."))
      # propagate HMM to target locus, ready to calculate matrix and call sprigs
      
      # calculate clade mat and sprigs
      if (fwd$l != bck$l){
        stop("forward and backward matrices are at different locations.")
      }

      M <- matrix(0,n,n)
      neigh <- CladeMat(fwd,bck,M,unit.dist = -log(pars$pars$mu),
                        thresh = thresh, max1var = max1var, nthreads = nthreads)
      # calculate neigh object, for effective number of clades 
      eff.num.clades <- neigh[[3]]
      
      # simulate pseudo-phenotype
      g <- c(locater:::Haps2Genotypes(QueryCache(loci.idx = causal.idx),ploidy = ploidy))
      g <- g/sqrt(sum(g^2))
      y <- rnorm(n)
      y <- matrix(y + g * (sqrt(80 + 2*log(eff.num.clades)) - sum(y*g)),ncol=1)
      # this is giving a lot of signal, also normalized based on number of effective clades.
    
      ########## start of calculate p-values and residuals for SD ################
      
      start1 <- proc.time()[3]
      sprigs <- Sprigs(neigh[[1]],old.sprigs = FALSE)
      PruneCladeMat(M,neigh,sprigs,prune="singleton.info")
      PruneCladeMat(M,neigh,sprigs,prune="sprigs")
      print(paste("Call and prune sprigs at target cost ",
                  signif(proc.time()[3] - start1,
                         digits = 3), "seconds."))

      gc()
      M <-  Matrix::symmpart(M)
      gc()
      ro.res <- TestSprigs(y, sprigs, ortho = TRUE,
                           Q = Q, use.forking = use.forking)
      ########## end of calculate p-values and residuals for SD ################
      
      ######## start of calculate p-values for QForm ############################
      
      start1 <- proc.time()[3]
      qf.res <- TestCladeMat(y = ro.res$y,M = M,Q = Q,
                             k = 0, use.forking = use.forking,
                             nthreads = 1L)
      # we recomment using nthreads = 1 for QForm
      
      print(paste("TestCladeMat at target cost ",
                  signif(proc.time()[3] - start1,
                         digits = 3), "seconds."))
      
      ######## end of calculate p-values for QForm ############################
      
      # Store results
      row.idx <- which(res$rep==i & res$par_id==j)
      smt.res[[row.idx]] <- TestCachedMarkers(y,A=A)
      res$eff.num.clades[row.idx] <- eff.num.clades
      res$num.sprigs[row.idx] <- sprigs$num.sprigs
      res$sd[row.idx] <- -log10(ro.res$p.value)
      res$qform[row.idx] <- qf.res$qform
      res$signal[row.idx] <- locater:::fishv(c(-log10(ro.res$p.value),qf.res$qform))
      # use Fisher's method to combine SD and QForm result
      
      print(paste(j,"of",nrow(res),"tuning runs done"))
    }, error = function(e){error.list <<- e},
    warning=function(w){warning.list <<- w})
  }
}



saveRDS(list(
  "sim" = sim,
  "res" = res, 
  "HMM_pars" = HMM_pars,
  "smt.res" = smt.res,
  "error.list"=error.list,
  "warning.list"=warning.list), paste0(res.storage.path,"tuning_res_",job.index,".rds"))

print(paste("FULL JOB DONE:",proc.time()[3] - start0,"SECONDS TOTAL"))

