##################
# this is a example code to run putative investigation
##################
rm(list=ls())

#######################
# pipe in parameter
#######################
nthreads <- as.integer(commandArgs(TRUE)[1])
chr <- as.integer(commandArgs(TRUE)[2])
start <- as.integer(commandArgs(TRUE)[3])
end <- as.integer(commandArgs(TRUE)[4])

##########################
# Specify Directories    #
##########################

run_group <- "WashU_CCDG"
run_name <- "2-wg-screening"

grand.data.dir <- "/gpfs/gibbs/pi/ycgh/xw445/projects/"

res.storage.path <- paste0(grand.data.dir,run_group,"/screening/to_share/",run_name,"/")
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
# Specify parameters
####################################################################

ploidy <- 2L
# diploid
num.ckpts <- 10L
# checkpoints for speed

use.forking <- FALSE
use.speidel <- TRUE
# if using speidel et al version of LS model

segment.length.cM <- 10
# segment length in cM
#################################
# generate sim object, an object storing information for each segment
#################################

segments.file <- paste0(grand.data.dir,run_group,"/chr_segments_4000segs_new_pass.txt")
segments.table <- read.table(segments.file,header=TRUE)
# this is to specify the boundaries for all small segments 

burnin.file <- paste0(grand.data.dir,run_group,"/phasing-vcfs.compute1.MANIFEST.burnin.table.txt")
burnin.table <- read.table(burnin.file,header=TRUE)

orig.seg.index <- which(burnin.table$chr==chr & 
                          burnin.table$core.start <= start & 
                          burnin.table$core.end >= end)
# find the original segment that this small segment is from

file_prefix <- paste(burnin.table[orig.seg.index,1:3],collapse = "-")

# read in new segment file
#segment.info <- segments.table[job.index,]


start0 <- proc.time()[3]
start1 <- proc.time()[3]


sim <- list()
# create haplotype file path
haps.path <- paste0(input.genome.dir,file_prefix,".phased_multi_merged_unique_ID_ancestral_split_sample_AC_var_metsim6806_nosingleton.hap.gz")
sim$haps <- haps.path

# read in map and position information for whole large segment
map.dir <- paste0(grand.data.dir,run_group,"/recomb_map_new_pass_6806metsim/")
map.file <- read.table(paste0(map.dir,file_prefix,"_AF_weighted_cM.txt"),header=TRUE)
map <- map.file$cM

pos <- map.file$POS
pass <- map.file$REF != "N" & nchar(map.file$REF) + nchar(map.file$ALT) <= 2
# need to have high quality ancestral coding, and need to be SNPs.

caching.index <- which(pos>= start & pos <= end)
# first get all variants in the core region
idx.middle <- floor(mean(caching.index))
map.middle <- map[idx.middle]
boundary.cM <- c(map.middle-segment.length.cM/2,map.middle+segment.length.cM/2)
# get the cM of boundaries 
caching.index <- which(map >= boundary.cM[1] & map <= boundary.cM[2])
# sometimes boundary cM is smaller than the real boundary of this big segment. This could will automatically take care of this.
caching.index <- intersect(caching.index,which(pass==TRUE))
# remove anything that did not pass

if (length(caching.index)<2000){stop("not enough pass var in this smaller segment.")}
print(paste0("length of this segment is ",diff(map[range(caching.index)]),"cM"))

sim$pos <- pos[caching.index]
sim$map <- map[caching.index]

sim$caching.index <- caching.index

print("sim object created!")

#################################################################
# subset genome and A based on interested individuals, cache haplotypes
#################################################################

# since we have decided on samples to run, I load the names of these samples and subset genome and A matrix based on that
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

# read in background covariate matrix 
pc.file <- paste0(pc.dir,pc.prefix,".tsv")
A <- read.table(pc.file,stringsAsFactors = FALSE,header = TRUE)

# sort A based on sample with phenotypes
sample.in.pc <- match(samples.in.experiment,A$ID)
sample.in.pc <- na.omit(sample.in.pc)
A <- A[sample.in.pc,]

A <- cbind(rep(1,nrow(A)),A[,2:11])
# as required by LOCATER, add a column of "1" before the background covariate matrices
A <- as.matrix(A)
Q <- qr.Q(qr(A))

# Load haplotypes and recombination map
CacheHaplotypes(sim$haps,
                loci.idx = sim$caching.index,
                hap.idx = hap.idx,
                warn.singletons = TRUE)
print(paste0("before subset, N=",N(),", L=",L()))

# since we subsetted the haplotypes. 
# need to see if we have singletons or monomorphics.
MAC <- rep(NA,L())
for (i in 1:L()){
  ac <- sum(QueryCache(loci.idx = i)) 
  MAC[i] <- min(ac,N()-ac)
}
sim$MAC <- MAC

if (sum(MAC==1 | MAC==0)>0){
  # if subsetting individuals created more singletons and monomorphics,
  ClearHaplotypeCache()
  new.loci.idx <- sim$caching.index[which(MAC>=2)]
  # get new loci index 
  CacheHaplotypes(sim$haps,
                  loci.idx = new.loci.idx,
                  hap.idx = hap.idx,
                  warn.singletons = TRUE)
  sim$map <- map[new.loci.idx]
  sim$pos <- pos[new.loci.idx]
  sim$pass <- pass[new.loci.idx]
  sim$new.loci.idx <- new.loci.idx
  sim$MAC <- MAC[which(MAC >= 2)]
  # update map, pos, pass, new.loci.idx and MAC to match up with newly cached haplotypes
  
} 

print(paste("region created and cached:",proc.time()[3] - start1,"seconds."))


#################################################################
# specify HMM parameters, read in rank matched y
#################################################################

n <- N() / ploidy

HMM_pars <- as.data.table(expand.grid("neglog10Ne" = c(6),"neglog10mu" = c(8),
                                      # pars$pars$mu need to be larger than 1
                                      "eff.num.clades" = NA_real_,
                                      "num.sprigs" = NA_integer_,
                                      "rd" = NA_real_,"qform"=NA_real_,"signal" = NA_real_,
                                      "DAC" = NA_real_,"causal.idx"=NA_real_,"target.idx" = NA_real_))

HMM_pars$par_id <- 1:nrow(HMM_pars)

y <- readRDS(paste0(grand.data.dir,run_group,
                    "/screening/data/whole_genome_yale/best_jitter_and_rank_normalized_METSIM_6795samples.rds"))
# load in phenotypes that contain the best rank matched and rank normalized

n.phenos <- ncol(y)

#################################################################
# run SMT
#################################################################
res <- replicate(n.phenos,
                 list("smt.max" = NA_real_,
                      "smt.results"=list()),
                 simplify = FALSE)

# initialize res to store SMT results

print(paste0("Start SMT, N=",N(),", L=",L()))
start1 <- proc.time()[3]

smt.res <- TestCachedMarkers(y, A = A)
# run SMT for all phenos

jj <- 0
for(j in 1:n.phenos){
  jj <- jj+1
  res[[j]]$smt.results <- c(smt.res[,jj])
  res[[j]]$smt.max <- max(res[[j]]$smt.results,na.rm=T)
  # store SMT results and the max value
}

print(paste("phenotypes loaded and tested with SMT:",proc.time()[3] - start1,"seconds."))

#################################################################
# specify target loci, run LOCATER
#################################################################

target.loci.smt <- FindTargetVars(map, min.cM = 0.1, initial.targets = smt.res, smt.thresh = 3)
# find target loci based on SMT results: only preserve variants that have at least one variant that is significant
# (-log10 p value > 3)
# also make sure the adjacent target variant are at most 0.1 cM away from each other

core.idx <- which(sim$pos >= start & sim$pos <= end)

print(paste0("There are ",L(), " variants in the cache."))

target.idx <- intersect(core.idx,target.loci.smt)
target.idx <- intersect(target.idx,which(sim$pass))
print(paste0("Upstream burnin: ",min(target.idx),"; downstream burnin: ",L()-max(target.idx)))
# make sure target index are in core and also pass the filter

out.list <- vector("list",nrow(HMM_pars))
for(i in 1:nrow(HMM_pars)){
  
  pars <- Parameters(CalcRho(diff(sim$map),s = 10^(-HMM_pars$neglog10Ne[i])),
                     mu = 10^(-HMM_pars$neglog10mu[i]),
                     use.speidel = use.speidel)
  # specify HMM parameters
  
  start_locater <- proc.time()[3]
  
  test.opts <- data.frame(
    "smt.noise" = c("raw"), # Clade-free testing options (eg: SMT, might be more complex)
    "thresh" = 0.2, # threshold for the distinction of noise and signal in distance matrix
    "old.sprigs" = c(FALSE), # Clade calling options
    "max1var" = c(TRUE),
    "max.k" = 512, # maximum number of eigen values
    "sw.thresh" = 6, # in QForm, we start with Satterthwaite approximation
    #          if approximate -log10pval > sw.thresh, then we do eigen decomposition 
    "eig.thresh" = 6,  # if -log10pval < eig.thresh, stop eigen decomposition early
    "calc.obs.T" = TRUE# Clade testing options
  )

  exact.res <- TestLoci(y, pars, target.loci = target.idx,
                          A = A,
                          test.opts = test.opts,
                          verbose = TRUE,
                          num.ckpts = num.ckpts,
                          ckpt.first.locus = FALSE,
                          use.forking = use.forking,
                          nthreads = nthreads)
  # y: phenotypes
  # pars: HMM parameters
  # target.loci: target loci that LOCATER will evaluate
  # A: background covariate matrix
  # test.opts: testing options
  # verbose: print out progress
  # num.ckpts: number of checkpoints
  # ckpt.first.locus: whether use the first locus as one checkpoint
  #  use.forking: whether use forking to parallelize
  # nthreads: number of threads
  
  if (exists("exact.res")){
    out.list[[i]] <- exact.res
  } else {
    stop("Job Failed.")
  }
}

print(paste("locater done:",proc.time()[3] - start_locater,"seconds."))

saveRDS(list(
  "sim" = sim,
  "out.list" = out.list,  
  "HMM_pars" = HMM_pars,
  "smt.res" = res), paste0(res.storage.path,"res_",chr,"_",start,"-",end,"_sprig_investigation.rds"))

print(paste("FULL JOB DONE:",proc.time()[3] - start0,"SECONDS TOTAL"))

