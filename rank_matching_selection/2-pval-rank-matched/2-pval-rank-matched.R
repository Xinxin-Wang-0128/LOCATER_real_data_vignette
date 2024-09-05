##################
# this is a example code to evaluate the p-value distribution for all rank matched versions of all phenotypes
##################
rm(list=ls())

#######################
# pipe in parameter
#######################
nthreads <- as.integer(commandArgs(TRUE)[1])
job.index <- commandArgs(TRUE)[2]
set.seed(job.index)
# this job index will be in the range of 1-number of small segments
# it is essential to give each small segment a different seed to make sure p-values from SD is not highly correlated.
# correlation will lead to inflation. 

##########################
# Specify Directories    #
##########################

run_group <- "WashU_CCDG"
run_name <- "2-pval-rank-matched"

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
num.ckpts <- 2L
# checkpoints for speed

use.forking <- FALSE
use.speidel <- TRUE
# if using speidel et al version of LS model

burnin.count <- 6000
# make sure N variants on each end of the segment is burnin region and will not be evaluated

total.targets <- 30000
total.segs <- 4000
# total targets: how many variants you want to evaluate whole genome
# total.segs: how many small segments you seperate the genome to
# in preprocess, we segmented the genome into 4000 segments

#################################
# generate sim object, an object storing information for each segment
#################################

segments.file <- paste0(grand.data.dir,run_group,"/chr_segments_4000segs_new_pass.txt")
segments.table <- read.table(segments.file,header=TRUE)
# this is to specify the boundaries for all small segments 

burnin.file <- paste0(grand.data.dir,run_group,"/phasing-vcfs.compute1.MANIFEST.burnin.table.txt")
burnin.table <- read.table(burnin.file,header=TRUE)

job.row <- segments.table[job.index,]
orig.seg.index <- which(burnin.table$chr==job.row$chr & 
                          burnin.table$core.start <= job.row$core.start.pos & 
                          burnin.table$core.end >= job.row$core.end.pos)
# find the original segment that this small segment is from

file_prefix <- paste(burnin.table[orig.seg.index,1:3],collapse = "-")

# read in new segment file
segment.info <- segments.table[job.index,]


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

# use position info to subset pos and map for the core region of this small segment
caching.index <- which(pos>= segment.info$core.start.pos & pos <= segment.info$core.end.pos)
# first get all variants in the core region
caching.index <- c((min(caching.index) - burnin.count):(max(caching.index) + burnin.count))
# expand this region to include some burnin
caching.index <- caching.index[caching.index >= 1 & caching.index <= length(pos)]
# make sure the burnin region is inside the big segment
caching.index <- intersect(caching.index,which(pass==TRUE))
# remove anything that did not pass

if (length(caching.index)<burnin.count){stop("not enough pass var in this smaller segment.")}
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


rank.matched.y.rds <- paste0(grand.data.dir,run_group,"/screening/data/whole_genome_yale/rank_matched_y_seed28_100cp.rds")
y <- readRDS(rank.matched.y.rds)
# load in all rank matched versions of all phenotypes

n.phenos <- ncol(y)


#################################################################
# specify target loci, run LOCATER
#################################################################

target.loci.smt <- 1:L()
# we don't subset the indices based on SMT value, so setting target.loci.smt to 1:L()

core.idx <- which(sim$pos >= segment.info$core.start.pos & sim$pos <= segment.info$core.end.pos)
# make sure only the core region is included in the target loci

print(paste0("There are ",L(), " variants in the cache."))

target.idx <- intersect(core.idx,target.loci.smt)
target.idx <- intersect(target.idx,which(sim$pass))
# make sure only pass variants are included

var.per.seg <- ceiling(total.targets/total.segs)
d <- length(target.idx) / var.per.seg
# calculate how many variants each segment should have
# and the distance between them

# Create a sequence of indices with equal distances
indices <- seq( 0.5 * d, length(target.idx) - 0.5 * d, by = d)

# Round the indices to the nearest integers
rounded_indices <- round(indices)

target.idx <- target.idx[rounded_indices]

print(paste0("Upstream burnin: ",min(target.idx),"; downstream burnin: ",L()-max(target.idx)))

out.list <- vector("list",nrow(HMM_pars))
for(i in 1:nrow(HMM_pars)){
  
  pars <- Parameters(CalcRho(diff(sim$map),s = 10^(-HMM_pars$neglog10Ne[i])),
                     mu = 10^(-HMM_pars$neglog10mu[i]),
                     use.speidel = use.speidel)
  # specify parameters
  
  start_locater <- proc.time()[3]
  
  test.opts <- data.frame(
    "smt.noise" = c("raw"), # Clade-free testing options (eg: SMT, might be more complex)
    "thresh" = 0.2, # threshold for the distinction of noise and signal in distance matrix
    "old.sprigs" = c(FALSE), # Clade calling options
    "max1var" = c(TRUE), 
    "max.k" = 512, # maximum number of eigen values
    "sw.thresh" = 0,  # this means we only use Satterthwaite approximation for p-value calculation in QForm
    "eig.thresh" = 0,  # setting it to zero means no eigen decomposition
    "calc.obs.T" = TRUE# Clade testing options
  )
  
  # see more details about options in LOCATER package
  
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
  "smt.res" = NULL), paste0(res.storage.path,"res_",job.index,".rds"))

print(paste("FULL JOB DONE:",proc.time()[3] - start0,"SECONDS TOTAL"))

