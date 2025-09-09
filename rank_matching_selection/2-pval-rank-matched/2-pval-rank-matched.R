rm(list=ls())

####################################################################
# Function 1: Create simulation object 
####################################################################
create_sim_object <- function(job.index, segments.table, burnin.table, input.genome.dir, map.dir, burnin.count) {
  sim <- list()
  
  # Get segment information for current job
  job.row <- segments.table[job.index, ]
  
  # Find original segment index
  orig.seg.index <- which(burnin.table$chr == job.row$chr & 
                            burnin.table$core.start <= job.row$core.start.pos & 
                            burnin.table$core.end >= job.row$core.end.pos)
  file_prefix <- paste(burnin.table[orig.seg.index, 1:3], collapse = "-")
  
  # Create haplotype file path
  sim$haps <- paste0(input.genome.dir, file_prefix, ".phased_multi_merged_unique_ID_ancestral_split_sample_AC_var_metsim6806_nosingleton.hap.gz")
  
  # Read map and position information
  map.file <- read.table(paste0(map.dir, file_prefix, "_AF_weighted_cM.txt"), header = TRUE)
  map <- map.file$cM
  pos <- map.file$POS
  pass <- map.file$REF != "N" & nchar(map.file$REF) + nchar(map.file$ALT) <= 2
  
  # Determine caching index with burnin regions
  caching.index <- which(pos >= job.row$core.start.pos & pos <= job.row$core.end.pos)
  caching.index <- c((min(caching.index) - burnin.count):(max(caching.index) + burnin.count))
  caching.index <- caching.index[caching.index >= 1 & caching.index <= length(pos)]
  caching.index <- intersect(caching.index, which(pass == TRUE))
  
  if (length(caching.index) < burnin.count) {
    stop("Not enough passing variants in this segment.")
  }
  
  # Update simulation object with position information
  sim$pos <- pos[caching.index]
  sim$map <- map[caching.index]
  sim$caching.index <- caching.index
  
  print("sim object created!")
  return(sim)
}

####################################################################
# Function 2: Subset genome and cache haplotypes 
####################################################################
subset_genome_and_cache_haplotypes <- function(sim, samples.in.experiment, sample.file, pc.file, ploidy) {
  # Process sample information
  sample.table <- read.table(sample.file, header = TRUE)
  sample.idx <- match(samples.in.experiment, sample.table$sample)
  sample.idx <- na.omit(sample.idx)
  sample.idx <- as.vector(sample.idx)
  
  # Validate sample order and create haplotype indices
  if (!all.equal(sort(sample.idx), sample.idx)) {
    stop("Phenotype file needs sorting for proper haplotype subsetting.")
  }
  hap.idx <- sort(c(2 * sample.idx, 2 * sample.idx - 1))
  hap.idx <- as.integer(hap.idx)
  
  # Process background covariate matrix
  A <- read.table(pc.file, stringsAsFactors = FALSE, header = TRUE)
  sample.in.pc <- match(samples.in.experiment, A$ID)
  sample.in.pc <- na.omit(sample.in.pc)
  A <- A[sample.in.pc, ]
  A <- cbind(rep(1, nrow(A)), A[, 2:11])  # Add intercept as required by LOCATER
  A <- as.matrix(A)
  Q <- qr.Q(qr(A))
  
  # Cache haplotypes and check for singletons/monomorphics
  CacheHaplotypes(sim$haps,
                  loci.idx = sim$caching.index,
                  hap.idx = hap.idx,
                  warn.singletons = TRUE)
  print(paste0("before subset, N=", N(), ", L=", L()))
  
  # Filter out singletons and monomorphics
  MAC <- sapply(1:L(), function(i) {
    ac <- sum(QueryCache(loci.idx = i))
    min(ac, N() - ac)
  })
  sim$MAC <- MAC
  
  if (sum(MAC == 1 | MAC == 0) > 0) {
    ClearHaplotypeCache()
    new.loci.idx <- sim$caching.index[which(MAC >= 2)]
    CacheHaplotypes(sim$haps,
                    loci.idx = new.loci.idx,
                    hap.idx = hap.idx,
                    warn.singletons = TRUE)
    
    # Update simulation object with filtered loci
    sim$map <- sim$map[which(MAC >= 2)]
    sim$pos <- sim$pos[which(MAC >= 2)]
    sim$pass <- pass[new.loci.idx]
    sim$new.loci.idx <- new.loci.idx
    sim$MAC <- MAC[which(MAC >= 2)]
  }
  
  print(paste("region created and cached:", proc.time()[3] - start1, "seconds."))
  return(list(sim = sim, A = A, Q = Q))
}

####################################################################
# Function 3: Setup LOCATER parameters and target loci
####################################################################
setup_locater_parameters <- function(sim, segment.info, total.targets, total.segs, rank.matched.y.rds) {
  # Setup HMM parameters
  HMM_pars <- as.data.table(expand.grid(
    "neglog10Ne" = c(6),
    "neglog10mu" = c(8),
    "eff.num.clades" = NA_real_,
    "num.sprigs" = NA_integer_,
    "rd" = NA_real_,
    "qform" = NA_real_,
    "signal" = NA_real_,
    "DAC" = NA_real_,
    "causal.idx" = NA_real_,
    "target.idx" = NA_real_
  ))
  HMM_pars$par_id <- 1:nrow(HMM_pars)
  
  # Load rank-matched phenotypes
  y <- readRDS(rank.matched.y.rds)
  n.phenos <- ncol(y)
  
  # Determine target loci
  target.loci.smt <- 1:L()  # Include all loci
  core.idx <- which(sim$pos >= segment.info$core.start.pos & sim$pos <= segment.info$core.end.pos)
  target.idx <- intersect(core.idx, target.loci.smt)
  target.idx <- intersect(target.idx, which(sim$pass))
  
  # Subset to appropriate number of targets per segment
  var.per.seg <- ceiling(total.targets / total.segs)
  d <- length(target.idx) / var.per.seg
  indices <- seq(0.5 * d, length(target.idx) - 0.5 * d, by = d)
  rounded_indices <- round(indices)
  target.idx <- target.idx[rounded_indices]
  
  print(paste0("Upstream burnin: ", min(target.idx), "; downstream burnin: ", L() - max(target.idx)))
  
  return(list(HMM_pars = HMM_pars, y = y, target.idx = target.idx))
}




#######################
# pipe in parameter
#######################
nthreads <- as.integer(commandArgs(TRUE)[1])
job.index <- commandArgs(TRUE)[2]
set.seed(job.index)

##########################
# Specify Directories    #
##########################
run_group <- "WashU_CCDG"
run_name <- "2-pval-rank-matched"
grand.data.dir <- "/gpfs/gibbs/pi/ycgh/xw445/projects/"

res.storage.path <- paste0(grand.data.dir, run_group, "/screening/to_share/", run_name, "/")
if(!dir.exists(res.storage.path)){dir.create(res.storage.path, recursive = TRUE)}

input.genome.dir <- paste0(grand.data.dir, run_group, "/haplotypes/segmented_phased_merged_reannotated_ancestral_split_new_pass_singleton_filtered_6806metsim_hap/")
pc.dir <- paste0(grand.data.dir, run_group, "/pc_pheno/pc_finn_only/")
pc.prefix <- "finnish_SD8_iteration4_to_keep_table"

#####################
# Load Dependencies #
#####################
.libPaths(c("/rchrist/lib/R.4", .libPaths()))
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
num.ckpts <- 2L
use.forking <- FALSE
use.speidel <- TRUE
burnin.count <- 6000
total.targets <- 30000
total.segs <- 4000

####################################################################
# Main execution
####################################################################
start0 <- proc.time()[3]

# Read segment and burnin tables
segments.file <- paste0(grand.data.dir, run_group, "/chr_segments_4000segs_new_pass.txt")
segments.table <- read.table(segments.file, header = TRUE)
burnin.file <- paste0(grand.data.dir, run_group, "/phasing-vcfs.compute1.MANIFEST.burnin.table.txt")
burnin.table <- read.table(burnin.file, header = TRUE)

# Create simulation object
sim <- create_sim_object(
  job.index = job.index,
  segments.table = segments.table,
  burnin.table = burnin.table,
  input.genome.dir = input.genome.dir,
  map.dir = paste0(grand.data.dir, run_group, "/recomb_map_new_pass_6806metsim/"),
  burnin.count = burnin.count
)

# Subset genome and cache haplotypes
start1 <- proc.time()[3]
samples.in.experiment <- readRDS("/gpfs/gibbs/pi/ycgh/xw445/projects/WashU_CCDG/screening/data/whole_genome_yale/wg_unrelated_METSIMonlyData_101phenos_local_10cM/samples_in_experiment.RDS")
sample.file <- paste0(input.genome.dir, 
                      paste(burnin.table[orig.seg.index, 1:3], collapse = "-"),
                      ".phased_multi_merged_unique_ID_ancestral_split_sample_AC_var_metsim6806_nosingleton.samples")

genome_subset_result <- subset_genome_and_cache_haplotypes(
  sim = sim,
  samples.in.experiment = samples.in.experiment,
  sample.file = sample.file,
  pc.file = paste0(pc.dir, pc.prefix, ".tsv"),
  ploidy = ploidy
)
sim <- genome_subset_result$sim
A <- genome_subset_result$A
Q <- genome_subset_result$Q

# Setup LOCATER parameters
locater_setup <- setup_locater_parameters(
  sim = sim,
  segment.info = segments.table[job.index, ],
  total.targets = total.targets,
  total.segs = total.segs,
  rank.matched.y.rds = paste0(grand.data.dir, run_group, "/screening/data/whole_genome_yale/rank_matched_y_seed28_100cp.rds")
)
HMM_pars <- locater_setup$HMM_pars
y <- locater_setup$y
target.idx <- locater_setup$target.idx

#################################################################
# Run LOCATER
#################################################################
n <- N() / ploidy
out.list <- vector("list", nrow(HMM_pars))

for(i in 1:nrow(HMM_pars)){
  pars <- Parameters(
    CalcRho(diff(sim$map), s = 10^(-HMM_pars$neglog10Ne[i])),
    mu = 10^(-HMM_pars$neglog10mu[i]),
    use.speidel = use.speidel
  )
  
  start_locater <- proc.time()[3]
  
  test.opts <- data.frame(
    "smt.noise" = c("raw"),
    "thresh" = 0.2,
    "old.sprigs" = c(FALSE),
    "max1var" = c(TRUE),
    "max.k" = 512,
    "sw.thresh" = 0,
    "eig.thresh" = 0,
    "calc.obs.T" = TRUE
  )
  
  exact.res <- TestLoci(
    y, pars, target.loci = target.idx,
    A = A,
    test.opts = test.opts,
    verbose = TRUE,
    num.ckpts = num.ckpts,
    ckpt.first.locus = FALSE,
    use.forking = use.forking,
    nthreads = nthreads
  )
  
  if (exists("exact.res")) {
    out.list[[i]] <- exact.res
  } else {
    stop("Job Failed.")
  }
}

print(paste("locater done:", proc.time()[3] - start_locater, "seconds."))

# Save results
saveRDS(
  list(
    "sim" = sim,
    "out.list" = out.list,
    "HMM_pars" = HMM_pars,
    "smt.res" = NULL
  ),
  paste0(res.storage.path, "res_", job.index, ".rds")
)

print(paste("FULL JOB DONE:", proc.time()[3] - start0, "SECONDS TOTAL"))
