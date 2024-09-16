##################
# this is a example code to run whole genome screening for best rank matched phenotypes
##################
rm(list=ls())

#######################
# pipe in parameter
#######################
args <- commandArgs(trailingOnly = TRUE)

nthreads <- as.integer(args[1])
chr <- as.integer(args[2])
start <- as.integer(args[3])
end <- as.integer(args[4])

interested.pos <- args[5]

if (grepl(",", interested.pos)) {
  interested.pos <- as.integer(unlist(strsplit(interested.pos, ",")))
} else {
  interested.pos <- as.integer(interested.pos)
}
print(interested.pos)


##########################
# Specify Directories    #
##########################

run_group <- "WashU_CCDG"
run_name <- "sprig_object_generation"

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

diagnostics <- TRUE
# if using diagnostics mode

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

smt.res.external <- res
rm(res)
rm(smt.res)
# this is to not confuse with later res and smt.res inside TestLoci (which is expanded in this script in diagnostic mode)

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

  
  if (diagnostics){
    start_inside_0 <- proc.time()[3]
    target.loci <- sort(target.idx)
    h0 <- locater:::FitNull(y, A)
    m <- ncol(y)
    #test.opts <- list()
    new.test.opts <- locater:::make.call.clade(test.opts)
    test.opts <- new.test.opts[[1]]
    call.clade <- new.test.opts[[2]]
    sw.approx <- TRUE
    nthreads <- as.integer(nthreads)
    fwd <- MakeForwardTable(pars)
    bck <- MakeBackwardTable(pars)
    if (length(target.loci) == 1) {
      num.ckpts <- 0L
    }
    ckpt.first.locus <- FALSE
    if (length(target.loci) > 1 & num.ckpts) {
      if (ckpt.first.locus) {
        fwd.baseline <- MakeForwardTable(pars)
        Forward(fwd.baseline, pars, target.loci[1], nthreads)
        suppressMessages(Iter <- ForwardIterator(pars, num.ckpts -
                                                   1, target.loci, fwd.baseline, force.unif = TRUE))
      }
      else {
        suppressMessages(Iter <- ForwardIterator(pars, num.ckpts,
                                                 target.loci, force.unif = TRUE))
      }
    }
    M <- matrix(0, N()/2, ncol(fwd$alpha)/2)
    template.res <- test.opts
    template.res[, c("num.sprigs", "k","exit.status")] <- NA_integer_
    template.res[, c("precise")] <- FALSE
    template.res[, c("obs.qform", "obs.qform.T", "smt", "rd","qform")] <- NA_real_    
    template.res <- tidyr::expand_grid(template.res, 
                                       phenotype = if (is.null(colnames(y))) {1:m}else {colnames(y)})
    res <- replicate(length(target.loci), template.res, simplify = FALSE)
    clade.details <- as.list(seq_len(length(target.loci)))
    verbose <- TRUE
    if (verbose) {
      print(paste("Starting loop over", length(target.loci),
                  "target loci..."))
    }
    for (t in length(target.loci):1) {
      start1 <- proc.time()[3]
      if (num.ckpts) {
        Iter(fwd, pars, target.loci[t], nthreads)
      }else {
        if (fwd$l > target.loci[t]) {
          ResetTable(fwd)
        }
        Forward(fwd, pars, target.loci[t], nthreads = nthreads)
      }
      Backward(bck, pars, target.loci[t], nthreads = nthreads)
      if (verbose) {
        print(paste("Propagating HMM to target", length(target.loci) - t + 1L, 
                    "took", signif(proc.time()[3] - start1,
                                   digits = 3), "seconds."))
      }
      
      if (t %in% interested.local.idx){
        hapM <- kalis::DistMat(fwd, bck, type = "minus.min", nthreads = nthreads)
        saveRDS(hapM, file = paste0(res.storage.path,"hapM_",t,"_chr",chr,".rds"))
        rm(hapM)
      }
      
      start1 <- proc.time()[3]
      for (ii in 1:length(call.clade)) {
        start2 <- proc.time()[3]
        neigh <- CladeMat(fwd, bck, M, unit.dist = -log(pars$pars$mu),
                          thresh = call.clade[[i]]$opts$thresh, max1var = call.clade[[i]]$opts$max1var,
                          nthreads = nthreads)
        
        if (verbose) {
          print(paste("Calling CladeMat @ target", length(target.loci) -
                        t + 1L, "took", signif(proc.time()[3] - start2,
                                               digits = 3), "seconds."))
        }
        
        if (t %in% interested.local.idx){
          saveRDS(neigh, file = paste0(res.storage.path,"neigh_",t,"_chr",chr,".rds"))
          M_sym <- Matrix::symmpart(M)
          saveRDS(M_sym, file = paste0(res.storage.path,"M_sym_",t,"_chr",chr,".rds"))
          rm(M_sym)
        }
        
        start2 <- proc.time()[3]
        sprigs <- Sprigs(neigh[[1]], old.sprigs = call.clade[[i]]$opts$old.sprigs)
        PruneCladeMat(M, neigh, sprigs, prune = "singleton.info")
        PruneCladeMat(M, neigh, sprigs, prune = "sprigs")
        ############ ended here ######
        if (verbose) {
          print(paste("Calling and Pruning sprigs at target", length(target.loci) -
                        t + 1L, "out of", length(target.loci), "took",
                      signif(proc.time()[3] - start2, digits = 3),
                      "seconds."))
        }
        
        
        start2 <- proc.time()[3]
        gc()
        M <- Matrix::symmpart(M)
        gc()
        if (verbose) {
          print(paste("Symmetrizing M at target", length(target.loci) -
                        t + 1L, "out of", length(target.loci), "took",
                      signif(proc.time()[3] - start2, digits = 3),
                      "seconds."))
        }
        pre.clade <- call.clade[[ii]]$pre.clade
        for (j in 1:length(pre.clade)) {
          g <- t(locater:::Haps2Genotypes(QueryCache(target.loci[t]),
                                          ploidy = 2L, method = "additive"))
          
          smt.res <- locater:::TestMarker(h0, g, add.noise = pre.clade[[j]]$opts$smt.noise)
          test.clade <- pre.clade[[j]]$test.clade
          for (k in 1:length(test.clade)) {

            start2 <- proc.time()[3]
   
            p <- sprigs$num.sprigs
            x <- sprigs$assignments
            x[is.na(x)] <- p + 1L
            n.inside <- length(x)/2
            X <- Matrix::sparseMatrix(i = rep(1:n.inside, each = 2), j = x, 
                                      x = 1L, dims = c(n.inside, p + 1))
            X <- X[, -(p + 1)]
            
            res.inside <- rdistill::rdistill_pivot_par(y = smt.res$y, x = X, Q = smt.res$Q, 
                                                       max_num_causal = 16)
            ro.res <-list(p.value = res.inside$p_value, y = res.inside$y, 
                          num.layers = nrow(res.inside$u) - 
                            c(Matrix::colSums(is.na(res.inside$u))))
            
            if (verbose) {
              print(paste("Call TestSprigs at target",
                          length(target.loci) - t + 1L, "took", signif(proc.time()[3] -
                                                                         start2, digits = 3), "seconds."))
            }
            
            
            if (t %in% interested.local.idx){
              saveRDS(sprigs, file = paste0(res.storage.path,"sprigs_",t,"_chr",chr,".rds"))
              saveRDS(A, file = paste0(res.storage.path,"A_",t,"_chr",chr,".rds"))
              saveRDS(y,file = paste0(res.storage.path,"y_",t,"_chr",chr,".rds"))
              saveRDS(res.inside,file = paste0(res.storage.path,"res_inside_",t,"_chr",chr,".rds"))
              saveRDS(g,file = paste0(res.storage.path,"g_",t,"_chr",chr,".rds"))
            }
            
            start2 <- proc.time()[3]
            
            qf.res <- TestCladeMat(ro.res$y, M, smt.res$Q,
                                   k = locater:::calc_k_sched(test.clade[[k]]$opts$max.k),
                                   stop.eval.func = function(x, prop.var) {
                                     if (prop.var == 0) {
                                       all(msse.test(-log10(smt.res$p.value),
                                                     -log10(ro.res$p.value), -log10(x),
                                                     test.1.solo = TRUE) < test.clade[[k]]$opts$sw.thresh)
                                     }
                                     else {
                                       all(msse.test(-log10(smt.res$p.value),
                                                     -log10(ro.res$p.value), -log10(x),
                                                     test.1.solo = TRUE) < test.clade[[k]]$opts$eig.thresh)
                                     }
                                   }, calc.obs.T = test.clade[[k]]$opts$calc.obs.T,
                                   use.forking = use.forking, nthreads = 1L)
            clade.details[[t]] <- attr(qf.res, "details")
            
            
            if (verbose) {
              print(paste("Run TestCladeMat @ target",
                          length(target.loci) - t + 1L, "took", signif(proc.time()[3] -
                                                                         start2, digits = 3), "seconds."))
            }
            temp <- 1:m + m * (test.clade[[k]]$test.config -
                                 1L)
            
            res[[t]][temp, c("num.sprigs", "num.layers",
                             "obs.qform", "obs.qform.T", "exit.status",
                             "precise", "smt", "rd", "qform")] <- cbind(sprigs$num.sprigs,
                                                                        ro.res$num.layers, qf.res$obs, qf.res$obs.T,
                                                                        qf.res$exit.status, qf.res$precise, -log10(smt.res$p.value),
                                                                        -log10(ro.res$p.value), qf.res$qform)
            
            
            # }
            
          }
        }
      }
      if (verbose) {
        print(paste("Running diagnostic tests at target", length(target.loci) -
                      t + 1L, "out of", length(target.loci), "took",
                    signif(proc.time()[3] - start1, digits = 3),
                    "seconds."))
      }
    }
    if (verbose) {
      print(paste("Iterating over all", length(target.loci),
                  "target loci took", signif(proc.time()[3] - start_inside_0,
                                             digits = 3), "seconds."))
    }
    names(res) <- as.character(target.loci)
    res <- data.table::rbindlist(res, idcol = "locus.idx")
    res[, `:=`(locus.idx, as.integer(locus.idx))]
    res[, `:=`(tot, msse.test(smt, rd, qform, test.1.solo = TRUE))]
    names(clade.details) <- as.character(target.loci)
    attr(res, "details") <- clade.details
    
    exact.res <- res
    
  } else {
    
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
  }
  
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
  "smt.res" = smt.res.external), paste0(res.storage.path,"res_",chr,"_",start,"-",end,"_sprig_investigation.rds"))

print(paste("FULL JOB DONE:",proc.time()[3] - start0,"SECONDS TOTAL"))

