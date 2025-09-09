################################
# this script is to generate rank matched phenotypes
###############################
rm(list=ls())
set.seed(28)

# core function to generate rank matched phenotypes
# y_orig: Original phenotype matrix
# n_copies: Number of copies to generate
# pheno_names: Vector of phenotype names
# sample_ids: Vector of sample IDs
#return generated rank-matched phenotype matrix
generate_rank_matched_phenotypes <- function(y.orig, n.copies, interested.phenos, sample.ids) {
  # Initialize result matrix
  rank.matched.y <- matrix(0, 
                           nrow = nrow(y.orig), 
                           ncol = length(interested.phenos) * n.copies)
  
  # Set names
  colnames(rank.matched.y) <- paste0(rep(interested.phenos, each = n.copies), 1:n.copies)
  rownames(rank.matched.y) <- sample.ids
  
  # Generate values for each phenotype and copy
  for (i in 1:length(interested.phenos)) {
    # Calculate column indices for current phenotype
    idx <- (n.copies * (i - 1) + 1):(n.copies * i)
    
    # Calculate ranks for original phenotype
    ranked <- rank(y.orig[, i], na.last = "keep")
    
    # Generate each copy
    for (j in idx) {
      # Original rank-matching formula preserved
      rank.matched.y[, j] <- qnorm(
        sort(c(0.5, runif(n - 2) * (n - 1) + 0.5, n - 0.5) / n)[ranked]
      )
    }
  }
  
  return(rank.matched.y)
}



# load the data
run_group <- "WashU_CCDG"
grand.data.dir <- "/gpfs/gibbs/pi/ycgh/xw445/projects/"


pheno.file<- paste0(grand.data.dir,run_group,"/pc_pheno/QT/qt_finnseq_20161017.5k_ids.ped")
pheno.table <- read.table(pheno.file,header=TRUE,sep="\t",comment.char = "")
# load phenotype file for all individuals and all phenotypes

interested.phenos  <- c("eGFR_MDRD_combined", "S_krea_combined", "S_totalc_combined", "ln_S_tottg_combined",
                        "ln_P_CRP_combined", "height_combined", "bmi_combined", "weight_combined",
                        "systbp_combined", "diastbp_combined", "pulsepress_combined", "ln_P_adipon_combined",
                        "waist_combined", "S_ApoA1_combined", "S_ApoB_combined", "S_hdlc_combined",
                        "S_hdlc_males", "hip_combined", "whr_combined", "S_ldlc_combined",
                        "fatmass_combined", "IDL_C_combined", "XXL_VLDL_P_combined", "ln_IDL_TG_combined",
                        "L_HDL_C_combined", "L_HDL_P_combined", "ln_L_HDL_TG_combined", "L_LDL_C_combined",
                        "L_LDL_P_combined", "ln_L_LDL_TG_combined", "L_VLDL_C_combined", "L_VLDL_P_combined",
                        "ln_L_VLDL_TG_combined", "M_HDL_C_combined", "M_HDL_P_combined", "ln_M_HDL_TG_combined",
                        "M_LDL_C_combined", "M_LDL_P_combined", "ln_M_LDL_TG_combined", "M_VLDL_C_combined",
                        "M_VLDL_P_combined", "ln_M_VLDL_TG_combined", "S_HDL_C_combined", "S_HDL_P_combined",
                        "ln_S_HDL_TG_combined", "S_LDL_C_combined", "S_LDL_P_combined", "ln_S_LDL_TG_combined",
                        "S_VLDL_C_combined", "S_VLDL_P_combined", "ln_S_VLDL_TG_combined", "XL_HDL_C_combined",
                        "XL_HDL_P_combined", "ln_XL_HDL_TG_combined", "XL_VLDL_C_combined", "XL_VLDL_P_combined",
                        "ln_XL_VLDL_TG_combined", "XS_VLDL_C_combined", "XS_VLDL_P_combined",
                        "ln_XS_VLDL_TG_combined", "XXL_VLDL_C_combined", "ln_XXL_VLDL_TG_combined",
                        "HDL2_C_combined", "HDL3_C_combined", "VLDL_C_combined", "Remnant_C_combined",
                        "IDL_P_combined", "ln_HDL_TG_combined", "ln_LDL_TG_combined", "ln_VLDL_TG_combined",
                        "Ala_combined", "Phe_combined", "ln_Ace_combined", "ln_AcAce_combined",
                        "His_combined", "Gp_combined", "Tyr_combined", "Val_combined", "DHA_combined",
                        "FAw3_combined", "FAw6_combined", "LA_combined", "MUFA_combined", "PUFA_combined",
                        "DHA_FA_combined", "FAw3_FA_combined", "FAw6_FA_combined", "MUFA_FA_combined",
                        "PUFA_FA_combined", "TotPG_combined", "PC_combined", "SM_combined", "TotCho_combined",
                        "Leu_combined", "SFA_combined", "SFA_FA_combined", "Cit_combined", "Gln_combined",
                        "Ile_combined", "Alb_combined", "Pyr_combined")

# specify interested phenotypes
samples.in.experiment <- readRDS("/gpfs/gibbs/pi/ycgh/xw445/projects/WashU_CCDG/screening/data/whole_genome_yale/wg_unrelated_METSIMonlyData_101phenos_local_10cM/samples_in_experiment.RDS")
# get the list of interested individuals
pheno.table <- pheno.table[pheno.table$IND_ID %in% samples.in.experiment,]
pheno.table <- pheno.table[match(samples.in.experiment, pheno.table$IND_ID),]

# load phenotype
y.orig <- pheno.table[,match(interested.phenos,names(pheno.table))]
y.orig <- as.matrix(y.orig)
# subset phenotype 

n.copies <- 100
# this is the number of rank matching versions that you want to generate
n <- length(samples.in.experiment)

#rank.matched.y <-  matrix(0, nrow = nrow(y.orig), ncol = length(interested.phenos)*n.copies)  
# Initialize matrix 

#colnames(rank.matched.y) <- paste0(rep(interested.phenos, each = n.copies), 1:n.copies)
#rownames(rank.matched.y) <- samples.in.experiment
# specify column and row names

# Generate values for each column
#for (i in 1:length(interested.phenos)) {
#  idx <- (n.copies * (i-1)+1):(n.copies*i)
# specify the columns that this phenotype will fill in
#  ranked <- rank(y.orig[,i],na.last="keep")
#  for (j in idx){
#   rank.matched.y[, j] <- qnorm(sort(c(0.5, runif(n - 2) * (n - 1) + 0.5, n - 0.5) / n)[ranked])
# generate rank matched phenotype with small truncation on both ends of the Gaussian distribution
#  }

#}
# Generate the rank-matched phenotypes using core function
rank.matched.y <- generate_rank_matched_phenotypes(
  y.orig = y.orig,
  n.copies = n.copies,
  interested.phenos = interested.phenos,
  sample.ids = samples.in.experiment
)

rank.matched.y.rds <- paste0(grand.data.dir,run_group,"/screening/data/whole_genome_yale/rank_matched_y_seed28_100cp.rds")
saveRDS(rank.matched.y,file = rank.matched.y.rds)
# save the rank matched phenotype


