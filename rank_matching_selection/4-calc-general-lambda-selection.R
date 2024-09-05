# this is to calculate a more general version of lambda (slope and intercept)
# and based on this lambda, select the best rank matched version of the phenotype

########################################################################
# collect data from all segments.
########################################################################

rm(list=ls())

require(locater)

experiment.name <- "WashU_CCDG"
test.name <- "2-pval-rank-matched"

data.dir <- paste0("/gpfs/gibbs/pi/ycgh/xw445/projects/",experiment.name,"/screening/to_share/locater_only/",test.name,"/")
grand.data.dir <- "/gpfs/gibbs/pi/ycgh/xw445/projects/"
# specify input directories

chr.segments <- read.table("/gpfs/gibbs/pi/ycgh/xw445/projects/WashU_CCDG/chr_segments_4000segs_new_pass.txt",header=TRUE)

rank.matched.y.rds <- paste0(grand.data.dir,experiment.name,"/screening/data/whole_genome_yale/rank_matched_y_seed28_100cp.rds")
y <- readRDS(rank.matched.y.rds)
phenos <- colnames(y)
# specify phenotype names

iters <- 1:4000

all.final.locater <- vector(mode = "list",length = length(iters))
# create a list to store all locater results

for (iter in iters){
  file <- paste0(data.dir,"locater_res_",as.character(iter),".txt")
  if (!file.exists(file)){
    print(paste0("file for iter ",iter," doesn't exist."))
    next
  } else{ 
    
    locater.res <- data.table::fread(file)
    all.final.locater[[match(iter,iters)]]  <-  locater.res

    rm(locater.res)
    gc()
  }
  
}


all.final.locater <- data.table::rbindlist(all.final.locater)

names(all.final.locater)[which(names(all.final.locater)=="rd")] <- "sd"
# fix the test name for sprig testing column

# Create a list with elements for each phenotype
all.final.locater <- split(all.final.locater, all.final.locater$phenotype)

# Order the list elements based on the order of interested phenotypes
all.final.locater <- all.final.locater[phenos]

saveRDS(all.final.locater,paste0(data.dir,"all_final_locater.rds"))

########################################################################
# calculate slope and intercept from data
########################################################################

NovelLambdaCalc <- function(x,lower,upper){
  
  # function to calculate general version of lambda
  
  # x: vector of -log10 p-values
  # lower: lower bound of the data used for lambda calculation
  # upper: upper bound of the data used for lambda calculation
  
  # Subset of data
  expected <- qexp(ppoints(length(x)),log(10))
  x <- sort(x)
  
  subset_expected <- expected[expected >= lower & expected <= upper]
  subset_x <- x[expected >= lower & expected <= upper]
  
  if (length(subset_x) < 2){
    return(list(slope=NA,intercept=NA,rsquared=NA))
  }
  
  # Linear regression
  lm_model <- lm(subset_x ~ subset_expected)
  
  # Slope
  slope <- coef(lm_model)[2]
  
  # Intercept
  intercept <- coef(lm_model)[1]
  
  # R-squared value
  rsquared <- summary(lm_model)$r.squared
  return(list(slope=slope,intercept=intercept,rsquared=rsquared))
}

viz.dir <- paste0(grand.data.dir,experiment.name,"/screening/to_share/inflation_QQ/",test.name,"/")
if(!dir.exists(viz.dir)){dir.create(viz.dir,recursive=TRUE)}


interested.phenos <- c("eGFR_MDRD_combined", "S_krea_combined", "S_totalc_combined", "ln_S_tottg_combined",
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

multicopy <- 100
# this is the number of versions in the rank matching process

lambda.table <- data.frame(matrix(ncol = 6, nrow = ncol(y)))
names(lambda.table) <- c("qform_slope","qform_intercept","qform_rsquared",
                         "sd_slope","sd_intercept","sd_rsquared")
row.names(lambda.table) <- colnames(y)
# initialize lambda table

for (pheno.idx in 1:length(interested.phenos)){
  gc()
  all.copies.index <- (multicopy * (pheno.idx-1)+1):(multicopy*pheno.idx)
  x <- all.final.locater[all.copies.index]
  
  for (test in c("qform","sd")){
    matrix <- sapply(x, function(table) table[[test]])
    
    tmp <- apply(matrix,MARGIN = 2, NovelLambdaCalc,lower=2,upper=2.5)
    # used the tail of the distirbution to calculate more general version of lambda
    
    tmp <- data.table:::rbindlist(tmp)
    
    names(tmp) <- paste0(test,"_",names(tmp))
    lambda.table[all.copies.index,names(tmp)] <- tmp
    
  }
  
}

write.table(lambda.table,file=paste0(viz.dir,"lambda_table_all_rankmatched_lm.txt"),
            quote = FALSE,row.names = TRUE,col.names = TRUE,sep = "\t")


########################################################################
# select best rank matched version based on slope and intercept of SD and QForm
########################################################################


lambda.table  <- read.table(paste0(viz.dir,"lambda_table_all_rankmatched_lm.txt"),
                             header = TRUE,sep = "\t")

findGoodLambdaName3 <- function(x,test,good.slope.range){
  
  # function to select rank matched versions that meet criteria
  
  # x: slope, intercept and r^2 information for all rank matched versions for this current phenotype
  # test: qform or sd
  # good.slope.range: range of slope values that are considered good
  
  slope <- x[,paste0(test,"_slope"),drop=FALSE]
  intercept <- x[,paste0(test,"_intercept")]
  r2 <- x[,paste0(test,"_rsquared")]
  # Extract the row names and values from the 'qform' column
  vector_names <- rownames(slope)
  slope_values <- as.vector(slope[[1]])

  # Create a named vector
  slope <- setNames(slope_values, vector_names)
  intercept <- setNames(intercept, vector_names)
  r2 <- setNames(r2, vector_names)
  
  if (test == "qform"){
    index <- which(slope >= good.slope.range[1] & slope <= good.slope.range[2] & intercept >= -0.1 & intercept <= 0.1 & r2 >= 0.8)
    # close to x=y
    index2 <- which(slope <= 0.8 &  slope >= 0.6 & intercept >= 0 & intercept <= 0.4 & r2 >= 0.8)
    # for qform, we want to additionally select versions that have small slope and positive intercept
  } else {
    index <- which(slope >= good.slope.range[1] & slope <= good.slope.range[2] & intercept >= -0.1 & intercept <= 0.1 & r2 >= 0.8)
    index2 <- NULL
    # close to x=y
    
  }

  slope[unique(c(index,index2))]
  
}

best.rank.matched.names <- rep(NA,length(interested.phenos))


best.rank.matched.lambda.table <- data.frame(matrix(ncol = 6, nrow = length(interested.phenos)))
names(best.rank.matched.lambda.table) <- c("qform_slope","qform_intercept","qform_rsquared",
                                     "sd_slope","sd_intercept","sd_rsquared")
# initialize best.rank.matched.names and their lambda table

for (pheno.idx in 1:length(interested.phenos)){ # iterate through all phenotypes

  gc()
  all.copies.index <- (multicopy * (pheno.idx-1)+1):(multicopy*pheno.idx)
  # get the index of all rank matched versions for this phenotype
  pheno <- interested.phenos[pheno.idx]
  
  slope.qform <- findGoodLambdaName3(lambda.table[all.copies.index,],"qform",c(0.7,1.2))
  slope.sd <- findGoodLambdaName3(lambda.table[all.copies.index,],"sd",c(0.8,1.1))
  # subset best rank matched version fo qform and sd
  # preserving the name and slope information
  
  shared_names <- intersect(names(slope.sd), names(slope.qform))
  # find the shared names between qform and sd
  
  if (length(shared_names)==0){
    next
  }
  sub.slope.qform <- slope.qform[shared_names]
  sub.slope.sd <- slope.sd[shared_names]
  # only preserve shard names and the slope information
  
  best.rank.matched.index <- which.max(pmin(sub.slope.qform,sub.slope.sd))
  # calculate the smaller lambda of qform and sd for each shared names
  # and select the name with the largest value
  
  best.rank.matched.name <- names(sub.slope.qform)[best.rank.matched.index]
  best.rank.matched.names[pheno.idx] <- best.rank.matched.name
  best.rank.matched.lambda.table[pheno.idx,]  <- lambda.table[match(best.rank.matched.name,row.names(lambda.table)),]
  
  # preserve the best rank matched phenotype, 
  # and their slope and intercept information
}

null.index <- which(is.na(best.rank.matched.names))
# these are the index for phenotypes that does not have a best rank matched phenotype
# we plan to use the rank normalize version of them

best.rank.matched.lambda.table <- best.rank.matched.lambda.table[complete.cases(best.rank.matched.lambda.table), ]

best.rank.matched.names <- best.rank.matched.names[!is.na(best.rank.matched.names)]
row.names(best.rank.matched.lambda.table) <- best.rank.matched.names
# subset best.rank.matched.lambda.table and best.rank.matched.names

saveRDS(best.rank.matched.names,paste0(data.dir,"best_rank_matched_seed28_lm.rds"))
write.table(best.rank.matched.lambda.table,file = paste0(data.dir,"best_rank_matched_lambda_table_seed28_lm.txt"),
            quote = FALSE,sep = "\t",row.names = TRUE,col.names = TRUE)

# after this, the user will need to find the phenotypes that cannot choose a best version of rank matched phenotypes
# then use the rank normalized version of them
# LOCATER could test multiple phenotypes together,
# so we recommend combine the phenotypes into matrices
# rows: individuals; columns: phenotypes



