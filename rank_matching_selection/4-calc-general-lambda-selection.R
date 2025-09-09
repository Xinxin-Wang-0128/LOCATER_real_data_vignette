rm(list=ls())
require(locater)
require(data.table)

########################################################################
# Utility Functions
########################################################################

# Create directory if it doesn't exist
create_directory <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
}

# Collect locater results from all segments
collect_locater_data <- function(data_dir, iters) {
  all_final_locater <- vector(mode = "list", length = length(iters))
  
  for (iter in iters) {
    file <- paste0(data_dir, "locater_res_", as.character(iter), ".txt")
    if (!file.exists(file)) {
      print(paste0("file for iter ", iter, " doesn't exist."))
      next
    }
    
    locater_res <- fread(file)
    all_final_locater[[match(iter, iters)]] <- locater_res
    rm(locater_res)
    gc()
  }
  
  combined_data <- rbindlist(all_final_locater)
  names(combined_data)[which(names(combined_data) == "rd")] <- "sd"  # Fix column name
  return(combined_data)
}

# Process locater data by phenotype
process_by_phenotype <- function(combined_data, phenos) {
  split_data <- split(combined_data, combined_data$phenotype)
  ordered_data <- split_data[phenos]  # Order by phenotype interest
  return(ordered_data)
}

# Calculate lambda table for all phenotypes
calculate_lambda_table <- function(locater_data, interested_phenos, multicopy, novel_lambda_func) {
  lambda_table <- data.frame(
    matrix(ncol = 6, nrow = length(interested_phenos) * multicopy)
  )
  names(lambda_table) <- c(
    "qform_slope", "qform_intercept", "qform_rsquared",
    "sd_slope", "sd_intercept", "sd_rsquared"
  )
  row.names(lambda_table) <- unlist(lapply(interested_phenos, function(p) 
    paste0(p, "_", 1:multicopy)))
  
  for (pheno_idx in seq_along(interested_phenos)) {
    gc()
    copy_indices <- (multicopy * (pheno_idx - 1) + 1):(multicopy * pheno_idx)
    x <- locater_data[copy_indices]
    
    for (test in c("qform", "sd")) {
      # Extract test values across all copies
      test_matrix <- sapply(x, function(table) table[[test]])
      
      # Calculate lambda parameters for each copy
      tmp_results <- apply(test_matrix, MARGIN = 2, novel_lambda_func, lower = 2, upper = 2.5)
      tmp_df <- rbindlist(tmp_results)
      names(tmp_df) <- paste0(test, "_", names(tmp_df))
      
      # Update lambda table
      lambda_table[copy_indices, names(tmp_df)] <- tmp_df
    }
  }
  
  return(lambda_table)
}

# Select best rank matched version
select_best_versions <- function(lambda_table, interested_phenos, multicopy) {
  best_names <- rep(NA, length(interested_phenos))
  best_lambda_table <- data.frame(
    matrix(ncol = 6, nrow = length(interested_phenos))
  )
  names(best_lambda_table) <- c(
    "qform_slope", "qform_intercept", "qform_rsquared",
    "sd_slope", "sd_intercept", "sd_rsquared"
  )
  
  for (pheno_idx in seq_along(interested_phenos)) {
    gc()
    copy_indices <- (multicopy * (pheno_idx - 1) + 1):(multicopy * pheno_idx)
    pheno <- interested_phenos[pheno_idx]
    
    # Get valid slope ranges for each test
    slope_qform <- findGoodLambdaName3(
      lambda_table[copy_indices, ], 
      "qform", 
      good.slope.range = c(0.7, 1.2)
    )
    slope_sd <- findGoodLambdaName3(
      lambda_table[copy_indices, ], 
      "sd", 
      good.slope.range = c(0.8, 1.1)
    )
    
    # Find overlapping good versions
    shared_names <- intersect(names(slope_sd), names(slope_qform))
    if (length(shared_names) == 0) next
    
    # Select best version based on combined criteria
    sub_slope_qform <- slope_qform[shared_names]
    sub_slope_sd <- slope_sd[shared_names]
    best_idx <- which.max(pmin(sub_slope_qform, sub_slope_sd))
    best_name <- names(sub_slope_qform)[best_idx]
    
    # Update results
    best_names[pheno_idx] <- best_name
    best_lambda_table[pheno_idx, ] <- lambda_table[match(best_name, row.names(lambda_table)), ]
  }
  
  # Clean up results
  valid_idx <- !is.na(best_names)
  best_names <- best_names[valid_idx]
  best_lambda_table <- best_lambda_table[valid_idx, ]
  row.names(best_lambda_table) <- best_names
  
  return(list(
    best_names = best_names,
    best_lambda_table = best_lambda_table
  ))
}

########################################################################
# Core Analysis Functions
########################################################################

# Calculate general version of lambda (slope, intercept, R-squared)
NovelLambdaCalc <- function(x, lower, upper) {
  expected <- qexp(ppoints(length(x)), log(10))
  x <- sort(x)
  
  # Subset data based on bounds
  subset_mask <- expected >= lower & expected <= upper
  subset_expected <- expected[subset_mask]
  subset_x <- x[subset_mask]
  
  if (length(subset_x) < 2) {
    return(list(slope = NA, intercept = NA, rsquared = NA))
  }
  
  # Linear regression
  lm_model <- lm(subset_x ~ subset_expected)
  
  return(list(
    slope = coef(lm_model)[2],
    intercept = coef(lm_model)[1],
    rsquared = summary(lm_model)$r.squared
  ))
}

# Identify good lambda candidates
findGoodLambdaName3 <- function(x, test, good.slope.range) {
  slope_col <- paste0(test, "_slope")
  intercept_col <- paste0(test, "_intercept")
  r2_col <- paste0(test, "_rsquared")
  
  slope <- x[[slope_col]]
  intercept <- x[[intercept_col]]
  r2 <- x[[r2_col]]
  
  # Create named vectors
  slope <- setNames(slope, rownames(x))
  intercept <- setNames(intercept, rownames(x))
  r2 <- setNames(r2, rownames(x))
  
  # Define selection criteria
  base_idx <- which(
    slope >= good.slope.range[1] & slope <= good.slope.range[2] &
      intercept >= -0.1 & intercept <= 0.1 &
      r2 >= 0.8
  )
  
  # Additional criteria for qform
  if (test == "qform") {
    extra_idx <- which(
      slope <= 0.8 & slope >= 0.6 &
        intercept >= 0 & intercept <= 0.4 &
        r2 >= 0.8
    )
    selected_idx <- unique(c(base_idx, extra_idx))
  } else {
    selected_idx <- base_idx
  }
  
  return(slope[selected_idx])
}

########################################################################
# Main Execution
########################################################################

# Configuration
experiment_name <- "WashU_CCDG"
test_name <- "2-pval-rank-matched"
grand_data_dir <- "/gpfs/gibbs/pi/ycgh/xw445/projects/"
data_dir <- paste0(grand_data_dir, experiment_name, "/screening/to_share/locater_only/", test_name, "/")
chr_segments_file <- "/gpfs/gibbs/pi/ycgh/xw445/projects/WashU_CCDG/chr_segments_4000segs_new_pass.txt"
rank_matched_y_rds <- paste0(grand_data_dir, experiment_name, "/screening/data/whole_genome_yale/rank_matched_y_seed28_100cp.rds")
iters <- 1:4000
multicopy <- 100  # Number of rank matching versions

# Create visualization directory
viz_dir <- paste0(grand_data_dir, experiment_name, "/screening/to_share/inflation_QQ/", test_name, "/")
create_directory(viz_dir)

# Load input data
chr_segments <- read.table(chr_segments_file, header = TRUE)
y <- readRDS(rank_matched_y_rds)
phenos <- colnames(y)

# 1. Collect and process locater data
all_final_locater_combined <- collect_locater_data(data_dir, iters)
all_final_locater <- process_by_phenotype(all_final_locater_combined, phenos)
saveRDS(all_final_locater, paste0(data_dir, "all_final_locater.rds"))

# 2. Define interested phenotypes
interested_phenos <- c(
  "eGFR_MDRD_combined", "S_krea_combined", "S_totalc_combined", "ln_S_tottg_combined",
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
  "Ile_combined", "Alb_combined", "Pyr_combined"
)

# 3. Calculate lambda parameters
lambda_table <- calculate_lambda_table(
  locater_data = all_final_locater,
  interested_phenos = interested_phenos,
  multicopy = multicopy,
  novel_lambda_func = NovelLambdaCalc
)
write.table(
  lambda_table,
  file = paste0(viz_dir, "lambda_table_all_rankmatched_lm.txt"),
  quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t"
)

# 4. Select best rank matched versions
lambda_table <- read.table(
  paste0(viz_dir, "lambda_table_all_rankmatched_lm.txt"),
  header = TRUE, sep = "\t"
)

best_matches <- select_best_versions(
  lambda_table = lambda_table,
  interested_phenos = interested_phenos,
  multicopy = multicopy
)

# Save final results
saveRDS(
  best_matches$best_names,
  paste0(data_dir, "best_rank_matched_seed28_lm.rds")
)
write.table(
  best_matches$best_lambda_table,
  file = paste0(data_dir, "best_rank_matched_lambda_table_seed28_lm.txt"),
  quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE
)
