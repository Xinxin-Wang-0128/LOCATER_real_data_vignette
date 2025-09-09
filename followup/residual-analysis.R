rm(list=ls())

require(locater)
require(scattermore)
library(ggplot2)
require(dplyr)
require(data.table)

experiment.name <- "WashU_CCDG"
test.name <- "wg_unrelated_METSIMonlyData_101phenos_Bjitter_local_10cM_clean"
grand.data.dir <- "/gpfs/gibbs/pi/ycgh/xw445/projects/"


smt.dir <- paste0("/gpfs/gibbs/pi/ycgh/xw445/projects/",experiment.name,
                  "/screening/data/whole_genome_yale/smt_only/",test.name,"/")
res.storage.path <- paste0(grand.data.dir,experiment.name,"/screening/data/whole_genome_yale/",test.name,"/")

smt.threshold <- 8.144726
conditional.threshold <- 8.144726-2

y <- readRDS(paste0(grand.data.dir,experiment.name,
                    "/screening/data/whole_genome_yale/best_jitter_and_rank_normalized_METSIM_6795samples.rds"))
phenos <- colnames(y)


pattern <- "^res_(\\d+)_(\\d+)-(\\d+)_rores_y_at_(\\d+)\\.rds$"
file_list <- list.files(res.storage.path, pattern = pattern)

chr_values <- numeric()
start_values <- numeric()
end_values <- numeric()
rores_index <- numeric()

# Loop through each file and parse the filename
for (file_name in file_list) {
  match <- regmatches(file_name, regexec(pattern, file_name))
  
  if (length(match) > 0 && all(sapply(match, function(x) !any(x == -1)))) {
    chr_values <- c(chr_values, as.numeric(match[[1]][2]))
    start_values <- c(start_values, as.numeric(match[[1]][3]))
    end_values <- c(end_values, as.numeric(match[[1]][4]))
    rores_index <- c(rores_index, as.numeric(match[[1]][5]))
  }
}

# Create a data frame with parsed values
parsed_data <- data.frame(
  chr = chr_values,
  start = start_values,
  end = end_values,
  rores_index = rores_index
)
parsed_data$interested.phenos <- c("ln_M_VLDL_TG_combined13","HDL2_C_combined92",
                                   "ln_M_HDL_TG_combined83","S_ApoA1_combined47","Remnant_C_combined35",                                  
                                   "ln_S_HDL_TG_combined25","HDL3_C_combined79",
                                   "MUFA_combined51","ln_VLDL_TG_combined3",
                                   "L_VLDL_P_combined4")

parsed_data$readable.phenos <- c("Triglycerides in Medium VLDL","HDL2 Cholesterol",
                                 "Triglycerides in Medium HDL","Apolipoprotein A1","Remnant Cholesterol",                                  
                                 "Triglycerides in Small HDL","HDL3 Choleterol",
                                 "MUFA","Triglycerides in VLDL",
                                 "Concentration of Large VLDL Particles")


parsed_data$plot.start <- c(104727888,61243278,
                            49038347,49217040,44308684,
                            45306012,15718536,
                            72840219,72832414,
                            73043687) 
parsed_data$plot.end <- c(105927888,62443278,
                          50253146,50417040,45809149,
                          46561659,16918536,
                          74243687,74243687,
                          74243687)

parsed_data$lead.marker <- c(105327888,61843278,49653146, 49817040,44922203,
                             45906012,16318536,73440219,73482065,73643687)

options(scipen = 999)

for (i in 1:nrow(parsed_data)){
  
  chr <- parsed_data$chr[i]
  interested.phenotype <- parsed_data$interested.phenos[i]
  pheno.index <- match(interested.phenotype,phenos)
  readable.phenotype <- parsed_data$readable.phenos[i]
  lead.marker <- parsed_data$lead.marker[i]
  rores.index <- parsed_data$rores_index[i]
  start <- parsed_data$start[i]
  end <- parsed_data$end[i]
  
  viz.dir <- paste0("/gpfs/gibbs/pi/ycgh/xw445/projects/",experiment.name,"/screening/data/whole_genome_yale/locater_screening_local_plots/",test.name,"/gc_controlled_alt/chr",chr,"/")
  if (!dir.exists(viz.dir)) {
    dir.create(viz.dir, recursive = TRUE)
  }
  
  original_pattern <- paste0("res_",chr,"_",start,"-",end,".rds")
  condition_on_SMT_pattern <- paste0("res_",chr,"_",start,"-",end,"_conditioned_on_",lead.marker,".rds")
  condition_on_SMT_RD_pattern <- paste0("res_",chr,"_",start,"-",end,"_rores_y_at_",rores.index,".rds")
  
  original.res <- readRDS(file = paste0(res.storage.path,original_pattern))
  first_res <- readRDS(paste0(res.storage.path,condition_on_SMT_pattern))
  second_res <- readRDS(paste0(res.storage.path,condition_on_SMT_RD_pattern))
  
  # Process first condition results
  first.in.plot <- cbind(first_res$sim$pos,first_res$smt.res[[pheno.index]]$smt.results)
  colnames(first.in.plot) <- c("smt.pos","smt.p")
  first.in.plot <- as.data.frame(first.in.plot)
  
  # Process second condition results
  tmp.second <- cbind(second_res$sim$pos,second_res$smt.res[[pheno.index]]$smt.results)
  colnames(tmp.second) <- c("smt.pos","smt.p")
  tmp.second <- as.data.frame(tmp.second)
  
  # Calculate differences and identify significant points
  diff <- first.in.plot$smt.p - tmp.second$smt.p
  index <- which(diff > 1 & first.in.plot$smt.p > 3)
  
  # Separate points into regular dots and triangles
  first.in.plot.dot <- first.in.plot[-index,]
  first.in.plot.triangle <- first.in.plot[index,]
  
  tmp.second.dot <- tmp.second[-index,]
  tmp.second.triangle <- tmp.second[index,]
  
  # Determine plot range based on triangle points
  if (nrow(tmp.second.triangle) > 0) {
    plot.start <- min(tmp.second.triangle$smt.pos) - 1e4
    plot.end <- max(tmp.second.triangle$smt.pos) + 1e4
  } else {
    # Fallback range if no triangle points
    plot.start <- min(first.in.plot$smt.pos) - 1e4
    plot.end <- max(first.in.plot$smt.pos) + 1e4
  }
  
  # Filter data to plot range
  first.in.plot.dot <- first.in.plot.dot %>%
    dplyr::filter(smt.pos >= plot.start & smt.pos <= plot.end)
  
  tmp.second.dot <- tmp.second.dot %>%
    dplyr::filter(smt.pos >= plot.start & smt.pos <= plot.end)
  
  # Combine data for Manhattan plot
  combined_plot_data_dot <- rbind(
    transform(first.in.plot.dot, group_type = "PS"),
    transform(tmp.second.dot, group_type = "PD")
  )
  
  # Calculate plot maximum y-value
  plot.max <- max(c(combined_plot_data_dot$smt.p, 
                    first.in.plot.triangle$smt.p, 
                    tmp.second.triangle$smt.p)) + 2
  
  # Generate Manhattan plot
  cairo_pdf(paste0(viz.dir,"Manhattan_",interested.phenotype,"_",plot.start,"_",plot.end,"_sequential_conditional_1st&2nd.pdf"), 
            width = 16, height = 8)
  
  color_map <- c("PS" = "#F0E442", "PD" = "#009E73")
  colors <- color_map[combined_plot_data_dot$group_type]
  
  plot(combined_plot_data_dot$smt.pos/1e6, combined_plot_data_dot$smt.p, 
       col = colors, pch = "*", cex = 0.6 * 1.3,
       xlim = plot_lim/1e6, ylim = c(0, plot.max), 
       xlab = "Position (Mb)", ylab = "-log10(p-value)",
       main = paste0(readable.phenotype," (Sequential Conditional)"),
       cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
  
  # Add triangle points
  points(first.in.plot.triangle$smt.pos/1e6, first.in.plot.triangle$smt.p, 
         pch = 24, col = "black", bg = color_map["PS"], cex = 1.2 * 1.3)
  points(tmp.second.triangle$smt.pos/1e6, tmp.second.triangle$smt.p, 
         pch = 24, col = "black", bg = color_map["PD"], cex = 1.2 * 1.3)
  
  # Add legend
  legend("topright", legend = names(color_map), pch = 24, 
         col = "black", pt.bg = color_map, cex = 1.2 * 1.3)
  
  dev.off()
  
  # Generate scatter plot
  cairo_pdf(paste0(viz.dir,"scatter_",interested.phenotype,"_",plot.start,"_",plot.end,"_serial_conditional.pdf"), 
            width = 10, height = 10)
  
  plot(first.in.plot.dot$smt.p, tmp.second.dot$smt.p, 
       xlab = "PS -log10(p-value)", ylab = "PD -log10(p-value)",
       main = paste0(interested.phenotype," (Serial Conditional)"),
       xlim = c(0, plot.max - 1), ylim = c(0, plot.max - 1), pch = "*",
       cex.main = 1.8, cex.axis = 1.8, cex.lab = 1.8)
  
  points(first.in.plot.triangle$smt.p, tmp.second.triangle$smt.p, 
         col = "red", pch = 17)
  abline(0, 1, col = "red")
  
  dev.off()
}
