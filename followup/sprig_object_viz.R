########################
# the following is the investigation of sprigs object at the lead marker
########################
rm(list=ls())

library(kalis)
library(patchwork)
library(ggplot2)
library(data.table)

run_group <- "WashU_CCDG"
run_name <- "whole_genome_yale/wg_unrelated_METSIMonlyData_101phenos_Bjitter_local_10cM_clean"

grand.data.dir <- "/gpfs/gibbs/pi/ycgh/xw445/projects/"


res.storage.path <- paste0(grand.data.dir,run_group,"/screening/data/",run_name,"/")

draw_confidence_interval <- function(n, conf.points = 1000, conf.col = "gray", conf.alpha = 0.05) {
  # Limit the number of confidence points to the size of the data minus one
  conf.points <- min(conf.points, n - 1)
  
  # Initialize a matrix to store the x and y coordinates of the confidence interval polygon
  mpts <- matrix(nrow = conf.points * 2, ncol = 2)
  
  # Loop through each point to calculate the x and y coordinates
  for (i in seq(from = 1, to = conf.points)) {
    # Calculate the x-coordinate as the negative log10 of the quantile position
    mpts[i, 1] <- -log10((i - 0.5) / n)
    
    # Calculate the upper bound of the y-coordinate using the beta quantile function
    mpts[i, 2] <- -log10(qbeta(1 - conf.alpha / 2, i, n - i))
    
    # Mirror the x-coordinate for the lower bound
    mpts[conf.points * 2 + 1 - i, 1] <- -log10((i - 0.5) / n)
    
    # Calculate the lower bound of the y-coordinate using the beta quantile function
    mpts[conf.points * 2 + 1 - i, 2] <- -log10(qbeta(conf.alpha / 2, i, n - i))
  }
  
  # Draw the confidence interval polygon using the calculated coordinates
  polygon(mpts[, 1], mpts[, 2], col = conf.col, border = NA)
}

QQwInterval <- function(x, ...) {
  n <- length(x)
  
  # Initialize an empty plot
  plot.new()
  plot.window(xlim = range(qexp(ppoints(n), log(10))), ylim = range(sort(x)))
  
  # Draw the confidence interval
  draw_confidence_interval(n, conf.col = "gray", conf.alpha = 0.05)
  
  # Add the QQ plot points
  points(qexp(ppoints(n), log(10)), sort(x), ...)
  
  # Add axes and labels
  box()
  axis(1)
  axis(2)
  title(xlab = "Expected Quantiles", ylab = "Observed Quantiles")
  
  # Add a reference line
  abline(0, 1)
}



PlotSprigPhenosSep <- function(y, index, leading_sprig_pheno_list,indices) {
  
  data_list <- leading_sprig_pheno_list[indices]
  
  # Determine common x-axis limits based on data range
  common_x_limits <- c(-4.4,4.4)
  
  # Prepare data for plotting
  y_df <- data.frame(value = y[, index])
  data_points <- lapply(seq_along(data_list), function(i) {
    data.frame(x = data_list[[i]], y = rep(i, length(data_list[[i]])))
  }) %>% bind_rows()
  
  # Create segments data
  segments_data <- data_points %>% 
    group_by(y) %>%
    summarise(xmin = min(x), xmax = max(x)) %>%
    ungroup()
  
  # Histogram plot
  p1 <- ggplot(y_df, aes(x = value)) +
    geom_histogram(binwidth = 0.1, fill = "grey", color = "white") +
    theme_minimal() +
    labs(x = "", y = "Frequency") +
    #coord_cartesian(xlim = common_x_limits) +
    xlim(common_x_limits) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text = element_text(size = 14),  # Increase axis number size
          legend.text = element_text(size = 14), axis.title.x = element_text(size = 14),  # Increase x-axis label size
          axis.title.y = element_text(size = 14))
  
  p2 <- ggplot(data_points, aes(x = x, y = y)) +
    # Plot segments first
    geom_segment(data = segments_data, aes(x = xmin, xend = xmax, y = y, yend = y), color = "grey", linewidth = 1) +
    # Then plot points on top of the segments
    geom_point(aes(color = x), size = 3) +
    scale_color_gradient2(low = "black", high = "black", mid = "black", midpoint = 0) +
    xlim(common_x_limits) +
    theme_minimal() +
    labs(x = "Phenotype value", y = "", color = "") +  # Change 'color' to your desired legend title
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none",
          axis.text = element_text(size = 14),  # Increase axis number size
          legend.text = element_text(size = 14), axis.title.x = element_text(size = 14),  # Increase x-axis label size
          axis.title.y = element_text(size = 14))
  
  
  # p2 <- ggplot(data_points, aes(x = x, y = y)) +
  #   # Plot segments first
  #   geom_segment(data = segments_data, aes(x = xmin, xend = xmax, y = y, yend = y), color = "grey", linewidth = 1) +
  #   # Then plot points on top of the segments
  #   geom_point(aes(color = x), size = 3) +
  #   scale_color_gradient2(low = "#D55E00", high = "#D55E00", mid = "#56B4E9", midpoint = 0,limits=common_x_limits,guide = "colourbar") +
  #   xlim(common_x_limits) +
  #   theme_minimal() +
  #   labs(x = "Phenotype value", y = "") +
  #   theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "bottom")
  
  
  # Arrange the plots
  top.ratio <- 4/ (4 + 0.3*length(data_list) + 1)
  top.panel.height <- 1.8 * top.ratio
  bottom.panel.height <- 1.8 - top.panel.height
  gridExtra::grid.arrange(p1, p2, ncol = 1, heights = c(top.panel.height, bottom.panel.height))
  
}

plotYLabels <- function(y_labels) {
  # Create a data frame for plotting
  label_data <- data.frame(y = seq_along(y_labels), labels = y_labels)
  
  # Plot the labels
  p <- ggplot(label_data, aes(y = y, label = labels)) +
    # geom_text(aes(x = 1), hjust = 0, vjust = 0.5, size = 5) +
    scale_y_continuous(breaks = 1:length(y_labels), labels = y_labels) +
    theme_void() +
    theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"), # Adjust right margin to create space
          axis.text.y = element_text(size = 14),
          axis.ticks.y = element_blank()) +
    labs(y = "", x = "") +
    xlim(c(0.5, 1.5)) # Narrow xlim to keep the plot focused on labels
  
  print(p)
}


PlotSprigPCs <- function(A, interesting.individuals) {
  # Predefined set of distinguishable colors
  # colors <- c("#E6194B", "#3CB44B", "#4363D8", "#F58231", 
  #             "#911EB4", "#42D4F4", "#F032E6", "#BFEF45", 
  #             "#FABEBE", "#469990", "#DCBEFF", "#9A6324", 
  #             "#FFFAC8", "#800000", "#AAF0D1", "#000075")
  
  colors <- c("#000000", "#D55E00","#56B4E9","#F0E442","#009E73",
              "#CC79A7","#E69F00","#0072B2")
  
  # Ensure there are enough colors for the groups
  # if (length(interesting.individuals) > length(colors)) {
  #   stop("Not enough colors for the number of groups")
  # }
  
  # Convert A to a data frame if it's not already
  if (!is.data.frame(A)) {
    A <- as.data.frame(A)
  }
  
  # Define the layout for the plots
  layout(matrix(1:4, nrow = 2))
  
  # Loop through the PCs and plot
  for (i in 1:4) {
    pc_x <- A[[paste0("pc", 2 * i - 1)]]
    pc_y <- A[[paste0("pc", 2 * i)]]
    
    plot(pc_x, pc_y, pch = "*", col = "#AAAAAA", xlab = paste0("PC", 2 * i - 1),ylab= paste0("PC", 2 * i))
    
    
    if (length(interesting.individuals)> length(colors)){
      for (j in seq_along(colors)) {
        points(pc_x[interesting.individuals[[j]]], pc_y[interesting.individuals[[j]]], col = colors[j], pch = "*")
      }
      for (j in (length(colors)+1):length(interesting.individuals)) {
        points(pc_x[interesting.individuals[[j]]], pc_y[interesting.individuals[[j]]], col = colors[j-length(colors)], pch = 17)
      }
    } else {
      for (j in seq_along(interesting.individuals)) {
        points(pc_x[interesting.individuals[[j]]], pc_y[interesting.individuals[[j]]], col = colors[j], pch = "*")
      }
    }
    
  }
  
}




PlotSprigPCsLegend <- function(interesting.individuals,sprig_names) {
  colors <- c("#000000", "#D55E00", "#56B4E9", "#F0E442", "#009E73",
              "#CC79A7", "#E69F00", "#0072B2")
  
  # Adjust the length of colors to match the number of groups, if necessary
  extended_colors <- if (length(interesting.individuals) > length(colors)) {
    rep(colors, length.out = length(interesting.individuals))
  } else {
    colors
  }
  
  # Define the point shapes
  pch_values <- if (length(interesting.individuals) > length(colors)) {
    c(rep(8, length(colors)), rep(17, length(interesting.individuals) - length(colors)))
  } else {
    rep(8, length(interesting.individuals))
  }
  legend_labels <- paste("sprig", sprig_names)
  
  plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=c(0, 1), xaxt='n', yaxt='n')
  legend("center", legend = legend_labels, col = extended_colors, pch = pch_values, cex = 1.2)
}

# Example usage
# Assuming interesting.individuals is defined
# PlotSprigPCsLegend(interesting.individuals)

############# change these #############

significant.sprig.info.file <- fread(paste0(res.storage.path,"significant_sprig_info_highlights.txt"))
significant.sprig.info.file$significant.sprigs <- rep(NA,nrow(significant.sprig.info.file))
      # create a blank column to store the index of significant sprigs                                                         
viz.dir <- paste0(grand.data.dir,run_group,"/screening/to_share/locater_screening_local_plots/",run_name,"/misc/")
if(!dir.exists(viz.dir)){ dir.create(viz.dir,recursive = TRUE)}

for (row in 1:nrow(significant.sprig.info.file)){
  chr <- as.integer(significant.sprig.info.file[row,]$chr)
  local.index <- as.integer(significant.sprig.info.file[row,]$local.index)
  #all.interested.index <- 200
  start <- as.integer(significant.sprig.info.file[row,]$start)
  end <- as.integer(significant.sprig.info.file[row,]$end)
  
  interested.pheno <- significant.sprig.info.file[row,]$interested.pheno
 # index.of.significant.sprigs <- significant.sprig.info.file[row,]$index.of.significant.sprigs
  #index.of.significant.sprigs <- eval(parse(text=index.of.significant.sprigs))
  
  ############# ############# #############
  
  
  
  sprigs <- readRDS(paste0(res.storage.path,"sprigs_",local.index,"_chr",chr,".rds"))
  A <- readRDS(paste0(res.storage.path,"A_",local.index,"_chr",chr,".rds"))
  y <- readRDS(paste0(res.storage.path,"y_",local.index,"_chr",chr,".rds"))
  res.inside <- readRDS(paste0(res.storage.path,"res_inside_",local.index,"_chr",chr,".rds"))
  g <- readRDS(paste0(res.storage.path,"g_",local.index,"_chr",chr,".rds"))
 
  
  index <- which(colnames(y) == interested.pheno)
  h0 <- locater:::FitNull(y, A)
  
  # 
  smt.res <- locater:::TestMarker(h0, g, add.noise = "raw")
  # reconstruct the phenotype residuals provided for SD test.
  
  # Find how many of the top sprigs we care about
  cairo_pdf(paste0(viz.dir,"QQ_test_chr",chr,"_",local.index,"_inside.pdf"))
  QQwInterval(x=-log10(res.inside$u[,index]),pch = 4)
  dev.off()

  
  index.of.significant.sprigs <- ...
  # here we will need to find the significant sprigs based manually looking at the QQ plot
  
  # for all top sprigs, lets get the underlying phenotypes
  leading_sprig_pheno_list<- sapply(head(order(res.inside$u[,index]),16),
                                    function(x){ smt.res$y[ceiling(which(sprigs$assignments == x)/2),index]},simplify = F)
  
  
  
  #plot the whole distribution VS the sprigs
  sprig_names <- order(res.inside$u[,index])[index.of.significant.sprigs]
  
  index.of.significant.sprigs$significant.sprigs[row] <- paste(sprig_names,collapse = ",")
  
  
  plot.height <- 4 + 0.3*length(sprig_names) + 1
  cairo_pdf(paste0(viz.dir,"pheno_hist_new_chr",chr,"_",local.index,"_inside.pdf"),
             width=7, height=plot.height)
  PlotSprigPhenosSep(y=smt.res$y, index=index,
                  leading_sprig_pheno_list=leading_sprig_pheno_list,
                  indices=index.of.significant.sprigs)

  dev.off()

  # plot a seperate legend file
  cairo_pdf(paste0(viz.dir,"pheno_labels_chr",chr,"_",local.index,"_inside.pdf"),width=2, height=0.3*length(sprig_names)+1)
  plotYLabels(sprig_names)
  dev.off()
  

  #plot the PC values for individuals carrying the sprig
  interested.sprig.indices <- order(res.inside$u[,index])[index.of.significant.sprigs]

  interesting.individuals <- vector(mode="list",length=length(interested.sprig.indices))
  for (i in 1:length(interested.sprig.indices)){
    interesting.individuals[[i]] <- ceiling(which(sprigs$assignments == interested.sprig.indices[i])/2)
  }
  

  cairo_pdf(paste0(viz.dir,"PC_interesting_individuals_chr",chr,"_",local.index,"_inside_additional.pdf"),width = 8, height = 8)
  PlotSprigPCs(A, interesting.individuals) 
  dev.off()
  
  
  
  # Save the legend plot to a file
  # Specify the desired file path for the legend image
  legend_file_path <- paste0(viz.dir,"PC_interesting_individuals_chr",chr,"_",local.index,"_inside_legend.pdf")
  
  # Open a png device, plot the legend, then close the device
  cairo_pdf(filename = legend_file_path)
  PlotSprigPCsLegend(interesting.individuals,sprig_names=order(res.inside$u[,index])[index.of.significant.sprigs])
  dev.off()
  
}

# update the significant.sprig.info.file
write.table(significant.sprig.info.file, file=paste0(res.storage.path,"significant_sprig_info_highlights.txt"),
            sep="\t", quote=FALSE, row.names=FALSE)
