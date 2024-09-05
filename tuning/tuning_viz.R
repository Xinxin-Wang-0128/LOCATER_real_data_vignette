##########  visualization ####################
####### we recommend run this process in R studio
# since saving a widget with selfcontained = TRUE requires pandoc
# usually it is installed in Rstudio

rm(list=ls())

splot <- function(x, y, z,z.name,z.range){
  # function for plot 3D plot for tuning result
  # x: x-axis data
  # y: y-axis data
  # z: z-axis data
  # z.name: z-axis label
  # z.range: z-axis range
  
  # set axis label
  axx <- list(title = "-log(Ne)")
  axy <- list(title = "-log(µ)")
  axz <- list(title = z.name,range = z.range)
  
  A <- akima::interp(x= x,y = y,z = z,
                     xo = seq(min(x), max(x), length = 50),
                     yo = seq(min(y), max(y), length = 50),
                     linear = TRUE, extrap = TRUE,duplicate ="mean")
  # interpolate data for a smooth surface
  
  A$z <- t(A$z)
  names(A) <- c("int","slope","neglog10pval")
  
  p <-  plotly::plot_ly() %>%
    plotly::add_surface(x= A$int , y= A$slope, z= A$neglog10pval) %>%
    plotly::add_markers(x= x , y= y, z= z) %>%
    # mark data points
    plotly::layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz, # set axis lable
                                # set the aspectmode to fix the ratio of 3 axes 
                                aspectmode='cube'))
  
  p
}



library(dplyr)
library(tibble)
require(data.table)
require(parallel)
require(akima)
require(plotly)


run_group <- "WashU_CCDG"
run_name <- "METSIM1"

res.storage.path <- paste0("/Users/xinxin.wang/Library/CloudStorage/Dropbox/",run_group,"/tuning/",run_name,"/")
all.res <- read.table(paste0(res.storage.path,"all_res_remove_nosignal.txt"),header=TRUE)
# load data


result  <- as.data.table(expand.grid("neglog10Ne" = c(-1,2,6,8),"neglog10mu" = c(2,5,8,10,20),
                                                "eff.num.clades" = NA_real_,
                                                "num.sprigs" = NA_integer_,
                                                "sd" = NA_real_,"qform"=NA_real_,"signal" = NA_real_))
# define result table
# neglog10Ne and neglog10mu are parameters we tuned

all.res <- as_tibble(all.res)

plot.dir <- paste0("/Users/xinxin.wang/Library/CloudStorage/Dropbox/",run_group,"/tuning/",run_name,"/plot/normalize_remove_nosignal/")
if (!dir.exists(plot.dir)){dir.create(plot.dir,recursive = TRUE)}
# define output directory

n.individuals <- 6795
n.haps <- n.individuals*2
AC.boundaries <- matrix(data=c(2,3,4,7,0.001*n.haps, 0.01*n.haps,0.1*n.haps,2,3,6,0.001*n.haps,0.01*n.haps,0.1*n.haps,n.haps),ncol = 2)
# after this, we created a AC.boundaries matrix, 
# which is a 6*2 matrix, each row represent a boundary for allele count
# we will use this as cutoffs when we visualize surfaces for allele count bins (separate data to ultra rare, rare, common etc)

####### visualize the trimmed mean of normalized signal - all data, LOCATER ########
z.name <- "trimmed mean - signal / causal SMT"
temp <- all.res %>%
  group_by(neglog10Ne,neglog10mu)

test <- temp %>% summarise(avg = mean(normalized.signal,trim=0.1))
# trim top and bottom 10% of data to calcuate trimmed mean

test$avg <- ifelse(!is.finite(test$avg),0,test$avg) 
p<- splot(test$neglog10Ne,test$neglog10mu,test$avg,z.name,z.range=c(0,0.5))

invisible(mcparallel({
  htmlwidgets::saveWidget(widget = p,
                          file = paste0(plot.dir,
                                        "0_1_trim_mean_normalized_surface.html"))
},detached = TRUE))

####### visualize the trimmed mean of normalized signal - separated by allele count bins, LOCATER ########
  
# Note: it is normal to have limited data for common variant bin, since there are not many common variants in the genome
# since LOCATER have its best power for low frequency and rare variants, we will focus on those bins
# so it is fine to not have data for common variants
  
for (i in 1:nrow(AC.boundaries)){
  temp <- all.res[all.res$DAC >=AC.boundaries[i,1]&all.res$DAC <=AC.boundaries[i,2], ]
  # seperate the data based on AC bin boundaries
  z.name <- "trimmed mean - signal / causal SMT"
  temp <- temp %>% 
    group_by(neglog10Ne,neglog10mu) 
  
  test <- temp %>% summarise(avg = mean(normalized.signal,trim=0.1))
  # trim top and bottom 10% of data to calcuate trimmed mean
  
  test$avg <- ifelse(!is.finite(test$avg),0,test$avg)
  p<- splot(test$neglog10Ne,test$neglog10mu,test$avg,z.name,z.range=c(0,0.5))
  
  invisible(mcparallel({
    htmlwidgets::saveWidget(widget = p,
                            file = paste0(plot.dir,
                                          "0_1_trim_mean_surface_AC_",i,".html"))
  },detached = TRUE))
}

####### visualize the trimmed mean of normalized signal - separated by allele count bins, SD ########
# Note: it is normal to have limited data for common variant bin, since there are not many common variants in the genome
# since LOCATER have its best power for low frequency and rare variants, we will focus on those bins
# so it is fine to not have data for common variants

all.res$normalized.sd <- all.res$sd / all.res$smt.causal
all.res$normalized.qform <- all.res$qform / all.res$smt.causal

for (i in 1:nrow(AC.boundaries)){

  temp <- all.res[all.res$DAC >=AC.boundaries[i,1]&all.res$DAC <=AC.boundaries[i,2], ]
  z.name <- "trimmed mean - sd / causal SMT"
  temp <- temp %>% 
    group_by(neglog10Ne,neglog10mu)  
  
  test <- temp %>% summarise(avg = mean(normalized.sd,trim=0.1))
  # trim top and bottom 10% of data to calcuate trimmed mean
  
  test$avg <- ifelse(!is.finite(test$avg),0,test$avg)
  p<- splot(test$neglog10Ne,test$neglog10mu,test$avg,z.name,z.range=c(0,0.6))
  
  invisible(mcparallel({
    htmlwidgets::saveWidget(widget = p,
                            file = paste0(plot.dir,
                                          "0_1_trim_mean_sd_surface_AC_",i,".html"))
  },detached = TRUE))
}

####### visualize the trimmed mean of normalized signal - separated by allele count bins, QForm ########
# Note: it is normal to have limited data for common variant bin, since there are not many common variants in the genome
# since LOCATER have its best power for low frequency and rare variants, we will focus on those bins
# so it is fine to not have data for common variants

for (i in 1:nrow(AC.boundaries)){
  temp <- all.res[all.res$DAC >=AC.boundaries[i,1]&all.res$DAC <=AC.boundaries[i,2], ]
  z.name <- "trimmed mean - qform / causal SMT"
  temp <- temp %>% 
    group_by(neglog10Ne,neglog10mu) 
  
  test <- temp %>% summarise(avg = mean(normalized.qform,trim=0.1))
  # trim top and bottom 10% of data to calcuate trimmed mean
  
  test$avg <- ifelse(!is.finite(test$avg),0,test$avg)
  p<- splot(test$neglog10Ne,test$neglog10mu,test$avg,z.name,z.range=c(0,0.5))
  
  invisible(mcparallel({
    htmlwidgets::saveWidget(widget = p,
                            file = paste0(plot.dir,
                                          "0_1_trim_mean_qform_surface_AC_",i,".html"))
  },detached = TRUE))
}


########## end of visualization ####################

