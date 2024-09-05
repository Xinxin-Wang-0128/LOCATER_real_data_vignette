# Divide Genome Into Segments for Analysis
# for efficiency, usually we need to seperate the genome into small segments
# so that LOCATER could run multiple parallel jobs on each segment

rm(list=ls())

require(purrr)

grand.data.dir <- "/gpfs/gibbs/pi/ycgh/xw445/projects/WashU_CCDG/"
map.dir<- paste0(grand.data.dir,"recomb_map_new_pass_6806metsim")
# this is the path to recombination map for original segments 
# (usually we already segment up the genome for phasing)
# which contain positions for variants (POS column)
# we want to use this to count number of variants 

# Although we're probably not going to keep all of these loci, 
# they should give us a rather robust way of evenly breaking
# up the genome into segments a priori.  

original.segment.count <- 155
# after phasing, we have 155 segments whole genome
total.segments <- 4000
# this is the target segment count


loci.pos  <- as.list(1:original.segment.count)

burnin.file <- paste0(grand.data.dir,"phasing-vcfs.compute1.MANIFEST.burnin.table.txt")
burnin.table <- read.table(burnin.file,header=TRUE)
# burnin.table need to have information to query the map file (in this case, chr, filename.start and filename.end)
# it also need to have information about core region (non overlapping region) 
# since phasing usually is done with segments with overlap
# in this case, core.start and core.end are the start and end of core region

for (iter in 1:length(loci.pos)){
  
  file_prefix <- paste(burnin.table[iter,1:3],collapse = "-")
  # query the map file with prefix
  table <- read.table(paste0(map.dir,"/",file_prefix,"_AF_weighted_cM.txt"),header=TRUE)
  tmp <- table$POS[table$REF!="N"]
  # this is to only preserve variants with ancestral allele annotation
  # so that the estimation for small segment length in the next step is more accurate
  loci.pos[[iter]] <- tmp
  print(iter)
}

# loci pos should be a list with 155 elements
# each element is the loci pos of all variants in each original segment

orig.seg.lengths <- sapply(loci.pos,length)
# number of variants in each original segment

num.segments <- rep(1,original.segment.count)
# this is to store number of small segments for each current segment

if(total.segments <= original.segment.count){stop("total.segments desired must be greater than original segment count.")}

core.loci.pos  <- as.list(1:original.segment.count)

for(iter in 1:original.segment.count){
  orig.segment.core <- c(NA,NA)
  orig.segment.core[1] <- burnin.table$core.start[iter]
  orig.segment.core[2] <- burnin.table$core.end[iter]
  core.loci.pos[[iter]] <- loci.pos[[iter]][which(loci.pos[[iter]] >=orig.segment.core[1] &  loci.pos[[iter]] <=orig.segment.core[2])]
}
# core.loci.pos should be a list with 155 elements
# each element is the loci pos of all variants in each original segment (core region only)

core.lengths <- sapply(core.loci.pos,length)

for(i in (original.segment.count+1):total.segments){
  num.segments[which.max(floor(core.lengths/num.segments))] <- num.segments[which.max(floor(core.lengths/num.segments))] + 1
}

# this segmentation is trying to make every small segment have relatively the same number of variants
# this aim to make each LOCATER run to have relatively same run time

# we still need to add burnin region on both ends during tuning, screening etc.
# to make sure local ancestry inference in stable in core region
# to make that more versatile, I will not add them now
# will add that in tuning and screening scripts

chr.segments <- as.list(1:original.segment.count)
# initialize new segments

for(iter in 1:original.segment.count){
  # find grid of segment boundaries that lie at observed loci
  seg.borders <- core.loci.pos[[iter]][floor(seq(1,core.lengths[[iter]],len=num.segments[[iter]]+1))]

  chr <- burnin.table[iter,1]
  chr.segments[[iter]] <- data.frame("chr" = chr, 
                                     "core.start.pos" = seg.borders[1], "core.end.pos" = seg.borders[2],
                                     "core.loci" = NA
  )
  if(length(seg.borders)==2){
    # if there's only one small segment in this original segment
    core.start.pos <- seg.borders[1]
    core.end.pos <- seg.borders[2]

    
    core.loci <- sum(loci.pos[[iter]] <= core.end.pos & loci.pos[[iter]] >= core.start.pos)
    chr.segments[[iter]]$core.loci <- core.loci

    next
  }else{
    chr.segments[[iter]] <- chr.segments[[iter]][-1,]
    # if this chr have multiple regions, delete it to recalculate.
    }
  j <- 2 
  while(j <= length(seg.borders)){
    # till the last segment border
    core.start.pos <- seg.borders[j-1]
    core.end.pos <- seg.borders[j]
  
    core.loci <- sum(loci.pos[[iter]] <= core.end.pos & loci.pos[[iter]] >= core.start.pos)

    chr.segments[[iter]] <- rbind(chr.segments[[iter]],
                                  data.frame("chr" = chr, 
                                             "core.start.pos" = core.start.pos, "core.end.pos" = core.end.pos,
                                             "core.loci" = core.loci
                                  )
    )
    j <- j+1
  }
  print(iter)
}


chr.segments.table <- purrr::reduce(chr.segments,rbind)
# combine all segments into one dataframe

write.table(chr.segments.table, file = paste0(grand.data.dir,"/chr_segments_",total.segments,"segs_new_pass.txt"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)
