#! ~/miniconda3/bin/Rscript
# 
# dedup by R2
# @author zliu
#
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("usage: Rscript CHARMdedup.R 10 30 cellname R2.bam.bed R2.dedup.bed", call.=FALSE)
}

library(dplyr)
library(readr)
library(valr)

#config
dedup_dist = as.numeric(args[1])
min_mapq = as.numeric(args[2])
cellname = args[3]
path_in = args[4]
path_out = args[5]

R2bed <- read_tsv(path_in,col_names = F)

# when there is no such tag
if(dim(R2bed)[1]==0){
    R2bed %>% write_tsv(path_out,col_names = F)
} else {
    bed5 <- read_tsv(path_in,col_names = F)
    names(bed5) <- c("chrom","start","end","readid","mapq","strand")
    # dedup by 5 prime end
    bed5_processing <- bed5 %>% arrange(chrom,start,end) %>% filter(mapq > min_mapq,chrom %in% paste0("chr",c(seq(1,19),"X","Y"))) %>%
        mutate(start = ifelse(strand == "+",start,end),end = start + 1) %>% bed_cluster(max_dist = dedup_dist)
    bed5_processing %>% group_by(.id) %>% arrange(desc(mapq)) %>% slice(1) -> bed5_processing

    # extention
    bed5_processing %>% ungroup() %>% select(1:3) %>% mutate(start = start - 75,end=end+75,cellname = cellname,count = 1) %>% write_tsv(path_out,col_names = FALSE) 
}

