#!/usr/bin/env Rscript
#Written by Yanjun, very slightly modified by Tilman.
#make genotype matrix and filter for all bins that have fewer than 10 markers.
## this script contains quite a few sample specific regular expressions & will likely break at some point.
## Heres to hoping that i will replace it with something generic before it does.
source("./scripts/functions.R")
library("data.table")
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop(" Arguments must be supplied (input file).n", call.=FALSE)
  }

## input:
input_folder <- args[1]
tiger_folder <- args[2]
chr_ranks_file <- args[3]

## parameters:
 binsize <- as.numeric(args[4]) #size of bin in bp
 cutoff <- as.numeric(args[5]) # min number of markers per bin

# output:
outfile <- args[6]

all_genotypes <- setdiff(list.files(path = input_folder),list.dirs(path = input_folder,full.names = F)) #  extract the name of output
all_id <- gsub(pattern = "(\\d+)\\.genotype$",replacement = "\\1",x = all_genotypes) # extract the ID

# determine number of bins:
chr_ranks_table <- fread(chr_ranks_file)
setnames(chr_ranks_table, old=c("chr_name","chr_size","chr_sizerank"), new=c("name","size","rank" ))

chr_ranks_table <- chr_ranks_table[order(chr_ranks_table$rank)]
chr_ranks_table <- data.frame(chr_ranks_table)
chroms <- chr_ranks_table$name
chroms_len <- chr_ranks_table$size
info_obj <- get.info(chroms=chroms,chroms.len=chroms_len,bin.size=binsize)
num_bins <- data.frame(info_obj["num.bin"])
#print(num_bins)
num_bins <- sum(num_bins$num.bin)

out.put <- data.frame(array(NA,dim=c(length(all_genotypes),num_bins))) # generate a holder ##
rownames(out.put) <- all_id
total <- numeric(length(all_id))


for( i in 1:nrow(out.put)){
    input_now <- paste0(input_folder,all_genotypes[i])
    total[i] <- nrow(fread(input_now))
    num_now<- get_density_input(inputfile = input_now,binsize = binsize,cutoff = cutoff, chr.ranks = chr_ranks_table)
    out.put[i,] <- wrap_get_density(test =num_now,chr.ranks = chr_ranks_table, binsize=binsize )
    #cat(i, "\n")
}

density <- apply(out.put,1,mean)
id.keep <- names(density)[density >5]
length(density)-length(id.keep)

################# get Tiger output ########################
all_vcf <- list.dirs(path = tiger_folder)
reg_expr <- paste0(tiger_folder,"/", "(\\d+)\\..*") ##gotcha
id_all <- gsub(pattern = reg_expr,replacement = "\\1",x = all_vcf) #
index.keep <- which(id_all %in% id.keep)
all <- list.files(all_vcf,pattern = "\\d+\\.genotype\\.(\\d+)\\.rough_COs\\.refined\\.breaks.txt")
chr <- sort(as.numeric(unique(gsub(pattern = "\\d+\\.genotype\\.(\\d+)\\.rough_COs\\.refined\\.breaks.txt",replacement = "\\1",x = all))))
out_3 <- Extract_all(chromosome = chr,id_all = id_all[index.keep],all_vcf = all_vcf[index.keep],gap=3e6,filter = T )
#print(colnames(out_3))
chr.match <- chr_ranks_table
length.chr <- ceiling(chr.match$size/binsize)
#print(length.chr)
names(length.chr) <- chr.match$rank
mat_3 <- Mat2geno(out = out_3,id = id_all[index.keep],length.chr = length.chr, binsize=binsize)
#print(colnames(mat_3))
write.csv(mat_3, outfile)
