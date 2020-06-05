#!/bin/bash -l
#SBATCH -A snic2017-7-53
#SBATCH -t 8:00:00 
#SBATCH -p core -n 2
#SBATCH -J mektrex
#SBATCH -o mektrex_%j.out
#SBATCH -e mektrex_%j.error
#SBATCH --mail-user tilman.ronneburg@imbim.uu.se
#SBATCH --mail-type=FAIL
#SBATCH --get-user-env
##command underneath this##
source activate v3


module load R

main_path="/home/tilman/storage3/process_stripes/"
gen_folder="stripes_optimise_F7/"


chr_ranks="/home/tilman/storage3/stripes_optimise/aux/Index_fastq_contig_chr.txt"

gt_folder="with.fam.f2.call2.Genotype/"
tg_folder="with.fam.f2.call2.TIGER_OUT/"
binsize=1000000
cutoff=10 
max_num=30
windowsize=200

output=${main_path}${gen_folder}binsize${binsize}_cutoff${cutoff}_wsize${windowsize}.csv

Rscript --vanilla scripts/format_zy.R  ${main_path}${gen_folder}${gt_folder} \
                                       ${main_path}${gen_folder}${tg_folder} \
                                       ${chr_ranks}\
                                       ${binsize} \
                                       ${cutoff} \
                                       ${max_num} \
                                       ${windowsize} \
                                       ${output}



windowsize=50
output=${main_path}${gen_folder}binsize${binsize}_cutoff${cutoff}_wsize${windowsize}.csv

Rscript --vanilla scripts/format_zy.R  ${main_path}${gen_folder}${gt_folder} \
                                       ${main_path}${gen_folder}${tg_folder} \
                                       ${chr_ranks}\
                                       ${binsize} \
                                       ${cutoff} \
                                       ${max_num} \
                                       ${windowsize} \
                                       ${output}
