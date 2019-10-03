#!/bin/bash -l
#SBATCH -A snic2018-3-170
#SBATCH -t 10:0:0 
#SBATCH -p core -n 15
#SBATCH -J tiger_snek
#SBATCH -o tiger_snek_chunk3_%j.out
#SBATCH -e tiger_snek_chunk3_%j.error
#SBATCH --mail-user tilman.ronneburg@imbim.uu.se
#SBATCH --mail-type=ALL
#SBATCH --get-user-env
##command underneath this##
module load java/sun_jdk1.8.0_151
source activate v3
snakemake --unlock -s ./tiger.snek --configfile config/config_rackham_all_samples.yaml  --rerun-incomplete 
snakemake --keep-going -j 13 -s ./tiger.snek  --configfile config/config_rackham_all_samples.yaml  --rerun-incomplete --use-conda


