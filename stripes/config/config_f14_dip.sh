# The reference genome, used for the mapping
reference:
  /home/tilman/nas/ref_gg6a/ref_unzipped/GCF_000002315.5_GRCg6a_genomic.fna

# If we want to have an output that contains only the sites selected between 2 lines #TODO needed? FIXME
selectBetweenLines:
  /home/tilman/storage2/stripes_2_data/Founder_fixed_sites_chr_sites.txt

# Folder containing the fastq (considered as paired and named {sample}_forward.fastq.gz and {sample}_reverse.fastq.gz)
Raw_Folder:
  /home/tilman/nas/stripes_AIL_gg6a/F14/stripes_input/raw

# file containing all sample wildcards, in case no individual input files are there anymore
# TODO have this derived from the VCF samples ?
sample_file:
  /home/tilman/nas/stripes_AIL_gg6a/F14/stripes_input/new_ids_F14.txt


conda:
  ./config/env.yaml

python:
  /home/tilman/miniconda3/envs/v3/bin/python

# The VCF file of founders (created separatly with a more classic approach, in our case GATK and a joined call)
founders:
  /home/tilman/nas/all_GQ20_MQ50_noIndels_callrate0.3.vcf

TIGER_WD:
  /home/tilman/nas/stripes_AIL_gg6a/F14/stripes_data

CALL_WD:
  VCF

cohort:
  /home/tilman/nas/stripes_AIL_gg6a/F14/stripes_input/all_F14_only_founder_snps_NoDuplicates_NewNames.vcf


# Was used to test the difference of quality in pipeline using markers fixed in family and in line
indiv:
  fam:
    with.fam.f2.call2.
#  lines: with.lines.f2.call2.

# To create a folder of the VCFs
CALLING_SUFFIX:
  Call

# To create a folder for individual genotypes per sample
GENOTYPE_SUFFIX:
  Genotype


# The results of TIGER, one subfolder per sample, splited by chromosomes
TIGER_SUFFIX:
  TIGER_OUT


# A file with the chromosomes and their size in base (see example)
genomeSizeFile:
  /home/tilman/nas/stripes_AIL_gg6a/F14/stripes_data/data/size.genome

# java options for rule caller and rule estimate:
java_options:
  -Xmx4G

# Number of SNP per windows
winsize:
  50

# Number of scaffolds to treat , the methods is not working well for smaller scaffolds
scaffoldNumber:
  24

# The pedigree file (see example for format)
pedigree:
  /home/tilman/nas/pedigree/AIL_pedigree_20190813.tsv
  
  #/home/tilman/nas/stripes_AIL_gg6a/F14/stripes_data/input_aux/combined_pedigree_20181217.tsv

# File containing the full names (from the reference, a rank (col 3) and the size of the chromosomes
contigNames:
  /home/tilman/nas/stripes_AIL_gg6a/F14/stripes_data/input_aux/Index_fastq_contig_chr.txt


# If not in the script folder use an absolute path
individualGenotypes2breaksPath:
  scripts/individualGenotypes2breaks

high:
  /home/tilman/nas/stripes_AIL_gg6a/F14/stripes_data/input_aux/high_full.txt
low:
  /home/tilman/nas/stripes_AIL_gg6a/F14/stripes_data/input_aux/low_full.txt
generation:
  14
