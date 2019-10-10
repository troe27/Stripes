"""Create index file for the tiger pipeline"""
import subprocess
import pandas as pd
import os

def make_index_file(fai, index_file, genome_size_file):
    """
    Strip the first two columns ( name, length) from .fai
    add ranking by size.
    create new file with these three columns ( name, length, )
    """
    with open(fai, "r") as handle:
        l = [[str(i.rstrip().split()[0]), int(i.rstrip().split()[1])] for i in handle.readlines()]
    df = pd.DataFrame(l)
    df['rank'] = df[1].rank(method='first', ascending=False)
    df["rank"] = df["rank"].astype(int)
    df.columns = ["chr_name", "chr_size", "chr_sizerank"]
    df.to_csv(index_file,header=True,index=False, sep="\t")
    df_subset = df[["chr_sizerank", "chr_size"]]
    df_subset.to_csv(genome_size_file, header=False, index=False, sep="\t")
    return None

def check_and_create_fai(ref):
    """ Check if .fai file exists, else create it."""
    if not ref.rsplit("/",1)[1]+".fai" in os.listdir(ref.rsplit("/",1)[0]):
        subprocess.call(["samtools", "faidx", ref])
    return None

def main():
    ref=snakemake.input[0]
    #ref="/home/tilman/storage/ref_gg6a/ref_unzipped/GCF_000002315.5_GRCg6a_genomic.fna"
    index_file=snakemake.output[0]
    #index_file = "test_index.tsv"
    #genome_size_file="test_size_file.tsv"
    genome_size_file=snakemake.output[1]
    fai=ref+".fai"

    check_and_create_fai(ref=ref)
    make_index_file(fai=fai, index_file=index_file ,genome_size_file=genome_size_file)
    return None

if __name__ == "__main__":
    main()
