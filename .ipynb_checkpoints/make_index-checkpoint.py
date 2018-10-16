"""Create index file for the tiger pipeline"""
import subprocess
import pandas as pd
import os

def make_index_file(fai, index_file):
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
    df.to_csv(index_file,header=False,index=True, sep="\t")
    columnsTitles=[0,"rank",1]
    df=df.reindex(columns=columnsTitles)
    return None
    
def check_and_create_fai(ref):
    """ Check if .fai file exists, else create it."""
    if not ref.rsplit("/",1)[1]+".fai" in os.listdir(ref.rsplit.("/",1)[0]):
        subprocess.call(["samtools", "faidx", ref])
    return None
    
def main():
    ref=snakemake.input
    index_file=snakemake.output
    fai=ref+".fai"
    
    check_and_create_fai(ref=ref)
    make_index_file(fai=fai, index_file=index_file)
    return None

if __name__ == "__main__":
    main()

