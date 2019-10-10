"""
Format Stripes output for r/qtl
Author:Tilman.Ronneburg@imbim.uu.se
20190718

"""
import pandas as pd
import numpy as np
from datetime import datetime

def make_rqtl_input(df, phe_out, gen_out, phenotype,generation="02", old_id_pheno=True, jitter=0.2, phe_id_col="ID", phe_trait_col="BW8"):
    """ makes input files for r/qtl  """
    if isinstance(phenotype, str) ==True:
         phenotypes_raw = pd.read_csv(phenotype, sep=",")
    else:
        phenotypes_raw = phenotype
    if old_id_pheno == True:
        # some files still follow the old naming convention,
        # which can lead to confusion since the IDs are not unique across generations.
        # we therefore add the generation to the end of the ID.
        # make sure that you have only individuals from one generation in here.
        phenotypes_raw[phe_id_col] = [str(i)+str(generation) for i in np.array(phenotypes_raw[phe_id_col])]
    phenotypes_raw[phe_id_col] = phenotypes_raw[phe_id_col].astype(str)
    df.index = df.index.astype(str)
    # r/qtl requires intersection between ids, apparently.
    set_gen = set(list(df.index))
    set_phe = set(list(phenotypes_raw[phe_id_col]))
    intersect_id = set.intersection(set_gen, set_phe)
    phe_intersect = phenotypes_raw[phenotypes_raw[phe_id_col].isin(list(intersect_id))]
    gen_intersect = df.loc[list(intersect_id)]
    ## phenotypes are done, lets export them
    phe_intersect.to_csv(phe_out, index=False, sep=",",)
    make_rqtl_geno_input(df=gen_intersect, outfile=gen_out, jitter=jitter)


def make_rqtl_geno_input(df, outfile, jitter=0.2):
    """ make a genotype input_file that r/qtl can understand."""
    df = df.fillna("NA")  # fill NaNs with "NA"

    def remove_amb(x):
        """ remove ambiguous calls and format to 1/0/-1/NA only """
        if x=="NA":
            return x
        elif x >= 1.0-jitter:
            return int(1)
        elif 1.0-jitter>= x >0.0+jitter:
            return "NA"
        elif 0.0+jitter >= x >0.0-jitter:
            return int(0)
        elif 0.0-jitter >= x > -1.0+jitter:
            return "NA"
        elif x<= -1.0+jitter:
            return int(-1)

    df = df.applymap(remove_amb)
    df.to_csv(".make_input_tmp.csv", sep=",")
    t = open(".make_input_tmp.csv").read().split("\n")
    first_line = t.pop(0).split(",")
    first_line[0]="id"
    chrom = []
    for i in first_line[1:]:
        chrom.append(i.split("-")[0])
    chrom = ","+",".join(chrom)
    first_line = ",".join(first_line)
    with open(outfile, "w") as handle:
        handle.write(first_line+"\n")
        handle.write(chrom+"\n")
        handle.write("\n".join(t))
