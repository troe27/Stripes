"""
Script to find all sites fixed between ancestors
that contribute to each individual in a specific generation.

Takes:
    - Pedigree
    - VCF of F0
    - Generation of interest ( e.g. 8)
    - Ancestor state. (here: high, low)
Returns:
    - .VCF file annotated with all sites fixed between unique sets of ancestors.
    - .CSV file with ID, annotating which set of ancestors an individual has.

Author: tilman.ronneburg@imbim.uu.se
"""

import pandas as pd
from collections import Counter
from cyvcf2 import VCF
import numpy as np
import argparse
from argparse import RawTextHelpFormatter


def cli_parser():
    """ Parse command line input."""
    parser_main = argparse.ArgumentParser(prog='find_fixed_TGR.py',
                                          description='''
    Takes:
        - Pedigree
        - VCF of F0
        - Generation of interest ( e.g. 8)
        - Ancestor state. (here: high, low)
    Returns:
        - .VCF file annotated with all sites fixed between unique sets of ancestors.
        - .TSV file with ID, annotating which set of ancestors an individual has.
    ''',formatter_class=RawTextHelpFormatter)
    
    parser_main.add_argument("-i","--infile",
                             help="/path/to/all.vcf",
                             required = True)
    parser_main.add_argument("-o", "--outfile",
                             help="/path/to/output.vcf",
                             required = True)
    parser_main.add_argument("-p", "--pedigree",
                             help="path/to/pedigree.csv",
                             required = True)
    parser_main.add_argument("-t", "--threads",
                             help="number of threads to use for vcf parsing.",
                             default=1)
    parser_main.add_argument("--high",
                             help="File detailing which F0 are from high-line.",
                             required = True)
    parser_main.add_argument("--low",
                             help="File detailing which F0 are from low-line.",
                             required = True)
    parser_main.add_argument("-g","--generation",
                             help="generation of interest, numeric.",
                             required = True)

    args = parser_main.parse_args()
    return args

def load_pedigree_15(pedigree_file):
    """ Load a pedigree file and return as a pandas dataframe. """
    data = pd.read_csv(pedigree_file, sep="\t", dtype={"ID": str, "Sire":str, "Dam":str, "Sex": str})
    pedigree = dict()
    for index, row in data.iterrows():
        #print(row["ID"])
        if str(row["Sex"]) == "1":
            sex = "female"
        elif str(row["Sex"]) == "2":
            sex = "male"
        pedigree[row["ID"]] = {"gen":int(str(row["ID"][-2:])),
                           "sire":str(row["Sire"]),
                           "dam":str(row["Dam"]),
                           "sex":sex}
    return pedigree

def get_fam(pedigree, ID, target="00"):
    """
    Return a list(set()) of all related individuals from target
    for an individual [ID]
    """
    ancestors = [pedigree[ID]["sire"],pedigree[ID]["dam"]]
    all_ancestors = [pedigree[ID]["sire"],pedigree[ID]["dam"]]
    for i in reversed(range(int(pedigree[ID]["gen"])-1)):
        new_ancestors = []
        while ancestors:
            ancestor = ancestors.pop(0)
            new_ancestors.append(pedigree[ancestor]["sire"])
            new_ancestors.append(pedigree[ancestor]["dam"])
        ancestors = new_ancestors
        if i == int(target):
            break
    return list(set(ancestors))
    
def get_sets_for_generation(pedigree, generation, target=0):
    """
    Take pedigree dict and generation.
    Returns dictionary with tuple of founders as key, and all IDs of generation that have these as ancestors.
    """
    gen = [key for key, item in pedigree.items() if item["gen"]==int(generation)]
    founders = {}
    for ID in gen:
        fam_ancestors = tuple(get_fam(pedigree=pedigree, ID=ID, target=int(target)))
        founders.setdefault(fam_ancestors, [])
        founders[fam_ancestors].append(ID)
    return founders

def load_pop(pfile):
    """ Load a one individual per line file as keys of a dictionary. """
    with open(pfile, "rt") as handle:
        pdict =  {i.split("_")[0].rsplit("-",1)[1]+"00":i for i in handle.read().rsplit("\n") if not i ==""} # specific to our naming-convention
    return pdict

def split_founders_high_low(founders, high_file, low_file):
    """
    Take founder dict, file with high-line individuals, low-line individuals.
    return two dicts:
     - one linking sample ID to setname
     - two linking setname to sets of ancestors, split by high-low (e.g. d["setname"]["high"])
    """
    high = load_pop(high_file)
    low = load_pop(low_file)
    id_to_set = dict()
    sets = dict()
    for i, (se, IDs) in enumerate(founders.items()):
        sets["set_"+str(i)]=dict()
        sets["set_"+str(i)]["high"]=[]
        sets["set_"+str(i)]["low"]=[]
        for founder in list(se):
            if str(founder) in high:
                sets["set_"+str(i)]["high"].append(str(founder))
            elif str(founder) in low:
                sets["set_"+str(i)]["low"].append(str(founder))
            else:
                print("WARNING: Founder {} was neither found in population {} nor in population {}. skipping...".format(founder, high_file, low_file))
        for ID in IDs:
            id_to_set[ID]="set_"+str(i)
    return id_to_set, sets


def find_fixed_sites(vcf_file, sets, outfile, threads=1, full=True, tsv=False):
    """
    Find sites fixed between the High- and Low-line founders of all sets.
    writes an annotated VCF file to [outfile].
    depending on wether full==True, writes the whole vcf, or only sites fixed between any of the sets.
    if tsv == True, writes fixed sites as tab separated table.
    """
    popcomb = {}
    vcf = VCF(vcf_file, threads=int(threads), gts012=True)
    ##### high - low specific, generate popcomb.
    ##### swap out if generalising again
    for setid, setitem in sets.items():
        pop1 = dict.fromkeys(setitem["high"])
        p1_indices = [i for i,sample in enumerate([i.split("_")[0].rsplit("-",1)[1]+"00" for i in vcf.samples]) if sample in pop1] #sample naming specific
        pop2 = dict.fromkeys(setitem["low"])
        p2_indices = [i for i,sample in enumerate([i.split("_")[0].rsplit("-",1)[1]+"00" for i in vcf.samples]) if sample in pop2] #sample naming specific
        popcomb[setid]=[p1_indices, p2_indices]
    #####
    vcf.add_info_to_header({'ID': 'fixed',
                            'Description': 'fixed_between_pops',
                            'Type': 'Character',
                            'Number': '1'
                            })
    if tsv == False:
        w = Writer(outfile, vcf)
    else:
        w = open(outfile, "wt")
    for v in vcf:
        fixed = []
        for key, (p1_indices, p2_indices) in popcomb.items():
            pop1_state = "Empty"
            pop2_state = "Empty"
            if 0 in set(np.take(v.gt_types, p1_indices)) and not 1 in set(np.take(v.gt_types, p1_indices)) and not 2 in set(np.take(v.gt_types,p1_indices)):
                pop1_state = "0"
            elif 2 in set(np.take(v.gt_types, p1_indices)) and not 1 in set(np.take(v.gt_types, p1_indices)) and not 0 in set(np.take(v.gt_types,p1_indices)):
                pop1_state = "2"
            if 0 in set(np.take(v.gt_types, p2_indices)) and not 1 in set(np.take(v.gt_types, p2_indices)) and not 2 in set(np.take(v.gt_types,p2_indices)):
                pop2_state = "0"
            elif 2 in set(np.take(v.gt_types, p2_indices)) and not 1 in set(np.take(v.gt_types, p2_indices)) and not 0 in set(np.take(v.gt_types,p2_indices)):
                pop2_state = "2"
            if not pop1_state == "Empty" and not pop2_state == "Empty":
                if not pop1_state == pop2_state:
                    if pop1_state == "2":
                        fixed.append(str(key)+"=HIGH") # if ALT is fixed in Pop1 #FIXME terminology is not generic enough.
                    else:
                        fixed.append(str(key)+"=LOW") # if ALT is fixed in Pop2
        if tsv == False:
            if fixed:
                v.INFO["fixed"] = ";".join(fixed)
            if full == True:
                w.write_record(v)
            else:
                if fixed:
                    w.write_record(v)
        else:
            if fixed:
                l = [str(v.CHROM), str(v.POS)]
                l.append(";".join(fixed))
                w.write("\t".join(l)+"\n")
                # writes table such as:
                # CHROM   POS    set1=HIGH;set2=HIGH;set3=LOW
                # chr2    2345   set23=HIGH;set41=HIGH
    w.close();vcf.close()


def main():
    """ Execute workflow. """
    args = cli_parser()
    pedigree = load_pedigree_15(pedigree_file=args.pedigree) #studied = generation ?
    founders = get_sets_for_generation(pedigree=pedigree,generation=args.generation)
    id_to_set, sets = split_founders_high_low(founders=founders,
                                              high_file=args.high,
                                              low_file=args.low)
    # write id_to_set to file:
    id_to_set_out = args.outfile.rsplit("/",1)[0]+"/id_to_set_F{}.csv".format(args.generation)
    with open(id_to_set_out, "w") as handle:
        for key, item in id_to_set.items():
            handle.write(key+";"+item+"\n")
    # write annotated vcf/tsv to file:
    find_fixed_sites(vcf_file=args.infile,
                     sets=sets,
                     outfile=args.outfile,
                     threads=args.threads,
                     full=False, tsv=True) # exclude all never-fixed sites & make tsv output instead of vcf


if __name__ == "__main__":
    main()
