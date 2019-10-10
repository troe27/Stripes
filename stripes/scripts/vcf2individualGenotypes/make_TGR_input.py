"""
Script to take the files derived from the founder-populations
about which sites are fixed between which sets of ancestors, and in combination
with the variant data of the individual creates the input for TIGER.

Author: tilman.ronneburg@imbim.uu.se
"""
from cyvcf2 import VCF
import argparse
from argparse import RawTextHelpFormatter

def cli_parser():
    """ Parse command line input."""
    parser_main = argparse.ArgumentParser(prog='make_TGR_input.py',
                                          description='''
    Takes:
        - sample_ID
        - vcf of offspring generation / samples of interest
        - file that details all sites fixed between populations.
        - id_to_set file that details which sample has which ancestors.

    Returns:
        - Input for TIGER.
          Format is:
         [CHROM    POS    HIGH    OFHIGH    LOW    OFLOW]
    e.g.  Chr33    44621  REF     0         ALT    1
    e.g.  Chr2     3217   ALT     0         REF    1
          ...
    ''',
                                         formatter_class=RawTextHelpFormatter
                                         )
    parser_main.add_argument("-i","--infile",
                             help="/path/to/vcf",
                             required = True)
    parser_main.add_argument("-o", "--outfile",
                             help="/path/to/output.tsv",
                             required = True)
    parser_main.add_argument("-t", "--threads",
                             help="number of threads to use for vcf parsing.",
                             default=1)
    parser_main.add_argument("-a","--annot",
                             help="table detailing fixed sites",
                             required = True)
    parser_main.add_argument("-s","--set_to_id",
                             help="File detailing which IDs have which ancestors",
                             required = True)
    parser_main.add_argument("-d","--sample_id",
                                 help="name of sample for which the output shall be generated.",
                                 required = True)

    args = parser_main.parse_args()
    return args

def get_id_set(sample_ID, set_to_id_file):
    """ Extract which Set of ancestors belong to the sample_ID. """
    set_ID = None
    with open(set_to_id_file, "rt") as handle:
        for line in handle.readlines():
            if line.split(";")[0] ==sample_ID:
                set_ID = line.split(";")[1].rstrip()
                break
    if set_ID == None:
        raise KeyError("the Sample_ID {} does not seem to be in the ID_to_set_file {}".format(sample_ID, set_to_id_file))
    return set_ID

def get_annotation_dict(annot_file, set_ID):
    """Parse the POS-CHR-SET file."""
    adict = dict()
    with open(annot_file, "rt") as handle:
        for line in handle.readlines():
            l = line.split("\t")
            fdict = {i.split("=")[0]:i.split("=")[1] for i in l[2].split(";")}
            if set_ID in fdict:
                adict[str(l[0])+"_"+str(l[1])]=fdict
    return adict

def make_tiger_input(sample_ID, vcf_file, adict, outfile, set_id, threads=1):
    """
     - Take the dictionary with annotations, sample and set IDs
     - Create the genotype input for TIGER.
    """
    vcf = VCF(vcf_file, threads=int(threads) , gts012=True)
    index = vcf.samples.index(sample_ID)
    with open(outfile, "wt") as handle:
        for v in vcf:
            if str(v.CHROM)+"_"+str(v.POS) in adict:
                gt = v.gt_types[index]
                if gt == 3:
                    continue
                # create line:
                l = [str(v.CHROM), str(v.POS), "HIGH", "OFHIGH", "LOW", "OFLOW"]
                if adict[str(v.CHROM)+"_"+str(v.POS)][set_id] == "HIGH":
                    l[2] = "ALT"
                    l[4] = "REF"
                    if gt == 2:
                        l[3] = "1"
                        l[5] = "0"
                    elif gt == 0:
                        l[3] = "0"
                        l[5] = "1"
                else:
                    l[2] = "REF"
                    l[4] = "ALT"
                    if gt == 2:
                        l[3] = "0"
                        l[5] = "1"
                    elif gt == 0:
                        l[3] = "1"
                        l[5] = "0"
                handle.write("\t".join(l)+"\n")

def main():
    """ Execute workflow. """
    args = cli_parser()
    set_ID = get_id_set(sample_ID=args.sample_id,
                        set_to_id_file=args.set_to_id)
    adict = get_annotation_dict(annot_file=args.annot,
                                set_ID=set_ID)
    make_tiger_input(sample_ID=args.sample_id,
                     vcf_file=args.infile,
                     adict=adict,
                     outfile=args.outfile,
                     set_id=set_ID,
                     threads=args.threads)

if __name__ == "__main__":
    main()
