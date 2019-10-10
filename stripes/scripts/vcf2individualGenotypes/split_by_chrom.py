"""script to
    a) split the genotype files per chromosome, and
    b) rename them to a numerical nomenclature. (e.g. NT_XFYFFA123a > 1  ; NT_XFYFFX233b > 2  ,etc.)
    Author: tilman.ronneburg@imbim.uu.se
    """
# imports
import argparse
from argparse import RawTextHelpFormatter

def cli_parser():
    """ Parse command line input."""
    parser_main = argparse.ArgumentParser(prog='split_by_chrom.py',
                                          description='''
    does the following:
    a ) retain the n largest scaffolds/chromosomes
    b ) Split the genotype files per chromosome
    c ) rename the chromosomes numerically, based on sizerank
    ''',
                                         formatter_class=RawTextHelpFormatter
                                         )
    parser_main.add_argument("-i","--infile",
                             help="/path/to/sample.genotype.sorted",
                             required = True)
    parser_main.add_argument("-o", "--outfolder",
                             help="/path/to/output_folder",
                             required = True)
    parser_main.add_argument("-s","--sizerankingfile",
                             help="File detailing the site and size-rank of each chromosomes.",
                             required = True)
    parser_main.add_argument("-m","--maxchr",
                             help="number of scaffolds/chromosomes to process.",
                             required = True)
    parser_main.add_argument("-d","--sampleid",
                                 help="name of sample for which the output shall be generated.",
                                 required = True)

    args = parser_main.parse_args()
    return args

def load_name_conversion(size_genome_file, max_chr):
    """Load the Genome size file and give a dictionary of name:number. retain only [max_chr] number of chromosomes. """
    chr_dict = {}
    with open(size_genome_file,"r") as handle:
        ll = [line.rstrip().split("\t") for line in handle.readlines()]
    for line in ll:
        if not line[0].startswith("chr"):
            if int(line[2]) <= int(max_chr):
                chr_dict[line[0]]=line[2]
    return chr_dict

def iter_over_file(genotype_name, chr_dict, gt_file, outfolder):
    """ iter over genotype file and split it into different chromosomes."""
    try:
        #open all relevant output files:
        files = {}
        for key, item in chr_dict.items():
            files[key] = open("{outfolder}/{genotype}.chr.{chr_nr}".format(outfolder=outfolder, genotype=genotype_name, chr_nr=item),"w")
        # iter over genotype file:
        with open(gt_file, "r") as handle:
            for line in handle.readlines():
                l = line.split()
                if l[0] in chr_dict:
                    fkey=l[0]
                    new_line = chr_dict[l[0]]+"\t"+"\t".join(l[1:])
                    files[fkey].write(new_line+"\n")
    finally:
        # close all files
        for key,item in files.items():
            item.close()
    return None

def main():
    """Execute workflow."""
    args = cli_parser()
    chr_dict = load_name_conversion(size_genome_file=args.sizerankingfile, max_chr=int(args.maxchr))
    iter_over_file(genotype_name=args.sampleid, chr_dict=chr_dict, gt_file=args.infile, outfolder=args.outfolder)

if __name__ == "__main__":
    main()
