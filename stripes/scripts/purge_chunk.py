import os
import argparse
from argparse import RawTextHelpFormatter

def cli_parser():
    '''Parse command line input.'''
    parser_main = argparse.ArgumentParser(prog='purge_chunk.py',
                                          description='''WHAT THIS SCRIPT DOES:
    ######################################################
    remove all files associated with a list of samples from output.
    optionally, remove snakemake metadata
    #####################################################''',
                                        formatter_class=RawTextHelpFormatter
                                         )
    parser_main.add_argument("-c","--chunkfile",
                             help="path/to/chunk_file.txt",
                             required = True)
    parser_main.add_argument("-d", "--dryrun",
                             help="dont purge anything, just list files.",
                             action = "store_true")
    parser_main.add_argument("-m", "--remove_metadata",
                             help="remove all snakemake metadata, including conda environments",
                             action = "store_true")

    args = parser_main.parse_args()
    return args

input_folder = "/home/tilman/storage2/stripes_2_data/data/"
gt_prefix = "with.fam.f2.call2.Genotype/"
tg_prefix = "with.fam.f2.call2.TIGER_OUT/"
chunkfile = "config/samples_new_chunk2.txt"



def get_files(chunkfile, input_folder,gt_prefix, tg_prefix):
    """Find all files associated with samples in chunkfile"""
    samples_to_purge = [line.rstrip() for line in open(chunkfile,"r").readlines()] 
    gt_folder = os.path.join(input_folder, gt_prefix)
    tg_folder = os.path.join(input_folder,tg_prefix)
    samples_to_purge = [line.rstrip() for line in open(chunkfile,"r").readlines()] 
    gt_files = [f for f in os.listdir(gt_folder) if f.split(".")[0] in samples_to_purge ]
    gt_files = [os.path.join(gt_folder,f) for f in gt_files]
    tg_files = [f for f in os.listdir(tg_folder) if f.split(".")[0] in samples_to_purge ]
    tg_files = [os.path.join(tg_folder,f) for f in tg_files]
    all_files = gt_files + tg_files
    return all_files

def purge_files(file_list):
    """Remove all files and folders in file_list."""
    for file in file_list:
        os.system("rm -r {file}".format(file=file))

def list_files(file_list):
    """Print all files and folders in file_list."""
    print("DRYRUN:the following files and folders with all their content will be removed:")
    for f in file_list:
        print(f)


def main():
    """Execute workflow."""
    args = cli_parser()
    ## Hardcode input_folder & prefixes:
    input_folder = "/home/tilman/storage2/stripes_2_data/data/"
    gt_prefix = "with.fam.f2.call2.Genotype/"
    tg_prefix = "with.fam.f2.call2.TIGER_OUT/"
    file_list = get_files(chunkfile=args.chunkfile,
                          input_folder=input_folder,
                          gt_prefix=gt_prefix,
                          tg_prefix=tg_prefix)
    if args.remove_metadata == True:
        file_list.append("./.snakemake")
    if args.dryrun == True:
        list_files(file_list=file_list)
    else:
        #print("this is where the remove function would come in after debugging")
        purge_files(file_list=file_list)


if __name__ == "__main__":
    main()
