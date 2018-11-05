import os
import pickle
#import pysam
#import vcf

from SomaticHaplotype import *

################################################################################
# functions
################################################################################
#(functions here)


################################################################################
# main
################################################################################   

def main(args):

  # path to input pickle file
  pickle_path = args.ps1
  phase_set_dict = pickle.load(pickle_path)
  variants_dict = pickle.load(pickle_path)

  # parse the genomic range argument
  if args.range is None:
    chrom = None
    start = None
    end = None
  else:  
    chrom = args.range.split(":")[0]
    try:
      start = int(args.range.split(":")[1].split("-")[0])
    except:
      start = None
    try:
      end = int(args.range.split(":")[1].split("-")[1])
    except:
      end = None

if __name__ == '__main__':
  main(args)
