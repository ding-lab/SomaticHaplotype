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

  # parse the genomic range argument
  chrom = args.range.split(":")[0]
  try:
    start = int(args.range.split(":")[1].split("-")[0])
  except:
    start = None
  try:
    end = int(args.range.split(":")[1].split("-")[1])
  except:
    end = None

  # path to input pickle file
  pickle_path = args.ps1
  phase_set_dict = pickle.load(pickle_path)
  variants_dict = pickle.load(pickle_path)

if __name__ == '__main__':
  main(args)
