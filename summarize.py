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
  bam_phase_set_dictionary = pickle.load(pickle_path)
  vcf_variants_dictionary = pickle.load(pickle_path)

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

  # Lan uses functions here

  # write output and close files
  os.makedirs(args.output_directory, exist_ok = True)
  output_file_path = os.path.join(args.output_directory, args.output_prefix + "YOUR SUFFIX HERE.tsv")
  output_file = open(output_file_path, 'w')
  # write results to output file here
  output_file.close()

if __name__ == '__main__':
  main(args)
