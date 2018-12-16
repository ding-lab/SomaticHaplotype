import os
import pickle
#import pysam
#import vcf

from SomaticHaplotype import *

################################################################################
# functions
################################################################################

def extend_phase_set_dictionary(ps_dict1, ps_dict2):
  return(extended_ps_dict1)


################################################################################
# main
################################################################################   

def main(args):

  # path to input pickle file of first sample
  pickle_path_ps1 = args.ps1
  bam_phase_set_dictionary_ps1 = pickle.load(pickle_path_ps1)
  vcf_variants_dictionary_ps1 = pickle.load(pickle_path_ps1)

  # path to input pickle file of second sample
  pickle_path_ps2 = args.ps2
  bam_phase_set_dictionary_ps2 = pickle.load(pickle_path_ps2)
  vcf_variants_dictionary_ps2 = pickle.load(pickle_path_ps2)

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

  # Steven write functions here

  extended_phase_set_dictionary = extend_phase_sets(
    bam_phase_set_dictionary_ps1,
    bam_phase_set_dictionary_ps2)

  # write output and close files
  os.makedirs(args.output_directory, exist_ok = True)
  output_file_path = os.path.join(args.output_directory, args.output_prefix + ".extended_phase_sets.tsv")
  output_file = open(output_file_path, 'w')
  # write results to output file here
  output_file.close()

if __name__ == '__main__':
  main(args)
