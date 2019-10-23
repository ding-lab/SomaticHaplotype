import os
import pickle
#import pysam
#import vcf

from SomaticHaplotype import *

################################################################################
# functions
################################################################################

def extend_phase_sets(ps_dict1, ps_dict2, chrom, start, end):

  # for each phase set in ps1 check if ps overlaps range
  for ps in [ x for x in list(ps_dict['phase_sets'].keys()) if x.startswith(chrom) ]:
    if ranges_overlap(ps.return_Chromosome(), ps.return_FirstVariantPosition(), ps.return_LastVariantPosition(), chrom, start, end): # ps range overlaps given range
      # DO WORK HERE
    else: # ps range does not overlap given range

  ps_dict['phase_sets']['chr1:109933436'].return_FirstVariantPosition()
  ps_dict['phase_sets']['chr1:109933436'].return_LastVariantPosition()
  ps_dict['phase_sets']['chr1:109933436'].return_Chromosome()

  return(extended_ps_dict1)

def ranges_overlap(chr1, start1, end1, chr2, start2, end2):
  if start2 is None: # if no start2, assume it is start1
    start2 = start1
  if end2 is None: # if no end2, assume it is end1
    end2 = end1

  if start1 > end1 or start2 > end2:
    sys.exit("Error in extend module ranges_overlap: start1 > end1 or start2 > end2.")

  if chr2 is None: # no range given
    return(True)
  elif chr1 != chr2: # ranges on different chromosomes
    return(False)
  elif len(set(range(start1, end1)).intersection(range(start2, end2))) > 0:
      return(True)
  else:
    return(False)

################################################################################
# main
################################################################################   

def main(args):

  # path to input pickle file of first sample
  pickle_path_ps1 = args.ps1
  with open(pickle_path_ps1, 'rb') as pickle_file_ps1:
    bam_phase_set_dictionary_ps1 = pickle.load(pickle_file_ps1)
    vcf_variants_dictionary_ps1 = pickle.load(pickle_file_ps1)

  # path to input pickle file of second sample
  pickle_path_ps2 = args.ps2
  with open(pickle_path_ps2, 'rb') as pickle_file_ps2:
    bam_phase_set_dictionary_ps2 = pickle.load(pickle_file_ps2)
    vcf_variants_dictionary_ps2 = pickle.load(pickle_file_ps2)

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
    bam_phase_set_dictionary_ps2,
    chrom, start, end)

  # write output and close files
  os.makedirs(args.output_directory, exist_ok = True)
  output_file_path = os.path.join(args.output_directory, args.output_prefix + ".extended_phase_sets.tsv")
  output_file = open(output_file_path, 'w')
  # write results to output file here

  output_file.close()

if __name__ == '__main__':
  main(args)
