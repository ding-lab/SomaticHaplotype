import os
import pickle
#import pysam
#import vcf
import pandas as pd

from SomaticHaplotype import *

################################################################################
# functions
################################################################################

# Extract phase sets within the genomic region
def extract_phase_sets(bam_phase_set_dictionary, chrom, start, end):
  target_phase_sets = {}
  phase_sets = bam_phase_set_dictionary['phase_sets']
  if chrom == None: # get all phase sets
    for psid in phase_sets.keys():
      if psid != "variant_not_phased_heterozygote":
        target_phase_sets[psid] = phase_sets[psid]
  else:
    if start == None: 
      start = 0
    if end == None:
      end = float('inf')
    for psid in phase_sets.keys():
      if psid != "variant_not_phased_heterozygote" and phase_sets[psid]._chr == chrom and phase_sets[psid].psStart() >= start and phase_sets[psid].psEnd() <= end:
        target_phase_sets[psid] = phase_sets[psid]
  return target_phase_sets

# Find the distribution of phase set lengths, summary stats
def summarize_phase_set_lengths(target_phase_sets):
  phase_set_lengths = []
  for psid in target_phase_sets.keys():
    phase_set_lengths.append(target_phase_sets[psid].length())
  summary_stats = pd.Series(phase_set_lengths).describe()
  print("--Distribution of phase set lengths--")
  print(summary_stats)
  return phase_set_lengths

# Find the distribution of phase set gaps, summary stats
def summarize_phase_set_gaps(target_phase_sets):
  phase_set_gaps = []
  for psid in sorted(target_phase_sets.keys()):
    print(target_phase_sets[psid])

# Find N50 of phase set lengths 
# (ref: http://seqanswers.com/forums/showthread.php?t=2332)
def compute_N50(phase_set_lengths):
  # the Broad way
  l = sorted(phase_set_lengths)
  l_prime = [[length]*length for length in l]
  flat_l_prime = [item for sublist in l_prime for item in sublist] #ref: https://stackoverflow.com/questions/952914/making-a-flat-list-out-of-list-of-lists-in-python
  N50 = pd.Series(flat_l_prime).median()
  print("Broad-style:", N50)
  # Miller definition (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2874646/)
  half_total_length = sum(phase_set_lengths)/2
  l = sorted(phase_set_lengths, reverse = True)
  running_total = 0
  for length in l:
    if (running_total + length) >= half_total_length: 
      N50 = length
      break
    running_total += length
  print("Miller-style:", N50)

# Find number of variants on H1, H2, unphased, and total


################################################################################
# main
################################################################################   

def main(args):

  # path to input pickle file
  pickle_path = args.ps1
  # bam_phase_set_dictionary = pickle.load(pickle_path)
  # vcf_variants_dictionary = pickle.load(pickle_path)
  with open(pickle_path, 'rb') as pickle_file:
    bam_phase_set_dictionary = pickle.load(pickle_file)
    vcf_variants_dictionary = pickle.load(pickle_file)

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
  target_phase_sets = extract_phase_sets(bam_phase_set_dictionary, chrom, start, end)
  phase_set_lengths = summarize_phase_set_lengths(target_phase_sets)
  # summarize_phase_set_gaps(target_phase_sets)
  compute_N50(phase_set_lengths)
'''
  # write output and close files
  os.makedirs(args.output_directory, exist_ok = True)
  output_file_path = os.path.join(args.output_directory, args.output_prefix + "YOUR SUFFIX HERE.tsv")
  output_file = open(output_file_path, 'w')
  # write results to output file here
  output_file.close()
'''
if __name__ == '__main__':
  main(args)
