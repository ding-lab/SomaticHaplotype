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
      if psid not in ["variant_not_phased_heterozygote", "variant_phase_set_not_in_bam"]:
        target_phase_sets[psid] = phase_sets[psid]
  else:
    if start == None: 
      start = 0
    if end == None:
      end = float('inf')
    for psid in phase_sets.keys():
      if psid not in ["variant_not_phased_heterozygote", "variant_phase_set_not_in_bam"] \
        and phase_sets[psid]._chr == chrom and phase_sets[psid].psStart() >= start and phase_sets[psid].psEnd() <= end:
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
  return phase_set_lengths, summary_stats

# Find the distribution of phase set gaps, summary stats
def summarize_phase_set_gaps(target_phase_sets):
  phase_set_gaps = []
  for psid in sorted(target_phase_sets.keys()):
    # print(target_phase_sets[psid]._key, target_phase_sets[psid]._chr, target_phase_sets[psid].getFirstVariantPosition(), target_phase_sets[psid].getLastVariantPosition())
    pass

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
  return N50

# Find number of variants on H1, H2, unphased, and total
# def summarize_phase_set_variants(target_phase_sets):
#   num_H1 = 0
#   num_H2 = 0
#   num_unphased = 0
#   for psid in target_phase_sets.keys():

# Write output file with summary stats
def write_output_summary(output_summary_file, phase_set_lengths_summary_stats, N50):
  output_summary_file.write("\t".join(["", "count", "mean", "std", "min", "25%", "50%", "75%", "max"])+"\n")
  output_summary_file.write("phase_set_lengths\t"+"\t".join([str(x) for x in list(phase_set_lengths_summary_stats)])+"\n")
  # output_summary_file.write("phase_set_gaps\t"+"\t".join(list(phase_set_gaps_summary_stats))+"\n")
  output_summary_file.write("N50\t" + "NA\t" * 5 + str(N50) + "\tNA\tNA\n")

# Write output file with info of each phase set
def write_output_ps(output_ps_file, target_phase_sets):
  output_ps_file.write("\t".join(["ps_id", "chr", "start", "end", "length", "first_variant_pos", "last_variant_pos", \
    "n_variants_H1", "n_variants_H2", "n_variants_unphased", "n_variants_total"])+"\n")
  for psid in target_phase_sets.keys():
    ps = target_phase_sets[psid]
    output_ps_file.write("\t".join([psid, ps._chr, str(ps._start), str(ps._end), str(ps.length()), str(ps.getFirstVariantPosition()), str(ps.getLastVariantPosition()), \
      "NA", "NA", "NA", str(len(ps.getVariants().keys()))])+"\n")


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
  phase_set_lengths, phase_set_lengths_summary_stats = summarize_phase_set_lengths(target_phase_sets)
  summarize_phase_set_gaps(target_phase_sets)
  N50 = compute_N50(phase_set_lengths)

  # write output and close files
  os.makedirs(args.output_directory, exist_ok = True)
  output_summary_file_path = os.path.join(args.output_directory, args.output_prefix + ".summary.tsv")
  output_ps_file_path = os.path.join(args.output_directory, args.output_prefix + ".ps.tsv")
  output_summary_file = open(output_summary_file_path, 'w')
  output_ps_file = open(output_ps_file_path, 'w')
  # write results to output file here
  write_output_summary(output_summary_file, phase_set_lengths_summary_stats, N50)
  write_output_ps(output_ps_file, target_phase_sets)
  output_summary_file.close()
  output_ps_file.close()

if __name__ == '__main__':
  main(args)
