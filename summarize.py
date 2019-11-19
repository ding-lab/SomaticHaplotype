import os
import pickle
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
        and phase_sets[psid].return_Chromosome() == chrom and phase_sets[psid].return_Start() >= start and phase_sets[psid].return_End() <= end:
        target_phase_sets[psid] = phase_sets[psid]
  return target_phase_sets

# Find the distribution of phase set lengths (defined by reads), summary stats
def summarize_phase_set_lengths(target_phase_sets):
  phase_set_lengths = []
  for psid in target_phase_sets.keys():
    phase_set_lengths.append(target_phase_sets[psid].return_LengthReads())
  summary_stats = pd.Series(phase_set_lengths).describe()
  return phase_set_lengths, summary_stats

# Find N50 of phase set lengths 
# (ref: http://seqanswers.com/forums/showthread.php?t=2332)
def compute_N50(phase_set_lengths):
  # the Broad way
  l = sorted(phase_set_lengths)
  l_prime = [[length]*length for length in l]
  # ref: https://stackoverflow.com/questions/952914/making-a-flat-list-out-of-list-of-lists-in-python
  flat_l_prime = [item for sublist in l_prime for item in sublist]
  N50_broad = pd.Series(flat_l_prime).median()
  # Miller definition (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2874646/)
  half_total_length = sum(phase_set_lengths)/2
  l = sorted(phase_set_lengths, reverse = True)
  running_total = 0
  for length in l:
    if (running_total + length) >= half_total_length: 
      N50_miller = length
      break
    running_total += length
  return N50_broad, N50_miller

# Write output file with summary stats
def write_output_summary(output_summary_file, phase_set_lengths_summary_stats, N50_broad, N50_miller):
  output_summary_file.write("\t".join(["", "count", "mean", "std", "min", "25%", "50%", "75%", "max"])+"\n")
  output_summary_file.write("phase_set_lengths\t"+"\t".join([str(x) for x in list(phase_set_lengths_summary_stats)])+"\n")
  output_summary_file.write("N50_Broad\t" + "NA\t" * 5 + str(N50_broad) + "\tNA\tNA\n")
  output_summary_file.write("N50_Miller\t" + "NA\t" * 5 + str(N50_miller) + "\tNA\tNA\n")

# Write output file with info of each phase set
def write_output_ps(output_ps_file, target_phase_sets):
  output_ps_file.write("\t".join(["ps_id", "chr", "start", "end", "length_reads", "first_variant_pos", "last_variant_pos", \
    "length_variants", "n_variants_H1", "n_variants_H2", "n_variants_total"])+"\n")
  for psid in target_phase_sets.keys():
    ps = target_phase_sets[psid]
    output_ps_file.write("\t".join([psid, ps.return_Chromosome(), str(ps.return_Start()), str(ps.return_End()), str(ps.return_LengthReads()), str(ps.return_FirstVariantPosition()), str(ps.return_LastVariantPosition()), \
      str(ps.return_LengthVariants()), str(ps.return_nVariantsH1()), str(ps.return_nVariantsH2()), str(len(ps.return_Variants().keys()))])+"\n")


################################################################################
# main
################################################################################   

def main(args):

  # path to input pickle file
  pickle_path = args.ps1
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

  target_phase_sets = extract_phase_sets(bam_phase_set_dictionary, chrom, start, end)
  phase_set_lengths, phase_set_lengths_summary_stats = summarize_phase_set_lengths(target_phase_sets)
  N50_broad, N50_miller = compute_N50(phase_set_lengths)

  # write output and close files
  os.makedirs(args.output_directory, exist_ok = True)
  output_summary_file_path = os.path.join(args.output_directory, args.output_prefix + ".phase_set_summary.tsv")
  output_ps_file_path = os.path.join(args.output_directory, args.output_prefix + ".phase_sets.tsv")
  output_summary_file = open(output_summary_file_path, 'w')
  output_ps_file = open(output_ps_file_path, 'w')
  write_output_summary(output_summary_file, phase_set_lengths_summary_stats, N50)
  write_output_ps(output_ps_file, target_phase_sets)
  output_summary_file.close()
  output_ps_file.close()

if __name__ == '__main__':
  main(args)
