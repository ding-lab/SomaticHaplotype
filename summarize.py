import os
import pickle
import pandas as pd

from SomaticHaplotype import *

################################################################################
# functions
################################################################################

# Extract phase blocks within the genomic region
def extract_phase_blocks(bam_phase_block_dictionary, chrom, start, end):
  target_phase_blocks = {}
  phase_blocks = bam_phase_block_dictionary['phase_blocks']
  if chrom == None: # get all phase blocks
    for pbid in phase_blocks.keys():
      if pbid not in ["variant_not_phased_heterozygote", "variant_phase_block_not_in_bam"]:
        target_phase_blocks[pbid] = phase_blocks[pbid]
  else:
    if start == None: 
      start = 0
    if end == None:
      end = float('inf')
    for pbid in phase_blocks.keys():
      if pbid not in ["variant_not_phased_heterozygote", "variant_phase_block_not_in_bam"] \
        and phase_blocks[pbid].return_Chromosome() == chrom and phase_blocks[pbid].return_Start() >= start and phase_blocks[pbid].return_End() <= end:
        target_phase_blocks[pbid] = phase_blocks[pbid]
  return target_phase_blocks

# Find the distribution of phase block lengths (defined by reads), summary stats
def summarize_phase_block_lengths(target_phase_blocks):
  phase_block_lengths = []
  for pbid in target_phase_blocks.keys():
    phase_block_lengths.append(target_phase_blocks[pbid].return_LengthReads())
  summary_stats = pd.Series(phase_block_lengths).describe()
  return phase_block_lengths, summary_stats

# Find N50 of phase block lengths 
# (ref: http://seqanswers.com/forums/showthread.php?t=2332)
def compute_N50(phase_block_lengths):
  # the Broad way
  l = sorted(phase_block_lengths)
  l_prime = [[length]*length for length in l]
  # ref: https://stackoverflow.com/questions/952914/making-a-flat-list-out-of-list-of-lists-in-python
  flat_l_prime = [item for sublist in l_prime for item in sublist]
  N50_broad = pd.Series(flat_l_prime).median()
  # Miller definition (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2874646/)
  half_total_length = sum(phase_block_lengths)/2
  l = sorted(phase_block_lengths, reverse = True)
  running_total = 0
  for length in l:
    if (running_total + length) >= half_total_length: 
      N50_miller = float(length)
      break
    running_total += length
  return N50_broad, N50_miller

# Write output file with summary stats
def write_output_summary(output_summary_file, phase_block_lengths_summary_stats, N50_broad, N50_miller):
  output_summary_file.write("\t".join(["count", "mean", "std", "min", "25%", "50%", "75%", "max", "N50_broad", "N50_miller"])+"\n")
  output_summary_file.write("\t".join([str(x) for x in list(phase_block_lengths_summary_stats)])+"\t"+str(N50_broad)+"\t"+str(N50_miller)+"\n")

# Write output file with info of each phase block
def write_output_pb(output_pb_file, target_phase_blocks):
  output_pb_file.write("\t".join(["pb_id", "chr", "start", "end", "length_reads", "first_variant_pos", "last_variant_pos", \
    "length_variants", "n_variants_H1", "n_variants_H2", "n_variants_total"])+"\n")
  for pbid in target_phase_blocks.keys():
    pb = target_phase_blocks[pbid]
    output_pb_file.write("\t".join([pbid, pb.return_Chromosome(), str(pb.return_Start()), str(pb.return_End()), str(pb.return_LengthReads()), str(pb.return_FirstVariantPosition()), str(pb.return_LastVariantPosition()), \
      str(pb.return_LengthVariants()), str(pb.return_nVariantsH1()), str(pb.return_nVariantsH2()), str(len(pb.return_Variants().keys()))])+"\n")


################################################################################
# main
################################################################################   

def main(args):

  # path to input pickle file
  pickle_path = args.pb1
  with open(pickle_path, 'rb') as pickle_file:
    bam_phase_block_dictionary = pickle.load(pickle_file)
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

  target_phase_blocks = extract_phase_blocks(bam_phase_block_dictionary, chrom, start, end)
  phase_block_lengths, phase_block_lengths_summary_stats = summarize_phase_block_lengths(target_phase_blocks)
  N50_broad, N50_miller = compute_N50(phase_block_lengths)

  # write output and close files
  os.makedirs(args.output_directory, exist_ok = True)
  output_summary_file_path = os.path.join(args.output_directory, args.output_prefix + ".phase_block_summary.tsv")
  output_pb_file_path = os.path.join(args.output_directory, args.output_prefix + ".phase_blocks.tsv")
  output_summary_file = open(output_summary_file_path, 'w')
  output_pb_file = open(output_pb_file_path, 'w')
  write_output_summary(output_summary_file, phase_block_lengths_summary_stats, N50_broad, N50_miller)
  write_output_pb(output_pb_file, target_phase_blocks)
  output_summary_file.close()
  output_pb_file.close()

if __name__ == '__main__':
  main(args)
