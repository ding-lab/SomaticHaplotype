import os
import pickle
import pysam
import vcf

from SomaticHaplotype import *

################################################################################
# VCF functions
################################################################################

def extract_variants_from_VCF(vcf_filename, sample_id, chr = None, start_bp = None, end_bp = None):
  # build a dictionary of phase sets present in VCF
  # only includes variants that are phased heterozygotes (useful for phase-related activities)
  this_vcf = vcf.Reader( filename = vcf_filename )
  variant_dict = {} # dictionary to hold all variants

  for record in this_vcf.fetch( str(chr) , start_bp, end_bp ): # loop over each record in VCF
    this_variant = Variant(record, sample_id)
    if this_variant.return_VariantKey() in variant_dict: # check if variant already in dictionary
      variant_dict[this_variant.return_VariantKey()].append(this_variant)
    else:
      variant_dict[this_variant.return_VariantKey()] = [this_variant]

  return(variant_dict)

################################################################################
# Summary file functions
################################################################################

def assign_variants_to_each_phaseset(summary_file_path, vcf_variants_dictionary, chr = None, start_bp = None, end_bp = None):
  # read in summary file and discard header line
  summary_file = open(summary_file_path, "r")
  summary_file.readline()

  ps_id_variants_dict = {}

  for line in summary_file:
    ps_id, this_chrom, start, end, length_reads, first_variant_pos, last_variant_pos, length_variants, n_variants_H1, n_variants_H2, n_variants_total = line.strip().split()
    if ps_id in ps_id_variants_dict:
      sys.exit("Phase set ID " + ps_id + " already in ps_id_variants_dict.")
    else:
      ps_id_variants_dict[ps_id] = []
      if first_variant_pos == "NA":
        continue
      else:
        first_variant_pos = int(first_variant_pos)
        last_variant_pos = int(last_variant_pos)
        if this_chrom == chr and first_variant_pos > start_bp and last_variant_pos < end_bp:
          for k,variants in vcf_variants_dictionary.items():
            for v in variants:
              if v.return_Chromosome() == chr and v.return_Position() >= first_variant_pos and v.return_Position() <= last_variant_pos:
                keep_variant = all([v.return_IsHeterozygote(), v.return_IsSNP(), v.return_Filter() in [[], ['10X_PHASING_INCONSISTENT']]])
                if keep_variant:
                  ps_id_variants_dict[ps_id].append(v.return_VariantKey)
  summary_file.close()
  return(ps_id_variants_dict)

def create_phaseset_summary_dict(summary_file_path, chr = None, start_bp = None, end_bp = None):
  # read in summary file and discard header line
  summary_file = open(summary_file_path, "r")
  summary_file.readline()

  ps_id_dict = {}

  for line in summary_file:
    ps_id = line.strip().split()[0]
    if ps_id in ps_id_dict:
      sys.exit("Phase set ID " + ps_id + " already in ps_id_variants_dict.")
    else:
      ps_id_dict[ps_id] = line.strip().split()
        
  summary_file.close()
  return(ps_id_dict)

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
  
  # vcf variant dictionary
  vcf_variants_dictionary = extract_variants_from_VCF(
    args.vcf,
    args.vcf_id,
    chr = chrom,
    start_bp = start,
    end_bp = end)

  # create dictionary with phase set keys and list of variants as values
  phasesets_with_variants_assigned_dict = assign_variants_to_each_phaseset(args.sum, vcf_variants_dictionary, chr = chrom, start_bp = start, end_bp = end)

  phaseset_summary_dict = create_phaseset_summary_dict(args.sum, chr = chrom, start_bp = start, end_bp = end)

  for k,v in phaseset_summary_dict.items():
    print(v)

  for k,v in phasesets_with_variants_assigned_dict.items():
    print(k, len(v)) 

  # add variants to bam phase set dictionary
  #add_variants_to_phase_sets(bam_phase_set_dictionary, vcf_variants_dictionary)

  os.makedirs(args.output_directory, exist_ok = True)
  output_file_path = os.path.join(args.output_directory, args.output_prefix + ".phaseset.pkl")
  #output_file = open(output_file_path, 'w')
  #pickle.dump(bam_phase_set_dictionary, output_file, pickle.HIGHEST_PROTOCOL)
  #pickle.dump(vcf_variants_dictionary, output_file, pickle.HIGHEST_PROTOCOL)
  #output_file.close()

if __name__ == '__main__':
  main(args)
