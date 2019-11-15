import numpy
import os
import pickle
import vcf

from SomaticHaplotype import *

################################################################################
# functions
################################################################################

def detect_run(flip_list):
  if len(flip_list) == 1:
    return("NA")
  else:
    n_shared = len(flip_list)
    n_flip = sum(flip_list)
    runs = 1
    for i in range(1, n_shared):
      if flip_list[i] != flip_list[i-1]:
        runs += 1
    background_runs = []
    for j in range(1000):
      x = numpy.random.binomial(1, float(n_flip)/n_shared, n_shared)
      x_runs = 1
      for i in range(1, n_shared):
        if x[i] != x[i-1]:
          x_runs += 1
      background_runs.append(x_runs)
    probability_extreme_flip = numpy.mean(runs <= numpy.array(background_runs))
    return(probability_extreme_flip)

def extract_variants_from_VCF(vcf_filename, sample_id, chr = None, start_bp = None, end_bp = None):
  # build a dictionary of phase sets present in VCF
  this_vcf = vcf.Reader( filename = vcf_filename )
  variant_dict = {} # dictionary to hold all variants

  for record in this_vcf.fetch( str(chr) , start_bp, end_bp ): # loop over each record in VCF
    this_variant = Variant(record, sample_id)
    if this_variant.return_VariantKey() in variant_dict: # check if variant already in dictionary
      variant_dict[this_variant.return_VariantKey()].append(this_variant)
    else:
      variant_dict[this_variant.return_VariantKey()] = [this_variant]

  return(variant_dict)

def overlap_phase_sets_with_phased_variants(ps_dict1, vcf_dict2, chrom, start, end):

  # value that will be returned:
  overlapping_variants_dict = {}

  # for each phase set in ps1 check if ps overlaps range
  ps1_phase_sets_in_range = []
  for ps in [ x for x in list(ps_dict1['phase_sets'].keys()) if x.startswith(chrom + ":") ]:
    phase_set_1 = ps_dict1['phase_sets'][ps]

    if phase_set_1.return_FirstVariantPosition() == "NA" or phase_set_1.return_LastVariantPosition() == "NA":
      next
    elif ranges_overlap(phase_set_1.return_Chromosome(), phase_set_1.return_FirstVariantPosition(), phase_set_1.return_LastVariantPosition(), chrom, start, end): # ps range overlaps given range
      ps1_phase_sets_in_range.append(ps)

  # for each ps1 phase set in range, iterate over each variant in vcf_dict2
  for ps1 in ps1_phase_sets_in_range:
    overlapping_variants_dict[ps1] = []
    n_variants_overlap = 0
    n_variants_flip_to_match = 0
    total_variants_this_phase_set = 0
    run_of_flips = []
    phase_set_1 = ps_dict1['phase_sets'][ps1]
    for variant_key in phase_set_1.return_Variants().keys():
      total_variants_this_phase_set += 1
      variant1 = phase_set_1.return_Variants()[variant_key][0]
      if variant_key in vcf_dict2:
        variant2 = vcf_dict2[variant_key][0]
        n_variants_overlap += 1
        GT1 = variant1.return_Genotype()
        GT2 = variant2.return_Genotype()
        if GT1 == GT2[::-1]:
          n_variants_flip_to_match += 1
          run_of_flips.append(1)
        elif GT1 == GT2:
          run_of_flips.append(0)

    run_probability = detect_run(run_of_flips)

    # define null as random fair coin flip
    # much more conservative, especially for low n_variants_overlap
    p_value_switch = numpy.mean(numpy.random.binomial(n_variants_overlap, .5, 1000) >= n_variants_flip_to_match)
    p_value_no_switch = numpy.mean(numpy.random.binomial(n_variants_overlap, .5, 1000) <= n_variants_flip_to_match)

    min_p_value = min(p_value_switch, p_value_no_switch)

    if p_value_switch < p_value_no_switch and p_value_switch < 0.01:
      recommendation = "Switch"
      graph_weight = 1
    elif p_value_no_switch < p_value_switch and p_value_no_switch < 0.01:
      recommendation = "No Switch"
      graph_weight = 2
    else:
      recommendation = "No Recommendation"
      graph_weight = 0

    overlapping_variants_dict[ps1] = [ps1, phase_set_1.return_Chromosome(), phase_set_1.return_FirstVariantPosition(), phase_set_1.return_LastVariantPosition(), total_variants_this_phase_set, n_variants_overlap, n_variants_flip_to_match, run_probability, p_value_switch, p_value_no_switch, min_p_value, recommendation, graph_weight] # graph_weight must be final element of list to work with create_graph()
  
  return(overlapping_variants_dict)

def ranges_overlap(chr1, start1, end1, chr2, start2, end2):
  if start2 is None: # if no start2, assume it is start1
    start2 = start1
  if end2 is None: # if no end2, assume it is end1
    end2 = end1

  if chr2 is None: # no range given
    return(True)
  elif chr1 != chr2: # ranges on different chromosomes
    return(False)
  elif(end2 < start1 or start2 > end1):
    return(False)
  elif len(set(range(start1, end1)).intersection(range(start2, end2))) > 0:
      return(True)
  else:
    return(False)

################################################################################
# main
################################################################################   

def main(args):

  # open up summary file
  #summary_file = open(args.sum, 'r')

  # path to input pickle file of first sample
  pickle_path_ps1 = args.ps1
  with open(pickle_path_ps1, 'rb') as pickle_file_ps1:
    bam_phase_set_dictionary_ps1 = pickle.load(pickle_file_ps1)
    vcf_variants_dictionary_ps1 = pickle.load(pickle_file_ps1)

  # parse the genomic range argument (which is required)
  chrom = str(args.range.split(":")[0])
  try:
    start = int(args.range.split(":")[1].split("-")[0])
  except:
    start = None
  try:
    end = int(args.range.split(":")[1].split("-")[1])
  except:
    end = None

  # run functions
  vcf_variants_dictionary = extract_variants_from_VCF(
    args.vcf,
    args.vcf_id,
    chr = chrom,
    start_bp = start,
    end_bp = end)

  overlapping_variants_dict = overlap_phase_sets_with_phased_variants(bam_phase_set_dictionary_ps1, vcf_variants_dictionary, chrom, start, end)

  for k,v in overlapping_variants_dict.items():
    print(k, v)

  # write outputs
  #s.makedirs(args.output_directory, exist_ok = True)
  #output_file_path = os.path.join(args.output_directory, args.output_prefix + "")
  #output_file = open(output_file_path, 'w')
  #output_file.close()

if __name__ == '__main__':
  main(args)
