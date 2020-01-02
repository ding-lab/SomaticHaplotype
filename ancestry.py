import gzip
import numpy
import os
import pickle
import vcf

from SomaticHaplotype import *

################################################################################
# functions
################################################################################

def ancestry_overlap(ibd_dictionary, demographic_dictionary):
  # Gather ancestry information associated with each IBD segment
  print_dictionary = {}

  for index in ibd_dictionary.keys():
    other_id = ibd_dictionary[index]["other_id"]
    print_list = [index, ibd_dictionary[index]["chr"], ibd_dictionary[index]["start"], ibd_dictionary[index]["end"], 
    ibd_dictionary[index]["this_hap"], ibd_dictionary[index]["hbd"], 
    other_id, demographic_dictionary[other_id]["super_pop"], demographic_dictionary[other_id]["pop"]]
    print_string = "\t".join([str(x) for x in print_list]) + "\n"
    print_dictionary[index] = print_string

  return(print_dictionary)

def overlap_IBD_segments_with_phase_sets(ibd_dictionary, ps_dict, ibd_variants_dictionary):
  
  print_dict = {}
  print_dict["column_names"] = ["IBD_segment_index", "IBD_or_HBD", "Sample1_ID", "Sample1_haplotype", "Sample2_ID", "Sample2_haplotype", "Chromosome", "IBD_start_position", "IBD_end_positions", "IBD_LOD_score", "IBD_cM_length"]
  print_dict["column_names"].extend(["Phase_set_key", "Phase_set_start_position", "Phase_set_end_position"])
  print_dict["column_names"].extend(["Number_of_overlapping_phased_heterozygous_variants", "Proportion_of_exact_genotype_agreement", "Proportion_of_IBD_variants_overlapping_H1", "Proportion_of_IBD_variants_overlapping_H2", "Overlapping_variant_positions", "Overlapping_variant_keys", "Phase_set_H1_alleles", "Phase_set_H2_alleles", "IBD_alleles"])
 

  for ibd_segment_index, ibd_segment_dict in ibd_dictionary.items():

    if ibd_segment_dict["hbd"]:
      
      print_dict[ibd_segment_index] = [ibd_segment_index, "HBD"]
      print_dict[ibd_segment_index].extend(["NA"]*(len(print_dict["column_names"]) - 2))
    
    else:
    
      for phase_set_key in ps_dict["phase_sets"].keys():

        if phase_set_key.startswith(ibd_segment_dict["chr"]):
        
          phase_set = ps_dict["phase_sets"][phase_set_key]

          if phase_set.return_FirstVariantPosition() == "NA" or phase_set.return_LastVariantPosition() == "NA":
        
            pass
        
          elif ranges_overlap(ibd_segment_dict["chr"], ibd_segment_dict["start"], ibd_segment_dict["end"], phase_set.return_Chromosome(), phase_set.return_FirstVariantPosition(), phase_set.return_LastVariantPosition()):
        
            print_list = [ibd_segment_index, "IBD", ibd_segment_dict["this_id"], ibd_segment_dict["this_hap"], ibd_segment_dict["other_id"], ibd_segment_dict["other_hap"], ibd_segment_dict["chr"], ibd_segment_dict["start"], ibd_segment_dict["end"], ibd_segment_dict["LOD"], ibd_segment_dict["cM"]]
            print_list.extend([phase_set_key, phase_set.return_FirstVariantPosition(), phase_set.return_LastVariantPosition()]) # modified print_list        

            overlapping_phased_het_variants = 0

            ps_variants_dictionary = phase_set.return_Variants()
            ps_variant_keys = ps_variants_dictionary.keys()

            ibd_hap = int(ibd_segment_dict["this_hap"])
            if ibd_hap == 1:
              ibd_index = 0
    
            elif ibd_hap == 2:
              ibd_index = 2

            ibd_alleles = []
            H1_alleles = []
            H2_alleles = []
            variant_positions = []
            overlapping_variant_keys = []
            overlaps_H1 = 0
            overlaps_H2 = 0
            exact_agreement = 0

            for variant_key in ps_variant_keys:
    
              if variant_key in ibd_variants_dictionary.keys():
                ps_variant = ps_variants_dictionary[variant_key][0]
                ibd_variant = ibd_variants_dictionary[variant_key][0]
    
                if ibd_variant.return_IsPhasedHeterozygote() and ps_variant.return_IsPhasedHeterozygote():
                  overlapping_variant_keys.append(variant_key)
                  variant_positions.append(ps_variant.return_Position())
                  H1_alleles.append(ps_variant.return_Genotype()[0])
                  H2_alleles.append(ps_variant.return_Genotype()[2])
                    
                  overlapping_phased_het_variants += 1
                  ibd_allele = ibd_variant.return_Genotype()[ibd_index] #0|0 0|1 1|0 1|1 #ibd_index (0 or 2) corresponds to ibd_hap (1 or 2)
                  ibd_alleles.append(ibd_allele)
                    
                  if ps_variant.return_Genotype() == ibd_variant.return_Genotype():
                    exact_agreement += 1

                  if ibd_allele == ps_variant.return_Genotype()[0]:
                    overlaps_H1 += 1
    
                  elif ibd_allele == ps_variant.return_Genotype()[2]:
                    overlaps_H2 += 1
    
                  else:
                    continue
    
                else:
                  continue
    
              else:
                continue

            if overlapping_phased_het_variants == 0:
              exact_agreement_proportion = "NA"
              proportion_ibd_overlaps_H1 = "NA"
              proportion_ibd_overlaps_H2 = "NA"
            else:
              exact_agreement_proportion = exact_agreement/overlapping_phased_het_variants
              proportion_ibd_overlaps_H1 = overlaps_H1/overlapping_phased_het_variants
              proportion_ibd_overlaps_H2 = overlaps_H2/overlapping_phased_het_variants
            print_list.extend([overlapping_phased_het_variants, exact_agreement_proportion, proportion_ibd_overlaps_H1, proportion_ibd_overlaps_H2, ",".join([str(x) for x in variant_positions]), ",".join([str(x) for x in overlapping_variant_keys]), ",".join([str(x) for x in H1_alleles]), ",".join([str(x) for x in H2_alleles]), ",".join([str(x) for x in ibd_alleles])])
            print_dict[ibd_segment_index] = print_list

          else:
            continue

  return(print_dict)

# def detect_run(flip_list):
#   if len(flip_list) == 1:
#     return("NA")
#   else:
#     n_shared = len(flip_list)
#     n_flip = sum(flip_list)
#     runs = 1
#     for i in range(1, n_shared):
#       if flip_list[i] != flip_list[i-1]:
#         runs += 1
#     background_runs = []
#     for j in range(1000):
#       x = numpy.random.binomial(1, float(n_flip)/n_shared, n_shared)
#       x_runs = 1
#       for i in range(1, n_shared):
#         if x[i] != x[i-1]:
#           x_runs += 1
#       background_runs.append(x_runs)
#     probability_extreme_flip = numpy.mean(runs <= numpy.array(background_runs))
#     return(probability_extreme_flip)

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

def ranges_overlap(chr1, start1, end1, chr2, start2, end2):

  if start1 is None: # if no start1, assume it is start2
    start1 = int(start2)
  if end1 is None: # if no end1, assume it is end2
    end1 = int(end2)

  if start2 is None: # if no start2, assume it is start1
    start2 = int(start1)
  if end2 is None: # if no end2, assume it is end1
    end2 = int(end1)

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

  # hbd_dictionary
  hbd_file = gzip.open(args.hbd, 'rb')
  hbd_dict = {}
  hbd_segment = 0
  for line in hbd_file:
    line_strip_split = [x.decode("utf-8") for x in line.strip().split()]

    if args.vcf_id == line_strip_split[0] or args.vcf_id == line_strip_split[2]:
      # HBD segment involves sample of interest  
      hbd_chr = line_strip_split[4]
      hbd_start = int(line_strip_split[5])
      hbd_end = int(line_strip_split[6])

      #if ranges_overlap(chrom, start, end, hbd_chr, hbd_start, hbd_end):
        #HBD segment overlap range of interest
      if chrom == hbd_chr:
        hbd_segment += 1
        hbd_dict[hbd_segment] = {k:v for k,v in zip(["chr", "start", "end"], [hbd_chr, hbd_start, hbd_end])}
  
  hbd_file.close()

  # ibd_dictionary
  ibd_file = gzip.open(args.ibd, 'r')
  ibd_dict = {}
  ibd_segment = 0
  vcf_id_b = args.vcf_id.encode('ascii')
  for line in ibd_file:
    line_strip_split_bytes = line.strip().split()

    if vcf_id_b == line_strip_split_bytes[0] or vcf_id_b == line_strip_split_bytes[2]:
      # IBD segment involves sample of interest
      line_strip_split = [x.decode("utf-8")  for x in line.strip().split()]

      ibd_chr = line_strip_split[4]
      ibd_start = int(line_strip_split[5])
      ibd_end = int(line_strip_split[6])

      if ranges_overlap(chrom, start, end, ibd_chr, ibd_start, ibd_end):
        #IBD segment overlap range of interest
  
        hbd_status = False
        for k in hbd_dict.keys():
          if ranges_overlap(ibd_chr, ibd_start, ibd_end, hbd_dict[k]["chr"], hbd_dict[k]["start"], hbd_dict[k]["end"]):
            hbd_status = True
            break

        if args.vcf_id == line_strip_split[0]:
          ibd_segment += 1
          ibd_dict[ibd_segment] = {k:v for k,v in zip(["this_id", "this_hap", "other_id", "other_hap", "chr", "start", "end", "LOD", "cM"], line_strip_split)}
        elif args.vcf_id == line_strip_split[2]:
          ibd_segment += 1
          ibd_dict[ibd_segment] = {k:v for k,v in zip(["other_id", "other_hap", "this_id", "this_hap", "chr", "start", "end", "LOD", "cM"], line_strip_split)}
          
        ibd_dict[ibd_segment]["hbd"] = hbd_status
        ibd_dict[ibd_segment]["start"] = int(ibd_dict[ibd_segment]["start"])
        ibd_dict[ibd_segment]["end"] = int(ibd_dict[ibd_segment]["end"])
        ibd_dict[ibd_segment]["LOD"] = float(ibd_dict[ibd_segment]["LOD"])
        ibd_dict[ibd_segment]["cM"] = float(ibd_dict[ibd_segment]["cM"])
    
    else: # neither sample is sample of interest
      continue

  ibd_file.close()

  # demographic dictionary
  dem_file = open(args.dem, 'r')
  dem_file.readline()
  dem_dict = {}
  for line in dem_file:
    line_strip_split = line.strip().split()
    dem_dict[line_strip_split[0]] = {k:v for k,v in zip(["sample", "pop", "super_pop", "gender"], line_strip_split)}

  # run functions
  ibd_variants_dictionary = extract_variants_from_VCF(
    args.vcf,
    args.vcf_id,
    chr = chrom,
    start_bp = start,
    end_bp = end)

  overlap_dict = overlap_IBD_segments_with_phase_sets(ibd_dict, bam_phase_set_dictionary_ps1, ibd_variants_dictionary)

  # write outputs
  os.makedirs(args.output_directory, exist_ok = True)

  ibd_overlap_output_file_path = os.path.join(args.output_directory, args.output_prefix + ".IBD_segment_overlap.tsv")
  ibd_overlap_output_file = open(ibd_overlap_output_file_path, 'w')
  ibd_overlap_output_file.write('\t'.join([str(x) for x in overlap_dict["column_names"]]) + '\n')
  for k,v in overlap_dict.items():
    if k != "column_names":
      ibd_overlap_output_file.write('\t'.join([str(x) for x in v]) + '\n')
  ibd_overlap_output_file.close()

  ancestry_output_file_path = os.path.join(args.output_directory, args.output_prefix + ".IBD_segment_ancestry.tsv")
  ancestry_output_file = open(ancestry_output_file_path, 'w')
  ancestry_overlap_dict = ancestry_overlap(ibd_dict, dem_dict)
  ancestry_output_file.write('\t'.join(["IBD_segment_index", "Chromosome", "IBD_start_position", "IBD_end_position", "Sample1_haplotype", "HBD_status", "Sample2_ID", "Sample2_super_population", "Sample2_population"]) + "\n")
  for k,v in ancestry_overlap_dict.items():
    ancestry_output_file.write(v)
  ancestry_output_file.close()

if __name__ == '__main__':
  main(args)