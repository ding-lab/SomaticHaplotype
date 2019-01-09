import os
import pickle

from SomaticHaplotype import *

################################################################################
# functions
################################################################################

def return_barcodes_supporting_variant(variant_key, vcf_variants_dictionary):

  if variant_key not in vcf_variants_dictionary:
    sys.exit("Variant " + variant_key + " not in vcf_variants_dictionary.")
  elif len(vcf_variants_dictionary[variant_key]) > 1:
    sys.exit("Variant " + variant_key + " has more than one VCF record.")
  else:
    this_variant = vcf_variants_dictionary[variant_key][0]
    if this_variant.return_Genotype() in ["1/0", "1|0"]:
      variant_allele = 0
    elif this_variant.return_Genotype() in ["0/1", "0|1"]:
      variant_allele = 1
    else:
      sys.exit("Variant " + variant_key + " is not 0/1, 1/0, 0|1, or 1|0.")
    other_allele = (variant_allele + 1) % 2
    molecules_supporting_this_position = list(this_variant.return_Molecules()[variant_allele].keys())
    #molecules_supporting_this_position.extend(list(this_variant.return_Molecules()[other_allele].keys()))
    return(molecules_supporting_this_position)

def return_variants_covered_by_barcodes(barcode_list, base_variant_key, vcf_variants_dictionary):

  variant_list = []
  this_phase_set = vcf_variants_dictionary[base_variant_key][0].return_PhaseSetID()

  for this_variant_key in vcf_variants_dictionary:
    if len(vcf_variants_dictionary[this_variant_key]) > 1:
      print("Variant " + this_variant_key + " has more than one VCF record.")
    else:
      this_variant = vcf_variants_dictionary[this_variant_key][0]
      #if this_variant.return_VariantKey() == base_variant or (this_variant.determine_if_phased_heterozygote() and this_variant.return_Genotype() in ["0|1", "1|0"] and this_variant.return_IsSNP()):
      if this_variant.return_IsSNP() and this_variant.return_PhaseSetID() == this_phase_set:
        barcodes_allele0 = list(this_variant.return_Molecules()[0].keys())
        barcodes_allele1 = list(this_variant.return_Molecules()[1].keys())

        bx_overlap_allele0 = [ bx in barcodes_allele0 for bx in barcode_list]
        bx_overlap_allele1 = [ bx in barcodes_allele1 for bx in barcode_list]

        if any(bx_overlap_allele0) or any(bx_overlap_allele1):
          variant_list.append(this_variant)

  return(variant_list)

def return_coverage_vector_for_variant(variant_key, barcode_list, vcf_variants_dictionary):

  this_variant = vcf_variants_dictionary[variant_key]
  coverage_vector = []

  if this_variant.return_Genotype() == "0|1":
    variant_allele = 1
  elif this_variant.return_Genotype == "1|0":
    variant_allele = 0
  else:
    sys.exit("Variant allele of " + variant + " is not 0 or 1.")

  barcodes_allele0 = this_variant.return_Molecules()[0].keys()
  barcdoes_allele1 = this_variant.return_Molecules()[1].keys()

  for bx in barcode_list:
    if bx in barcodes_allele0 and variant_allele == 0:
      alt_supported = True
    if bx in barcodes_allele0 and variant_allele == 1:
      ref_supported = True
    if bx in barcodes_allele1 and variant_allele == 1:
      alt_supported = True
    if bx in barcodes_allele1 and variant_allele == 0:
      ref_supported = True
    if alt_supported and ref_supported:
      sys.exit("Both REF and ALT of " + variant_key + " supported by bx " + bx)
    elif alt_supported and not ref_supported:
      coverage_vector.append(1)
    elif ref_supported and not alt_supported:
      coverage_vector.append(0)
    else:
      coverage_vector.append(-1)

  if len(coverage_vector) == len(barcode_list):
    return(coverage_vector)
  else:
    sys.exit("Coverage vector is different length than barcode vector.")

def create_coverage_dictionary(variant_key, vcf_variants_dictionary):

  bx_supporting_variants = return_barcodes_supporting_variant(variant_key, vcf_variants_dictionary)
  variants_covered = return_variants_covered_by_barcodes(bx_supporting_variants, variant_key, vcf_variants_dictionary)
  #coverage_dict = {}
  #coverage_dict["barcodes"].extend(bx_supporting_variants)
  
  #for variant in variants_covered:
  #  coverage_dict[variant] = return_coverage_vector_for_variant()

  #return(coverage_dict)

def write_coverage_dictionary(coverage_dictionary, output_file_path):

  variants = sorted([ x for x in coverage_dictionary.keys() if x != "barcodes" ])
  column_headers = ["barcodes"]
  column_headers.extend(variants)
  n_bx = len(coverage_dictionary["barcodes"])

  output_file = open(output_file_path, "w")

  #output_file.write("\t".join(column_headers) + "\n")
  print("\t".join(column_headers))

  for i in range(n_bx):
    write_line = []
    for var in column_headers:
      write_line.append(coverage_dictionary[var][i])
    #output_file.write("\t".join(write_line) + "\n")
    print("\t".join(write_line))

  output_file.close()


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

  # use this variant as basis of analysis
  variant_id = args.variant

  # uses functions here
  coverage_dictionary = create_coverage_dictionary(variant_id, vcf_variants_dictionary)

  # write output and close files
  os.makedirs(args.output_directory, exist_ok = True)
  output_file_path = os.path.join(args.output_directory, args.output_prefix + ".FILENAME.tsv")
  #write_coverage_dictionary(coverage_dictionary, output_file_path)  
  
  # write results to output file here
  
if __name__ == '__main__':
  main(args)
