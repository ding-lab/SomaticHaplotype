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
    if this_variant.return_Genotype() not in ["1/0", "1|0", "0/1", "0|1", "1/1", "1|1", "0/0", "0|0"]:
      sys.exit("Variant " + variant_key + " is not 0/1, 1/0, 0|1, 1|0, 0/0, 0|1, 1/1, or 1|1.")
    else:
      molecules_supporting_this_position = list(this_variant.return_Molecules()[0].keys()) # 0 refers to return_Molecules() dictionary {0:{ref bx}, 1:{alt bx}}
      molecules_supporting_this_position.extend(list(this_variant.return_Molecules()[1].keys())) # 1 refers to return_Molecules() dictionary {0:{ref bx}, 1:{alt bx}}
      return(molecules_supporting_this_position) # List of barcodes supporting either REF or ALT allele

def return_variants_covered_by_barcodes(barcode_list, base_variant_key, vcf_variants_dictionary):

  variant_list = []
  this_phase_set = vcf_variants_dictionary[base_variant_key][0].return_PhaseSetID()

  for this_variant_key in vcf_variants_dictionary:
    if len(vcf_variants_dictionary[this_variant_key]) > 1:
      print("Variant " + this_variant_key + " has more than one VCF record.")
    else:
      this_variant = vcf_variants_dictionary[this_variant_key][0]
      if this_variant.return_IsSNP() and this_variant.return_PhaseSetID() == this_phase_set:
        if this_variant.return_Genotype() not in ["1/0", "1|0", "0/1", "0|1", "1/1", "1|1", "0/0", "0|0"]:
          continue
        else:
          barcodes_allele0 = list(this_variant.return_Molecules()[0].keys())
          barcodes_allele1 = list(this_variant.return_Molecules()[1].keys())

        bx_overlap_allele0 = [ bx in barcodes_allele0 for bx in barcode_list ]
        bx_overlap_allele1 = [ bx in barcodes_allele1 for bx in barcode_list ]

        if any(bx_overlap_allele0) or any(bx_overlap_allele1):
          variant_list.append(this_variant)

  return(variant_list)

def return_allele_supported_by_barcode(barcode, variant_key, vcf_variants_dictionary):
  if variant_key not in vcf_variants_dictionary:
    sys.exit("Variant " + variant_key + " not in vcf_variants_dictionary.")
  elif len(vcf_variants_dictionary[variant_key]) > 1:
    sys.exit("Variant " + variant_key + " has more than one VCF record.")
  else:
    barcode_supports_this_allele = "No Coverage"
    this_variant = vcf_variants_dictionary[variant_key][0]
    n_alleles = len(this_variant.return_Molecules())
    for i in range(n_alleles):
      if barcode in this_variant.return_Molecules()[i]:
        if i == 0:
          barcode_supports_this_allele = "REF"
        if i == 1:
          barcode_supports_this_allele = "ALT"
    return(barcode_supports_this_allele)

def return_haplotype_supported_by_barcode(barcode, variant_key, vcf_variants_dictionary):
  if variant_key not in vcf_variants_dictionary:
    sys.exit("Variant " + variant_key + " not in vcf_variants_dictionary.")
  elif len(vcf_variants_dictionary[variant_key]) > 1:
    sys.exit("Variant " + variant_key + " has more than one VCF record.")
  else:
    this_variant = vcf_variants_dictionary[variant_key][0]
    if this_variant.return_IsPhasedHeterozygote():
      barcode_supports_this_allele = return_allele_supported_by_barcode(barcode, variant_key, vcf_variants_dictionary)
      if this_variant.return_Genotype()[0] == barcode_supports_this_allele:
        barcode_supports_this_haplotype = "H1"
      elif this_variant.return_Genotype()[2] == barcode_supports_this_allele:
        barcode_supports_this_haplotype = "H2"
      elif "No Coverage" == barcode_supports_this_allele:
        barcode_supports_this_haplotype = "No Coverage"
    else:
      barcode_supports_this_haplotype = "Not Phased Heterozygote"
    return(barcode_supports_this_haplotype)

def create_coverage_dictionary(variant_key, vcf_variants_dictionary):

  bx_supporting_variants = return_barcodes_supporting_variant(variant_key, vcf_variants_dictionary)
  variants_covered = return_variants_covered_by_barcodes(bx_supporting_variants, variant_key, vcf_variants_dictionary)
  
  coverage_dictionary = {}
  for bx in bx_supporting_variants:
    for var in variants_covered:
      this_variant_key = var.return_VariantKey()
      this_coverage_key = bx + "--" + this_variant_key
      coverage_dictionary[this_coverage_key] = [bx, this_variant_key]
      allele_supported_by_barcode = return_allele_supported_by_barcode(bx, this_variant_key, vcf_variants_dictionary)
      if not allele_supported_by_barcode == "NA":
        haplotype_supported_by_barcode = return_haplotype_supported_by_barcode(bx, this_variant_key, vcf_variants_dictionary)
      else:
        haplotype_supported_by_barcode = "NA"
      if var.return_Filter() == []:
        filter_string = "PASS"
      else:
        filter_string = ", ".join(var.return_Filter())
      coverage_dictionary[this_coverage_key].extend([
        allele_supported_by_barcode, 
        haplotype_supported_by_barcode, 
        var.return_IsPhasedHeterozygote(),
        var.return_Chromosome(), 
        var.return_Position(), 
        var.return_Genotype(), 
        filter_string])
  return(coverage_dictionary)

def write_coverage_dictionary(coverage_dictionary, output_file_path):

  output_file = open(output_file_path, "w")
  output_file.write('\t'.join(["Barcode", "Variant", "Allele", "Haplotype", 
    "Phased_Heterozygote", "Chromosome", "Position", "Genotype", "Filter"]) + '\n')
  for x in coverage_dictionary:
    output_file.write('\t'.join([str(x) for x in coverage_dictionary[x]]) + '\n')

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
  write_coverage_dictionary(coverage_dictionary, output_file_path)
  
  # write results to output file here
  
if __name__ == '__main__':
  main(args)
