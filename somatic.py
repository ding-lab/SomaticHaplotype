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
        barcode_supports_this_allele = str(i)
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

def create_coverage_dictionary(variant_key, vcf_variants_dictionary, phase_set_dictionary):

  if variant_key not in vcf_variants_dictionary:
    sys.exit("Variant " + variant_key + " not in vcf_variants_dictionary.")
  elif len(vcf_variants_dictionary[variant_key]) > 1:
    sys.exit("Variant " + variant_key + " has more than one VCF record.")

  bx_supporting_variants = return_barcodes_supporting_variant(variant_key, vcf_variants_dictionary)
  variants_covered = return_variants_covered_by_barcodes(bx_supporting_variants, variant_key, vcf_variants_dictionary)
  
  n_REF_H1, n_REF_H2, n_ALT_H1, n_ALT_H2, n_not_phased_heterozygote = 0, 0, 0, 0, 0
  phase_set_of_variant = vcf_variants_dictionary[variant_key][0].return_PhaseSetID()
  variant_phased_by_longranger = vcf_variants_dictionary[variant_key][0].return_IsPhasedHeterozygote()
  variant_GT = vcf_variants_dictionary[variant_key][0].return_Genotype()

  coverage_dictionary = {}
  for bx in bx_supporting_variants:
    this_bx_supports_somatic_01 = return_allele_supported_by_barcode(bx, variant_key, vcf_variants_dictionary)
    for var in variants_covered:
      this_variant_key = var.return_VariantKey()
      this_coverage_key = bx + "--" + this_variant_key
      coverage_dictionary[this_coverage_key] = [bx, this_variant_key]
      allele_supported_by_barcode = return_allele_supported_by_barcode(bx, this_variant_key, vcf_variants_dictionary)
      haplotype_supported_by_barcode = return_haplotype_supported_by_barcode(bx, this_variant_key, vcf_variants_dictionary)
      
      if var.return_Filter() == []:
        filter_string = "PASS"
      else:
        filter_string = ", ".join(var.return_Filter())

      if filter_string in ["PASS", "10X_PHASING_INCONSISTENT"] and var.return_Genotype() != "1|1":
        if this_bx_supports_somatic_01 == "0" and haplotype_supported_by_barcode == "H1":
          n_REF_H1 += 1
        elif this_bx_supports_somatic_01 == "0" and haplotype_supported_by_barcode == "H2":
          n_REF_H2 += 1
        elif this_bx_supports_somatic_01 == "1" and haplotype_supported_by_barcode == "H1":
          n_ALT_H1 += 1
        elif this_bx_supports_somatic_01 == "1" and haplotype_supported_by_barcode == "H2":
          n_ALT_H2 += 1
        elif haplotype_supported_by_barcode == "Not Phased Heterozygote":
          n_not_phased_heterozygote += 1

      if allele_supported_by_barcode != "No Coverage":
        coverage_dictionary[this_coverage_key].extend([
          phase_set_of_variant,
          allele_supported_by_barcode,
          haplotype_supported_by_barcode, 
          var.return_IsPhasedHeterozygote(),
          var.return_Chromosome(), 
          var.return_Position(), 
          var.return_Genotype(), 
          filter_string])

  pct_REF_on_H1 = float(n_REF_H1)/float(n_REF_H1 + n_REF_H2)
  pct_REF_on_H2 = 1 - pct_REF_on_H1
  pct_ALT_on_H1 = float(n_ALT_H1)/float(n_ALT_H1 + n_ALT_H2)
  pct_ALT_on_H2 = 1 - pct_ALT_on_H1
  chrom, pos, ref, alt = variant_key.split(":")
  ps_length = phase_set_dictionary[phase_set_of_variant]
  variant_phasing = [variant_key, chrom, pos, ref, alt, phase_set_of_variant, ps_length, variant_phased_by_longranger, variant_GT, 
  pct_REF_on_H1, pct_REF_on_H2, pct_ALT_on_H1, pct_ALT_on_H2, n_REF_H1, n_REF_H2, n_ALT_H1, n_ALT_H2, n_not_phased_heterozygote]
  
  return(coverage_dictionary, variant_phasing)

def compare_coverage_dictionaries(somatic_variants_dictionary, phasing_dictionary):

  variant_list = list(somatic_variants_dictionary.keys())
  n_variants = len(variant_list)
  
  pairs_dictionary = {}

  for i in range(n_variants - 1):
    var1 = variant_list[i]
    var1_ps = phasing_dictionary[var1][5]
    var1_bx_allele0 = []
    var1_bx_allele1 = []
    for bx_key1 in somatic_variants_dictionary[var1].keys():
      if len(somatic_variants_dictionary[var1][bx_key1]) == 2:
        continue # bx has no coverage at var1
      if somatic_variants_dictionary[var1][bx_key1][1] == var1:
        allele_supported_by_barcode = somatic_variants_dictionary[var1][bx_key1][3]
        if allele_supported_by_barcode == "0":
          var1_bx_allele0.append(bx_key1.split("--")[0])
        elif allele_supported_by_barcode == "1":
          var1_bx_allele1.append(bx_key1.split("--")[0])
        else:
          sys.exit("Whoa, variant1 bx key " + bx_key1 + " supports non-01 allele.")
    for j in range(i + 1, n_variants):
      var2 = variant_list[j]
      this_pair = var1 + "-" + var2
      var2_ps = phasing_dictionary[var2][5]
      if var1_ps != var2_ps:
        continue
      var2_bx_allele0 = []
      var2_bx_allele1 = []
      
      for bx_key2 in somatic_variants_dictionary[var2].keys():
        if len(somatic_variants_dictionary[var2][bx_key2]) == 2:
          continue # bx has no coverage at var2
        if somatic_variants_dictionary[var2][bx_key2][1] == var2:
          allele_supported_by_barcode = somatic_variants_dictionary[var2][bx_key2][3]
          if allele_supported_by_barcode == "0":
            var2_bx_allele0.append(bx_key2.split("--")[0])
          elif allele_supported_by_barcode == "1":
            var2_bx_allele1.append(bx_key2.split("--")[0])
          else:
            sys.exit("Whoa, variant2 bx key ", bx_key2, " supports non-01 allele.")

      #print(var1_bx_allele0)
      #print(var1_bx_allele1)
      #print(var2_bx_allele0)
      #print(var2_bx_allele1)

      bx_overlap_00 = list(set(var1_bx_allele0) & set(var2_bx_allele0))
      bx_overlap_01 = list(set(var1_bx_allele0) & set(var2_bx_allele1))
      bx_overlap_10 = list(set(var1_bx_allele1) & set(var2_bx_allele0))
      bx_overlap_11 = list(set(var1_bx_allele1) & set(var2_bx_allele1))

      bx_overlap_00_csv = ",".join(bx_overlap_00)
      bx_overlap_01_csv = ",".join(bx_overlap_01)
      bx_overlap_10_csv = ",".join(bx_overlap_10)
      bx_overlap_11_csv = ",".join(bx_overlap_11)

      n_bx_overlap_00 = len(bx_overlap_00)
      n_bx_overlap_01 = len(bx_overlap_01)
      n_bx_overlap_10 = len(bx_overlap_10)
      n_bx_overlap_11 = len(bx_overlap_11)

      pair_list = [var1, var2, var1_ps, 
      bx_overlap_00_csv, bx_overlap_01_csv, bx_overlap_10_csv, bx_overlap_11_csv,
      n_bx_overlap_00, n_bx_overlap_01, n_bx_overlap_10, n_bx_overlap_11]
      pairs_dictionary[this_pair] = pair_list
      
  return(pairs_dictionary)

def write_somatic_variants_dictionary(somatic_variants_dictionary, output_file_path):

  output_file = open(output_file_path, "w")
  output_file.write('\t'.join(["Barcode", "Variant", "Phase_Set", "Allele", "Haplotype", 
    "Phased_Heterozygote", "Chromosome", "Position", "Genotype", "Filter", "Somatic_Variant"]) + '\n')
  print_these_combinations = {}
  for var in somatic_variants_dictionary:
    for bx_pos in somatic_variants_dictionary[var]:
      if len(somatic_variants_dictionary[var][bx_pos]) == 2:
        continue
      elif bx_pos not in print_these_combinations:
        print_these_combinations[bx_pos] = '\t'.join([str(x) for x in somatic_variants_dictionary[var][bx_pos]]) + '\t' + var + '\n'
      else:
        continue
  for bx_pos in sorted(print_these_combinations.keys()):
    output_file.write(print_these_combinations[bx_pos])

  output_file.close()

def write_phasing_dictionary(phasing_dictionary, output_file_path):
  
  output_file = open(output_file_path, "w")
  output_file.write('\t'.join(["Variant", "Chromosome", "Position", "Reference", "Alternate", "Phase_Set", "Phase_Set_Length", "Variant_Phased_by_longranger", "Genotype", 
  "pct_REF_on_H1", "pct_REF_on_H2", "pct_ALT_on_H1", "pct_ALT_on_H2", "n_REF_H1", "n_REF_H2", "n_ALT_H1", "n_ALT_H", "n_not_phased_heterozygote"]) + '\n')

  for var in sorted(phasing_dictionary.keys()):
    output_file.write('\t'.join([str(x) for x in phasing_dictionary[var]]) + '\n')

  output_file.close()

def create_somatic_variants_dictionary(maf_filepath, variant_filepath):
  if maf_filepath is not None and variant_filepath is None:
    somatic_variants_dictionary = extract_maf_variants(maf_filepath)
  elif variant_filepath is not None and maf_filepath is None:
    somatic_variants_dictionary = extract_variant_variants(variant_filepath)
  else:
    sys.error("Whoa, either maf_filepath and variant_filepath are both given or neither maf_filepath nor variant_filepath is given.")
  return(somatic_variants_dictionary)

def extract_maf_variants(maf_filepath):
  maf_file = open(maf_filepath, 'r')
  somatic_variants_dictionary = {}
  for line in maf_file:
    if line.startswith("#"):
      continue
    else:
      if line.startswith("Hugo_Symbol"):
        column_names = line.strip().split("\t")
      else:
        variant_values = line.strip().split("\t")
        this_variant_dict = {k:v for k,v in zip(column_names, variant_values)}
        this_variant_CHROM = this_variant_dict["Chromosome"]
        this_variant_POS = this_variant_dict["Start_Position"]
        this_variant_REF = this_variant_dict["Reference_Allele"]
        this_variant_ALT = this_variant_dict["Tumor_Seq_Allele2"]
        this_variant_type = this_variant_dict["Variant_Type"]
        if this_variant_type != "SNP":
          continue
        else:
          this_variant_key = ":".join([this_variant_CHROM, this_variant_POS, this_variant_REF, this_variant_ALT])
          somatic_variants_dictionary[this_variant_key] = None
  maf_file.close()
  return(somatic_variants_dictionary)

def extract_variant_variants(variant_filepath):
  variant_file = open(variant_filepath, 'r')
  somatic_variants_dictionary = {}
  for line in variant_file:
    somatic_variants_dictionary[line.strip()] = None
  return(somatic_variants_dictionary)

def create_phase_set_dictionary(phase_set_filepath):
  ps_file = open(phase_set_filepath, 'r')
  ps_file.readline() # for the header
  ps_dict = {}
  for line in ps_file:
    ps_id, chrom, start, end, length_reads, first_variant_pos, last_variant_pos, length_variants, n_variants_H1, n_variants_H2, n_variants_total = line.strip().split()
    ps_dict[ps_id] = length_variants
  ps_file.close()
  return(ps_dict)

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

  # phase set summary dictionary
  phase_set_dictionary = create_phase_set_dictionary(phase_set_filepath = args.sum)

  # use these somatic variants as basis of analysis
  somatic_variants_dictionary = create_somatic_variants_dictionary(maf_filepath = args.maf, variant_filepath = args.variant)

  # create coverage dictionary for each somatic variant
  phasing_dictionary = {}
  for variant_key in somatic_variants_dictionary.keys():
    if variant_key in vcf_variants_dictionary:
      somatic_variants_dictionary[variant_key], phasing_dictionary[variant_key] = create_coverage_dictionary(variant_key, vcf_variants_dictionary, phase_set_dictionary)
    else:
      sys.exit("Variant " + variant_key + "not called by longranger.")

  # determine relationship of each pair of variants in same phase set
  variant_comparison = compare_coverage_dictionaries(somatic_variants_dictionary, phasing_dictionary)
  print(variant_comparison)

  sys.exit()

  # write output and close files
  os.makedirs(args.output_directory, exist_ok = True)
  output_file_path_somatic_variants_dictionary = os.path.join(args.output_directory, args.output_prefix + ".barcodes_variants.tsv")
  output_file_path_phasing_dictionary = os.path.join(args.output_directory, args.output_prefix + ".phasing_variants.tsv")
  write_somatic_variants_dictionary(somatic_variants_dictionary, output_file_path_somatic_variants_dictionary)
  write_phasing_dictionary(phasing_dictionary, output_file_path_phasing_dictionary)
  
  # write results to output file here
  
if __name__ == '__main__':
  main(args)
