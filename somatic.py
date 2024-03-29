import os
import pickle

from SomaticHaplotype import *

################################################################################
# functions
################################################################################

def compare_coverage_dictionaries(somatic_variants_dictionary, phasing_dictionary):

  variant_list = list(somatic_variants_dictionary.keys())
  n_variants = len(variant_list)
  
  pairs_dictionary = {}

  for i in range(n_variants - 1):
    var1 = variant_list[i]
    
    if phasing_dictionary[var1] is None:
      continue

    var1_pb = phasing_dictionary[var1][5]
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
      
      if phasing_dictionary[var2] is None:
        continue

      this_pair = var1 + "-" + var2
      var2_pb = phasing_dictionary[var2][5]
      
      if var1_pb != var2_pb:
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

      bx_overlap_00 = list(set(var1_bx_allele0) & set(var2_bx_allele0))
      bx_overlap_01 = list(set(var1_bx_allele0) & set(var2_bx_allele1))
      bx_overlap_10 = list(set(var1_bx_allele1) & set(var2_bx_allele0))
      bx_overlap_11 = list(set(var1_bx_allele1) & set(var2_bx_allele1))

      n_bx_overlap_00 = len(bx_overlap_00)
      n_bx_overlap_01 = len(bx_overlap_01)
      n_bx_overlap_10 = len(bx_overlap_10)
      n_bx_overlap_11 = len(bx_overlap_11)

      if bx_overlap_00 == []:
        bx_overlap_00_csv = 'No_shared_barcodes'
      else:
        bx_overlap_00_csv = ",".join(bx_overlap_00)

      if bx_overlap_01 == []:
        bx_overlap_01_csv = 'No_shared_barcodes'
      else:
        bx_overlap_01_csv = ",".join(bx_overlap_01)

      if bx_overlap_10 == []:
        bx_overlap_10_csv = 'No_shared_barcodes'
      else:
        bx_overlap_10_csv = ",".join(bx_overlap_10)

      if bx_overlap_11 == []:
        bx_overlap_11_csv = 'No_shared_barcodes'
      else:
        bx_overlap_11_csv = ",".join(bx_overlap_11)

      pair_list = [var1, var2, var1_pb, 
      bx_overlap_00_csv, bx_overlap_01_csv, bx_overlap_10_csv, bx_overlap_11_csv,
      n_bx_overlap_00, n_bx_overlap_01, n_bx_overlap_10, n_bx_overlap_11]
      pairs_dictionary[this_pair] = pair_list
      
  return(pairs_dictionary)

def create_coverage_dictionary(variant_key, vcf_variants_dictionary, phase_block_dictionary, somatic_barcodes_dictionary_by_haplotype, phase_block_of_variant):

  if variant_key in vcf_variants_dictionary and len(vcf_variants_dictionary[variant_key]) > 1:
    sys.exit("Variant " + variant_key + " has more than one VCF record.")

  bx_supporting_variants1 = []
  bx_supporting_variants2 = []
  
  if variant_key in vcf_variants_dictionary:
    # get barcode information as reported by longranger in VCF
    bx_supporting_variants1 = return_barcodes_supporting_variant_vcf(variant_key, vcf_variants_dictionary)
  
  if variant_key in somatic_barcodes_dictionary_by_haplotype:
    # if variant is not in VCF, use somatic barcodes reported in --sombx input
    bx_supporting_variants_by_haplotype = return_barcodes_supporting_variant_bam(variant_key, somatic_barcodes_dictionary_by_haplotype)
    bx_supporting_variants2 = []
    for k,v in bx_supporting_variants_by_haplotype.items():
      bx_supporting_variants2.extend(v)

  bx_supporting_variants = list(set(bx_supporting_variants1 + bx_supporting_variants2))

  variants_covered_in_vcf = return_variants_covered_by_barcodes(bx_supporting_variants, phase_block_of_variant, vcf_variants_dictionary)

  if variant_key in vcf_variants_dictionary:
    variant_phased_by_longranger = vcf_variants_dictionary[variant_key][0].return_IsPhasedHeterozygote()
    variant_GT = vcf_variants_dictionary[variant_key][0].return_Genotype()
  else:
    variant_phased_by_longranger = "Not called by longranger"
    variant_GT = "NA"

  coverage_dictionary = {}
  n_REF_H1, n_REF_H2, n_ALT_H1, n_ALT_H2, n_not_phased = 0, 0, 0, 0, 0
  just_barcode_REF_H1, just_barcode_REF_H2, just_barcode_ALT_H1, just_barcode_ALT_H2, just_barcode_not_phased = 0, 0, 0, 0, 0
  for bx in bx_supporting_variants:
    this_bx_supports_somatic_01 = return_allele_supported_by_barcode(bx, variant_key, vcf_variants_dictionary, somatic_barcodes_dictionary_by_haplotype)
    for var in variants_covered_in_vcf:
      this_variant_key = var.return_VariantKey()
      this_coverage_key = bx + "--" + this_variant_key
      coverage_dictionary[this_coverage_key] = [bx, this_variant_key]
      allele_supported_by_barcode = return_allele_supported_by_barcode(bx, this_variant_key, vcf_variants_dictionary, somatic_barcodes_dictionary_by_haplotype)
      haplotype_supported_by_barcode = return_haplotype_supported_by_barcode(bx, this_variant_key, vcf_variants_dictionary, somatic_barcodes_dictionary_by_haplotype)
      
      if var.return_Filter() == []:
        filter_string = "PASS"
      else:
        filter_string = ",".join(var.return_Filter())

      if filter_string in ["PASS", "10X_PHASING_INCONSISTENT"] and var.return_Genotype() != "1|1":
        # not including variants with 10X_ALLELE_FRACTION_FILTER excludes variants with VAF < 0.15
        if this_bx_supports_somatic_01 == "0" and haplotype_supported_by_barcode == "H1":
          n_REF_H1 += 1
        elif this_bx_supports_somatic_01 == "0" and haplotype_supported_by_barcode == "H2":
          n_REF_H2 += 1
        elif this_bx_supports_somatic_01 == "1" and haplotype_supported_by_barcode == "H1":
          n_ALT_H1 += 1
        elif this_bx_supports_somatic_01 == "1" and haplotype_supported_by_barcode == "H2":
          n_ALT_H2 += 1
        elif haplotype_supported_by_barcode == "Not_Phased/Not_Heterozygote/Allele_Not_Matching":
          n_not_phased += 1

      if allele_supported_by_barcode != "No Coverage":
        coverage_dictionary[this_coverage_key].extend([
          phase_block_of_variant,
          allele_supported_by_barcode,
          haplotype_supported_by_barcode, 
          var.return_IsPhasedHeterozygote(),
          var.return_Chromosome(), 
          var.return_Position(), 
          var.return_Genotype(), 
          filter_string])

    if bx + "--" + variant_key not in coverage_dictionary:
      coverage_dictionary[bx + "--" + variant_key] = [bx, variant_key]
      allele_supported_by_barcode = return_allele_supported_by_barcode(bx, variant_key, vcf_variants_dictionary, somatic_barcodes_dictionary_by_haplotype)
      haplotype_supported_by_barcode = return_haplotype_supported_by_barcode(bx, variant_key, vcf_variants_dictionary, somatic_barcodes_dictionary_by_haplotype)
      if bx in somatic_barcodes_dictionary_by_haplotype[variant_key]['ref_H1']:
        n_REF_H1 += 1
      elif bx in somatic_barcodes_dictionary_by_haplotype[variant_key]['ref_H2']:
        n_REF_H2 += 1
      elif bx in somatic_barcodes_dictionary_by_haplotype[variant_key]['alt_H1']:
        n_ALT_H1 += 1
      elif bx in somatic_barcodes_dictionary_by_haplotype[variant_key]['alt_H2']:
        n_ALT_H2 += 1
      elif haplotype_supported_by_barcode == "Not_Phased/Not_Heterozygote/Allele_Not_Matching":
          n_not_phased += 1

      if allele_supported_by_barcode != "No Coverage":
        coverage_dictionary[bx + "--" + variant_key].extend([
          phase_block_of_variant, 
          allele_supported_by_barcode, 
          haplotype_supported_by_barcode, 
          variant_phased_by_longranger, 
          variant_key.split(":")[0],
          variant_key.split(":")[1], 
          variant_GT, 
          "PASS"])
    
    # phase just based on haplotype assigned to a barcode
    if variant_key in somatic_barcodes_dictionary_by_haplotype:
      if bx in somatic_barcodes_dictionary_by_haplotype[variant_key]['ref_H1']:
        just_barcode_REF_H1 += 1
      if bx in somatic_barcodes_dictionary_by_haplotype[variant_key]['ref_H2']:
        just_barcode_REF_H2 += 1
      if bx in somatic_barcodes_dictionary_by_haplotype[variant_key]['alt_H1']:
        just_barcode_ALT_H1 += 1
      if bx in somatic_barcodes_dictionary_by_haplotype[variant_key]['alt_H2']:
        just_barcode_ALT_H2 += 1
      if bx in somatic_barcodes_dictionary_by_haplotype[variant_key]['ref_None']:
        just_barcode_not_phased += 1
      if bx in somatic_barcodes_dictionary_by_haplotype[variant_key]['alt_None']:
        just_barcode_not_phased += 1

  if n_REF_H1 > 0 or n_REF_H2 > 0:
    pct_REF_on_H1 = float(n_REF_H1)/float(n_REF_H1 + n_REF_H2)
    pct_REF_on_H2 = 1 - pct_REF_on_H1
  else:
    pct_REF_on_H1 = 'NA'
    pct_REF_on_H2 = 'NA'

  if n_ALT_H1 > 0 or n_ALT_H2 > 0 :
    pct_ALT_on_H1 = float(n_ALT_H1)/float(n_ALT_H1 + n_ALT_H2)
    pct_ALT_on_H2 = 1 - pct_ALT_on_H1
  else:
    pct_ALT_on_H1 = 'NA'
    pct_ALT_on_H2 = 'NA'

  if just_barcode_REF_H1 > 0 or just_barcode_REF_H2 > 0:
    pct_barcode_REF_on_H1 = float(just_barcode_REF_H1)/float(just_barcode_REF_H1 + just_barcode_REF_H2)
    pct_barcode_REF_on_H2 = 1 - pct_barcode_REF_on_H1
  else:
    pct_barcode_REF_on_H1 = 'NA'
    pct_barcode_REF_on_H2 = 'NA'

  if just_barcode_ALT_H1 > 0 or just_barcode_ALT_H2 > 0:
    pct_barcode_ALT_on_H1 = float(just_barcode_ALT_H1)/float(just_barcode_ALT_H1 + just_barcode_ALT_H2)
    pct_barcode_ALT_on_H2 = 1 - pct_barcode_ALT_on_H1
  else:
    pct_barcode_ALT_on_H1 = 'NA'
    pct_barcode_ALT_on_H2 = 'NA'

  chrom, pos, ref, alt = variant_key.split(":")
  pb_length = phase_block_dictionary[phase_block_of_variant][6]
  variant_phasing = [variant_key, chrom, pos, ref, alt, phase_block_of_variant, pb_length, variant_phased_by_longranger, variant_GT, 
  pct_REF_on_H1, pct_REF_on_H2, pct_ALT_on_H1, pct_ALT_on_H2, n_REF_H1, n_REF_H2, n_ALT_H1, n_ALT_H2, n_not_phased,
  pct_barcode_REF_on_H1, pct_barcode_REF_on_H2, pct_barcode_ALT_on_H1, pct_barcode_ALT_on_H2, just_barcode_REF_H1, just_barcode_REF_H2, just_barcode_ALT_H1, just_barcode_ALT_H2, just_barcode_not_phased]
  
  return([coverage_dictionary, variant_phasing])

def create_phase_block_dictionary(phase_block_filepath):
  pb_file = open(phase_block_filepath, 'r')
  pb_file.readline() # for the header
  pb_dict = {}
  for line in pb_file:
    pb_id, chrom, start, end, length_reads, first_variant_pos, last_variant_pos, length_variants, n_variants_H1, n_variants_H2, n_variants_total = line.strip().split()
    pb_dict[pb_id] = [chrom, start, end, length_reads, first_variant_pos, last_variant_pos, length_variants, n_variants_H1, n_variants_H2, n_variants_total]
  pb_file.close()
  return(pb_dict)

def create_somatic_barcodes_dictionary(somatic_barcodes_filepath):
  sombx_file = open(somatic_barcodes_filepath, 'r')
  sombx_file.readline()
  sombx_dict = {}
  sombx_dict_by_haplotype = {}
  for line in sombx_file:
    variant_key, chromosome, position, ref, alt, ref_barcodes_H1, ref_barcodes_H2, ref_barcodes_None, alt_barcodes_H1, alt_barcodes_H2, alt_barcodes_None = line.strip().split()
    sombx_dict[variant_key] = []
    sombx_dict_by_haplotype[variant_key] = {'ref_H1':[], 'ref_H2':[], 'ref_None':[], 'alt_H1':[], 'alt_H2':[], 'alt_None':[]}
    if ref_barcodes_H1 != 'NA':
      sombx_dict[variant_key].extend(ref_barcodes_H1.split(";"))
      sombx_dict_by_haplotype[variant_key]['ref_H1'] = ref_barcodes_H1.split(";")
    if ref_barcodes_H2 != 'NA':
      sombx_dict[variant_key].extend(ref_barcodes_H2.split(";"))
      sombx_dict_by_haplotype[variant_key]['ref_H2'] = ref_barcodes_H2.split(";")
    if ref_barcodes_None != 'NA':
      sombx_dict[variant_key].extend(ref_barcodes_None.split(";"))
      sombx_dict_by_haplotype[variant_key]['ref_None'] = ref_barcodes_None.split(";")  
    if alt_barcodes_H1 != 'NA':
      sombx_dict[variant_key].extend(alt_barcodes_H1.split(";"))
      sombx_dict_by_haplotype[variant_key]['alt_H1'] = alt_barcodes_H1.split(";")
    if alt_barcodes_H2 != 'NA':
      sombx_dict[variant_key].extend(alt_barcodes_H2.split(";"))
      sombx_dict_by_haplotype[variant_key]['alt_H2'] = alt_barcodes_H2.split(";")
    if alt_barcodes_None != 'NA':
      sombx_dict[variant_key].extend(alt_barcodes_None.split(";"))
      sombx_dict_by_haplotype[variant_key]['alt_None'] = alt_barcodes_None.split(";")

    sombx_dict[variant_key] = list(set(sombx_dict[variant_key]))
    sombx_dict_by_haplotype[variant_key]['ref_H1'] = list(set(sombx_dict_by_haplotype[variant_key]['ref_H1']))
    sombx_dict_by_haplotype[variant_key]['ref_H2'] = list(set(sombx_dict_by_haplotype[variant_key]['ref_H2']))
    sombx_dict_by_haplotype[variant_key]['ref_None'] = list(set(sombx_dict_by_haplotype[variant_key]['ref_None']))
    sombx_dict_by_haplotype[variant_key]['alt_H1'] = list(set(sombx_dict_by_haplotype[variant_key]['alt_H1']))
    sombx_dict_by_haplotype[variant_key]['alt_H2'] = list(set(sombx_dict_by_haplotype[variant_key]['alt_H2']))
    sombx_dict_by_haplotype[variant_key]['alt_None'] = list(set(sombx_dict_by_haplotype[variant_key]['alt_None']))

  return([sombx_dict, sombx_dict_by_haplotype])

def create_somatic_variants_dictionary(maf_filepath, variant_filepath, chrom, start_bp, end_bp):
  if maf_filepath is not None and variant_filepath is None:
    somatic_variants_dictionary = extract_maf_variants(maf_filepath, chrom, start_bp, end_bp)
  elif variant_filepath is not None and maf_filepath is None:
    somatic_variants_dictionary = extract_variant_variants(variant_filepath, chrom, start_bp, end_bp)
  else:
    sys.error("Whoa, somatic module requires one and only one of maf_filepath and variant_filepath to be given.")

  # CHECK THERE IS AT LEAST ONE VARIANT
  if len(somatic_variants_dictionary) > 0:
    return(somatic_variants_dictionary)
  else:
    sys.exit("No somatic variants found. No output generated.")

def extract_maf_variants(maf_filepath, chrom, start_bp, end_bp):
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
        elif this_variant_CHROM == chrom and (start_bp is None or int(this_variant_POS) >= int(start_bp)) and (end_bp is None or int(this_variant_POS) <= int(end_bp)):
          this_variant_key = ":".join([this_variant_CHROM, this_variant_POS, this_variant_REF, this_variant_ALT])
          somatic_variants_dictionary[this_variant_key] = None
        else:
          continue

  maf_file.close()
  return(somatic_variants_dictionary)

def extract_variant_variants(variant_filepath, chrom, start_bp, end_bp):
  variant_file = open(variant_filepath, 'r')
  somatic_variants_dictionary = {}
  for line in variant_file:
    this_variant_CHROM, this_variant_POS, this_variant_REF, this_variant_ALT = line.strip().split(":")
    if this_variant_REF in ["A", "T", "C", "G"] and this_variant_ALT in ["A", "T", "C", "G"]:
      if this_variant_CHROM == chrom and (start_bp is None or int(this_variant_POS) >= int(start_bp)) and (end_bp is None or int(this_variant_POS) <= int(end_bp)):
        somatic_variants_dictionary[line.strip()] = None
      else:
        continue
    else:
      continue
  
  return(somatic_variants_dictionary)

def return_allele_supported_by_barcode(barcode, variant_key, vcf_variants_dictionary, somatic_barcodes_dictionary_by_haplotype):
  if variant_key in vcf_variants_dictionary and len(vcf_variants_dictionary[variant_key]) > 1:
    sys.exit("Variant " + variant_key + " has more than one VCF record.")
  
  if variant_key in vcf_variants_dictionary:
    barcode_supports_this_allele = "No Coverage"
    this_variant = vcf_variants_dictionary[variant_key][0]
    n_alleles = len(this_variant.return_Molecules())
    if barcode in this_variant.return_Molecules()[0]:
      barcode_supports_this_allele = str(0)
    else:
      for i in range(1, n_alleles):
        if barcode in this_variant.return_Molecules()[i]:
          barcode_supports_this_allele = str(1) # treats multiple ALT alleles as same
    
  elif variant_key in somatic_barcodes_dictionary_by_haplotype:
    
    ref_barcodes = somatic_barcodes_dictionary_by_haplotype[variant_key]['ref_H1'] + somatic_barcodes_dictionary_by_haplotype[variant_key]['ref_H2'] + somatic_barcodes_dictionary_by_haplotype[variant_key]['ref_None']
    alt_barcodes = somatic_barcodes_dictionary_by_haplotype[variant_key]['alt_H1'] + somatic_barcodes_dictionary_by_haplotype[variant_key]['alt_H2'] + somatic_barcodes_dictionary_by_haplotype[variant_key]['alt_None']

    if barcode in ref_barcodes:
      barcode_supports_this_allele = str(0)
    elif barcode in alt_barcodes:
      barcode_supports_this_allele = str(1)
    else:
      barcode_supports_this_allele = "No Coverage"
  else: 
    sys.exit("Variant " + variant_key + " not in VCF or somatic barcodes.")
  return(barcode_supports_this_allele)


def return_barcodes_supporting_variant_vcf(variant_key, vcf_variants_dictionary):

  this_variant = vcf_variants_dictionary[variant_key][0]
  if this_variant.return_Genotype() not in ["1/0", "1|0", "0/1", "0|1", "1/1", "1|1", "0/0", "0|0"]:
    sys.exit("Variant " + variant_key + " is not 0/1, 1/0, 0|1, 1|0, 0/0, 0|1, 1/1, or 1|1.")
  else:
    molecules_supporting_this_position = list(this_variant.return_Molecules()[0].keys()) # 0 refers to return_Molecules() dictionary {0:{ref bx}, 1:{alt bx}}
    molecules_supporting_this_position.extend(list(this_variant.return_Molecules()[1].keys())) # 1 refers to return_Molecules() dictionary {0:{ref bx}, 1:{alt bx}}
  return(list(set(molecules_supporting_this_position))) # List of barcodes supporting either REF or ALT allele

def return_barcodes_supporting_variant_bam(variant_key, somatic_barcodes_dictionary):
  barcodes_list = somatic_barcodes_dictionary[variant_key]
  return(barcodes_list)

def return_haplotype_supported_by_barcode(barcode, variant_key, vcf_variants_dictionary, somatic_barcodes_dictionary_by_haplotype):
  if variant_key in vcf_variants_dictionary and len(vcf_variants_dictionary[variant_key]) > 1:
    sys.exit("Variant " + variant_key + " has more than one VCF record.")
  
  elif variant_key in vcf_variants_dictionary:
    this_variant = vcf_variants_dictionary[variant_key][0]
    if this_variant.return_IsPhasedHeterozygote():
      barcode_supports_this_allele = return_allele_supported_by_barcode(barcode, variant_key, vcf_variants_dictionary, somatic_barcodes_dictionary_by_haplotype)
      if this_variant.return_Genotype()[0] == barcode_supports_this_allele:
        barcode_supports_this_haplotype = "H1"
      elif this_variant.return_Genotype()[2] == barcode_supports_this_allele:
        barcode_supports_this_haplotype = "H2"
      elif "No Coverage" == barcode_supports_this_allele:
        barcode_supports_this_haplotype = "No Phased Coverage"
      else:
        barcode_supports_this_haplotype = "Barcode variant not in genotype"
    else:
      barcode_supports_this_haplotype = "Not Phased Heterozygote"
  
  else:
    
    h1_barcodes = somatic_barcodes_dictionary_by_haplotype[variant_key]['ref_H1'] + somatic_barcodes_dictionary_by_haplotype[variant_key]['alt_H1']
    h2_barcodes = somatic_barcodes_dictionary_by_haplotype[variant_key]['ref_H2'] + somatic_barcodes_dictionary_by_haplotype[variant_key]['alt_H2']

    if barcode in h1_barcodes:
      barcode_supports_this_haplotype = "H1"
    elif barcode in h2_barcodes:
      barcode_supports_this_haplotype = "H2"
    else:
      barcode_supports_this_haplotype = "No Phased Coverage"

  if barcode_supports_this_haplotype not in ["H1", "H2"]:
    barcode_supports_this_haplotype = "Not_Phased/Not_Heterozygote/Allele_Not_Matching"

  return(barcode_supports_this_haplotype)

def return_variant_phase_block(base_variant_key, vcf_variants_dictionary, phase_block_dictionary):
  if base_variant_key in vcf_variants_dictionary:
    # if variant was called by longranger and is in longranger VCF
    this_variant_phase_block = vcf_variants_dictionary[base_variant_key][0].return_PhaseBlockID()
    if this_variant_phase_block not in phase_block_dictionary:
      # if the variant's phase block wasn't included in the bam
      this_variant_phase_block = None
  else:
    # if the variant was NOT called by longranger (not in longranger VCF)
    # see if variant falls within range of any existing phase blocks and report that phase block
    chrom, pos, ref, alt = base_variant_key.split(":")
    this_variant_phase_block = None
    for pb_id in phase_block_dictionary.keys():
      if chrom == phase_block_dictionary[pb_id][0] and phase_block_dictionary[pb_id][4] != 'NA' and phase_block_dictionary[pb_id][5] != 'NA':
        # if same chromosome, and phase block has a start and end position that is not missing
        if int(pos) >= int(phase_block_dictionary[pb_id][4]) and int(pos) <= int(phase_block_dictionary[pb_id][5]):
          # if position of variant falls within chr:start-stop range
          this_variant_phase_block = pb_id
          break # break out of current loop
  return(this_variant_phase_block)

def return_variants_covered_by_barcodes(barcode_list, variant_phase_block, vcf_variants_dictionary):

  variant_list = []

  for this_variant_key in vcf_variants_dictionary:
    if len(vcf_variants_dictionary[this_variant_key]) > 1:
      # more than one variant at this key
      pass
    else:
      this_variant = vcf_variants_dictionary[this_variant_key][0]
      if this_variant.return_IsSNP() and this_variant.return_PhaseBlockID() == variant_phase_block:
        if this_variant.return_Genotype() in ["1/0", "1|0", "0/1", "0|1", "1/1", "1|1", "0/0", "0|0"]:
          barcodes_allele0 = list(this_variant.return_Molecules()[0].keys())
          barcodes_allele1 = list(this_variant.return_Molecules()[1].keys())

          bx_overlap_allele0 = [ bx in barcodes_allele0 for bx in barcode_list ]
          bx_overlap_allele1 = [ bx in barcodes_allele1 for bx in barcode_list ]

          if any(bx_overlap_allele0) or any(bx_overlap_allele1):
            variant_list.append(this_variant)

  return(variant_list)

def somatic_variants_per_phase_block(phase_block_dictionary, phasing_dictionary):
  pb_id_list = list(phase_block_dictionary.keys())
  for pb_id in pb_id_list:
    phase_block_dictionary[pb_id].append(0)

  for var in phasing_dictionary.keys():
    if phasing_dictionary[var] is not None:
      this_var_pb_id = phasing_dictionary[var][5]
      phase_block_dictionary[this_var_pb_id][-1] += 1

  return(phase_block_dictionary)

def write_phasing_dictionary(phasing_dictionary, output_file_path):
  
  output_file = open(output_file_path, "w")
  output_file.write('\t'.join(["Variant", "Chromosome", "Position", "Reference", "Alternate", "Phase_Block", "Phase_Block_Length", "Variant_Phased_by_longranger", "Genotype", 
  "pct_REF_on_H1", "pct_REF_on_H2", "pct_ALT_on_H1", "pct_ALT_on_H2", "n_REF_H1", "n_REF_H2", "n_ALT_H1", "n_ALT_H2", "n_not_phased",
  "pct_barcode_REF_on_H1", "pct_barcode_REF_on_H2", "pct_barcode_ALT_on_H1", "pct_barcode_ALT_on_H2", "barcode_REF_H1", "barcode_REF_H2", "barcode_ALT_H1", "barcode_ALT_H2", "just_barcode_not_phased"]) + '\n')

  for var in sorted(phasing_dictionary.keys()):
    if phasing_dictionary[var] is not None:
      output_file.write('\t'.join([str(x) for x in phasing_dictionary[var]]) + '\n')

  output_file.close()

def write_somatic_variants_dictionary(somatic_variants_dictionary, output_file_path):

  output_file = open(output_file_path, "w")
  output_file.write('\t'.join(["Barcode", "Variant", "Phase_Block", "Allele", "Haplotype", 
    "Phased_Heterozygote", "Chromosome", "Position", "Genotype", "Filter", "Somatic_Variant"]) + '\n')
  for var in somatic_variants_dictionary:
    if somatic_variants_dictionary[var] is not None:
      for bx_pos in somatic_variants_dictionary[var]:
        if len(somatic_variants_dictionary[var][bx_pos]) == 2:
          continue
        else:
          print_these_combinations = '\t'.join([str(x) for x in somatic_variants_dictionary[var][bx_pos]]) + '\t' + var + '\n'
          output_file.write(print_these_combinations)
  output_file.close()

def write_somatic_variants_per_phase_block(phase_block_dictionary, output_file_path):
  output_file = open(output_file_path, 'w')
  output_file.write('\t'.join(["pb_id", "chrom", "start", "end", "length_reads", "first_variant_pos", "last_variant_pos", "length_variants", "n_variants_H1", "n_variants_H2", "n_variants_total", "n_somatic_variants"]) + '\n')
  for pb_id in phase_block_dictionary.keys():
    output_file.write(pb_id + '\t' + '\t'.join([str(x) for x in phase_block_dictionary[pb_id]]) + '\n')
  output_file.close()

def write_variant_pairs_dictionary(variant_comparison, output_file_path):

  output_file = open(output_file_path, "w")
  output_file.write('\t'.join(["Variant1", "Variant2", "Phase_Block", "bx_overlap_00", "bx_overlap_01", "bx_overlap_10", "bx_overlap_11", "n_bx_overlap_00", "n_bx_overlap_01", "n_bx_overlap_10", "n_bx_overlap_11"]) + '\n')
  
  for pair in variant_comparison:
    output_file.write('\t'.join([str(x) for x in variant_comparison[pair]]) + '\n')

  output_file.close()

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

  # phase block summary dictionary
  phase_block_dictionary = create_phase_block_dictionary(phase_block_filepath = args.sum)

  # somatic barcodes dictionary
  if args.sombx is not None:
    somatic_barcodes_dictionary, somatic_barcodes_dictionary_by_haplotype = create_somatic_barcodes_dictionary(somatic_barcodes_filepath = args.sombx)
  else:
    somatic_barcodes_dictionary = {}
    somatic_barcodes_dictionary_by_haplotype = {}

  # use these somatic variants as basis of analysis
  # somatic_variants_dictionary[variant_key] = []
  # keys: variant_key
  # values: empty list
  somatic_variants_dictionary = create_somatic_variants_dictionary(maf_filepath = args.maf, variant_filepath = args.variant, chrom = chrom, start_bp = start, end_bp = end)

  # create coverage dictionary for each somatic variant
  phasing_dictionary = {}
  for variant_key in somatic_variants_dictionary.keys():
    this_pb_id = return_variant_phase_block(variant_key, vcf_variants_dictionary, phase_block_dictionary)
    if this_pb_id is None: # variant does not fall within an existing phase block
      somatic_variants_dictionary[variant_key], phasing_dictionary[variant_key] = [None, None]
    else: # variant is within the range of an existing phase block
      somatic_variants_dictionary[variant_key], phasing_dictionary[variant_key] = create_coverage_dictionary(variant_key, vcf_variants_dictionary, phase_block_dictionary, somatic_barcodes_dictionary_by_haplotype, this_pb_id)

  # determine relationship of each pair of variants in same phase block
  # produces 00, 01, 10, 11 shared barcode patterns for pairs of variants
  variant_comparison = compare_coverage_dictionaries(somatic_variants_dictionary, phasing_dictionary)

  # somatic variants per phase block
  phase_block_dictionary_with_n_somatic = somatic_variants_per_phase_block(phase_block_dictionary, phasing_dictionary)

  # write output and close files
  os.makedirs(args.output_directory, exist_ok = True)
  output_file_path_somatic_variants_dictionary = os.path.join(args.output_directory, args.output_prefix + ".barcodes_variants.tsv")
  output_file_path_phasing_dictionary = os.path.join(args.output_directory, args.output_prefix + ".phasing_variants.tsv")
  output_file_path_variant_pairs = os.path.join(args.output_directory, args.output_prefix + ".variant_pairs.tsv")
  output_file_path_somatic_per_phase_block = os.path.join(args.output_directory, args.output_prefix + ".somatic_per_phase_block.tsv")
  write_somatic_variants_dictionary(somatic_variants_dictionary, output_file_path_somatic_variants_dictionary)
  write_phasing_dictionary(phasing_dictionary, output_file_path_phasing_dictionary)
  write_variant_pairs_dictionary(variant_comparison, output_file_path_variant_pairs)
  write_somatic_variants_per_phase_block(phase_block_dictionary_with_n_somatic, output_file_path_somatic_per_phase_block)
  
  # write results to output file here
  
if __name__ == '__main__':
  main(args)