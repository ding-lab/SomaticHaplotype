import os
import pickle
import pysam
import vcf

from SomaticHaplotype import *

################################################################################
# bam functions
################################################################################

def extract_read_info(read):
  # Given a read from a bam, extract the phase block if it exists, else NA
  # Also report the position of the read and read quality metrics
  if read.is_duplicate or read.is_qcfail or read.is_secondary or not read.is_proper_pair:
    return("read is bad quality")
  elif read.has_tag("PB"): # read is part of phase block
    tags_dict = {x:y for (x,y) in read.get_tags()}
    if "MI" not in tags_dict:
      tags_dict["MI"] = None
    return(tags_dict)
  else: # read is not part of a phase block
    return("read is not part of phase block")

def extract_phase_blocks_from_bam(bam_filename, chr = None, start_bp = None, end_bp = None):
  
  samfile = pysam.AlignmentFile(bam_filename, "rb")
  
  phase_block_dict = {"n_total_reads" : 0, "n_reads_bad_quality" : 0,
    "n_reads_good_quality" : 0, "n_reads_phased" : 0, "n_reads_not_phased" : 0,
    "phase_blocks" : {} }
  
  for read in samfile.fetch(chr, start_bp, end_bp):
    phase_block_dict["n_total_reads"] += 1
    read_info = extract_read_info(read)
    if read_info == "read is bad quality":
      phase_block_dict["n_reads_bad_quality"] += 1
    elif read_info == "read is not part of phase block":
      phase_block_dict["n_reads_not_phased"] += 1
      phase_block_dict["n_reads_good_quality"] += 1
    else:
      
      phase_block_dict["n_reads_phased"] += 1
      phase_block_dict["n_reads_good_quality"] += 1
      
      pb_id = read.reference_name + ":" + str(read_info["PB"])

      if pb_id in phase_block_dict["phase_blocks"]:
        
        if int(read.reference_start) < phase_block_dict["phase_blocks"][pb_id].return_Start():
          phase_block_dict["phase_blocks"][pb_id].update_Start(read.reference_start)
        if int(read.reference_end) > phase_block_dict["phase_blocks"][pb_id].return_End():
          phase_block_dict["phase_blocks"][pb_id].update_End(read.reference_end)
        
        phase_block_dict["phase_blocks"][pb_id].add_SingleEndRead()

        if read_info["HP"] == 1:
          phase_block_dict["phase_blocks"][pb_id].add_SupportH1()
          phase_block_dict["phase_blocks"][pb_id].add_MoleculeH1(read_info["MI"])
        elif read_info["HP"] == 2:
          phase_block_dict["phase_blocks"][pb_id].add_SupportH2()
          phase_block_dict["phase_blocks"][pb_id].add_MoleculeH2(read_info["MI"])
        else:
          sys.exit("Whoa no H1 or H2 support for\n" + str(read))

      else:
      
        phase_block_dict["phase_blocks"][pb_id] = PhaseBlock(pb_id = pb_id, 
          chromosome = read.reference_name, 
          start_bp = int(read.reference_start), 
          end_bp = int(read.reference_end))

        phase_block_dict["phase_blocks"][pb_id].add_SingleEndRead()

        if read_info["HP"] == 1:
          phase_block_dict["phase_blocks"][pb_id].add_SupportH1()
          phase_block_dict["phase_blocks"][pb_id].add_MoleculeH1(read_info["MI"])
        elif read_info["HP"] == 2:
          phase_block_dict["phase_blocks"][pb_id].add_SupportH2()
          phase_block_dict["phase_blocks"][pb_id].add_MoleculeH2(read_info["MI"])
        else:
          sys.exit("Whoa no H1 or H2 support for\n" + str(read))
        
  samfile.close()
  return(phase_block_dict)

################################################################################
# VCF functions
################################################################################


def extract_variants_from_VCF(vcf_filename, sample_id, chr = None, start_bp = None, end_bp = None):
  # build a dictionary of phase blocks present in VCF
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
# Combined bam phase block and VCF variants functions
################################################################################

def add_variants_to_phase_blocks(bam_phase_block_dictionary, vcf_variants_dictionary):
  bam_phase_block_dictionary["phase_blocks"]["variant_not_phased_heterozygote"] = {}
  bam_phase_block_dictionary["phase_blocks"]["variant_phase_block_not_in_bam"] = {}
  for variant_key in vcf_variants_dictionary:
    for variant in vcf_variants_dictionary[variant_key]:
      variant_pbid = variant.return_VariantPhaseBlock()
      if not variant.return_IsPhasedHeterozygote():
        if variant_key in bam_phase_block_dictionary["phase_blocks"]["variant_not_phased_heterozygote"]:
          bam_phase_block_dictionary["phase_blocks"]["variant_not_phased_heterozygote"][variant_key].append(variant)
        else:
          bam_phase_block_dictionary["phase_blocks"]["variant_not_phased_heterozygote"][variant_key] = [variant]
      elif variant_pbid in bam_phase_block_dictionary["phase_blocks"]:
        bam_phase_block_dictionary["phase_blocks"][variant_pbid].add_Variant(variant)
      else:
        if variant_key in bam_phase_block_dictionary["phase_blocks"]["variant_phase_block_not_in_bam"]:
          bam_phase_block_dictionary["phase_blocks"]["variant_phase_block_not_in_bam"][variant_key].append(variant)
        else:
          bam_phase_block_dictionary["phase_blocks"]["variant_phase_block_not_in_bam"][variant_key] = [variant]
  for pbid in bam_phase_block_dictionary["phase_blocks"]:
    if pbid not in ["variant_not_phased_heterozygote", "variant_phase_block_not_in_bam"]:
      bam_phase_block_dictionary["phase_blocks"][pbid].add_FirstVariantPosition()
      bam_phase_block_dictionary["phase_blocks"][pbid].add_LastVariantPosition()

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

  # bam phase block dictionary
  bam_phase_block_dictionary = extract_phase_blocks_from_bam(
    args.bam,
    chr = chrom, 
    start_bp = start, 
    end_bp = end)
  
  # vcf variant dictionary
  vcf_variants_dictionary = extract_variants_from_VCF(
    args.vcf,
    args.vcf_id,
    chr = chrom,
    start_bp = start,
    end_bp = end)

  # add variants to bam phase block dictionary
  add_variants_to_phase_blocks(bam_phase_block_dictionary, vcf_variants_dictionary)

  os.makedirs(args.output_directory, exist_ok = True)
  output_file_path = os.path.join(args.output_directory, args.output_prefix + ".phaseblock.pkl")
  output_file = open(output_file_path, 'wb')
  pickle.dump(bam_phase_block_dictionary, output_file, pickle.HIGHEST_PROTOCOL)
  pickle.dump(vcf_variants_dictionary, output_file, pickle.HIGHEST_PROTOCOL)
  output_file.close()

if __name__ == '__main__':
  main(args)
