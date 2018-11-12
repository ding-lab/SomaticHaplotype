import os
import pickle
import pysam
import vcf

from SomaticHaplotype import *

################################################################################
# bam functions
################################################################################

def extract_read_info(read):
  # Given a read from a bam, extract the phase set if it exists, else NA
  # Also report the position of the read and read quality metrics
  if read.is_duplicate or read.is_qcfail or read.is_secondary or not read.is_proper_pair:
    return("read is bad quality")
  elif read.has_tag("PS"): # read is part of phase set
    tags_dict = {x:y for (x,y) in read.get_tags()}
    if "MI" not in tags_dict:
      tags_dict["MI"] = None
    return(tags_dict)
  else: # read is not part of a phase set
    return("read is not part of phase set")

def extract_phase_sets_from_bam(bam_filename, chr = None, start_bp = None, end_bp = None):
  
  samfile = pysam.AlignmentFile(bam_filename, "rb")
  
  phase_set_dict = {"n_total_reads" : 0, "n_reads_bad_quality" : 0,
    "n_reads_good_quality" : 0, "n_reads_phased" : 0, "n_reads_not_phased" : 0,
    "phase_sets" : {} }
  
  for read in samfile.fetch(chr, start_bp, end_bp):
    phase_set_dict["n_total_reads"] += 1
    read_info = extract_read_info(read)
    if read_info == "read is bad quality":
      phase_set_dict["n_reads_bad_quality"] += 1
    elif read_info == "read is not part of phase set":
      phase_set_dict["n_reads_not_phased"] += 1
      phase_set_dict["n_reads_good_quality"] += 1
    else:
      
      phase_set_dict["n_reads_phased"] += 1
      phase_set_dict["n_reads_good_quality"] += 1
      
      ps_id = read.reference_name + ":" + str(read_info["PS"])

      if ps_id in phase_set_dict["phase_sets"]:
        
        if int(read.reference_start) < phase_set_dict["phase_sets"][ps_id].psStart():
          phase_set_dict["phase_sets"][ps_id].updateStart(read.reference_start)
        if int(read.reference_end) > phase_set_dict["phase_sets"][ps_id].psEnd():
          phase_set_dict["phase_sets"][ps_id].updateEnd(read.reference_end)
        
        phase_set_dict["phase_sets"][ps_id].addSingleEndRead()

        if read_info["HP"] == 1:
          phase_set_dict["phase_sets"][ps_id].addSupportH1()
          phase_set_dict["phase_sets"][ps_id].addMoleculeH1(read_info["MI"])
        elif read_info["HP"] == 2:
          phase_set_dict["phase_sets"][ps_id].addSupportH2()
          phase_set_dict["phase_sets"][ps_id].addMoleculeH2(read_info["MI"])
        else:
          sys.exit("Whoa no H1 or H2 support for\n" + str(read))

      else:
      
        phase_set_dict["phase_sets"][ps_id] = PhaseSet(ps_id = ps_id, 
          chromosome = read.reference_name, 
          start_bp = int(read.reference_start), 
          end_bp = int(read.reference_end))

        phase_set_dict["phase_sets"][ps_id].addSingleEndRead()

        if read_info["HP"] == 1:
          phase_set_dict["phase_sets"][ps_id].addSupportH1()
          phase_set_dict["phase_sets"][ps_id].addMoleculeH1(read_info["MI"])
        elif read_info["HP"] == 2:
          phase_set_dict["phase_sets"][ps_id].addSupportH2()
          phase_set_dict["phase_sets"][ps_id].addMoleculeH2(read_info["MI"])
        else:
          sys.exit("Whoa no H1 or H2 support for\n" + str(read))
        
  samfile.close()
  return(phase_set_dict)

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
    if this_variant.getVariantKey() in variant_dict: # check if variant already in dictionary
      variant_dict[this_variant.getVariantKey()].append(this_variant)
      #sys.exit("Whoa variant already exists in variant dictionary\n" + str(record))
    else:
      variant_dict[this_variant.getVariantKey()] = [this_variant]

  return(variant_dict)


################################################################################
# Combined bam phase set and VCF variants functions
################################################################################

def add_variants_to_phase_sets(bam_phase_set_dictionary, vcf_variants_dictionary):
  bam_phase_set_dictionary["phase_sets"]["variant_not_phased_heterozygote"] = {}
  for variant_key in vcf_variants_dictionary:
    for variant in vcf_variants_dictionary[variant_key]:
      variant_psid = variant.getVariantPhaseSet()
      if variant_psid in bam_phase_set_dictionary["phase_sets"]:
        bam_phase_set_dictionary["phase_sets"][variant_psid].addVariant(variant)
      elif not variant.getPhasedHeterozygoteStatus():
        if variant_key in bam_phase_set_dictionary["phase_sets"]["variant_not_phased_heterozygote"]:
          bam_phase_set_dictionary["phase_sets"]["variant_not_phased_heterozygote"][variant_key].append(variant)
        else:
          bam_phase_set_dictionary["phase_sets"]["variant_not_phased_heterozygote"][variant_key] = [variant]
      else:
        sys.exit("Whoa, variant phase set not in bam phase set dictionary or not not phased heterozygote.")

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

  # bam phase set dictionary
  bam_phase_set_dictionary = extract_phase_sets_from_bam(
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

  # add variants to bam phase set dictionary
  add_variants_to_phase_sets(bam_phase_set_dictionary, vcf_variants_dictionary)

  os.makedirs(args.output_directory, exist_ok = True)
  output_file_path = os.path.join(args.output_directory, args.output_prefix + ".phaseset.pkl")
  output_file = open(output_file_path, 'wb')
  pickle.dump(bam_phase_set_dictionary, output_file, pickle.HIGHEST_PROTOCOL)
  pickle.dump(vcf_variants_dictionary, output_file, pickle.HIGHEST_PROTOCOL)
  output_file.close()

if __name__ == '__main__':
  main(args)
