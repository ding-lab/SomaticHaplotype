import pysam
import vcf

################################################################################
# bam functions
################################################################################

def extract_phase_set_from_read(read):
  # Given a read from a bam, extract the phase set if it exists, else NA
  # Also report the position of the read and read quality metrics
  if read.is_duplicate or read.is_qcfail or read.is_secondary or not read.is_proper_pair:
    return(1)
  elif read.has_tag("PS"): # read is part of phase set
    tags_dict = {x:y for (x,y) in read.get_tags()}
    tags_dict["mapping_quality"] = read.mapping_quality
    tags_dict["reference_start"] = read.reference_start
    tags_dict["reference_end"] = read.reference_end
    return(tags_dict)
  else: # read is not part of a phase set
    return(2)

def extract_phase_sets_from_bam(bam_filename, 
  chr = None, start_bp = None, end_bp = None):
  samfile = pysam.AlignmentFile(bam_filename, "rb")
  phase_set_dict = {"n_total_reads" : 0, "n_reads_bad_quality" : 0,
    "n_reads_good_quality" : 0, "n_reads_phased" : 0, "n_reads_not_phased" : 0,
    "phase_sets" : {} }
  for read in samfile.fetch(chr, start_bp, end_bp):
    phase_set_dict["n_total_reads"] += 1
    read_info = extract_phase_set_from_read(read)
    if read_info == 1:
      phase_set_dict["n_reads_bad_quality"] += 1
    elif read_info == 2:
      phase_set_dict["n_reads_not_phased"] += 1
      phase_set_dict["n_reads_good_quality"] += 1
    else:
      phase_set_dict["n_reads_phased"] += 1
      phase_set_dict["n_reads_good_quality"] += 1
      ps_id = read.reference_name + ":" + str(read_info["PS"])
      if ps_id in phase_set_dict["phase_sets"]:
        if int(read_info["reference_start"]) < phase_set_dict["phase_sets"][ps_id]["min_position"]:
          phase_set_dict["phase_sets"][ps_id]["min_position"] = int(read_info["reference_start"])
        if int(read_info["reference_start"]) > phase_set_dict["phase_sets"][ps_id]["max_position"]:
          phase_set_dict["phase_sets"][ps_id]["max_position"] = int(read_info["reference_end"])
        phase_set_dict["phase_sets"][ps_id]["n_single_end_reads"] += 1
        if read_info["HP"] == 1:
          phase_set_dict["phase_sets"][ps_id]["n_support_H1"] += 1
        else:
          phase_set_dict["phase_sets"][ps_id]["n_support_H2"] += 1
        if "MI" not in read_info:
          read_info["MI"] = None
        if read_info["MI"] in phase_set_dict["phase_sets"][ps_id]["molecules"]:
          phase_set_dict["phase_sets"][ps_id]["molecules"][read_info["MI"]] += 1
        else:
          phase_set_dict["phase_sets"][ps_id]["molecules"][read_info["MI"]] = 1
      else:
        phase_set_dict["phase_sets"][ps_id] = {"chromosome" : read.reference_name, "min_position" : int(read_info["reference_start"]),
          "max_position" : int(read_info["reference_end"]),
          "n_single_end_reads" : 1, "n_support_H1" : 0, "n_support_H2" : 0,
          "molecules" : {read_info["MI"] : 1} }
        if read_info["HP"] == 1:
          phase_set_dict["phase_sets"][ps_id]["n_support_H1"] += 1
        else:
          phase_set_dict["phase_sets"][ps_id]["n_support_H2"] += 1
  samfile.close()
  return(phase_set_dict)

################################################################################
# VCF functions
################################################################################

def extract_key_from_record(record):
  # extract the key information from pyVCF record
  key_list = [ record.CHROM, record.start, record.end, record.REF, ",".join([ str(x) for x in record.ALT ]) ]
  key = ':'.join( [ str(x) for x in key_list ] )
  return(key)

def extract_phase_set_from_record(record, sample_id):
  # extract the phase set from a pyVCF record
  chrom = record.CHROM
  ps = record.genotype(sample_id).data.PS
  return(str(chrom) + ":" + str(ps))

def extract_phase_sets_from_vcf(vcf_filename, sample_id, 
  chr = None, start_bp = None, end_bp = None):
  # build a dictionary of phase sets present in VCF
  # only includes variants that are phased heterozygotes (useful for phase-related activities)
  this_vcf = vcf.Reader( filename = vcf_filename )

  phase_set_dict = {} # dictionary to hold phase set information

  for record in this_vcf.fetch( str(chr) , start_bp, end_bp ): # loop over each record in VCF
    if phased_heterozygote(record, sample_id):
      this_key = extract_key_from_record(record)
      this_key += ":" + str(record.genotype(sample_id).data.GT)
      this_ps = extract_phase_set_from_record(record, sample_id)
      if this_ps not in phase_set_dict:
        # initialize phase set dict as [ [list of keys], chrom, min pos, max pos, number of variants, phase_set ]
        phase_set_dict[this_ps] = [ [this_key], record.CHROM, record.POS, record.POS, 1, this_ps]
      else:
        phase_set_dict[this_ps][0].append(this_key) # append variant key to list of variant keys
        phase_set_dict[this_ps][4] += 1 # add one to number of variants associated with this phase set
        if record.POS < phase_set_dict[this_ps][2]: # update the min position 
          phase_set_dict[this_ps][2] = record.POS
        elif record.POS > phase_set_dict[this_ps][3]: # update the max position
          phase_set_dict[this_ps][3] = record.POS
        else:
          continue
    else:
      continue

  return(phase_set_dict)

def phased_heterozygote(record, sample_id):
  # determine if pyVCF record refers to a phased heterozygote (return True)
  # otherwise, variants are not useful for phasing
  phased = record.genotype(sample_id).phased
  heterozygote = record.genotype(sample_id).is_het
  if phased and heterozygote:
    return(True)
  else:
    return(False)


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

  # initialize phase set dictionary
  bam_phase_set_dictionary = extract_phase_sets_from_bam(
    args.bam,
    chr = chrom, 
    start_bp = start, 
    end_bp = end)

  vcf_phase_set_dictionary = extract_phase_sets_from_vcf(
    args.vcf,
    args.vcf_id,
    chr = chrom,
    start_bp = start,
    end_bp = end)

  return(bam_phase_set_dictionary)

if __name__ == '__main__':
  main(args)