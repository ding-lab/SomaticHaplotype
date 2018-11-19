import argparse
import sys

# Import modules from this repo
import phaseset
import summarize
#import visualize
#import extend
#import somatic

class PhaseSet:
  def __init__(self, ps_id, chromosome, start_bp, end_bp):
    self._psid = ps_id
    self._chr = chromosome
    self._start = int(start_bp)
    self._end = int(end_bp)
    self._firstVariantPosition = None
    self._lastVariantPosition = None
    self._key = chromosome + ":" + ps_id
    self._n_single_end_reads = 0
    self._n_support_H1 = 0
    self._n_support_H2 = 0
    self._n_variants_H1 = 0
    self._n_variants_H2 = 0
    self._molecules_H1 = {}
    self._molecules_H2 = {}
    self._variants = {}

  def psStart(self):
    return(int(self._start))

  def updateStart(self, new_start_bp):
    self._start = int(new_start_bp)

  def psEnd(self):
    return(int(self._end))

  def updateEnd(self, new_end_bp):
    self._end = int(new_end_bp)

  def addSingleEndRead(self):
    self._n_single_end_reads += 1

  def addSupportH1(self):
    self._n_support_H1 += 1

  def addSupportH2(self):
    self._n_support_H2 += 1

  def addMoleculeH1(self, molecule):
    if molecule in self._molecules_H1:
      self._molecules_H1[molecule] += 1
    else:
      self._molecules_H1[molecule] = 1

  def addMoleculeH2(self, molecule):
    if molecule in self._molecules_H2:
      self._molecules_H2[molecule] += 1
    else:
      self._molecules_H2[molecule] = 1

  def length(self, start_bp, end_bp):
    phase_set_length = end_bp - start_bp
    return(phase_set_length)

  def addVariantH1(self):
    self._n_variants_H1 += 1

  def addVariantH2(self):
    self._n_variants_H2 += 1

  def addVariant(self, variant):
    if variant.getVariantKey() in self._variants:
      self._variants[variant.getVariantKey()].append(variant)
    else:
      self._variants[variant.getVariantKey()] = [variant]

    if variant.getPhasedHeterozygoteStatus():
      if variant.getGenotype()[0] != 0:
        self._n_variants_H1 += 1
      elif variant.getGenotype()[2] != 0:
        self._n_variants_H2 += 1

  def getVariants(self):
    return(self._variants)

  def addFirstVariantPosition(self):
    min_variant_position = float("inf")
    for variant_key in self.getVariants():
      for variant in self.getVariants()[variant_key]:
        if int(variant.getStartPosition()) < min_variant_position:
          min_variant_position = int(variant.getStartPosition())
    if min_variant_position == float("inf"):
      min_variant_position = "NA"
    self._firstVariantPosition = min_variant_position

  def addLastVariantPosition(self):
    max_variant_position = float("-inf")
    for variant_key in self.getVariants():
      for variant in self.getVariants()[variant_key]:
        if int(variant.getStartPosition()) > max_variant_position:
          max_variant_position = int(variant.getStartPosition())
    if max_variant_position == float("-inf"):
      max_variant_position = "NA"
    self._lastVariantPosition = max_variant_position

  def getFirstVariantPosition(self):
    return(self._firstVariantPosition)

  def getLastVariantPosition(self):
    return(self._lastVariantPosition)

  def __str__(self):
    print_result = self._psid + "\t" + self._chr + ":" + str(self._start) + "-" + str(self._end)
    return(print_result)

class Variant:
  def __init__(self, record, sample_id):
    self._sampleid = sample_id
    self._chr = record.CHROM
    self._start = record.start
    self._end = record.end
    self._ref = record.REF 
    self._alt = record.ALT 
    self._filter = record.FILTER
    self._genotype = record.genotype(sample_id).data.GT
    self._is_phased = record.genotype(sample_id).phased
    self._is_heterozygote = record.genotype(sample_id).is_het
    self._is_snp = record.is_snp
    self.extract_variant_phase_set_from_VCF_record(record, sample_id)
    self.extract_molecules_from_VCF_record(record, sample_id)
    self.determine_if_phased_heterozygote()

  def extract_molecules_from_VCF_record(self, record, sample_id):
    # add molecules supporting each allele
    molecule_dictionary = {}
    molecules_split_by_allele = record.genotype(sample_id).data.BX
    n_alleles = len(molecules_split_by_allele)
    for allele in range(n_alleles):
      if allele not in molecule_dictionary:
        molecule_dictionary[allele] = {}
      molecules_supporting_this_allele = molecules_split_by_allele[allele].split(";")
      for mol_qual in molecules_supporting_this_allele:
        molecule = mol_qual.split("_")[0]
        quality_list = mol_qual.split("_")[1:]
        if molecule in molecule_dictionary[allele]:
            sys.exit("Whoa, molecule already in molecule dictionary")
        else:
          molecule_dictionary[allele][molecule] = quality_list

    self._molecules = molecule_dictionary

  def extract_variant_phase_set_from_VCF_record(self, record, sample_id):
    # extract the phase set from a pyVCF record
    chrom = record.CHROM
    ps = record.genotype(sample_id).data.PS
    self._psid = str(chrom) + ":" + str(ps)

  def determine_if_phased_heterozygote(self):
  # determine if Varaint refers to a phased heterozygote (return True)
  # otherwise, variants are not useful for phasing
    if self._is_phased and self._is_heterozygote:
      self._is_phased_heterozygote = True
    else:
      self._is_phased_heterozygote = False

  def getPhasedHeterozygoteStatus(self):
    return(self._is_phased_heterozygote)

  def getVariantPhaseSet(self):
    return(self._psid)

  def getVariantKey(self):
    key_list = [ self._chr, self._start, self._end, self._ref, ",".join([ str(x) for x in self._alt ]) ]
    variant_key = ':'.join( [ str(x) for x in key_list ] )
    return(variant_key)

  def getStartPosition(self):
    return(self._start)

  def getGenotype(self):
    return(self._genotype)

  def __str__(self):
    print_result = self.getVariantKey()
    return(print_result)

def parse_input_arguments():
  
  parser = argparse.ArgumentParser()

  # Required positional arguments
  parser.add_argument("module", help = "Module the program should run. Could be phaseset, summarize, visualize, extend, somatic.")
  parser.add_argument("output_directory", help = "Absolute or relative path to output directory")
  parser.add_argument("output_prefix", help = "Prefix for file names in output directory")

  # Optional named arguments
  parser.add_argument('--bam', action = 'store', help = "Path to bam file")
  parser.add_argument('--vcf', action = 'store', help = "Path to VCF file")
  parser.add_argument('--vcf_id', action = 'store', help = "Sample ID from VCF file")
  parser.add_argument('--range', action = 'store', help = "Genomic range chr:start-stop, chr, chr:start, chr:-stop")
  parser.add_argument('--ps1', action = 'store', help = "Path to first phase set file")
  parser.add_argument('--ps2', action = 'store', help = "Path to second phase set file")
  parser.add_argument('--version', action = 'version', version = '%(prog)s 0.1')

  # Return arguments object
  return(parser.parse_args())


def main():
  args = parse_input_arguments()
  
  acceptable_modules = ["phaseset", "summarize", "visualize", "extend", "somatic"]
  if args.module in acceptable_modules:
    no_error = True
    error_message = []

    if args.module == "phaseset":
      if args.bam is None:
        no_error = False
        error_message.append("The phaseset module requires a --bam.")
      if args.vcf is None:
        no_error = False
        error_message.append("The phaseset module requires a --vcf.")
      if args.vcf_id is None:
        no_error = False
        error_message.append("The phaseset module requires a --vcf_id.")
      if args.range is None:
        no_error = False
        error_message.append("The phaseset module requires a --range.")
      if no_error:
        x = phaseset.main(args)
      else:
        sys.exit("\n".join(error_message))
    elif args.module == "summarize":
      if args.ps1 is None:
        no_error = False
        error_message.append("The summarize module requires a --ps1 (phase set file).")
      if no_error:
        x = summarize.main(args)
      else:
        sys.exit("\n".join(error_message))
    else:
      print("Not X")
  
  else:
    sys.exit("Module must be one of " + ', '.join(acceptable_modules) + ".")


if __name__ == '__main__':
  main()
