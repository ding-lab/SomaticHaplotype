import argparse
import sys

class PhaseSet:
  def __init__(self, ps_id, chromosome, start_bp, end_bp):
    self._chr = chromosome
    self._end = int(end_bp)
    self._firstVariantPosition = None
    self._key = chromosome + ":" + ps_id
    self._lastVariantPosition = None
    self._molecules_H1 = {}
    self._molecules_H2 = {}
    self._n_single_end_reads = 0
    self._n_support_H1 = 0
    self._n_support_H2 = 0
    self._n_variants_H1 = 0
    self._n_variants_H2 = 0
    self._psid = ps_id
    self._start = int(start_bp)
    self._variants = {}

  def add_FirstVariantPosition(self):
    min_variant_position = float("inf")
    for variant_key in self.return_Variants():
      for variant in self.return_Variants()[variant_key]:
        if int(variant.return_StartPosition()) < min_variant_position:
          min_variant_position = int(variant.return_StartPosition())
    if min_variant_position == float("inf"):
      min_variant_position = "NA"
    self._firstVariantPosition = min_variant_position

  def add_LastVariantPosition(self):
    max_variant_position = float("-inf")
    for variant_key in self.return_Variants():
      for variant in self.return_Variants()[variant_key]:
        if int(variant.return_StartPosition()) > max_variant_position:
          max_variant_position = int(variant.return_StartPosition())
    if max_variant_position == float("-inf"):
      max_variant_position = "NA"
    self._lastVariantPosition = max_variant_position

  def add_MoleculeH1(self, molecule):
    if molecule in self._molecules_H1:
      self._molecules_H1[molecule] += 1
    else:
      self._molecules_H1[molecule] = 1

  def add_MoleculeH2(self, molecule):
    if molecule in self._molecules_H2:
      self._molecules_H2[molecule] += 1
    else:
      self._molecules_H2[molecule] = 1

  def add_SingleEndRead(self):
    self._n_single_end_reads += 1

  def add_SupportH1(self):
    self._n_support_H1 += 1

  def add_SupportH2(self):
    self._n_support_H2 += 1

  def add_Variant(self, variant):
    if variant.return_VariantKey() in self._variants:
      self._variants[variant.return_VariantKey()].append(variant)
    else:
      self._variants[variant.return_VariantKey()] = [variant]

    if variant.return_IsPhasedHeterozygote():
      if variant.return_Genotype()[0] != "0":
        self._n_variants_H1 += 1
      elif variant.return_Genotype()[2] != "0":
        self._n_variants_H2 += 1  

  def add_VariantH1(self):
    self._n_variants_H1 += 1

  def add_VariantH2(self):
    self._n_variants_H2 += 1

  def update_End(self, new_end_bp):
    self._end = int(new_end_bp)

  def update_Start(self, new_start_bp):
    self._start = int(new_start_bp)

  # Functions to return values about the PhaseSet

  def return_Chromosome(self):
    return(self._chr)

  def return_End(self):
    return(self._end)

  def return_FirstVariantPosition(self):
    return(self._firstVariantPosition)

  def return_Key(self):
    return(self._key)

  def return_LastVariantPosition(self):
    return(self._lastVariantPosition)

  def return_LengthReads(self):
    return(self._end - self._start)

  def return_LengthVariants(self):
    # return(self._lastVariantPosition - self._firstVariantPosition)
    try:
      return(self._lastVariantPosition - self._firstVariantPosition)
    except:
      return("NA")

  def return_MoleculesH1(self):
    return(self._molecules_H1)

  def return_MoleculesH2(self):
    return(self._molecules_H2)

  def return_nSingleEndReads(self):
    return(self._n_single_end_reads)

  def return_nSupportH1(self):
    return(self._n_support_H1)

  def return_nSupportH2(self):
    return(self._n_support_H2)

  def return_nVariantsH1(self):
    return(self._n_variants_H1)

  def return_nVariantsH2(self):
    return(self._n_variants_H2)

  def return_PhaseSetID(self):
    return(self._psid)

  def return_Start(self):
    return(self._start)

  def return_Variants(self):
    return(self._variants)

  def __str__(self):
    print_result = self._psid + "\t" + self._chr + ":" + str(self._start) + "-" + str(self._end)
    return(print_result)

class Variant:
  def __init__(self, record, sample_id):
    self._alt = record.ALT 
    self._chr = record.CHROM
    self._end = record.end
    self._filter = record.FILTER
    self._genotype = record.genotype(sample_id).data.GT
    self._is_heterozygote = record.genotype(sample_id).is_het
    self._is_phased = record.genotype(sample_id).phased
    self._is_phased_heterozygote = self.determine_if_phased_heterozygote()
    self._is_snp = record.is_snp
    self._molecules = self.extract_molecules_from_VCF_record(record, sample_id)
    self._position = record.POS
    self._psid = self.extract_variant_phase_set_from_VCF_record(record, sample_id)
    self._ref = record.REF 
    self._sampleid = sample_id
    self._start = record.start

  def determine_if_phased_heterozygote(self):
  # determine if Varaint refers to a phased heterozygote (return True)
  # otherwise, variants are not useful for phasing
    if self._is_phased and self._is_heterozygote:
      return(True)
    else:
      return(False)

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
    return(molecule_dictionary)

  def extract_variant_phase_set_from_VCF_record(self, record, sample_id):
    # extract the phase set from a pyVCF record
    chrom = record.CHROM
    ps = record.genotype(sample_id).data.PS
    psid = str(chrom) + ":" + str(ps)
    return(psid)

  # Functions to return values about the Variant

  def return_AlternateAllele(self):
    return(self._alt)

  def return_Chromosome(self):
    return(self._chr)

  def return_EndPosition(self):
    return(self._end)

  def return_Filter(self):
    return(self._filter)

  def return_Genotype(self):
    return(self._genotype)

  def return_IsHeterozygote(self):
    return(self._is_heterozygote)

  def return_IsPhased(self):
    return(self._is_phased)

  def return_IsPhasedHeterozygote(self):
    return(self._is_phased_heterozygote)
 
  def return_IsSNP(self):
    return(self._is_snp)

  def return_Molecules(self):
    return(self._molecules)

  def return_PhaseSetID(self):
    return(self._psid)

  def return_Position(self):
    return(self._position)

  def return_ReferenceAllele(self):
    return(self._ref)

  def return_SampleID(self):
    return(self._sampleid)

  def return_StartPosition(self):
    return(self._start)

  def return_VariantKey(self):
    key_list = [ self._chr, self._position, self._ref, ",".join([ str(x) for x in self._alt ]) ]
    variant_key = ':'.join( [ str(x) for x in key_list ] )
    return(variant_key)

  def return_VariantPhaseSet(self):
    return(self._psid)

  def __str__(self):
    print_result = self.return_VariantKey()
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
  parser.add_argument('--sum', action = 'store', help = "Path to existing summary file")
  parser.add_argument('--maf', action = 'store', help = "Path to sample-specific somatic MAF (assumes all variants are associated with single sample)")
  parser.add_argument('--sombx', action = 'store', help = "Path to file containing barcodes supporting somatic MAF variants extracted from BAM")
  parser.add_argument('--variant', action = 'store', help = "Path to file containing newline-separated variant IDs, format CHROM:POS:REF:ALT (ALT is comma separated list of each ALT variant)")
  parser.add_argument('--plot', action = 'store', help = "If True, outputs necessary files for plotting.")
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
        import phaseset
        x = phaseset.main(args)
      else:
        sys.exit("\n".join(error_message))
    elif args.module == "extend":
      if args.ps1 is None:
        no_error = False
        error_message.append("The extend module requires a --ps1.")
      if args.ps2 is None:
        no_error = False
        error_message.append("The extend module requires a --ps2.")
      if args.range is None:
        no_error = False
        error_message.append("The extend module requires a --range.")
      if no_error:
        import extend
        x = extend.main(args)
      else:
        sys.exit("\n".join(error_message))
    elif args.module == "summarize":
      if args.ps1 is None:
        no_error = False
        error_message.append("The summarize module requires a --ps1 (phase set file).")
      if no_error:
        import summarize
        x = summarize.main(args)
      else:
        sys.exit("\n".join(error_message))
    elif args.module == "somatic":
      if args.ps1 is None:
        no_error = False
        error_message.append("The somatic module requires a --ps1 (phase set file).")
      if args.range is None:
        no_error = False
        error_message.append("The somatic module requires a --range (genomic range).")
      if args.maf is None and args.variant is None:
        no_error = False
        error_message.append("The somatic module requires a --maf (MAF) or --variant (variant IDs).")
      if args.maf is not None and args.variant is not None:
        no_error = False
        error_message.append("The somatic module can only have a --maf (MAF) or --variant (variant IDs), not both.")
      if args.sum is None:
        no_error = False
        error_message.append("The somatic module requires a --sum (phase set summary file).")
      if no_error:
        import somatic
        x = somatic.main(args)
      else:
        sys.exit("\n".join(error_message))
    else:
      print("Not X")  
  
  else:
    sys.exit("Module must be one of " + ', '.join(acceptable_modules) + ".")


if __name__ == '__main__':
  main()
