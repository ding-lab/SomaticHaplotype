import argparse
import sys

# Import modules from this repo
import phaseset
#import summarize
#import visualize
#import extend
#import somatic

class PhaseSet:
  def __init__(self, ps_id, chromosome, start_bp, end_bp):
    self._psid = ps_id
    self._chr = chromosome
    self._start = int(start_bp)
    self._end = int(end_bp)
    self._key = chromosome + ":" + ps_id
    self._n_single_end_reads = 0
    self._n_support_H1 = 0
    self._n_support_H2 = 0
    self._moleculesH1 = {}
    self._moleculesH2 = {}
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
    if molecule in self._moleculesH1:
      self._moleculesH1[molecule] += 1
    else:
      self._moleculesH1[molecule] = 1

  def addMoleculeH2(self, molecule):
    if molecule in self._moleculesH2:
      self._moleculesH2[molecule] += 1
    else:
      self._moleculesH2[molecule] = 1

  def length(self, start_bp, end_bp):
    phase_set_length = end_bp - start_bp
    return(phase_set_length)

  def addVariants(self, variant):
    if variant.getVariantKey in self._variants:
      sys.exit("Whoa, variant already in PhaseSet.")
    else:
      self._variants[variant.getVariantKey] = variant

  def __str__(self):
    print_result = self._psid + "\t" + self._chr + ":" + str(self._start) + "-" + str(self._end)
    return(print_result)

class Variant:
  def __init__(self, chromosome, position, REF, ALT, ps_id):
    self._chr = chromosome
    self._pos = position
    self._ref = REF
    self._alt = ALT
    self._psid = ps_id

  def getVariantKey(self):
    variant_key = ":".join([self._chr, self._pos, self._ref, self._alt, self._psid])
    return(variant_key)

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
        for y in x["phase_sets"]:
          print(x["phase_sets"][y])
      else:
        sys.exit("\n".join(error_message))
    
    else:
      print("Not X")
  
  else:
    sys.exit("Module must be one of " + ', '.join(acceptable_modules) + ".")


if __name__ == '__main__':
  main()