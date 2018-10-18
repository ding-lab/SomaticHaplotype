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
    self._start = start_bp
    self._end = end_bp
    self._key = chromosome + ":" + ps_id
    self._variants = {}

  def length(self, start_bp, end_bp):
    phase_set_length = end_bp - start_bp
    return(phase_set_length)

  def addVariants(self, variant):
    if variant.getVariantKey in self._variants:
      sys.exit("Whoa, variant already in PhaseSet.")
    else:
      self._variants[variant.getVariantKey] = variant

  def __str__(self):
    print_result = self._id + "\t" + self._chr + ":" + self._start + "-" + self._end
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
  parser.add_argument("out", help = "Path to output file")

  # Optional named arguments
  parser.add_argument('--bam', action = 'store', help = "Path to bam file")
  parser.add_argument('--vcf', action = 'store', help = "Path to VCF file")
  parser.add_argument('--vcf_id', action = 'store', help = "Sample ID from VCF file")
  parser.add_argument('--range', action = 'store', help = "Genomic range chr:start-stop")
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
        print(x)
      else:
        sys.exit("\n".join(error_message))
    
    else:
      print("Not X")
  
  else:
    sys.exit("Module must be one of " + ', '.join(acceptable_modules) + ".")


if __name__ == '__main__':
  main()