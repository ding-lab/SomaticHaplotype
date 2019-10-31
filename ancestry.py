import os
import pickle
import vcf

from SomaticHaplotype import *

################################################################################
# functions
################################################################################

def get_variant_positions(reference_vcf, args):
  os.makedirs(args.output_directory, exist_ok = True)
  output_file_path = os.path.join(args.output_directory, args.output_prefix + ".1000G_variant_positions.tsv")
  output_file = open(output_file_path, "w")
  # use vcf module to get print out variant positions
  #output_file.write(str(chr) + "\t" + str(pos) + "\n")
  output_file.close()

def run_beagle(reference_vcf, args):
  target_vcf = os.path.join(args.output_directory, args.output_prefix + ".1000G_variant_calls.vcf.gz") 
  output_vcf_prefix = os.path.join(args.output_directory, args.output_prefix + ".1000G_variant_calls.beagle_phased")
  subprocess.run(['beagle ref=', reference_vcf, 'gt=', target_vcf, 'out=', output_vcf_prefix], check = True, universal_newlines = True, shell = False)
  subprocess.run(['bcftools index -t', output_vcf_prefix+".vcf.gz"], check = True, universal_newlines = True, shell = False)

def run_call(bam, reference_fasta, args):
  os.makedirs(args.output_directory, exist_ok = True)
  variant_positions = os.path.join(args.output_directory, args.output_prefix + ".1000G_variant_positions.tsv")
  output_vcf = os.path.join(args.output_directory, args.output_prefix + ".1000G_variant_calls.vcf.gz") 
  subprocess.run(['bcftools mpileup -Ou -f', reference_fasta, bam, '-R', variant_positions, '| bcftools call -m -Oz -o', output_vcf], check = True, universal_newlines = True, shell = False)
  subprocess.run(['bcftools index -t', output_vcf], check = True, universal_newlines = True, shell = False)

def run_merge_vcf(reference_vcf, args):
  os.makedirs(args.output_directory, exist_ok = True)
  phased_vcf = os.path.join(args.output_directory, args.output_prefix + ".1000G_variant_calls.beagle_phased.vcf.gz")
  output_vcf = os.path.join(args.output_directory, args.output_prefix + ".1000G_variant_calls.beagle_phased.merged_samples.vcf.gz")
  subprocess.run(['bcftools merge -m both -Oz -o', output_vcf, reference_vcf, phased_vcf], check = True, universal_newlines = True, shell = False)

def run_refined_ibd(map_file, args):
  external_resources_directory = os.path.join(os.path.dirname(__file__), "external_resources")
  os.makedirs(args.output_directory, exist_ok = True)
  merged_vcf = os.path.join(args.output_directory, args.output_prefix + ".1000G_variant_calls.beagle_phased.merged_samples.vcf.gz")
  output_ibd = output_vcf = os.path.join(args.output_directory, args.output_prefix + ".1000G_variant_calls.beagle_phased.merged_samples.ibd_segments")
  subprocess.run(['java -jar', external_resources_directory, 'refined-ibd.jar gt=', merged_vcf, 'map=', map_file, 'out=', output_ibd], check = True, universal_newlines = True, shell = False)

def run_merge_ibd(map_file, args):
  external_resources_directory = os.path.join(os.path.dirname(__file__), "external_resources")
  os.makedirs(args.output_directory, exist_ok = True)
  merged_vcf = os.path.join(args.output_directory, args.output_prefix + ".1000G_variant_calls.beagle_phased.merged_samples.vcf.gz")
  output_ibd = output_vcf = os.path.join(args.output_directory, args.output_prefix + ".1000G_variant_calls.beagle_phased.merged_samples.merged_ibd_segments")
  subprocess.run(['java -jar', external_resources_directory, 'merge-ibd-segments.jar vcf=', merged_vcf, 'map=', map_file, 'gap=0.6 discord=1 > ', output_ibd], check = True, universal_newlines = True, shell = False)



################################################################################
# main
################################################################################   

def main(args):

  # open up summary file
  summary_file = open(args.sum, 'r')

  # path to input pickle file of first sample
  pickle_path_ps1 = args.ps1
  with open(pickle_path_ps1, 'rb') as pickle_file_ps1:
    bam_phase_set_dictionary_ps1 = pickle.load(pickle_file_ps1)
    vcf_variants_dictionary_ps1 = pickle.load(pickle_file_ps1)

  # path to input pickle file of second sample
  pickle_path_ps2 = args.ps2
  with open(pickle_path_ps2, 'rb') as pickle_file_ps2:
    bam_phase_set_dictionary_ps2 = pickle.load(pickle_file_ps2)
    vcf_variants_dictionary_ps2 = pickle.load(pickle_file_ps2)

  # parse the genomic range argument (which is required)
  chrom = str(args.range.split(":")[0])
  try:
    start = int(args.range.split(":")[1].split("-")[0])
  except:
    start = None
  try:
    end = int(args.range.split(":")[1].split("-")[1])
  except:
    end = None

  # run functions

  get_variant_positions(reference_vcf, args)
  run_call(bam, reference_fasta, args)
  run_beagle(reference_vcf, args)
  run_merge_vcf(reference_vcf, args)
  run_refined_ibd(map_file, args)
  run_merge_ibd(map_file, args)

  # write outputs
  #s.makedirs(args.output_directory, exist_ok = True)
  #output_file_path = os.path.join(args.output_directory, args.output_prefix + "")
  #output_file = open(output_file_path, 'w')
  #output_file.close()

if __name__ == '__main__':
  main(args)
