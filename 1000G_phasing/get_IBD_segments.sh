#!/bin/bash

# command line arguments
sample=$1 # sample name given in your bam (SM tag in bam header)
chromosome=$2 # which (one) chromosome you want to analyze (should match what is in bam/index)
bam=$3 # your bam file (remember to check if your bam has chr as part of chromosome name)
fa=$4 # fasta used to generate your bam

if [ -z $fa ]; then
  echo "Not enough command line arguments (need exactly 4). Make sure you have: sample chromosome bam fasta"
  exit 1
elif [ ! -z $5 ]; then
  echo "Too many command line arguments (need exactly 4). Make sure you have: sample chromosome bam fasta"
  exit 1
fi

output_dir="IBD/$sample/$chromosome"
if [ -d $output_dir ]; then
  echo Output directory already exists for this sample and chromsome: IBD/$sample/$chromosome\.
  exit 1
else
  data_dir=1000G_data/$chromosome
  results_dir=$output_dir/results
  tools_dir=tools
  temp_dir=$output_dir/temp
  mkdir -p $results_dir
  mkdir -p $temp_dir
fi

if [ ! -d $data_dir ]; then
  echo Data directory $data_dir does not exist
  exit 1
fi

if [ ! -d $tools_dir ]; then
  echo Tools directory $tools_dir does not exist
  exit 1
fi

echo Generating IBD segments for $sample on chromosome $chromosome
echo
echo Output directory: IBD/$sample/$chromosome
echo 
echo Start time: 
date
echo
echo BCFtools version
bcftools --version
echo
echo Java version
java -version
echo

#######################################
# Download beagle, Refined IBD, merge #
#######################################

# Get beagle.21Sep19.ec3.jar
beagle="java -jar $tools_dir/beagle.21Sep19.ec3.jar"

# Get Refined IBD
refined_ibd="java -jar $tools_dir/refined-ibd.16May19.ad5.jar"

# Get merge IBD
merge_ibd="java -jar $tools_dir/merge-ibd-segments.16May19.ad5.jar"

vcf_1000G=$data_dir/1000G.vcf.gz
targets_1000G=$data_dir/1000G.targets.tsv.gz
map=$data_dir/$chromosome\.map

######################################################
# Call variants on new sample, resolve VCF positions #
######################################################
vcf_sample=$temp_dir/$sample\.vcf.gz
vcf_sample_resolved=$temp_dir/$sample\.resolved.vcf.gz
vcf_beagle=$results_dir/beagle.vcf.gz
vcf_merged=$temp_dir/merged.vcf.gz

# run mpileup and call variants, index, and merge with 1000G VCF
bcftools mpileup -Ou -f $fa -r $chromosome -T $targets_1000G --annotate DP,AD $bam | 
  bcftools call --keep-alts -T $targets_1000G --consensus-caller | 
  bcftools view -i 'FORMAT/DP[0] >=10 && ((GT="hom" && AD[0:0]/FORMAT/DP[0] < 0.1) | (GT="het" && AD[0:0]/FORMAT/DP[0] > 0.4 && AD[0:0]/FORMAT/DP[0] < 0.6) | (GT="hom" && AD[0:0]/FORMAT/DP[0] > 0.9))' -Oz -o $vcf_sample
bcftools index -t $vcf_sample

bcftools merge -m all -Ou $vcf_1000G $vcf_sample | 
  bcftools norm --multiallelics - |
  bcftools norm --rm-dup all |
  bcftools view -v snps -s $sample -T $targets_1000G -Oz -o $vcf_sample_resolved
bcftools index -t $vcf_sample_resolved

# run beagle -- output is the sample of interest
$beagle nthreads=1 map=$map ref=$vcf_1000G gt=$vcf_sample_resolved out=$results_dir/beagle
bcftools index -t $vcf_beagle

# merge existing vcf files
bcftools merge -m none -Ou $vcf_1000G $vcf_beagle | bcftools norm --rm-dup all -Ou | bcftools view --genotype ^miss -Oz -o $vcf_merged 
bcftools index -t $vcf_merged

# run Refined IBD
$refined_ibd nthreads=1 length=0.5 gt=$vcf_merged map=$map out=$results_dir/refined_ibd

# merge IBD segements
gunzip -c $results_dir/refined_ibd.ibd.gz | $merge_ibd $vcf_merged $map 0.5 1 | bgzip -c > $results_dir/refined_ibd.ibd.merged.gz
                                                       
#############################
# Remove intermediate files #
#############################

rm -fr $temp_dir

###########################################
# Check if output exists and is non-empty #
###########################################
echo Finish time:
date
echo
if [ -z $(gzip -cd  $results_dir/refined_ibd.ibd.merged.gz | head -c1) ]; then
  echo $sample $chromosome: Output file $results_dir/refined_ibd.ibd.merged.gz does not exist or is empty.
  echo $sample $bam $fa $chromosome | tr ' ' '\t' >> rerun_input_file_list.tsv
else
  echo $sample $chromosome: passes basic QC
fi
