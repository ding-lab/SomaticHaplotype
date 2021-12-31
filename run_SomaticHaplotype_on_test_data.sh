#!/bin/bash

range="chr1:109000000-121000000"
somatic_haplotype="python SomaticHaplotype.py"
run_10Xmapping="perl 10Xmapping/10Xmapping.pl"
convert_10Xmapping="python convert_10Xmapping_output.py"

bam_27522_1="test_data/27522_1.chr1_109Mb_to_121Mb.test.bam"
vcf_27522_1="test_data/27522_1.chr1_109Mb_to_121Mb.test.vcf.gz"

bam_27522_2="test_data/27522_2.chr1_109Mb_to_121Mb.test.bam"
vcf_27522_2="test_data/27522_2.chr1_109Mb_to_121Mb.test.vcf.gz"

maf_27522_1="test_data/27522_1.test.maf"
variant_27522_1="test_data/27522_1.test.variant_list"

out_dir=test_output
mkdir -p $out_dir

date
echo Output stored in $out_dir directory...

echo Running SomaticHaplotype phaseblock module on 27522_1...
$somatic_haplotype phaseblock $out_dir 27522_1.test --bam $bam_27522_1 --vcf $vcf_27522_1 --range $range --vcf_id 27522_1
echo Running SomaticHaplotype phaseblock module on 27522_2...
$somatic_haplotype phaseblock $out_dir 27522_2.test --bam $bam_27522_2 --vcf $vcf_27522_2 --range $range --vcf_id 27522_2

echo Running SomaticHaplotype summarize module on 27522_1...
$somatic_haplotype summarize $out_dir 27522_1.test --pb1 $out_dir/27522_1.test.phaseblock.pkl
echo Running SomaticHaplotype summarize module on 27522_2...
$somatic_haplotype summarize $out_dir 27522_2.test --pb1 $out_dir/27522_2.test.phaseblock.pkl

echo Running SomaticHaplotype somatic module on 27522_1 somatic mutation calls using input MAF with and without 10Xmapping submodule...
$run_10Xmapping --mapq 20 --bam $bam_27522_1 --maf $maf_27522_1 --out $out_dir/27522_1.test.10Xmapping.txt
$convert_10Xmapping $out_dir/27522_1.test.10Xmapping.txt $out_dir/27522_1.test.sombx.tsv
$somatic_haplotype somatic $out_dir 27522_1.test.maf_with_sombx --pb1 $out_dir/27522_1.test.phaseblock.pkl --range $range --maf $maf_27522_1 --sum $out_dir/27522_1.test.phase_blocks.tsv --sombx $out_dir/27522_1.test.sombx.tsv
$somatic_haplotype somatic $out_dir 27522_1.test.maf_without_sombx --pb1 $out_dir/27522_1.test.phaseblock.pkl --range $range --maf $maf_27522_1 --sum $out_dir/27522_1.test.phase_blocks.tsv

echo Running SomaticHaplotype somatic module on 27522_1 somatic mutation calls using input variant list with and without 10Xmapping submodule...
$somatic_haplotype somatic $out_dir 27522_1.test.variant_with_sombx --pb1 $out_dir/27522_1.test.phaseblock.pkl --range $range --variant $variant_27522_1 --sum $out_dir/27522_1.test.phase_blocks.tsv --sombx $out_dir/27522_1.test.sombx.tsv
$somatic_haplotype somatic $out_dir 27522_1.test.variant_without_sombx --pb1 $out_dir/27522_1.test.phaseblock.pkl --range $range --variant $variant_27522_1 --sum $out_dir/27522_1.test.phase_blocks.tsv

echo Running SomaticHaplotype extend module on 27522_2 extended by 27522_1...
$somatic_haplotype extend $out_dir 27522_2.test.extended_by_27522_1 --pb1 $out_dir/27522_2.test.phaseblock.pkl --pb2 $out_dir/27522_1.test.phaseblock.pkl --range $range --sum $out_dir/27522_2.test.phase_blocks.tsv

echo Test complete!
date
