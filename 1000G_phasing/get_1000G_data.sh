#!/bin/bash

# command line arguments
chromosome=$1 # which (one) chromosome you want to download
add_chr=$2 # if you bam has chr as part of chromosome name, must be 'true' to add chr to 1000G data

if [ -z $add_chr ]; then
  echo "Not enough command line arguments (need exactly 2). Make sure you have: chromosome add_chr(true/false)"
  exit 1
elif [ ! -z $3 ]; then
  echo "Too many command line arguments (need exactly 2). Make sure you have: chromosome add_chr(true/false)"
  exit 1
fi

data_dir=1000G_data/$chromosome
mkdir -p $data_dir

echo Retrieving 1000G data for chromosome $chromosome
echo
echo Data directory: $data_dir
echo 
echo Start time: 
date
echo
echo BCFtools version
bcftools --version

#######################
# Retrieve 1000G data #
#######################

ftp_manifest_file_remote_path="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/20170504_GRCh38_positions_manifest.txt"
ftp_manifest_directory_remote=$(echo $ftp_manifest_file_remote_path | rev | cut -f2- -d'/' | rev)
ftp_file_local=$data_dir/$(echo $ftp_manifest_file_remote_path | rev | cut -f1 -d'/' | rev)
vcf_1000G=$data_dir/1000G.vcf.gz
targets_1000G=$data_dir/1000G.targets.tsv.gz
map=$data_dir/$chromosome\.map

wget $ftp_manifest_file_remote_path -P $data_dir
wget http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip -P $data_dir

echo 1 chr1 > $data_dir/chr_name_conv.txt
for i in {2..22}; do
  echo $i chr$i >> $data_dir/chr_name_conv.txt
done

# unzip and update the chromosome centiMorgan map file
if $add_chr; then
  map_chr=$chromosome #chr already part of chromosome
else
  map_chr=chr$chromosome #chr not already part of chromosome
fi

cat $ftp_file_local | cut -f1 | cut -f2 -d'/' | grep -P \.$map_chr\_ | grep genotypes | grep -v tbi | while read file; do
  base_file_name=$(echo $file | cut -f1-4 -d'.')
  if $add_chr; then
    curl -s $ftp_manifest_directory_remote/$file |
      bcftools annotate --rename-chrs $data_dir/chr_name_conv.txt -Ou |
      bcftools norm --multiallelics - |
      bcftools norm --rm-dup all |
      bcftools view --genotype ^miss |
      bcftools plugin fill-tags |
      bcftools view -m2 -M2 -v snps -q 0.01 -Q 0.99 -i 'HWE > 0.00001' -Oz -o $vcf_1000G
  else
    curl -s $ftp_manifest_directory_remote/$file |
      #bcftools annotate --rename-chrs $data_dir/chr_name_conv.txt -Ou |
      bcftools norm --multiallelics - |
      bcftools norm --rm-dup all |
      bcftools view --genotype ^miss |
      bcftools plugin fill-tags |
      bcftools view -m2 -M2 -v snps -q 0.01 -Q 0.99 -i 'HWE > 0.00001' -Oz -o $vcf_1000G
  fi
  bcftools index -t $vcf_1000G
  bcftools query -f '%CHROM\t%POS\t%REF,%ALT\n' $vcf_1000G | bgzip -c > $targets_1000G && tabix -s1 -b2 -e2 $targets_1000G
done

unzip -o $data_dir/plink.GRCh38.map.zip -d $data_dir
if $add_chr; then
  sed 's/^/chr/g' $data_dir/plink.$map_chr\.GRCh38.map > $map
fi

rm -f $data_dir/plink* $data_dir/README.txt $data_dir/20170504_GRCh38_positions_manifest.txt $data_dir/chr_name_conv.txt
