# SomaticHaplotype

Somatic mutations occur on specific haplotypes. Short-read sequencing methods obscure the original haplotype relationship and important information about the cis/trans relationship of somatic events is lost. This work utilizes linked-read WGS data to leverage germline haplotype structures and infer the haplotypic context of somatic mutations.

## Installation

This repository contains submodules. Use this command when you clone: 

`git clone --recurse-submodules https://github.com/ding-lab/SomaticHaplotype`

If you have already cloned this repository without using `--recurse-submodules`, you can retroactively initialize and update the submodules recursively using

`git submodule update --init --recursive`

Source: https://git-scm.com/book/en/v2/Git-Tools-Submodules

## Conda environment

To create and activate the SomaticHaplotype conda environment with correct python and python modules, run

```
cd SomaticHaplotype
bash external_downloads/set_up_SomaticHaplotype_environment.sh
```


## How to run

```python SomaticHaplotype.py [module] [output directory] [output prefix]```

```
python SomaticHaplotype.py --help

usage: SomaticHaplotype.py [-h] [--bam BAM] [--vcf VCF] [--vcf_id VCF_ID]
                           [--range RANGE] [--ps1 PS1] [--ps2 PS2] [--sum SUM]
                           [--maf MAF] [--sombx SOMBX] [--variant VARIANT]
                           [--ibd IBD] [--hbd HBD] [--dem DEM] [--version]
                           module output_directory output_prefix

positional arguments:
  module             Module the program should run. Could be one of phaseset,
                     summarize, extend, somatic, or ancestry.
  output_directory   Absolute or relative path to output directory
  output_prefix      Prefix for file names in output directory. Warning:
                     existing files in output_directory with same prefix will
                     be overwritten.

optional arguments:
  -h, --help         show this help message and exit
  --bam BAM          Path to bam file
  --vcf VCF          Path to VCF file
  --vcf_id VCF_ID    Sample ID from VCF file
  --range RANGE      Genomic range chr:start-stop, chr, chr:start, chr:-stop
  --ps1 PS1          Path to first phase set file
  --ps2 PS2          Path to second phase set file
  --sum SUM          Path to existing summary file
  --maf MAF          Path to sample-specific somatic MAF (assumes all variants
                     are associated with single sample)
  --sombx SOMBX      Path to file containing barcodes supporting somatic MAF
                     variants extracted from BAM
  --variant VARIANT  Path to file containing newline-separated variant IDs,
                     format CHROM:POS:REF:ALT (ALT is comma separated list of
                     each ALT variant)
  --ibd IBD          Path to file reporting IBD (identical-by-descent)
                     segments, reported in Refined-IBD format
  --hbd HBD          Path to file reporting HBD (homozygous-by-descent)
                     segments, reported in Refined-IBD format
  --dem DEM          Demographic information about reference population used
                     in IBD analysis. Tab-separated columns:
                     sample/pop/super_pop/sex
  --version          show program's version number and exit
  ```
  
## Input requirements

| module | arguments |
| ------ | --------- |
| phaseset | `--bam`, `--vcf`, `--vcf_id`, `--range`|
| summarize | `--ps1` |
| extend | `--sum`, `--ps1`, `--ps2`, `--range` |
| somatic | `--ps1`, `--range`, `--maf` xor `--variant`, `--sum` |
| ancestry | `--ps1`, `--vcf`, `--vcf_id`, `--range`, `--ibd`, `--hbd`, `--dem` |



## Test data

Download test data files (test_data.tar.gz) from [figshare](doi.org/10.6084/m9.figshare.12295883) (size: 1.3 Gb) and untar using this command.

```
tar xvf test_data.tar.gz
```

The following script runs `SomaticHaplotype.py` on data in `test_data/` and creates output in `test_output/`.
```
bash run_SomaticHaplotype_on_test_data.sh
```
