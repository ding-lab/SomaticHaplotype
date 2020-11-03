# SomaticHaplotype

Somatic mutations occur on specific haplotypes. Short-read sequencing methods obscure the original haplotype relationship and important information about the cis/trans relationship of somatic events is lost. This work utilizes linked-read WGS data to leverage germline haplotype structures and infer the haplotypic context of somatic mutations.

## Installation

This repository contains submodules. Use this command when you clone: 

`git clone --recurse-submodules https://github.com/ding-lab/SomaticHaplotype`

If you have already cloned this repository without using --recurse-submodules, you can retroactively initialize and update the submodules recursively using

`git submodule update --init --recursive`

Source: https://git-scm.com/book/en/v2/Git-Tools-Submodules

## How to run

```python SomaticHaplotype.py [module] [output directory] [output prefix]```

or

```python SomaticHaplotype.py --help```
