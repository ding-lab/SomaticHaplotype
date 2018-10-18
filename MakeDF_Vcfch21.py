#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 01:33:24 2018

@author: jessikabaral
"""

import vcf

import pysam
import vcf.filters
import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

#in VCFs, the rows represent variants and the columns represent patient samples

#insert the pathway to the file in katmai or denali
vcf_reader = vcf.Reader(open('/diskmnt/Projects/SomaticHaplotype/Data/1000G/chr21.1kg.phase3.v5a.vcf.gz.tbi', 'r'))

variantPosition_list = []
variantID_list = []
variantRef_list = []
variantAlt_list = []
variantHet_list = []
variantHomRef_list = []
variantHomAlt_list = []

for record in vcf_reader:
    if(record.is_snp):
        variantPosition_list.append(record.POS)
        variantID_list.append(record.ID)
        variantRef_list.append(record.REF)
        variantAlt_list.append(record.ALT)
        variantHet_list.append(record.num_het)
        variantHomRef_list.append(record.num_hom_ref)
        variantHomAlt_list.append(record.num_hom_alt)


df = pd.DataFrame({'Position': variantPosition_list, 'Chrom': "21", 'ID': variantID_list, 'Reference': variantRef_list, \
                   'Alternate': variantAlt_list, 'numHet': variantHet_list, \
                   'numberHomozygousRef': variantHomRef_list, 'numberHomozygousAlt': variantHomAlt_list})
    
df['PatientCount'] = (df['numberHomozygousRef'] + df['numberHomozygousAlt'] + df['numHet'])

df['MutantAlleleFrequency'] = (df['numHet'] + (2*df['numberHomozygousAlt']))/(2*(df['numberHomozygousRef'] + df['numberHomozygousAlt'] + df['numHet']))

