#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tues Nov 6 1:36 PM

@author: jessikabaral
"""
#the purpose of this code is to make a generalized pipeline for 
#generating the data table with all the VCF relevant information 
#for each vcf file for each chromosome on 1000G

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
filename_list = ['chr0.1kg.phase3.v5a.vcf.gz', 'chr1.1kg.phase3.v5a.vcf.gz', 'chr2.1kg.phase3.v5a.vcf.gz',
                 'chr3.1kg.phase3.v5a.vcf.gz', 'chr4.1kg.phase3.v5a.vcf.gz', 'chr5.1kg.phase3.v5a.vcf.gz', 
                 'chr6.1kg.phase3.v5a.vcf.gz', 'chr7.1kg.phase3.v5a.vcf.gz', 'chr8.1kg.phase3.v5a.vcf.gz',
                 'chr9.1kg.phase3.v5a.vcf.gz', 'chr10.1kg.phase3.v5a.vcf.gz', 'chr11.1kg.phase3.v5a.vcf.gz',
                 'chr12.1kg.phase3.v5a.vcf.gz', 'chr13.1kg.phase3.v5a.vcf.gz', 'chr14.1kg.phase3.v5a.vcf.gz', 
                 'chr15.1kg.phase3.v5a.vcf.gz', 'chr16.1kg.phase3.v5a.vcf.gz', 'chr17.1kg.phase3.v5a.vcf.gz',
                 'chr18.1kg.phase3.v5a.vcf.gz', 'chr19.1kg.phase3.v5a.vcf.gz', 'chr20.1kg.phase3.v5a.vcf.gz',
                 'chr21.1kg.phase3.v5a.vcf.gz', 'chr22.1kg.phase3.v5a.vcf.gz', 'chrX.1kg.phase3.v5a.vcf.gz']

chromosomes_list = ['0',
'1',
'2',
'3',
'4',
'5',
'6',
'7',
'8',
'9',
'10',
'11',
'12',
'13',
'14',
'15',
'16',
'17',
'18',
'19',
'20',
'21',
'22', 'X']
                 
i = 0

for chromosome in chromosomes_list:
    path = 'chr%s.1kg.phase3.v5a.vcf.gz' % chromosome
    #print path
    fileName = '/diskmnt/Projects/SomaticHaplotype/Data/1000G/%s' % path
    #print fileName
    vcf_reader = vcf.Reader(filename = fileName)
    
    variantPosition_list = []
    variantID_list = []
    variantRef_list = []
    variantAlt_list = []
    variantHet_list = []
    variantHomRef_list = []
    variantHomAlt_list = []
    
# =============================================================================
#     for record in vcf_reader:
#         if(record.is_snp):
#             variantPosition_list.append(record.POS)
#             variantID_list.append(record.ID)
#             variantRef_list.append(record.REF)
#             variantAlt_list.append(record.ALT)
#             variantHet_list.append(record.num_het)
#             variantHomRef_list.append(record.num_hom_ref)
#             variantHomAlt_list.append(record.num_hom_alt)
# =============================================================================
    
    for record in vcf_reader:
        if (i < 100):
            if(record.is_snp):
                variantPosition_list.append(record.POS)
                variantID_list.append(record.ID)
                variantRef_list.append(record.REF)
                variantAlt_list.append(record.ALT)
                variantHet_list.append(record.num_het)
                variantHomRef_list.append(record.num_hom_ref)
                variantHomAlt_list.append(record.num_hom_alt)
                i += 1
        else:
            break

    df = pd.DataFrame({'Position': variantPosition_list, 'Chrom': chromosome, 'ID': variantID_list, 'Reference': variantRef_list, \
                       'Alternate': variantAlt_list, 'numHet': variantHet_list, \
                       'numberHomozygousRef': variantHomRef_list, 'numberHomozygousAlt': variantHomAlt_list})
        
    df['PatientCount'] = (df['numberHomozygousRef'] + df['numberHomozygousAlt'] + df['numHet'])
    
    df['MutantAlleleFrequency'] = (df['numHet'] + (2*df['numberHomozygousAlt']))/(2*(df['numberHomozygousRef'] + df['numberHomozygousAlt'] + df['numHet']))
    
    path_underscores = path.replace(".", "_")
    fileNameString  = '/diskmnt/Projects/Users/jbaral/SomaticHaplotypeFiles/testCSV_%s.csv' % path_underscores
    df.to_csv(fileNameString, sep=',')


