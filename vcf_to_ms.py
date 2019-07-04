import msprime
import numpy as np
import math
import os
import argparse
import time
import re
import random
import sys
from multiprocessing import Process,Manager

#RUN COMMANDpython
#python3 vcf_to_ms.py --vcf total_chroms.vcf --numb_ind 30

parser = argparse.ArgumentParser()

parser.add_argument('--vcf',nargs='+',type=str)
parser.add_argument('--numb_ind',nargs='+',type=int)


args = parser.parse_args()


VCF_file=open('{}'.format(args.vcf[0]),'r')
MS_FILE=open('MS_FORMAT_OUT','w')
NUM_INDIVID=args.numb_ind[0]


def transform_genotypes(genotypelist):
    new_genotypelist=[]
    for geno in genotypelist:
        zygocity=geno.count('0')
        if zygocity==0:
            new_genotypelist.append(1)
            new_genotypelist.append(1)
        if zygocity==1:
            new_genotypelist.append(0)
            new_genotypelist.append(1)
        if zygocity==2:
            new_genotypelist.append(0)
            new_genotypelist.append(0)
    return new_genotypelist


begin=0
counter=0
CHUNK=[]
for line in VCF_file:
    if line[0]!='#' and line[0]!='!' and line[0]!='\\':
    #only lines with genotypes
        line=line.strip().split()
        genotypes=line[9:]
        if len(genotypes)!=NUM_INDIVID:
        #Error
            print('Number of Individuals not matching file!!!')
            break
        newgenotypes=transform_genotypes(genotypes)
        CHUNK.append(newgenotypes)
        counter+=1
        end=float(line[1])
        
        if end-begin>=1000000:
        #Chunk is larger than 100000 base length ,write out to ms style
            CHUNK=np.asarray(CHUNK)
            print(CHUNK.shape)
            CHUNK=[]
            begin=end
        if counter>=10000:
            break