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
        missingness=geno.count('.')
        if zygocity==0 and missingness==0:
            new_genotypelist.append(1)
            new_genotypelist.append(1)
        if zygocity==1 and missingness==0:
            new_genotypelist.append(0)
            new_genotypelist.append(1)
        if zygocity==2 and missingness==0:
            new_genotypelist.append(0)
            new_genotypelist.append(0)
        if missingness!=0: #REQUIRES CORRECTION
            new_genotypelist.append(0)
            new_genotypelist.append(0)
    return new_genotypelist

#gather data,split to chunks ,transform them
begin=0
counter=0
CHUNK=[] #genotypes of each chun go here
ALL_CHUNKS=[] #each transformed chunk goes in here
POSITIONS=[] #positions for each chunk go in here
ALL_POSITIONS=[] #all positions go in here
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
        if end-begin>=1000000 or begin>=end:
        #Chunk is larger than 100000 base length or we switch chromosome,write out to ms style
            CHUNK=np.asarray(CHUNK)
            CHUNK=np.transpose(CHUNK)
            ALL_CHUNKS.append(CHUNK)
            CHUNK=[]
            begin=end
        if counter>=10000:
            break
#Write to ms file
            
for CHUNK in ALL_CHUNKS:
    for atomo in CHUNK:
        atomo=np.array2string(atomo)[1:-1]
        MS_FILE.write(atomo+'\n')
     MS_FILE.write()