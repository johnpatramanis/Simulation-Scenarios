import msprime
import numpy as np
import math
import os
import time
import re
import random
import sys
from multiprocessing import Process,Manager

# SCENARIO 1 : Ambracians come from Corinth , very small migration from Locals, normal migration from Corinth  -fast growth of colony-
#              Corinthians Locals small migration between them but cinverge 5k years ago
start_time = time.time()


reps=1
for REPS in range(0,reps):

##############################################################################################################################################
#Simulation Parameters
    
    parametersfile=open('PARAMETERS_{}'.format(REPS),'w')
    N_OG=2000
    N_OUT=2000
    N_AB=2000
    N_A0=2000
    N_B0=200

    r_A=0.000
    r_B=0.5

    generation_time = 25

    T_split_OUT_AB=5000/generation_time
    T_split_AB=700/generation_time




    N_A=N_A0 / math.exp(-r_A * T_split_AB)
    N_B=N_B0 / math.exp(-r_B * T_split_AB)



    population_configurations = [
        msprime.PopulationConfiguration(initial_size=N_OUT),
        msprime.PopulationConfiguration(initial_size=N_A, growth_rate=r_A),
        msprime.PopulationConfiguration(initial_size=N_B, growth_rate=r_B)
    ]



    migration_matrix = [
        [0,0.0001,0.0001],
        [0.0001,0,0.0001],
        [0.0001,0.0001,0]]

    N1=20
    N2=20
    N3=20
    samples=[msprime.Sample(0,0)]*N1 + [msprime.Sample(1,0)]*N2 + [msprime.Sample(2,0)] *N3


    demographic_events = [
        # A and B merge
        msprime.MassMigration(time=T_split_AB, source=2, destination=1, proportion=1.0),
        msprime.MigrationRateChange(time=T_split_AB, rate=0),
        msprime.MigrationRateChange(time=T_split_AB, rate=0.0001, matrix_index=(0, 1)),
        msprime.MigrationRateChange(time=T_split_AB, rate=0.0001, matrix_index=(1, 0)),
        msprime.PopulationParametersChange(time=T_split_AB, initial_size=N_AB, growth_rate=0, population_id=1),
        # Population AB merges into OUT
        msprime.MassMigration(time=T_split_OUT_AB, source=1, destination=0, proportion=1.0),
        # Size changes to N_A at T_AF
        msprime.PopulationParametersChange(time=T_split_OUT_AB, initial_size=N_OUT, population_id=0)
    ]

    
    

######################################################################################################################################################
#RUN the simulation and output genotypes in vcfs and ms format files, one for each chrom 

    
    def SIMULATE(L,argument,samples,population_configurations,migration_matrix,demographic_events):
        j=int(argument)
        recomb_map=msprime.RecombinationMap.read_hapmap('genetic_map_GRCh37_chr{}.txt'.format(j))
        dd = msprime.simulate(samples=samples,
            population_configurations=population_configurations,
            migration_matrix=migration_matrix,mutation_rate=1e-8,
            demographic_events=demographic_events,recombination_map=recomb_map)
        outfile=open('ms_prime_{}'.format(j),'w')   
        for var in dd.variants():
            L.append([j,var.index,var.position])
            for genotype in var.genotypes:
                outfile.write(str(genotype))
            outfile.write('\n')
        outfile.close()
        wow=open('mynewvcf{}.vcf'.format(j),'w')
        dd.write_vcf(wow,2,str(j))
        wow.close()
        
        population_labels= ["locals"]*int(N1/2) + ["metropolis"]*int(N2/2) + ["apoikia"]*int(N3/2)
        d=0
        newlabels=[]
        for i in range(0,len(population_labels)):
            newlabels.append(population_labels[i]+str(d))
            d+=1
            if i==len(population_labels)-2:
                newlabels.append(population_labels[i]+str(d))
                break
            if population_labels[i]!=population_labels[i+1]:
                d=0
        population_labels=newlabels
        wow=open('mynewvcf{}.vcf'.format(j))
        wowzers=open('myvcf{}.vcf'.format(j),'w')
        for line in wow:
            line=line.strip().split()
            if line[0]=='#CHROM':
                line[9:]=population_labels
            wowzers.write("\t".join(line))
            wowzers.write("\n")
        wow.close()
        
        return j,L
    
    
    
    L=[]
    if __name__ == '__main__':
        with Manager() as manager:
            L=manager.list(L)
            processes=[]
            for loop in range(1,23):
                p=Process(target=SIMULATE,args=(L,loop,samples,population_configurations,migration_matrix,demographic_events,))
                processes.append(p)
                
                p.start()
        
                
            for p in processes:
                p.join()
            #print(len(L),'1')
            sys.stdout.flush()
            variants=sorted(list(L))


    variantinfo=['{}\t{}\t{}\n'.format(x[0],x[1],x[2])for x in variants]
    print(len(variants),len(variantinfo))

    variantinformation=open('variants_info.txt','w')
    variantinformation.write('CHROM\tVARIANT\tPOSITION\n')
    for loop in variantinfo:
        variantinformation.write(loop)
    
    variantinformation.close()

    elapsed_time_1 = time.time() - start_time        
        
    print('Step 1 : {} '.format(elapsed_time_1/60))        

#Arrange the vcf files into one, fix labels , bed file , transform to eigen, calculate pca and f stats
        
    os.system('rm mynewvcf*.vcf')
    os.system('bcftools concat -o total_chroms.vcf myvcf*.vcf')
    os.system('rm myvcf*.vcf')


    VCF=open('total_chroms.vcf','r')
    newVCF=open('newtotal_chroms.vcf','w')

    snpcount=0
    
    variants=sorted(variants)
    for line in VCF:
        if line[0]!='#' and snpcount<len(variants):
            line=line.strip().split()
            line[2]='rs{}'.format(snpcount)
            line[1]=str(variants[snpcount][2])
            line.append('\n')
            line='\t'.join(line)
            snpcount+=1
        newVCF.write(line)

    VCF.close
    newVCF.close
    os.system('rm ms_prime*')
    os.system('mv newtotal_chroms.vcf total_chroms.vcf')

