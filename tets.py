import msprime
import numpy as np
import math



N1=100
N_2_original=100 #Drastika Megethi
N12=200 #MEGETHOS progonikou plhthismou 1-2
N3=100

generation_time = 20 #Metatroph Xronou se genies kai anapoda

T_Split_1_2=1000/generation_time #edw topotheto ta gegonota xronika
T_Split_2_3=5000/generation_time


r_2=0.1



N2= N_2_original/ math.exp(-r_2 * (T_Split_1_2/generation_time)) #tupos upologismou megethous apo prohgoumeno megethos * xrono * ruthmo auxishs

print(r_2,N_2_original,N2)

###############################################################################################################################################################

#Orismos Plhthuismwn
population_configurations = [
    msprime.PopulationConfiguration(initial_size=N1,growth_rate=0),
    msprime.PopulationConfiguration(initial_size=N2,growth_rate=r_2),
    msprime.PopulationConfiguration(initial_size=N3,growth_rate=0)
]


#migration matrix
#Sto matrix h thesh 1(grammh),2(kolwna) adiproswpeuei thn metanasteush Plhthismo 1- Plhthismo2 
#H metanasteush 1 - 2 shmainei to pososto tou Plhth 1 pou einai migrants apo ton plhth 2 
# auto to matrix isxuei apo to 'shmera' mexri na to allaxei o xrhsths sto femographic events. O monos allos tropos na diagrafoun stoixeia einai otan exafanizetai enas plhtismos
migration_matrix = [
[0.0,0.0,0.0],
[0.0,0.0,0.05],
[0.0,0.05,0.0],
]


#deigmatoleipsia
#orizoume posa samples theloume apo kathe plithismo, einai arithms aploidwn deigmatwn
#ara S1=20 ,20 aploideis h 10 diplohseis
S1=20
S2=20
S3=20




#To sample (0,0) shmainei ena atomo apo ton plhtismo 0 (1os opws ton orisame sto population configuration) se xrono 0 (paron) 
#to pollaplasiazw me twn arithmo twn smaples S1

#To sample (1,0) shmainei ena atomo apo ton plhtismo 1 (2os opws ton orisame sto population configuration) se xrono 0 (paron)
samples=[msprime.Sample(0,0)]*S1 + [msprime.Sample(1,0)]*S2 + [msprime.Sample(2,0)]*S3 






#To pio shmadiko, ola ta events elenxontai apo edw
#Sta demographic events xekiname apo to paron pros to parelthon
demographic_events = [


msprime.MassMigration(time=T_Split_1_2,source=0,destination=1,proportion = 1.0), #Otan enwnw enan plhtismo se enan allon xrhsimopoiw to mass migration me proportion 1.0, dhldi to 100% tou plhthismou 1 bainei ston 2
msprime.PopulationParametersChange(time=T_Split_1_2,initial_size=N12,population_id=1), #thn stigmh Tsplit thelw o plhthismos 2 na allaxei kai se megethos, gia na adikatoptrizei ton progoniko plhthismo
msprime.MigrationRateChange(time=T_Split_1_2 , rate=0, matrix_index=(1,2)), #episis thelw na stamatisei na exei migration me ton 3 (h python thewrei to 0 ws 1, opote to (1,2) einai h 2h grammh,3h kolwna)
msprime.MigrationRateChange(time=T_Split_1_2 , rate=0, matrix_index=(2,1)), #mhdenizw to migration stis theseis sto matrix 
msprime.PopulationParametersChange(time=T_Split_1_2,growth_rate=0,population_id=1), #kai mhdenizw kai to growth tou


msprime.MassMigration(time=T_Split_2_3,source=1,destination=2,proportion = 1.0), #otikaname apo panw
msprime.PopulationParametersChange(time=T_Split_2_3,initial_size=N12,population_id=2) #thn stigmh Tsplit thelw o plhthismos 2 na allaxei kai se megethos, gia na adikatoptrizei ton progoniko plhthismo

]




SIMULATION= msprime.simulate(samples=samples, population_configurations=population_configurations, migration_matrix=migration_matrix, mutation_rate=1.45e-8,length=1e6, recombination_rate=2e-8,demographic_events=demographic_events)

FILE=open('mynewvcf.vcf','w') #anoigw ena arxeio gi na grapsw tous gonotypous
chromosome='LOL' #dinw ena onoma xrwmoswmatos gia to vcf pou paragetai
SIMULATION.write_vcf(FILE,2,chromosome) #dinw to simulation pou etrexe, to tonoma tou arxeiou vcf pou thata apothikeusw + to 2 shmainei diplohdhs





