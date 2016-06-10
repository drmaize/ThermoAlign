### vcf processing: pre processes the vcf files so that it can be used by TA in an efficient manner
### Version-1.05: 06/10/2016
### Author: Felix Francis (felixfrancier@gmail.com); Under the guidance and financial support of Randall J Wisser (rjw@udel.edu) 

### Input requirements: all the input vcf files must be named as chr1.vcf, chr2.vcf, ...; vcf files should follow the standard vcf format specified here https://samtools.github.io/hts-specs/VCFv4.2.pdf


### Import functions
import datetime
import os
import pandas as pd


############################################################
### Time to run the code: start timer
############################################################
import time
t0 = time.time()


############################################################
### Functions
############################################################


input_path   =   "./hapmap_vcf/"


def vcf2coords(input_file, input_path):
    file_name = input_file.split(".")[0]
    df = pd.read_csv(input_path + input_file, sep='\t', comment='#', skiprows=0, usecols=[1, 4], header=None)
    df.columns = ['coordinate', 'alternate_allele']
    df.to_csv(input_path + 'processed_hapmap_files/' + file_name, sep='\t', encoding='utf-8', index=False)


def vcf2coords_all(input_path):
    if not os.path.exists(input_path + 'processed_hapmap_files/'):
        os.makedirs(input_path + 'processed_hapmap_files/')
    for file in os.listdir(input_path):
        if file.endswith(".vcf"):
            input_file = file
            vcf2coords(input_file, input_path)


############################################################
### Run
############################################################

if __name__ == '__main__':

    vcf2coords_all(input_path)


##########################################################
###Time to run the code: end timer
##########################################################
t1 = time.time()
total = t1-t0
total = ("{0:.2f}".format(round(total,2)))
print '\n', "Time to run code = ", total, " seconds", '\n'


    
    
### Change log
### v1.03 -> v1.04 vcf to coord preprocessing
### v1.04 -> v1.05 process any vcf files

    