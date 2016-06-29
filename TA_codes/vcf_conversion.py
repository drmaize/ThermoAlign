### vcf_conversion: converts a vcf file to ThermoAlign input format
### Version 1.0.0: 06/28/2016
### Authors: Felix Francis (felixfrancier@gmail.com); Randall J. Wisser (rjw@udel.edu) 

### Requirements
### vcf files should follow the format described here: https://samtools.github.io/hts-specs/VCFv4.2.pdf
### all vcf input files must be named according to the corresponding fasta input files, e.g. chr1.vcf, chr2.vcf, ... 



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


input_path   =   "../sample_vcf/"


def vcf2coords(input_file, input_path):
    file_name = input_file.split(".")[0]
    df = pd.read_csv(input_path + input_file, sep='\t', comment='#', skiprows=0, usecols=[1, 4], header=None)
    df.columns = ['coordinate', 'alternate_allele']
    df.to_csv(input_path + 'processed_vcf_files/' + file_name, sep='\t', encoding='utf-8', index=False)


def vcf2coords_all(input_path):
    if not os.path.exists(input_path + 'processed_vcf_files/'):
        os.makedirs(input_path + 'processed_vcf_files/')
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


### Version log (SemVer format)
### 0.1.0: Initial development release
### 0.1.1: Generalized to process any vcf files (v4.0 or advanced as per the 1000 genome vcf file format)
    
