### Target Region Selection module: selects targeted region and masks vcf polymorphic sites
### Version-1.05: 06/10/2016
### Authors: Felix Francis (felixfrancier@gmail.com); Randall J. Wisser (rjw@udel.edu) 

### Import functions
import shutil
import datetime
import os
from Bio import SeqIO
import csv
import pandas as pd
from parameters import *
from sequence_processing_functions import gc_content
############################################################
### Time to run the code: start timer
############################################################
import time
t0 = time.time()


############################################################
### Functions
############################################################

### get locus sequence from maize b73 ref gen V2
def get_Locus(start_pos, stop_pos, chr_no, refDB_path):                                                                                                            
    fasta_sequences = SeqIO.parse(open(refDB_path+"chr"+str(chr_no)+".fasta"),'fasta')
    for fasta in fasta_sequences:
        chr_seq    = str(fasta.seq)
        Locus            =    chr_seq[(start_pos-1):(stop_pos)]
        return Locus.upper()

### mask HapMap polymorphic sites in the extracted locus (testing more efficient way)
def hapmap_mask(chr_no, start_pos, stop_pos, hapmap_path, refDB_path, write_path, time_stamp, hapmap_mask_condition):
    indel_count = 0
    snp_count   = 0
    filename = "chr" + str(chr_no)                                                              
    output_poly_masked = open(write_path + str(time_stamp)+"_"+ str(chr_no) + "_"+ str(start_pos) + "_"+ str((stop_pos-start_pos)+1) + "_PolyMasked.fasta", "w")
    output_poly_masked.write(str(">"+"TA_"+ str(chr_no) + "_"+ str(start_pos) + "_"+ str((stop_pos-start_pos)+1) + "_"+"HapMapPolyMasked"+'\n'))
    Locus = get_Locus(start_pos, stop_pos, chr_no, refDB_path)
    gap_count = Locus.count('N'*100)
    locus_gc = gc_content(Locus)
    if hapmap_mask_condition == 1:
        df = pd.read_csv(hapmap_path+filename, sep='\t', header=0) 
        df = df[(df.coordinate >= int(start_pos)) & (df.coordinate <= int(stop_pos))]
        for index, row in df.iterrows():
            if '<INS>' in row['alternate_allele'] or '<DEL>' in row['alternate_allele']:
                poly_coord    =    int(row['coordinate'])-(start_pos-1)
                Locus = Locus[:(poly_coord-1)] + "n" + Locus[poly_coord:]                  # Indels annotated as "n"
                indel_count += 1
            else :
                poly_coord    =    int(row['coordinate'])-(start_pos-1)
                Locus = Locus[:(poly_coord-1)] + "N" + Locus[poly_coord:]                  # Polymorphisms other than indels annotated as "N"
                snp_count += 1
    output_poly_masked.write(str(Locus))
    output_poly_masked.close()
    return [indel_count, snp_count, gap_count, locus_gc]

############################################################
### Run
############################################################
if __name__ == '__main__':

    ### Flanking primer conditon
    if flanking_primers == 1:
        user_given_start_pos    = start_pos
        user_given_stop_pos     = stop_pos
        start_pos   = start_pos - flanking_size
        stop_pos    = stop_pos + flanking_size

    ### create a directory with timestamp name
    time_stamp  =   str("TA_"+ datetime.datetime.now().isoformat()).replace(":", "_").replace(".", "_")
    out_path    =   './'+str(time_stamp) +'/'
    os.makedirs(out_path)
    
    ### write timestamp to parameters.py
    record_timestamp = open('parameters.py', 'a')
    record_timestamp.write("time_stamp_used          =   "+ "'" + str(time_stamp)+ "'"+"\n")
    record_timestamp.close()
    
    ### create temp timestamp file
    timestamp_file = open("time_stamp.py", "w")
    timestamp_file.write("Time_stamp          =   "+ "'" + str(time_stamp)+ "'"+"\n")
    timestamp_file.close()
    
    ### write parameters used out put file
    parameters_used = open(out_path +str(time_stamp)+'_run_summary.txt', 'w')
    ### Print parameters used  
    parameters_used.write   (   "##########################################################"+"\n"+
                                "### TRS Locus selection & HapMap masking parameters"+"\n"+
                                "### Start date: " + str(time.strftime("%m/%d/%Y"))+";"+ "\t"+ "Start time: "+ str(time.strftime("%H:%M:%S")) +"\n"+
                                "##########################################################"+"\n"+
                                "Chromosome number      =  "+str(chr_no)+"\n")

    if flanking_primers == 1:
        parameters_used.write   (   "Locus start position   =  "+str(user_given_start_pos)+"\n"+
                                    "Locus stop position    =  "+str(user_given_stop_pos)+"\n"+
                                    "Flanking region size   =  "+str(flanking_size)+"\n")
    else:
        parameters_used.write   (   "Locus start position   =  "+str(start_pos)+"\n"+
                                    "Locus stop position    =  "+str(stop_pos)+"\n")
                                
    parameters_used.write   (   "HapMap mask condition  =  "+str(hapmap_mask_condition)+ " (1: Yes; 0: No)" + "\n"+
                                "HapMap masked output file =  "+str("TA_"+ str(chr_no) + "_"+ str(start_pos) + "_"+ str((stop_pos-start_pos)+1) + "_PolyMasked.fasta")+"\n"+"\n"+
                                "## Whole genomic assembly and HapMap files input path"+"\n"+
                                "Genome assembly file path  =  "+str(refDB_path)+"\n"+
                                "HapMap file path           =  "+str(hapmap_path)+"\n"+          
                                "##########################################################"+"\n")

    ### Locus extraction and HapMap masking
    indel_snp_gap_info = hapmap_mask(chr_no, start_pos, stop_pos, hapmap_path, refDB_path, out_path, time_stamp, hapmap_mask_condition) 
    indel_count, snp_count, gap_count, locus_gc = indel_snp_gap_info[0], indel_snp_gap_info[1], indel_snp_gap_info[2], indel_snp_gap_info[3]
    
    ############################################################
    #Time to run the code: end timer
    ############################################################
    t1 = time.time()
    total = t1-t0
    total = ("{0:.2f}".format(round(total,2)))  
    parameters_used.write   (  
                            "Number of indels       =  "+str(indel_count)+"\n"+
                            "Number of SNPs         =  "+str(snp_count)+"\n"+
                            "Number of assembly gaps=  "+str(gap_count)+"\n"+
                            "Locus GC               =  "+str(locus_gc)+"\n"+
                            "### TRS run duration   : " + str(total) + " seconds"+'\n'
                            "##########################################################"+"\n"+
                            "\n"+"\n"
                            )
    parameters_used.close()

    ### Copy parameters file into the timestamp directory
    shutil.copy2('./parameters.py', str('./'+time_stamp + '/parameters.py'))

    
