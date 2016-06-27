### UOD: Unique Oligonucleotide Design algorithm; a part of Thermo-Align tool for the design of template specific hybridization and priming oligonucleotides
### Version-1.05: 06/10/2016
### Authors: Felix Francis (felixfrancier@gmail.com); Randall J. Wisser (rjw@udel.edu) 

############################################################
#Time to run the code: start timer
############################################################
import time
t0 = time.time()

############################################################
#### IMPORT FUNCTIONS
############################################################
import math
import re
from Santalucia_NN_Tm import NN_Tm, complement, mM_monovalent
import itertools
import primer3
import subprocess as sp                                   
from operator import itemgetter
import shutil
import os.path
import glob
from time_stamp import Time_stamp
from parameters import *
from sequence_processing_functions import rev_complement, gc_content, fasta_to_seq

############################################################
#### FUNCTIONS
############################################################

### Function to check if a string is comprised of only DNA input_sequences (ACGT)
def check_compostion(string):
    allowed_chars    =    set('ACGT')
    if set(string.upper())    <=    allowed_chars:
        return True
    else:
        return False

### Function to check if >= 4 continuous mono/single base repeats (example: ATATAT/CCCC); returns 1 if it has repeats. else 0
def di_single_nucleo_repeat_filter(seq, filter_di_si_repeats = 1, di_si_repeats_threshold = 4):
    if filter_di_si_repeats == 1:
        filter_condition    =   0    
        seq = seq.upper()
        seq_non_rep =   ''.join(ch for ch, _ in itertools.groupby(seq))
        set_seq_bases   =   set(seq)
        for base in set_seq_bases:
            if (base*di_si_repeats_threshold) in seq:
                filter_condition    =    1
                break
        seq_len    =    len(seq)
        if seq_len%2==0:
            for i in range(0, seq_len, 2):                                
                if seq[i:i+2]*di_si_repeats_threshold in seq:
                    filter_condition    =    1
                    break        
        else:
            for i in range(0, (seq_len-1), 2):
                if seq[i:i+2]*di_si_repeats_threshold in seq:
                    filter_condition    =    1
                    break
    else :
        filter_condition    =    0
    return filter_condition                      

### FUNCTION TO CALCULATE HAIRPIN TM
def hairpin_Tm(primer_sequence, mv_cation=0,primer_conc=0): 
    Tm_hairpin =  (primer3.calcHairpin(primer_sequence,mv_conc=mv_cation, dv_conc=0, dntp_conc=0, dna_conc=primer_conc, temp_c=37, max_loop=30)).tm
    return ("{0:.2f}".format(round(Tm_hairpin,2)))    

### FUNCTION TO CALCULATE HOMODIMER TM
def homodimer_Tm(primer_sequence, mv_cation=0,primer_conc=0): 
    Tm_homodimer = (primer3.calcHomodimer(primer_sequence,mv_conc=mv_cation, dv_conc=0, dntp_conc=0, dna_conc=primer_conc, temp_c=37, max_loop=30)).tm
    return ("{0:.2f}".format(round(Tm_homodimer,2)))        

### FUNCTION TO CHECK IF BOTH THE ENDS OF PRIMERS HAS G/C (if not, they will be filtered) if filter based on A/T end condition is true  
def check_ATends(sequence, filter_AT_3prime ):  #filter_AT_3prime= 0(do not filter); filter_AT_3prime= 1 (filter)
    if filter_AT_3prime == 1:
        sequence        = sequence.upper()
        sequence        =    sequence[:-1]+sequence[-1].replace("C","G")
        if sequence [-1]  == "G":    
            return "noATend_filter"
        else:
            return "ATend_filter"
    elif filter_AT_3prime== 0:
        return "noATend_filter"

def check_GC_clamp(sequence, filter_GC_clamp ):
    if filter_GC_clamp == 1:
        sequence        = sequence.upper()
        sequence        = sequence[-5:].replace("C","G")
        if sequence.count('G')    > 3:    
            return "GC_clamp_filter"
        else:
            return "noGC_clamp_filter"
    elif filter_AT_3prime== 0:
        return "noGC_clamp_filter"

### Chopped primer filter function
def primer_filter (primer, minTm, maxTm,  GC_range_min,GC_range_max, self_Tmdiff ,filter_AT_3prime, primer_dict, i, strand,  filter_di_si_repeats,  filter_GC_clamp):
        if check_compostion(primer):

            primer_GC    =    gc_content(primer)        
            if primer_GC >= GC_range_min and primer_GC <= GC_range_max:
            
                primer_end_filter = check_ATends(primer, filter_AT_3prime)
                if primer_end_filter =="noATend_filter":            
                            
                    primer_GC_clamp_filter = check_GC_clamp(primer, filter_GC_clamp )
                    if primer_GC_clamp_filter =="noGC_clamp_filter":            
                    
                        primer_di_single_repeat_filter    =    di_single_nucleo_repeat_filter(primer, filter_di_si_repeats = filter_di_si_repeats, di_si_repeats_threshold = di_si_repeats_threshold)
                        if primer_di_single_repeat_filter    ==0:
                    
                            primer_Tm    = float(NN_Tm(seq=primer, compl_seq=complement(primer), primer_conc=primer_conc, Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs, ion_corr=True))        
                            if primer_Tm >=minTm and primer_Tm <=maxTm:                
                        
                                Tm_hairpin            =    float(hairpin_Tm(primer, monovalent_cation_eq, primer_conc))
                                Tm_homodimer    =    float(homodimer_Tm(primer, monovalent_cation_eq, primer_conc))
                                
                                if (primer_Tm - Tm_hairpin ) >=  float(self_Tmdiff) and (primer_Tm - Tm_homodimer ) >=  float(self_Tmdiff) :                                   
                                    if primer not in primer_dict:            
                                        primer_dict[primer] = [i, (primer_size), 1, primer_Tm, strand]
                                    else:
                                        primer_dict[primer][2] +=1

### Chop input_sequence into primers of a given size range
def chop_input_seq(input_seq, primer_size, minTm=52, maxTm = 58,  GC_range_min = 0,GC_range_max = 0, self_Tmdiff = 10,  filter_AT_3prime=filter_AT_3prime,  filter_GC_clamp= filter_GC_clamp, flanking_primers=flanking_primers, target_start_pos=0,target_stop_pos=0):  # the given values are the default, if no user input is given
    primer_dict    = {}
    for i in xrange(len(input_seq)-primer_size+1):
        if flanking_primers == 1:
            if i <= target_start_pos:
                primer          =   input_seq[i:i+primer_size].upper()
                primer_filter(primer, minTm, maxTm, GC_range_min,GC_range_max, self_Tmdiff, filter_AT_3prime, primer_dict, i=i, strand = "+", filter_di_si_repeats=filter_di_si_repeats, filter_GC_clamp= filter_GC_clamp) #+ strand    
            elif (i+primer_size) >= target_stop_pos:
                primer          =   input_seq[i:i+primer_size].upper()
                rev_comp_primer =   rev_complement(primer)              
                primer_filter (rev_comp_primer, minTm, maxTm,  GC_range_min,GC_range_max, self_Tmdiff, filter_AT_3prime, primer_dict, i=i, strand = "-",  filter_di_si_repeats=filter_di_si_repeats,  filter_GC_clamp= filter_GC_clamp) # -  strand
        else:
            primer          =   input_seq[i:i+primer_size].upper()
            rev_comp_primer =   rev_complement(primer)
            primer_filter (primer, minTm, maxTm,  GC_range_min,GC_range_max, self_Tmdiff, filter_AT_3prime, primer_dict, i=i, strand = "+",  filter_di_si_repeats=filter_di_si_repeats,  filter_GC_clamp= filter_GC_clamp)  # + strand                     
            primer_filter (rev_comp_primer, minTm, maxTm,  GC_range_min,GC_range_max, self_Tmdiff, filter_AT_3prime, primer_dict, i=i, strand = "-",  filter_di_si_repeats=filter_di_si_repeats,  filter_GC_clamp= filter_GC_clamp) # -  strand
    
    return primer_dict

############################################################
#### CODE
############################################################

if __name__ == '__main__':
    path    =   Time_stamp + "/"
    ### Flanking primer conditon
    if flanking_primers == 1:
        target_start_pos    = flanking_size
        target_stop_pos     = flanking_size+(stop_pos - start_pos)
        start_pos   = start_pos - flanking_size
        stop_pos    = stop_pos + flanking_size
    else:
        target_start_pos    = 0
        target_stop_pos     = 0
    
    ### Locus file name
    locus       =    path + Time_stamp + "_" + str(chr_no) + "_" + str(start_pos) + "_" + str((stop_pos - start_pos)+1) +"_VariantMasked.fasta"
    ### Query file path(SWGA primer path)
    query_path      = os.getcwd()   # after blastn db creation, need to get back to current working directory

    ### Compute the monovalent cation equivalent
    monovalent_cation_eq    =    mM_monovalent(Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs)
      

    ############################################################
    #### OUTPUT FILES
    ############################################################
    FRprimer        = path + Time_stamp + "_" +'UPD_out1_1.fasta'
    output1         = open (FRprimer, 'w')
    primerFR_list_Tm_len_filtered_with_info = path + Time_stamp + "_" +'UPD_out1_2.txt'   
    output1_1       = open( primerFR_list_Tm_len_filtered_with_info, 'w')
    primerF_list_Tm_len_exact_match_filtered = path + Time_stamp + "_" +'UPD_out2_1.fasta'     
    output2_1       = open(primerF_list_Tm_len_exact_match_filtered, 'w')
    primerR_list_Tm_len_exact_match_filtered =  path + Time_stamp + "_" +'UPD_out2_2.fasta'
    output2_2       = open(primerR_list_Tm_len_exact_match_filtered, 'w')
    primerF_list_Tm_len_exact_match_filtered_withTmCoord =  path + Time_stamp + "_" +'UPD_out3_1.fasta'     # only for reference
    output3_1       = open( primerF_list_Tm_len_exact_match_filtered_withTmCoord, 'w')                      # only for reference
    primerR_list_Tm_len_exact_match_filtered_withTmCoord =  path + Time_stamp + "_" +'UPD_out3_2.fasta'     # only for reference
    output3_2       = open(primerR_list_Tm_len_exact_match_filtered_withTmCoord, 'w')                       # only for reference
    pooled_f_r_primer =   path + Time_stamp + "_" +'UPD_out4_1.fasta'  
    output4         = open( pooled_f_r_primer, 'w')
    
    parameters_used = open( path + Time_stamp +'_run_summary.txt', 'a')

    #### If a chr all merged file (b73v2all.fasta) does not exist, merge all chr files            
    if not os.path.exists(refDB_path+"b73v2all.fasta"):
        ### Merge all chromosome files
        os.chdir(refDB_path)
        with open(refDB_path+"b73v2all.fasta", 'w') as outfile:
            infiles = glob.glob("chr*.fasta")
            for file in infiles:
                shutil.copyfileobj(open(file), outfile)

    #### Make a database of the whole genome (if all db files do not exist)
    os.chdir(refDB_path)
    file_set = {"zmv2all.nin", "zmv2all.nhr", "zmv2all.nsq", "zmv2all.nsi", "zmv2all.nsd", "zmv2all.nog"}
    if set(glob.glob("zmv2all.*")) < file_set:
        f0 = open(os.devnull, 'w')
        sp.call(["makeblastdb","-in","%sb73v2all.fasta" %refDB_path, "-dbtype","nucl","-parse_seqids","-out","%szmv2all" %refDB_path], stdout=f0,stderr=f0)

    os.chdir(query_path)
    ### Get genomic start coord from parameters file
    actual_locus_start_pos        =    start_pos

    ### Process input loci sequence file
    ## If input file is .fasta file:
    with open (locus) as sequence_data:
        line = sequence_data.read()
        lines = line.split("\n")
        input_seq = lines[1]

    ## GET GC RANGE MIN AND MAX
    GC_range_values =   [size.strip() for size in GC_range.split('-')]
    if len(GC_range_values) >2 or len(GC_range_values) <2 or GC_range_values[0] == '0' or GC_range_values[1] == '0' :
        raise ValueError('Check your input GC values. Should only have min and max, non zero values! Given GC_range = ', GC_range)
    minGC    =  float(min(GC_range_values) )
    maxGC    =  float(max(GC_range_values) )
        
    ## GET Tm RANGE MIN AND MAX
    Tm_range_values =   [size.strip() for size in Tm_range.split('-')]
    if len(Tm_range_values) >2 or len(Tm_range_values) <2 or Tm_range_values[0] == '0' or Tm_range_values[1] == '0' :
        raise ValueError('Check your input Tm values. Should only have min and max, non zero values! Given Tm_range = ', Tm_range)
    minTm   =   float(min(Tm_range_values) )
    maxTm   =   float(max(Tm_range_values) )
        
    ############################################################
    ### Chop primers
    ############################################################
    primer_size_range    =     [size.strip() for size in primer_size_range.split('-')]
    if len(primer_size_range) >2 or primer_size_range[0] == '0':
        raise ValueError('primer_size_range should only have a min and max value! Given primer_size_range = ', primer_size_range)
    if int(primer_size_range[0]) < 7:
        raise ValueError('Minimum primer_size_range should be >= 12! Given primer_size_range = ', primer_size_range)
    if int(primer_size_range[0]) >= int(primer_size_range[1]):
        raise ValueError('Provide the minimum and maximum primer size in the format "12-14" ! Given primer_size_range = ', primer_size_range)
    if len(primer_size_range) == 1:
        final_primer_dict        =    chop_input_seq(input_seq, primer_size, minTm, maxTm, GC_range_min = minGC, GC_range_max = maxGC, self_Tmdiff=self_Tmdiff, filter_AT_3prime=filter_AT_3prime, filter_GC_clamp= filter_GC_clamp, flanking_primers=flanking_primers,target_start_pos=target_start_pos,target_stop_pos=target_stop_pos)
    elif len(primer_size_range) == 2:
        min        =    int(primer_size_range[0])
        max        =    int(primer_size_range[1])
        final_primer_dict = {}
        for primer_size in xrange(min,max+1):
            primer_dict            =    chop_input_seq(input_seq, primer_size, minTm, maxTm, GC_range_min = minGC, GC_range_max = maxGC, self_Tmdiff = self_Tmdiff, filter_AT_3prime=filter_AT_3prime, filter_GC_clamp= filter_GC_clamp, flanking_primers=flanking_primers, target_start_pos=target_start_pos,target_stop_pos=target_stop_pos)
            final_primer_dict.update(primer_dict)

    no_primers_designed =   len(final_primer_dict)

    ### Chopped primers written to output1 in order, based on their start coordinates. This will be used as the input for blastn to find exact matches
    ### Write sorted output (based on primer start position) to file (primer_list_Tm_len_filtered.fasta) 
    sorted_final_primer_dict    =    sorted(final_primer_dict.iteritems(), key=itemgetter(1), reverse=False)
    output1_1.write("Primer_seq"+"\t"+"Primer_size"+"\t"+"Primer_occurrence"+"\t"+"Primer_Tm"+"\t"+"Strand"+"\n")
    for primer_info in sorted_final_primer_dict:
        primer            =    primer_info[0]
        primer_start    =    int(primer_info[1][0])
        primer_size        =    int(primer_info[1][1])
        primer_occurrence        =    primer_info[1][2]
        primer_Tm        =    primer_info[1][3]
        primer_strand    =    primer_info[1][4]
        output1.write(">TA_"+str(chr_no)+"_"+str(primer_start+actual_locus_start_pos)+"_"+str(primer_size)+"_"+str(primer_Tm)+"_"+str(primer_strand)+"\n")
        output1.write(primer +"\n")
        output1_1.write(str(primer)+"\t"+str(primer_size)+"\t"+str(primer_occurrence)+"\t"+str(primer_Tm)+"\t"+str(primer_strand)+"\n")
    output1.close()
    output1_1.close()


    ############################################################
    ### Identify exact matches
    ############################################################
    ### Get exact matches if any for the primers, in non target genomic region
    fasta_input_file    =    FRprimer    
    ### Use the smallest value in the primer size range as the word size
    word_size   =    min

    p1 = sp.Popen(["blastn","-task","blastn","-db","%szmv2all" %refDB_path,"-query","%s" %fasta_input_file,"-evalue","%s" %em_e_value,"-word_size","%s" %word_size,"-gapopen","%s" %em_gapopen,"-gapextend","%s" %em_gapextend,"-reward","%s" %em_reward,"-penalty","%s" %em_penalty,"-dust","no","-perc_identity","%s" %em_perc_identity,"-max_target_seqs","%s" %em_max_target_seqs,"-max_hsps","%s" %em_max_hsps,"-outfmt","10 qseq qlen qseqid sacc sstart send sstrand", "-num_threads","%s" %em_num_threads],stdout=sp.PIPE)

    exact_match_output, error = p1.communicate()

    ############################################################
    ### Print parameters used
    ############################################################
    parameters_used.write   (  "##########################################################"+"\n"+
                    "### Summary of UOD parameters"+"\n"+
                    "### Start date: " + str(time.strftime("%m/%d/%Y"))+";"+ "\t"+ "Start time: "+ str(time.strftime("%H:%M:%S"))+"\n"+
                    "##########################################################"+"\n"+
                    "Primer_size_range          =   "+str(primer_size_range[0]) + "-" + str(primer_size_range[1])+"\n"+
                    "GC_range                   =   "+str(GC_range)+"\n"+
                    "Tm_range                   =   "+str(Tm_range)+ " (C)" +"\n"+
                    "Filter_A/T_3prime          =   "+str(filter_AT_3prime)+ " (1:Yes; 0:No)"+"\n"+
                    "Di_si_repeats_threshold    =   "+str(di_si_repeats_threshold) +"\n"+
                    "Filter_di/si_repeats       =   "+str(filter_di_si_repeats)+ " (1:Yes; 0:No)"+"\n"+
                    "Filter_GC_clamp            =   "+str(filter_GC_clamp)+ " (1:Yes; 0:No)"+"\n"+
                    "Self_Tm difference         =   "+str(self_Tmdiff)+ " (C)"+"\n"+"\n"+
                    "## NN Tm parameters ##"+"\n"+
                    "Primer_conc                =   "+str(primer_conc)+"     (nM)"+"\n"+
                    "Na                         =   "+str(Na)+"      (mM)"+"\n"+            
                    "K                          =   "+str(K)+"      (mM)"+"\n"+    
                    "Tris                       =   "+str(Tris)+"      (mM)"+"\n"+    
                    "Mg                         =   "+str(Mg)+"     (mM)"+"\n"+    
                    "dNTPs                      =   "+str(dNTPs)+"     (mM)"+"\n"+"\n"+
                    "## blastn parameters for exact match search"+"\n"+
                    "Query_file                 =   "+str(fasta_input_file)+"\n"+
                    "e_value                    =   "+str(em_e_value)+"\n"+
                    "Word_size                  =   "+str(word_size)+"\n"+
                    "Gapopen                    =   "+str(em_gapopen)+"\n"+
                    "Gapextend                  =   "+str(em_gapextend)+"\n"+
                    "Reward                     =   "+str(em_reward)+"\n"+
                    "Penalty                    =   "+str(em_penalty)+"\n"
                    "%_identity                 =   "+str(em_perc_identity)+"\n"+
                    "Max_target_seqs            =   "+str(em_max_target_seqs)+"\n"+
                    "Max_hsps                   =   "+str(em_max_hsps)+"\n"+
                    "Number_threads             =   "+str(em_num_threads)+"\n"+
                    "##########################################################"+"\n")

    ### Parse the exact match output to get a set of primers with exact matches in the genome
    exact_match_set    =    set()
    for exact_match_output_line in exact_match_output.split('\n')[:-1]:
        exact_match_output_line = exact_match_output_line.strip(' ').split(',')
        Primer            =    exact_match_output_line[0]
        qseqid            =    exact_match_output_line[2].split('_')
        qseq_chr        =    int(qseqid[1])
        qseq_start    =    int(qseqid[2])
        qseq_strand    =    qseqid[5]
        qseq_stop    =    int(qseq_start) + int(qseqid[3])
        targetseq_chr    =    int(exact_match_output_line[3].split(':')[2])
        targetseq_start    =    int(exact_match_output_line[4])
        targetseq_stop    =    int(exact_match_output_line[5])
        alignment_length    =    len(Primer)
        query_length    =    int(exact_match_output_line[1])
        if alignment_length == query_length:
            if qseq_chr != targetseq_chr:
                exact_match_set.add(Primer)
            if qseq_chr == targetseq_chr:
                if qseq_strand == "+":    
                    if qseq_start != targetseq_start and qseq_stop != targetseq_stop:    
                        exact_match_set.add(Primer)
                if qseq_strand == "-":    
                    if qseq_start != targetseq_stop and qseq_stop != targetseq_start:    
                        exact_match_set.add(Primer)
                
    ### Remove from original dictionary, those primers with exact matches elsewhere in the genome
    for primer_exact_match in exact_match_set:
        if primer_exact_match in final_primer_dict:
            final_primer_dict.pop(primer_exact_match, None)
    no_primers_no_exact_match   =   len(final_primer_dict)

    ### Write output to file (primer_list_Tm_len_exact_match_filtered.fasta) This will be used as the input for blastn to find mis-matches.   # later on make sure the path of this output file is also specified!!!!!
    for primer, value in final_primer_dict.items():
        primer_start_pos    =   int(value[0])
        primer_len  =   (value[1])
        primer_Tm   =        value[3]
        primer_strand   =   value[4]
        if primer_strand    == "+":
            output2_1.write(">"+primer+"\n")
            output2_1.write(primer+"\n")        
            output4.write(">"+primer+"\n")
            output4.write(primer+"\n")
            output3_1.write(">TA_"+str(chr_no)+"_"+str(primer_start_pos+actual_locus_start_pos)+"_"+str(primer_len)+"_"+str(primer_Tm)+"\n")
            output3_1.write(primer+"\n")
        elif primer_strand  == "-":
            output2_2.write(">"+primer+"\n")
            output2_2.write(primer+"\n")
            output4.write(">"+primer+"\n")
            output4.write(primer+"\n")
            output3_2.write(">TA_"+str(chr_no)+"_"+str(primer_start_pos+actual_locus_start_pos)+"_"+str(primer_len)+"_"+str(primer_Tm)+"\n")
            output3_2.write(primer+"\n")    
    output2_1.close()
    output3_1.close()
    output2_2.close()
    output3_2.close()
    output4.close()

    ############################################################
    #Time to run the code: end timer
    ############################################################
    t1 = time.time()
    total = t1-t0
    total = ("{0:.2f}".format(round(total,2)))
    parameters_used.write   (   "no. of primers designed based on filter criteria : " + str(no_primers_designed)+'\n'
                                "no. of primers without exact match               : " + str(no_primers_no_exact_match)+'\n'
                                "### UOD run duration : " + str(total) + " seconds"+'\n'
                                "##########################################################"+"\n"+
                                "\n"+"\n")
    parameters_used.close()

### Change log
### v1.02 -> v1.03
    # Removed output files used during the test phase for review
### v1.03 -> v1.04
    #Use parameters from parameters.py file, including generating the locus name for the input file
