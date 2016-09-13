### PSE: Primer Specificity Evaluation (PSE) algorithm; a part of Thermo-Align tool for the design of template specific hybridization and priming oligonucleotides
### Version 1.0.0: 06/28/2016
### Authors: Felix Francis (felixfrancier@gmail.com); Randall J. Wisser (rjw@udel.edu) 

############################################################
#Time to run the code: start timer
############################################################
import time
t0 = time.time()

############################################################
#### IMPORT PYTHON FUNCTIONS
############################################################
import subprocess as sp
import os, errno
import shutil
import os.path
import glob
import csv
import pandas as pd
import primer3    
import sys
import numpy as np
import re
from Santalucia_NN_Tm import NN_Tm, complement, mM_monovalent
from Bio import SeqIO
from time_stamp import Time_stamp
from parameters import *
from sequence_processing_functions import fasta_to_seq, gc_content, rev_complement

monovalent_cation_eq    =    mM_monovalent(Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs)

############################################################
#### FUNCTIONS
############################################################      

### FUNCTION TO CALCULATE HAIRPIN TM
def hairpin_Tm(primer_sequence, mv_cation=0,primer_conc=0): 
    Tm_hairpin =  (primer3.calcHairpin(primer_sequence,mv_conc=mv_cation, dv_conc=0, dntp_conc=0, dna_conc=primer_conc, temp_c=37, max_loop=30)).tm
    return ("{0:.2f}".format(round(Tm_hairpin,2)))    
    
### FUNCTION TO CALCULATE HOMODIMER TM
def homodimer_Tm(primer_sequence, mv_cation=0,primer_conc=0): 
    Tm_homodimer = (primer3.calcHomodimer(primer_sequence,mv_conc=mv_cation, dv_conc=0, dntp_conc=0, dna_conc=primer_conc, temp_c=37, max_loop=30)).tm
    return ("{0:.2f}".format(round(Tm_homodimer,2)))

### function to take each chromosome into memory
def genome_read(no_chrs_in_genome):
    for i in xrange(no_chrs_in_genome):
        chr_no = int(i)+1
        fasta_sequences = SeqIO.parse(open(refDB_path+"chr"+ str(chr_no) +".fasta"),'fasta')
        for fasta in fasta_sequences:
            genome_dict['chr'+ str(chr_no)] = [str(fasta.seq), len(fasta.seq)]


### sequence extraction from the genome
def seq_extraction(fasta_id, coordinates_3primeextend):
    chr_seq     =    "chr"+ str(fasta_id.split(':')[2])
    coordinates =    [x.strip() for x in coordinates_3primeextend.split(':')]
    start_coord =    int(coordinates[0])
    stop_coord  =    int(coordinates[1])
    if start_coord < 1:
        start_n = start_coord
        start_coord = 1
        select_seq  =    genome_dict[str(chr_seq)][0][(start_coord-1):(stop_coord)]
        select_seq  =    ('N'*start_n) + select_seq
    elif stop_coord > genome_dict[str(chr_seq)][1]:
        stop_n = stop_coord
        start_coord = int(genome_dict[str(chr_seq)][1])
        select_seq  =    genome_dict[str(chr_seq)][0][(start_coord-1):(stop_coord)]
        select_seq  =    select_seq+('N'*stop_n)
    else:
        select_seq  =    genome_dict[str(chr_seq)][0][(start_coord-1):(stop_coord)]
    return select_seq

### function to calculate continuous gc stretch
def continuous_gc(sequence):
    sequence     = sequence.upper()
    if "G" in sequence or "C" in sequence :
        sequence     = sequence.replace("C","G")
        if "GG" in sequence:
            out_put        = len(max(re.compile("(G+G)").findall(sequence))) 
        else:
            out_put        = 1
    else:
        out_put        = 0
    return out_put

### function to calculate Tm (old approach, not using NN Tm)
def Tm(sequence):
    sequence = sequence.upper()
    A = sequence.count('A')
    C = sequence.count('C')
    G = sequence.count('G')
    T = sequence.count('T')    
    Tm = (2*(A+T))+ (4*(G+C))
    return Tm

### function to count the number of mismatches from 3' end of two strands
def count_3prime_mismatches(seq1, seq2):
    num_mismatch = 0
    seq1 = (seq1[::-1]).upper()
    seq2 = (seq2[::-1]).upper()
    list = [seq1,seq2]
    for i in xrange(len(min(list, key=len))):
        if seq1[i] == seq2[i]:
            break
        elif seq1[i] != seq2[i]:    
            num_mismatch +=1
    return num_mismatch
    
### function to count the start position of mismatches from 3' end of two strands
def start_pos3prime_mismatches(seq1, seq2):    
    mismatches_start_pos3prime = "Nill"
    seq1 = (seq1[::-1]).upper()
    seq2 = (seq2[::-1]).upper()
    list = [seq1,seq2]
    for i in xrange(len(min(list, key=len))):    
        if seq1[i] != seq2[i]:
            mismatches_start_pos3prime = i+1
            break
    return mismatches_start_pos3prime
    
### function to count the hamming distance between two sequences    
def HammingDistance(seq1, seq2):
    if len(seq1) == len(seq2):
        score = 0
        seq1 = (seq1).upper()
        seq2 = (seq2).upper()
        for i in xrange(len(seq1)):
            if seq1[i] != seq2[i]:
                score = score+1
    else:
        score = "different_lengths!"
    return score
    
### function to calculate the percentile of the mipriming tm s for a particular primer
def misprime_percentile(Tm_list, misprime_Tm_percentile_value):
    Tm_list_cleaned = [ x for x in Tm_list if isinstance(x, float)]
    if len(Tm_list_cleaned) == 0:
        p = 'No match/3_p_mismatch'
    else:    
        p = np.percentile(Tm_list_cleaned, misprime_Tm_percentile_value)
        p =(np.float64(p).item())
        p = round(p,2)
    return p    
    
### function to calculate the max mipriming tm for a particular primer
def max_misprime_Tm(Tm_list):
    Tm_list_cleaned = [ x for x in Tm_list if isinstance(x, float)]
    if len(Tm_list_cleaned) == 0:
        max_Tm = "NA"
    else:
        max_Tm    = max(Tm_list_cleaned)
    return max_Tm

#########################################################
### CODE
#########################################################

if __name__ == '__main__':

    ### Query file path
    query_path  = str(os.getcwd()) +"/"
    path    =   Time_stamp + "/"
    
    three_prime_region = 5              
    misprime_Tm_percentile_value = 90
    
    ### Flanking primer conditon
    if flanking_primers == 1:
        start_pos   = start_pos - flanking_size
        stop_pos    = stop_pos + flanking_size

    ### input files
    fasta_input_file_F    = path + Time_stamp + "_" +'UOD_out2_1.fasta'  
    fasta_input_file_R    = path + Time_stamp + "_" +'UOD_out2_2.fasta'
          
    ### Locus file name
    locus       =    path + Time_stamp + "_" + str(chr_no) + "_" + str(start_pos) + "_" + str((stop_pos - start_pos)+1) +"_VariantMasked.fasta"
       
    ### Parameter out put file
    parameters_used = open( path + Time_stamp +'_run_summary.txt', 'a')

    ###take each chromosome into memory
    genome_dict ={}
    genome_read(no_chrs_in_genome)
    actual_locus_start_pos  =   start_pos
    actual_locus_stop_pos   =   stop_pos

    #### if a chr all merged files does not exist, merge all chr files            
    if not os.path.exists(refDB_path+"b73v2all.fasta"):
        ### Merge all chromosome files
        os.chdir(refDB_path)
        with open(refDB_path+"b73v2all.fasta", 'w') as outfile:
            infiles = glob.glob("chr*.fasta")
            for file in infiles:
                shutil.copyfileobj(open(file), outfile)

    #### make a database of the whole genome (if all db files do not exist)
    os.chdir(refDB_path)
    file_set = {"zmv2all.nin", "zmv2all.nhr", "zmv2all.nsq", "zmv2all.nsi", "zmv2all.nsd", "zmv2all.nog"}
    if set(glob.glob("zmv2all.*")) < file_set:
        f0 = open(os.devnull, 'w')
        sp.call(["makeblastdb","-in","%sb73v2all.fasta" %refDB_path, "-dbtype","nucl","-parse_seqids","-out","%szmv2all" %refDB_path], stdout=f0,stderr=f0)

    #### create a dictionary of all primers from the input(query) sequence
    d = dict()
    d2 = dict()
 
    with open(query_path+fasta_input_file_F, 'r') as f:
        for line in f.readlines():
            if line.startswith('>'):
                d[line[1:-1]] =    ['No match', 0,'No match', 0 ]
                d2[line[1:-1]] =[]
                
    min_pimer_len    =    len(min(d,key=len))        
    min_pooled_primer_len   =   min_pimer_len 
    word_size = 7

    #### Condition to raise an error if the user defined 3' region is > min_pimer_len
    if three_prime_region    > min_pimer_len:
            raise ValueError('3 prime region should not be greater than the minimum primer length! Given three_prime_region = ', three_prime_region, 'min_pimer_len = ', min_pimer_len)
            
    #### BLASTn alignment . . . . to get other close matches (less than 100 pid)
    f0 = open(os.devnull, 'w')
    ###New blast parameters 082415 
    p1 = sp.Popen(["blastn","-task","blastn","-db","%szmv2all" %refDB_path,"-query","%s%s" %(query_path,fasta_input_file_F),"-evalue","%s" %mp_e_value,"-word_size","%s" %word_size,"-gapopen","%s" %mp_gapopen,"-gapextend","%s" %mp_gapextend,"-reward","%s" %mp_reward,"-penalty","%s" %mp_penalty,"-dust","no","-perc_identity","%s" %mp_perc_identity,"-max_target_seqs", "%s" %mp_max_target_seqs,"-max_hsps", "%s" %mp_max_hsps,"-outfmt","10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq qseq", "-num_threads","%s" %mp_num_threads],stdout=sp.PIPE,stderr=f0) 
    output, error = p1.communicate()

    ### Print parameters used
    parameters_used.write   (  "##########################################################"+"\n"+
                    "### Summary of PSE parameters"+"\n"+
                    "### Start date: " + str(time.strftime("%m/%d/%Y"))+";"+ "\t"+ "Start time: "+ str(time.strftime("%H:%M:%S"))+"\n"+
                    "##########################################################"+"\n"+
                    "## blastn parameters for mis-prime alignment"+"\n"+
                    "Query_file_F_primers       =   "+str(fasta_input_file_F)+"\n"+
                    "Query_file_R_primers       =   "+str(fasta_input_file_R)+"\n"+
                    "e_value                    =   "+str(mp_e_value)+"\n"+
                    "Word_size_F_primers        =   "+str(word_size)+"\n"
                            )

    list_lines = []        
    db = "%szmv2all" %refDB_path            
    for line in output.split('\n')[:-1]:
        line = line.strip(' ').split(',')
        query_len   =    len(line[0])
        qstart        =    int(line[6])
        qend        =    int(line[7])
        sstart        =    float(line[8])
        send        =    float(line[9])
        pident        =    float(line[2])
        match_length    =    float(line[3])

        ### skip self alignments
        if pident == 100 and query_len == match_length:
            pass
        else:
            three_prime_mismatches    =    query_len-qend
            five_prime_mismatches    =    qstart-1
            fasta_id = line[1]
            Match_seq   = line[12]
            Query_seq   = line[13]
            match_gaps  = Match_seq.count('-') 
            query_gaps  = Query_seq.count('-') 
            if send-sstart > 0:
                match_direction = "right"
                actual_match_start = int(sstart)
                actual_match_end = int(send)
                actual_match_end_3prime_extended = actual_match_end + int(three_prime_mismatches)
                actual_match_start_5prime_extended = actual_match_start - int(five_prime_mismatches)
                actual_match_start_5prime_extended = actual_match_start_5prime_extended + query_gaps
                if int(three_prime_mismatches) == 0 and int(five_prime_mismatches) == 0 and match_gaps == 0 and query_gaps == 0:
                    match_extend_3_5 = Match_seq
                else:
                    coordinates_3_5primeextend = str(actual_match_start_5prime_extended)+":"+str(actual_match_end_3prime_extended)
                    match_extend_3_5 = seq_extraction(fasta_id, coordinates_3_5primeextend)
            else:
                match_direction = "left"
                actual_match_start = int(send)
                actual_match_end = int(sstart)
                actual_match_end_3prime_extended = actual_match_start - int(three_prime_mismatches)
                actual_match_start_5prime_extended = actual_match_end + int(five_prime_mismatches)
                actual_match_start_5prime_extended = actual_match_start_5prime_extended + match_gaps
                actual_match_start_5prime_extended = actual_match_start_5prime_extended - query_gaps
                if int(three_prime_mismatches) == 0 and int(five_prime_mismatches) == 0 and match_gaps == 0 and query_gaps == 0:
                    match_extend_3_5 = rev_complement(Match_seq) ### 040416 corrected
                else:
                    coordinates_3_5primeextend    = str(actual_match_end_3prime_extended)+":"+str(actual_match_start_5prime_extended)     
                    match_extend_3_5 = seq_extraction(fasta_id, coordinates_3_5primeextend)            
            Primer        =      line[0]    
            line.extend([match_direction])
            line.extend([three_prime_mismatches])
            line.extend([five_prime_mismatches])
            PrimerTm = Tm(Primer)
            line.append(PrimerTm)
            match_melt_temp        =     Tm(Match_seq)
            line.append(match_melt_temp)
            line.append(PrimerTm-match_melt_temp)
            fasta_id = line[1]
            qstart        =    int(line[6])-1
            qend        =    int(line[7])
            Match_seq_blast     = line[12]
            line.append(match_extend_3_5)    

            ### get actual_3prime_mismatches
            if match_direction == "right":
                Actual_3prime_mismatches    = count_3prime_mismatches(Primer, match_extend_3_5)
                Mismatch_3prime_start_pos    = start_pos3prime_mismatches(Primer, match_extend_3_5)
            elif match_direction == "left":
                Actual_3prime_mismatches    = count_3prime_mismatches(Primer, rev_complement(match_extend_3_5))
                Mismatch_3prime_start_pos    = start_pos3prime_mismatches(Primer, rev_complement(match_extend_3_5))
            elif "-" in Match_seq_blast:
                Actual_3prime_mismatches     = "gap_in_match"
                Mismatch_3prime_start_pos    = "gap_in_match"
            else:
                Actual_3prime_mismatches     = "check_your_sequences!"
                Mismatch_3prime_start_pos    = "check_your_sequences!"
            line.append(Actual_3prime_mismatches)    
            line.append(Mismatch_3prime_start_pos)
            Primer_NN_Tm = NN_Tm(seq=Primer, compl_seq=complement(Primer), primer_conc=primer_conc, Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs, ion_corr=True)
            if match_direction == "right":
                complementary_match    =    complement(match_extend_3_5)
            else:
                complementary_match    =    match_extend_3_5[::-1]
            Gap_adjusted_end_filling_Tm = NN_Tm(seq=Primer, compl_seq=complementary_match, primer_conc=primer_conc, Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs, ion_corr=True)
            Local_alignment_Tm = NN_Tm(seq=Query_seq, compl_seq=complement(Match_seq), primer_conc=primer_conc, Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs, ion_corr=True)
            line.append(Primer_NN_Tm)
            line.append(Local_alignment_Tm)  
            line.append(Gap_adjusted_end_filling_Tm)
            list_lines.append(line) 

            # Identify the orientation of the match and then calculate the respective 3' mismatches if any    
            if Actual_3prime_mismatches > 0 and (d[line[0]][0] == "3prime_mismatch" or d[line[0]][0] == "No match"):    
                d[line[0]][0] =    "3prime_mismatch"
            if Actual_3prime_mismatches == 0:
                if (d[line[0]][0]=="No match" )or (d[line[0]][0]=="3prime_mismatch") or (float(Gap_adjusted_end_filling_Tm) > float(d[line[0]][0])) :    
                    d[line[0]][0] =    Gap_adjusted_end_filling_Tm    
                if float(Gap_adjusted_end_filling_Tm) > float(d[line[0]][0]) :    
                    d[line[0]][3] = 1
                    if Mismatch_3prime_start_pos > three_prime_region:
                        d[line[0]][1]    += 1                    
                if  float(Gap_adjusted_end_filling_Tm) == float(d[line[0]][0]) :    
                    d[line[0]][3] += 1
                    if Mismatch_3prime_start_pos > three_prime_region:
                        d[line[0]][1]    += 1
            d2[line[0]].append(float(Gap_adjusted_end_filling_Tm))
        for key, value in d2.iteritems():
            d[key][2]    =    misprime_percentile(d2[key], misprime_Tm_percentile_value)
            maximum_misprime_Tm    =    max_misprime_Tm(d2[key])
            if maximum_misprime_Tm != 'NA' and d[key][0] != "3prime_mismatch":
                    d[key][0]    =    maximum_misprime_Tm
    ### insert header names to the beginning of the list of list
    list_lines.insert( 0, ["Primer","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","e-value","bitscore","Match_seq", "Query_seq", "match_direction","3prime_overhang","5prime_overhang","PrimerTm","MatchTm","PrimerTm-MatchTm", "3&5_prime_match_extend",  "Actual_3prime_mismatches","First3Prime_mismatch" ,"Primer_NN_Tm","Local_alignment_Tm", "Gap_adjusted_end_filling_Tm"])
        
    ### convert the list of list to a data frame    
    headers = list_lines.pop(0) 
    df = pd.DataFrame(list_lines, columns=headers)


    ### write to output file
    output_headers = ["Primer","qstart","qend","sstart","send", "Query_seq","Match_seq","match_direction","3prime_overhang","5prime_overhang","3&5_prime_match_extend","PrimerTm","MatchTm","PrimerTm-MatchTm","First3Prime_mismatch","Primer_NN_Tm","Local_alignment_Tm", "Gap_adjusted_end_filling_Tm"]
    ### write selected columns of the df into a csv file
    df.to_csv(query_path + path + Time_stamp + "_" + 'PSE_out1_1.csv', columns = output_headers)

    ###
    import primer3
    ###
    primer_f_dict    =    dict()
    for key, value1 in d.iteritems():
        value = value1[0]
        three_prime_mismatch_alignments    = value1[3] - value1[1]
        max_Tm_alignments                = value1[3]                 
        three_primer_mismatch_alignments = str(three_prime_mismatch_alignments)+" -/- "+str(max_Tm_alignments)
        misprime_Tm_percentile = value1[2]
        key_Tm = NN_Tm(seq=key, compl_seq=complement(key), primer_conc=primer_conc, Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs, ion_corr=True)
        if value != "No match" and value != "3prime_mismatch":
            Tm_difference = str(float(key_Tm)-float(value))
        else:
            Tm_difference = "-"
        TmHairpin    =    hairpin_Tm(key, monovalent_cation_eq, primer_conc) 
        TmHomodimer    =    homodimer_Tm(key, monovalent_cation_eq, primer_conc)
        primer_f_dict[key.upper()]    =    [key_Tm, value, Tm_difference, three_primer_mismatch_alignments, misprime_Tm_percentile,TmHairpin,TmHomodimer]
    f.close() 

    ############################################################
    ### reverse primers
    ############################################################
    os.chdir(query_path)
    #### create a dictionary of all oligos from the input(query) sequence
    d = dict()
    d2 = dict()
    with open(query_path+fasta_input_file_R, 'r') as f3:
        for line3 in f3.readlines():
            if line3.startswith('>'):
                d[line3[1:-1]] =    ['No match', 0,'No match', 0 ]
                d2[line3[1:-1]] =[]      
    min_pimer_len    =    len(min(d,key=len))
    if min_pimer_len < min_pooled_primer_len: 
        min_pooled_primer_len   =   min_pimer_len 
    word_size = 7
    #### condition to raise an error if the user defined 3' region is > min_pimer_len
    if three_prime_region    > min_pimer_len:
            raise ValueError('3 prime region should not be greater than the minimum primer length! Given three_prime_region = ', three_prime_region, 'min_pimer_len = ', min_pimer_len)
    ### nearly identical match identification using blastn
    f0 = open(os.devnull, 'w')
    p3 = sp.Popen(["blastn","-task","blastn","-db","%szmv2all" %refDB_path,"-query","%s%s" %(query_path,fasta_input_file_R),"-evalue","%s" %mp_e_value,"-word_size","%s" %word_size,"-gapopen","%s" %mp_gapopen,"-gapextend","%s" %mp_gapextend,"-reward","%s" %mp_reward,"-penalty","%s" %mp_penalty,"-dust","no","-perc_identity","%s" %mp_perc_identity,"-max_target_seqs", "%s" %mp_max_target_seqs,"-max_hsps", "%s" %mp_max_hsps,"-outfmt","10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq qseq", "-num_threads","%s" %mp_num_threads],stdout=sp.PIPE,stderr=f0)
    output, error = p3.communicate()
    parameters_used.write(                
                    "Word_size_R_primers        =   "+str(word_size)+"\n"+
                    "Gapopen                    =   "+str(mp_gapopen)+"\n"+
                    "Gapextend                  =   "+str(mp_gapextend)+"\n"+
                    "Reward                     =   "+str(mp_reward)+"\n"+
                    "Penalty                    =   "+str(mp_penalty)+"\n"
                    "%_identity                 =   "+str(mp_perc_identity)+"\n"+
                    "Max_target_seqs            =   "+str(mp_max_target_seqs)+"\n"+
                    "Max_hsps                   =   "+str(mp_max_hsps)+"\n"+
                    "Number_threads             =   "+str(mp_num_threads)+"\n"+"\n"+
                    "## PSE parameters"+"\n"+
                    "3_prime_region             =   "+str(three_prime_region)+"\n"+
                    "Misprime_Tm_percentile     =   "+str(misprime_Tm_percentile_value )+"\n"+
                    "##########################################################"+"\n"
                        )

    list_lines = []        
    db = "%szmv2all" %refDB_path            
    for line in output.split('\n')[:-1]: 
        line = line.strip(' ').split(',')
        query_len    =    len(line[0])
        qstart        =    int(line[6])
        qend        =    int(line[7])
        sstart        =    float(line[8])
        send        =    float(line[9])
        pident        =    float(line[2])                           
        match_length    =    float(line[3])
        
        #### skip self alignments
        if pident == 100 and query_len == match_length:
            pass
        else:
            three_prime_mismatches    =    query_len-qend
            five_prime_mismatches    =    qstart-1
            fasta_id = line[1]
            Match_seq    =    line[12]
            Query_seq   = line[13]
            match_gaps  = int(Match_seq.count('-') )
            query_gaps  = int(Query_seq.count('-') )
            if send-sstart > 0:
                match_direction = "right"
                actual_match_start = int(sstart)
                actual_match_end = int(send)
                actual_match_end_3prime_extended = actual_match_end + int(three_prime_mismatches)
                actual_match_start_5prime_extended = actual_match_start - int(five_prime_mismatches)
                actual_match_start_5prime_extended = actual_match_start_5prime_extended - match_gaps
                actual_match_start_5prime_extended = actual_match_start_5prime_extended + query_gaps
                if int(three_prime_mismatches) == 0 and int(five_prime_mismatches) == 0 and match_gaps ==0 and query_gaps == 0:
                    match_extend_3_5 = Match_seq
                else:
                    coordinates_3_5primeextend = str(actual_match_start_5prime_extended)+":"+str(actual_match_end_3prime_extended)
                    match_extend_3_5 =    seq_extraction(fasta_id, coordinates_3_5primeextend) 
            else:
                match_direction = "left"
                actual_match_start = int(send)
                actual_match_end = int(sstart)
                actual_match_end_3prime_extended = actual_match_start - int(three_prime_mismatches)
                actual_match_start_5prime_extended = actual_match_end + int(five_prime_mismatches)
                actual_match_start_5prime_extended = actual_match_start_5prime_extended + match_gaps
                actual_match_start_5prime_extended = actual_match_start_5prime_extended - query_gaps
                if int(three_prime_mismatches) == 0 and int(five_prime_mismatches) == 0 and match_gaps == 0 and query_gaps == 0:
                    match_extend_3_5 = rev_complement(Match_seq) ### 040416 corrected
                else:
                    coordinates_3_5primeextend	= str(actual_match_end_3prime_extended)+":"+str(actual_match_start_5prime_extended) 
                    
                    
                    
                    match_extend_3_5 =    seq_extraction(fasta_id, coordinates_3_5primeextend) 
            Primer        =      line[0]    
            line.extend([match_direction])
            line.extend([three_prime_mismatches])
            line.extend([five_prime_mismatches])
            PrimerTm = Tm(Primer)
            line.append(PrimerTm)
            match_melt_temp        =     Tm(Match_seq)
            line.append(match_melt_temp)
            line.append(PrimerTm-match_melt_temp)
            qstart        =    int(line[6])-1
            qend        =    int(line[7])
            Match_seq_blast     = line[12]
            line.append(match_extend_3_5)    

            ### actual_3prime_mismatches
            if match_direction == "right":
                Actual_3prime_mismatches    = count_3prime_mismatches(Primer, match_extend_3_5)
                Mismatch_3prime_start_pos    = start_pos3prime_mismatches(Primer, match_extend_3_5)
            elif match_direction == "left":
                Actual_3prime_mismatches    = count_3prime_mismatches(Primer, rev_complement(match_extend_3_5))
                Mismatch_3prime_start_pos    = start_pos3prime_mismatches(Primer, rev_complement(match_extend_3_5))
            elif "-" in Match_seq_blast:
                Actual_3prime_mismatches     = "gap_in_match"
                Mismatch_3prime_start_pos    = "gap_in_match"
            else:
                Actual_3prime_mismatches     = "check_your_sequences!"
                Mismatch_3prime_start_pos    = "check_your_sequences!"

            line.append(Actual_3prime_mismatches)    
            line.append(Mismatch_3prime_start_pos)
            Primer_NN_Tm = NN_Tm(seq=Primer, compl_seq=complement(Primer), primer_conc=primer_conc, Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs, ion_corr=True)
            if match_direction == "right":
                complementary_match    =    complement(match_extend_3_5)
            else:
                complementary_match    =    match_extend_3_5[::-1]
            Gap_adjusted_end_filling_Tm = NN_Tm(seq=Primer, compl_seq=complementary_match, primer_conc=primer_conc, Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs, ion_corr=True)
            Local_alignment_Tm = NN_Tm(seq=Query_seq, compl_seq=complement(Match_seq), primer_conc=primer_conc, Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs, ion_corr=True)
            line.append(Primer_NN_Tm)
            line.append(Local_alignment_Tm)   
            line.append(Gap_adjusted_end_filling_Tm)            
            list_lines.append(line) 
            if Actual_3prime_mismatches > 0 and (d[line[0]][0] == "3prime_mismatch" or d[line[0]][0] == "No match"):    
                d[line[0]][0] =    "3prime_mismatch"
            if Actual_3prime_mismatches == 0: 
                if d[line[0]][0]=="No match"or d[line[0]][0]=="3prime_mismatch" or float(Gap_adjusted_end_filling_Tm) > float(d[line[0]][0]) :        
                    d[line[0]][0] =    Gap_adjusted_end_filling_Tm    
                if float(Gap_adjusted_end_filling_Tm) > float(d[line[0]][0]) :    
                    d[line[0]][3] = 1
                    if Mismatch_3prime_start_pos > three_prime_region:
                        d[line[0]][1]    += 1                    
                if  float(Gap_adjusted_end_filling_Tm) == float(d[line[0]][0]) :    
                    d[line[0]][3] += 1
                    if Mismatch_3prime_start_pos > three_prime_region:
                        d[line[0]][1]    += 1
            d2[line[0]].append(float(Gap_adjusted_end_filling_Tm))
    for key, value in d2.iteritems():
        d[key][2]    =    misprime_percentile(d2[key], misprime_Tm_percentile_value)
        maximum_misprime_Tm    =    max_misprime_Tm(d2[key])
        if maximum_misprime_Tm != 'NA' and d[key][0] != "3prime_mismatch":
                d[key][0]    =    maximum_misprime_Tm
    list_lines.insert( 0, ["Primer","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","e-value","bitscore","Match_seq","Query_seq","match_direction","3prime_overhang","5prime_overhang","PrimerTm","MatchTm","PrimerTm-MatchTm", "3&5_prime_match_extend",  "Actual_3prime_mismatches","First3Prime_mismatch" ,"Primer_NN_Tm","Local_alignment_Tm", "Gap_adjusted_end_filling_Tm"])

    ### convert the list of list to a data frame    
    headers = list_lines.pop(0)
    df = pd.DataFrame(list_lines, columns=headers)

    ### write to output file
    output_headers = ["Primer","qstart","qend","sstart","send","Query_seq","Match_seq","match_direction","3prime_overhang","5prime_overhang","3&5_prime_match_extend","PrimerTm","MatchTm","PrimerTm-MatchTm","First3Prime_mismatch","Primer_NN_Tm","Local_alignment_Tm", "Gap_adjusted_end_filling_Tm"]
    df.to_csv(query_path + path + Time_stamp + "_" + 'PSE_out2_1.csv', columns = output_headers)

    ###
    import primer3                                                
    primer_r_dict    =    dict()
    for key, value1 in d.iteritems():
        value = value1[0]
        three_prime_mismatch_alignments = value1[3] - value1[1]
        max_Tm_alignments               = value1[3]    
        three_primer_mismatch_alignments = str(three_prime_mismatch_alignments)+" -/- "+str(max_Tm_alignments)
        misprime_Tm_percentile = value1[2]
        key_Tm = NN_Tm(seq=key, compl_seq=complement(key), primer_conc=primer_conc, Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs, ion_corr=True)
        if value != "No match" and value != "3prime_mismatch":
            Tm_difference = str(float(key_Tm)-float(value))
        else:
            Tm_difference = "-"
        TmHairpin    =    hairpin_Tm(key, monovalent_cation_eq, primer_conc)                ## hairpin_Tm        091715
        TmHomodimer    =    homodimer_Tm(key, monovalent_cation_eq, primer_conc)            ## homodimer_Tm     091715
        primer_r_dict[key.upper()]    =    [key_Tm, value, Tm_difference, three_primer_mismatch_alignments, misprime_Tm_percentile,TmHairpin,TmHomodimer]        ## hairpin, homodomer Tm 091715
    f.close() 

    ### Combine two dictionaries (primer_f_dict & primer_r_dict)
    pooled_primer_f_r_dict    =    dict()
    pooled_primer_f_r_dict.update(primer_f_dict)
    pooled_primer_f_r_dict.update(primer_r_dict)

    ############################################################
    ### coordinate identification
    ############################################################

    ### Get exact matches for the primers
    fasta_input_file_FR    =    path + Time_stamp + "_" +'UOD_out4_1.fasta' 
    ### Use the smallest value in the primer size range as the word size
    word_size            =    min_pooled_primer_len
    os.chdir(query_path)

    f0 = open(os.devnull, 'w')
    sp.call(["makeblastdb","-in","%s" %locus, "-dbtype","nucl","-parse_seqids","-out","%sloci_db" %query_path], stdout=f0,stderr=f0)

    ### Exact match to loci to get coords
    p5 = sp.Popen(["blastn","-db","%sloci_db" %query_path,"-query","%s" %fasta_input_file_FR,"-evalue","0.1","-word_size","%s" %word_size,"-gapopen","0","-gapextend","2","-reward","1","-penalty","-3","-dust","no","-perc_identity","100","-max_target_seqs", "13","-outfmt","10 qseqid length sstart send", "-num_threads","%s" %mp_num_threads],stdout=sp.PIPE,stderr=f0)
    exact_match_output_pooled_primers, error = p5.communicate()

    ### final pooled oligo coordinate output file
    f5 = open(query_path + path + Time_stamp + "_" + 'PSE_out3_1.csv','w')
    f5.write("Primer"+','+"Loci_start"+','+"Loci_stop"+','+"Genome_start"+','+"Genome_stop"+','+"Strand"+','+"Primer_Tm"+','+"Max_misprime_Tm"+','+"Tm_difference"+','+'Misprime_Tm_'+str(misprime_Tm_percentile_value)+'th_percentile'+','+"Primer_GC"+','+"Continuous_GC"+','+"3'_region_mismatches"+','+"Hairpin_Tm"+','+"Homodimer_Tm"+'\n')                            ### changed header name to Max_misprime_Tm

    for exact_match_pooled_output_line in exact_match_output_pooled_primers.split('\n')[:-1]:
        exact_match_pooled_output_line = exact_match_pooled_output_line.strip(' ').split(',')
        Primer_pooled        =    exact_match_pooled_output_line[0]
        primer_gc            =    gc_content(Primer_pooled)
        primer_continuous_gc=    continuous_gc(Primer_pooled)
        query_len_pooled    =    len(Primer_pooled)
        match_len_pooled    =    int(exact_match_pooled_output_line[1])
        sstart                 =    int(exact_match_pooled_output_line[2])
        sstart                =    (actual_locus_start_pos + sstart) - 1
        send                =    int(exact_match_pooled_output_line[3])
        send                =    (actual_locus_start_pos + send) - 1
        if query_len_pooled == match_len_pooled:
            Primer_Tm_p            =    pooled_primer_f_r_dict[Primer_pooled][0]
            Max_misprime_Tm_p    =    pooled_primer_f_r_dict[Primer_pooled][1]
            Tm_difference_p        =    pooled_primer_f_r_dict[Primer_pooled][2]
            three_prime_region_mismatches        =    pooled_primer_f_r_dict[Primer_pooled][3]
            misprime_Tm_percentile = pooled_primer_f_r_dict[Primer_pooled][4]
            Hairpin_Tm            = pooled_primer_f_r_dict[Primer_pooled][5]
            Homodimer_Tm        = pooled_primer_f_r_dict[Primer_pooled][6]
            if send-sstart > 0:
                Primer_pooled_start    =    sstart
                Primer_pooled_stop    =    send
                strand                =    "+"
                f5.write(str(Primer_pooled) +','+str(Primer_pooled_start-actual_locus_start_pos)+','+str(Primer_pooled_stop-actual_locus_start_pos)+','+str(Primer_pooled_start)+','+str(Primer_pooled_stop)+','+str(strand)+','+str(Primer_Tm_p)+','+str(Max_misprime_Tm_p)+','+str(Tm_difference_p)+','+str(misprime_Tm_percentile)+','+str(primer_gc)+','+str(primer_continuous_gc)+','+str(three_prime_region_mismatches)+','+str(Hairpin_Tm)+','+str(Homodimer_Tm)+'\n')
            else:
                Primer_pooled_start    =    send
                Primer_pooled_stop    =    sstart
                strand                =    "-"    
                f5.write(str(Primer_pooled) +','+str(Primer_pooled_start-actual_locus_start_pos)+','+str(Primer_pooled_stop-actual_locus_start_pos)+','+str(Primer_pooled_start)+','+str(Primer_pooled_stop)+','+str(strand)+','+str(Primer_Tm_p)+','+str(Max_misprime_Tm_p)+','+str(Tm_difference_p)+','+str(misprime_Tm_percentile)+','+str(primer_gc)+','+str(primer_continuous_gc)+','+str(three_prime_region_mismatches)+','+str(Hairpin_Tm)+','+str(Homodimer_Tm)+'\n')
    f5.close()            

    ############################################################
    #Time to run the code: end timer
    ############################################################
    t1 = time.time()
    total = t1-t0
    total = ("{0:.2f}".format(round(total,2)))
    parameters_used.write(  
                                "### PSE run duration : " + str(total) + " seconds"+'\n'
                                "##########################################################"+"\n"+
                                "\n"+"\n"
                    )
        
    parameters_used.close()

### Version log (SemVer format)
### 0.1.0: Initial development release
### 0.1.1: Use parameters from parameters.py file; time stamp based naming of the directory and output files
### 0.1.2: TA_timestamp_PSE_out3_1.csv has new format for "3'_region_mismatches" column: 2 -/- 3





