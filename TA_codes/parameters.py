##########################################################
### User input parameters for Thermo-Align primer design
##########################################################

### path to required files
# refDB_path          =   "/mnt/data27/ffrancis/ZmDB_RefGenv3_renamed/"
# hapmap_path         =   "/mnt/data27/ffrancis/HapMap3/all_lines/processed_hapmap_files/"
# multiplx_path       =   "/home/ffrancis/softwares/multiplx/"
refDB_path          =   "../sample_genome/"
hapmap_path         =   "../hapmap_vcf_sample/processed_hapmap_files/"
multiplx_path       =   "../"

## target region selection 
chr_no              =   3           # fast test
start_pos           =   100    # fast test
stop_pos            =   1000    # fast test
no_chrs_in_genome   =   5           # specify the number of chromosomes in the genome used
hapmap_mask_condition = 1       #(1:Yes; 0:No)   
flanking_primers    =   0       #(1:Yes; 0:No)
flanking_size       =   2000
                              
### PCR conditions
primer_size_range   =   "24-25"
GC_range            =   "40-55"
Tm_range            =   "64-70" 
# Tm_range            =   "68-70" 

### NN Tm parameters NEB Wisser lab
primer_conc         =   100    #nM  
Na                  =   0      #mM
K                   =   50     #mM
Tris                =   10     #mM
Mg                  =   1.5    #mM  
dNTPs               =   0.2    #mM  

## primer features
self_Tmdiff         =   20
filter_AT_3prime    =   0       #(1:Yes; 0:No)
di_si_repeats_threshold =   4
filter_di_si_repeats=   1       #(1:Yes; 0:No)
filter_GC_clamp     =   0       #(1:Yes; 0:No)   
            
## blastn parameters for exact-match(em) search
em_e_value          =   30000
em_gapopen          =   2
em_gapextend        =   2
em_reward           =   1
em_penalty          =   -3
em_perc_identity    =   100
em_max_target_seqs  =   2
em_max_hsps         =   2
em_num_threads      =   5

## blastn parameters for mis-prime(mp) search
mp_e_value          =   30000
mp_gapopen          =   2                                                                  
mp_gapextend        =   2
mp_reward           =   1
mp_penalty          =   -1
mp_num_threads      =   5
mp_perc_identity    =   70
mp_max_target_seqs  =   13
mp_max_hsps         =   20

### primer pair picking parameters                               
amplicon_size_min   =   100
amplicon_size_max   =   800
amplicon_gap_filter =   1       #(1:Yes; 0:No)
pair_Tm_diff        =   10        
pair_GC_diff        =   20        
pair_misprimeTm_diff=   10
min_misprime_dg     =   -2000

### multiplex primer picking parameters
multiplex_primers   =   1       #(1:Yes; 0:No)
multiplex_Tm_min    =   64
multiplex_Tm_max    =   70

multiplx_calcscores =   "12345"     #(1->Primer-primer both 3-p ends; 2->Primer-primer 3-p end with any region; 3->Primer-primer any regions; 4->Primer-product 3-p end with any region; 5->Primer-product any regions)
grouping_stringency =   "normal"    #(low/normal/high)
maxgroups           =   1000       
maxitemsingroup     =   1000

##########################################################
### time_stamp
##########################################################

