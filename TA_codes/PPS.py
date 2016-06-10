### PPS: Primer Pair Selection (PPS); a part of Thermo-Align tool for the design of template specific hybridization and priming oligonucleotides
### Version-1.05: 06/10/2016
### Authors: Felix Francis (felixfrancier@gmail.com); Randall J. Wisser (rjw@udel.edu) 


############################################################
# Required dependencies
############################################################

# sudo apt-get install python-pip
# sudo pip install numpy
# sudo pip install pandas
# pip install networkx
# pip install primer3-py


############################################################
# Time to run the code: start timer
############################################################
import time
t0 = time.time()


############################################################
# Import
############################################################
import numpy as np
import pandas as pd
import sys
import networkx as nx
import csv
import primer3
from Santalucia_NN_Tm import mM_monovalent
import sequence_processing_functions as spf
from Bio import SeqIO
from time_stamp import Time_stamp
from parameters import *                                                          
import subprocess as sp    
import os

############################################################
# Functions
############################################################

### function to check hetero dimer dg
def heterodimer_dg(seq1, seq2, mv_cation=0,primer_conc=0):
        dg =  (primer3.calcHeterodimer(seq1, seq2,mv_conc=mv_cation, dv_conc=0, dntp_conc=0, dna_conc=primer_conc, temp_c=60, max_loop=30)).tm
        return float(("{0:.2f}".format(round(dg,2))))

### monovalent_cation_eq
monovalent_cation_eq    =    mM_monovalent(Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs)

### function to calculate max_cover_sets from all dijkstra's shortest paths
def dijkstra_max_cover_amplicon_sets(dijkstra_shortest_paths):
    max_cover_sets  = [] 
    for shortpath in dijkstra_shortest_paths:
        shortpath_min, shortpath_max = int(shortpath[0][0]), int(shortpath[-1][-1]) 
        if not len(max_cover_sets):
            max_cover_sets.append(shortpath)
        else:
            flagGoodToInsert = 1
            for selSet in max_cover_sets:
                selSetStart, selSetEnd = int(selSet[0][0]), int(selSet[-1][-1])
                if shortpath_min>selSetStart and shortpath_max<selSetEnd:                                            # ignore amplicon sets if the corresponding range has already been covered
                    flagGoodToInsert = 0
                    break
                if shortpath_min==selSetStart and shortpath_max==selSetEnd and (len(shortpath) >= len(selSet)):      # ignore amplicon sets if the corresponding range has already been covered; same start-stop,more # primers
                    flagGoodToInsert = 0
                    break
                if shortpath_min==selSetStart and shortpath_max<selSetEnd:                                           # ignore amplicon sets if the corresponding range has already been covered; same start
                    flagGoodToInsert = 0
                    break
                if shortpath_min>selSetStart and shortpath_max==selSetEnd:                                           # ignore amplicon sets if the corresponding range has already been covered; same stop
                    flagGoodToInsert = 0
                    break
            if flagGoodToInsert:
                pairsToBeRemoved = []
                for selSet in max_cover_sets:
                    selSetStart, selSetEnd = int(selSet[0][0]), int(selSet[-1][-1])
                    if shortpath_min <= selSetStart and shortpath_max >= selSetEnd:
                        pairsToBeRemoved.append(selSet)
                if len(pairsToBeRemoved):
                    for removePair in pairsToBeRemoved:
                        max_cover_sets.remove(removePair)
                max_cover_sets.append(shortpath)
    return sorted(max_cover_sets) 

### function to pick primer pairs
def pick_primer_pairs(input_file):
    G=nx.Graph()
    primer_list =   []
    d_coords_plus   = {}
    d_coords_minus  = {}
    counter         = 0
    f = open(input_file)
    csv_f = csv.DictReader(f, delimiter=',')
    for row in csv_f:
        primer, p_start_pos, p_stop_pos, strand, Tm, max_misprimeTm, GC  = row['Primer'], int(row['Genome_start']), int(row['Genome_stop']), row['Strand'], \
        float(row['Primer_Tm']), (row['Max_misprime_Tm']), float(row['Primer_GC'])
        if max_misprimeTm == "3prime_mismatch":
            max_misprimeTm=0
        else:
            max_misprimeTm=float(max_misprimeTm)
            
        if Tm - max_misprimeTm >= pair_misprimeTm_diff:
            primer_list.append(primer)
            G.add_node(primer, p_start_pos=p_start_pos, p_stop_pos=p_stop_pos, strand=strand, Tm=Tm, max_misprimeTm=max_misprimeTm, GC=GC)
    unique_primers          =   {}
    primer_pair_coords      =   []
    amplicon_len_list       =   []
    amplicon_coords_list    =   []
    mplex_ampl_coords_list  =   []
    ### pick primer pairs that meet the user defined conditions
    ### write bed files (for picked oligos)
    f8 = open(path + Time_stamp + "_" +'bed_separate_tracks_selected_oligos.bed', 'w')
    f8.write('browser position chr '+str(chr_no)+':'+ str(start_pos)+'-'+str(stop_pos)+'\n')
    f8.write('track name="Primers" description="Primers on separate tracks" visibility=2 colorByStrand="255,0,0 0,0,255"' + '\n')
    Gfor_subsetting=nx.Graph()
    primer_pair_counter = 0
    for (nodeId, data_f) in G.nodes(data=True):
        if data_f['strand']   == "+":
            forward_start   =   data_f['p_start_pos']
            forward_Tm      =   int(data_f['Tm'])
            forward_GC      =   int(data_f['GC'])
            forward_MaxMispTm      =   int(data_f['max_misprimeTm'])
            query_nodeId    =   nodeId
            for (nodeId, data_r) in G.nodes(data=True):
                if data_r['strand']   == "-":
                    amplicon_len = (data_r['p_stop_pos']- forward_start)+1
                    reverse_Tm      =   int(data_r['Tm'])
                    reverse_GC      =   int(data_r['GC'])
                    reverse_MaxMispTm      =   int(data_r['max_misprimeTm'])
                    ### check if primer pair meets user defined amplicon length conditon
                    if amplicon_len >= amplicon_size_min and amplicon_len <= amplicon_size_max:
                        ### check if primer pair meets user defined pair_misprimeTm_diff condition
                        if min(forward_Tm, reverse_Tm) - max(forward_MaxMispTm, reverse_MaxMispTm) >= pair_misprimeTm_diff:
                            ### check if primer pair meets user defined pair_Tm_diff and pair_GC_diff conditions
                            if abs(forward_Tm-reverse_Tm) <= pair_Tm_diff:
                                ### check if primer pair meets user defined between primer interaction_dg conditions
                                interaction_dg  =    heterodimer_dg(query_nodeId, nodeId, mv_cation=monovalent_cation_eq,primer_conc=primer_conc)
                                if interaction_dg    >=  min_misprime_dg:
                                    amplicon_len_list.append(amplicon_len)
                                    ### get the expected amplicon sequence
                                    amplicon_seq    =    spf.seq_extraction_loci(locus, start_pos, data_f['p_start_pos'], data_r['p_stop_pos'])                               
                                                         
                                    ### eliminate amplicons with gaps
                                    if amplicon_gap_filter ==   1:
                                        if 'N'*100 in amplicon_seq:
                                            break
                                    
                                    amplicon_gc    =    spf.gc_content(amplicon_seq.upper()) 
                                    
                                    ### check if amplicons have gaps and indels
                                    if 'N'*100 in amplicon_seq:
                                        gaps    =   'Yes'
                                    else:
                                        gaps    =   'No'
                                    if 'n' in amplicon_seq:
                                        indel    =   'Yes'
                                    else:
                                        indel    =   'No'
                                    f_primer_length =   len(query_nodeId)    
                                    r_primer_length =   len(nodeId)
                                    primer_name_f     =   "TA_" + str(chr_no) + "_" + str(data_f['p_start_pos']) + "_" + str(f_primer_length) + "_F"
                                    primer_name_r     =   "TA_" + str(chr_no) + "_" + str(data_r['p_stop_pos']) + "_" + str(r_primer_length) + "_R"          ### now naming of r primers uses 5' end coords
                                    primer_pair_coords.append((data_f['p_start_pos'], data_r['p_stop_pos']))
                                    if query_nodeId not in unique_primers:
                                        f8.write('chr'+ str(chr_no)+'\t'+ str(data_f['p_start_pos'])+'\t'+ str(data_f['p_stop_pos'])+'\t'+ str(primer_name_f)+'\t'+ str(0) +'\t'+ "+"+'\n') 
                                    elif nodeId not in unique_primers:
                                        f8.write('chr'+ str(chr_no)+'\t'+ str(data_r['p_start_pos'])+'\t'+ str(data_r['p_stop_pos'])+'\t'+ str(primer_name_r)+'\t'+ str(0) +'\t'+ "-"+'\n') 
                                    unique_primers[query_nodeId]    =   [str(chr_no), str(data_f['p_start_pos']), str(data_f['p_stop_pos']), str(data_f['Tm']), \
                                    str(data_f['max_misprimeTm']), str(data_f['GC']), str(data_f['strand'])]
                                    unique_primers[nodeId]          =   [str(chr_no), str(data_r['p_start_pos']), str(data_r['p_stop_pos']), str(data_r['Tm']), \
                                    str(data_r['max_misprimeTm']), str(data_r['GC']), str(data_r['strand'])]
                                    amplicon_coords    = (data_f['p_start_pos'], data_r['p_stop_pos'])                                                     # nodes for creation of new network for subsetting connected components
                                    amplicon_coords_list.append(amplicon_coords) 
                                    ### info for multiplex primer picking
                                    ### select primers within the given Tm range if select multiplex primer select option is given in the parameters file
                                    if multiplex_primers   == 1 and data_f['Tm'] >= multiplex_Tm_min and data_r['Tm'] >= multiplex_Tm_min and data_f['Tm'] <= multiplex_Tm_max and data_r['Tm'] <= multiplex_Tm_max:
                                            f_primer_info   =   {'primer_name':primer_name_f, 'primer_sequence':query_nodeId, 'strand':data_f['strand'], 'p_start_pos':data_f['p_start_pos'], 'p_stop_pos':data_f['p_stop_pos'], 'Tm':data_f['Tm'], 'max_misprimeTm':data_f['max_misprimeTm'], 'GC':data_f['GC']}
                                            r_primer_info   =   {'primer_name':primer_name_r, 'primer_sequence':nodeId, 'strand':data_r['strand'], 'p_start_pos':data_r['p_start_pos'], 'p_stop_pos':data_r['p_stop_pos'], 'Tm':data_r['Tm'], 'max_misprimeTm':data_r['max_misprimeTm'], 'GC':data_r['GC']}
                                            amplicon_info   =   {'interaction_dg':interaction_dg, 'amplicon_len':amplicon_len, 'amplicon_gc':amplicon_gc,'gaps':gaps, 'indel':indel, 'amplicon_seq':amplicon_seq}
                                            Gfor_subsetting.add_node(amplicon_coords, f_primer_info = f_primer_info, r_primer_info = r_primer_info, amplicon_info = amplicon_info)
                                            mplex_ampl_coords       = (data_f['p_start_pos'], data_r['p_stop_pos'])     
                                            mplex_ampl_coords_list.append(mplex_ampl_coords)
                                    primer_pair_counter += 1
                                    output_primer_pairs.write(str(primer_pair_counter) + '\t' + str(primer_name_f)+'\t'+str(query_nodeId)  +'\t'+ str(data_f['strand']) +'\t'+ str(data_f['p_start_pos']) +'\t'+ str(data_f['p_stop_pos']) +'\t'+ str(data_f['Tm'])+'\t'+ str(data_f['max_misprimeTm'])+'\t'+ str(data_f['GC']) +'\t'+ str(interaction_dg)+'\t'+ str(amplicon_len) +'\t'+ str(amplicon_gc)+'\t'+ str(gaps)+'\t'+ str(indel)+'\t'+ str(amplicon_seq)+'\n')
                                    output_primer_pairs.write(str(primer_pair_counter) + '\t' + str(primer_name_r )+'\t'+str(nodeId)+'\t'+ str(data_r['strand'])+'\t'+ str(data_r['p_stop_pos']) +'\t'+ str(data_r['p_start_pos']) +'\t'+ str(data_r['Tm'])+'\t'+ str(data_r['max_misprimeTm'])+'\t'+ str(data_r['GC']) +'\t'+ '-'+'\t'+ '-'  +'\t'+ '-' +'\t'+ '-' +'\t'+ '-' +'\t'+ '-'+'\n')
                                    order_primer_pairs.write(str(primer_name_f )+'\t'+ str(query_nodeId) +'\t'+ str(data_f['Tm'])+'\n')
                                    order_primer_pairs.write(str(primer_name_r)+'\t'+ str(nodeId) +'\t'+ str(data_r['Tm']) +'\n') 
    f8.close()
    output_primer_pairs.close()                               
    order_primer_pairs.close()                               
    no_unique_primers_picked    =   len(unique_primers)                               
    return {'mplex_ampl_coords_list': mplex_ampl_coords_list, 'Gfor_subsetting': Gfor_subsetting, 'amplicon_len_list': amplicon_len_list, 'amplicon_coords_list':amplicon_coords_list, 'primer_pair_coords': primer_pair_coords,\
    'no_unique_primers_picked': no_unique_primers_picked}

### function to retrieve amplicons from each multiplex set:
def multiplx_input_create(list_amplicons, output_file):
    for amplicon_coord_info in list_amplicons:
        f_primer_seq    = f_primer_info_all[amplicon_coord_info]['primer_sequence']   
        r_primer_seq    = r_primer_info_all[amplicon_coord_info]['primer_sequence'] 
        f_primer_name   = f_primer_info_all[amplicon_coord_info]['primer_name']   
        r_primer_name   = r_primer_info_all[amplicon_coord_info]['primer_name']
        amplicon_name   = f_primer_name + "-" + r_primer_name   
        amplicon_seq    = amplicon_info_all[amplicon_coord_info]['amplicon_seq']
        output_file.write(amplicon_name + '\t' + f_primer_seq + '\t' + r_primer_seq + '\t' + amplicon_seq +'\n')
    output_file.close()

### function to retrieve multiplex compatible groups :
def multiplex_group_output(input_file, set_no, primer_pair_counter):
    with open(input_file) as f:
        content = f.readlines()
        for line in content:
            if line.startswith('Group'):
                line = line.split()
                group = line[1]
                for amplicon in line[2:]:
                    f_primer, r_primer = amplicon.split('-')[0], amplicon.split('-')[1]
                    forward_primer_start = int(f_primer.split('_')[2])   
                    for key, value in r_primer_info_all.iteritems():
                        if value['primer_name'] == r_primer and int(key[0]) == forward_primer_start :
                            r_primer_info = value
                            selected_coords = key
                    primer_pair_counter += 1
                    f_primer_name, f_primer_seq, fStrand, fStart_pos, fStop_pos, fTm, fMax_misprimeTm, fGC, Primer_dimer_dG, Amplicon_size, Amplicon_gc, Gaps, Polymorphisms, Amplicon_seq \
                                    = f_primer_info_all[selected_coords]['primer_name'],f_primer_info_all[selected_coords]['primer_sequence'],f_primer_info_all[selected_coords]['strand'],\
                                    f_primer_info_all[selected_coords]['p_start_pos'],f_primer_info_all[selected_coords]['p_stop_pos'],f_primer_info_all[selected_coords]['Tm'],\
                                    f_primer_info_all[selected_coords]['max_misprimeTm'],f_primer_info_all[selected_coords]['GC'],amplicon_info_all[selected_coords]['interaction_dg'],\
                                    amplicon_info_all[selected_coords]['amplicon_len'],amplicon_info_all[selected_coords]['amplicon_gc'],amplicon_info_all[selected_coords]['gaps'],amplicon_info_all[selected_coords]['indel'],\
                                    amplicon_info_all[selected_coords]['amplicon_seq']
                    
                    multiplx_pooled_out.write(str(set_no) +'\t'+ str(group) +'\t'+ str(primer_pair_counter) +'\t'+ f_primer_name+'\t'+ f_primer_seq +'\t'+ fStrand +'\t'+ str(fStart_pos) +'\t'+ str(fStop_pos) +'\t'+ str(fTm) +'\t'+ str(fMax_misprimeTm) +'\t'+ str(fGC) +'\t'+ str(Primer_dimer_dG) +'\t'+ str(Amplicon_size) +'\t'+ str(Amplicon_gc) +'\t'+ str(Gaps) +'\t'+ Polymorphisms +'\t'+ Amplicon_seq +'\n')
                    
                    r_primer_name, r_primer_seq, rStrand, rStart_pos, rStop_pos, rTm, rMax_misprimeTm, rGC \
                                    = r_primer_info['primer_name'],r_primer_info['primer_sequence'],r_primer_info['strand'],\
                                    r_primer_info['p_start_pos'],r_primer_info['p_stop_pos'],r_primer_info['Tm'],\
                                    r_primer_info['max_misprimeTm'],r_primer_info['GC']
                    multiplx_pooled_out.write(str(set_no) +'\t'+ str(group) +'\t'+ str(primer_pair_counter) +'\t'+ r_primer_name+'\t'+ r_primer_seq+'\t'+ rStrand +'\t'+ str(rStop_pos) +'\t'+ str(rStart_pos) +'\t'+ str(rTm) +'\t'+ str(rMax_misprimeTm) +'\t'+ str(rGC) +'\t'+ '-' +'\t'+ '-' +'\t'+ '-' +'\t'+ '-' +'\t'+ '-' +'\t'+ '-' +'\n')         
    f.close()
    return int(primer_pair_counter)
    
monovalent_salts = Na + K + (Tris/2)

############################################################
# Code
############################################################

if __name__ == '__main__':

    ### Flanking primer conditon
    if flanking_primers == 1:
        start_pos   = start_pos - flanking_size
        stop_pos    = stop_pos + flanking_size

    ### input_files
    path        =   Time_stamp + "/"
    input_file  = path + Time_stamp + "_" + 'HSE_out3_1.csv'

    ### Locus file name
    locus_len   =   (stop_pos - start_pos)+1
    locus       =    path + Time_stamp + "_" + str(chr_no) + "_" + str(start_pos) + "_" + str((stop_pos - start_pos)+1) +"_PolyMasked.fasta"

    ############################################################
    # Output files
    ############################################################    
    ### primer pair output files
    output_primer_pairs = open(path + Time_stamp + "_" +"primer_pairs_info.txt", "w")
    output_primer_pairs.write('Primer_pair#'+'\t'+'Primer_name'+'\t'+'Primer_seq'+'\t'+ 'Strand'+'\t'+ '5prime_pos' +'\t'+ '3prime_pos' +'\t'+ 'Tm' +'\t'+ \
    'Max_misprimeTm'+'\t'+ 'GC' +'\t'+ 'Primer_dimer_dG' +'\t'+ 'Amplicon_size' +'\t'+ 'Amplicon_GC' +'\t'+ 'Gaps'+'\t'+ 'Polymorphisms'+'\t'+ 'Amplicon_seq'+'\n')
    order_primer_pairs = open(path + Time_stamp + "_" +"primer_pairs_order.txt", "w")
    order_primer_pairs .write('Primer_name'+'\t'+'Primer_seq'+'\t'+ 'Tm' +'\n')

    ### multiplex output files
    if multiplex_primers   == 1:     
            ### multiplx input files
            multiplx_input_set1_1 = open(path + Time_stamp + "_" +"multiplx_input_set1_1.txt", "w")
            multiplx_input_set2_1 = open(path + Time_stamp + "_" +"multiplx_input_set2_1.txt", "w")
            multiplx_pooled_out = open(path + Time_stamp + "_" +"multiplx_pooled_output.txt", "w")
            multiplx_pooled_out.write('Multiplex_set#'+'\t'+'Multiplex_group#'+'\t'+'Primer_pair#'+'\t'+'Primer_name'+'\t'+'Primer_sequence'+'\t'+ 'Strand'+\
            '\t'+ '5prime_pos' +'\t'+ '3prime_pos' +'\t'+ 'Tm' +'\t'+ 'Max_misprimeTm'+'\t'+ 'GC' +'\t'+ 'Primer_dimer_dG' +'\t'+ 'Amplicon_size' +'\t'+ 'Amplicon_GC' +'\t'+ 'Gaps'+\
            '\t'+ 'Polymorphisms'+'\t'+ 'Amplicon_seq'+'\n')
    ### Parameter out put file
    parameters_used     = open(path + Time_stamp +'_run_summary.txt', 'a')

    ### print parameters used
    parameters_used.write(  "##########################################################"+"\n"+
                    "### Summary of PPS parameters"+"\n"+
                    "### Start date: " + str(time.strftime("%m/%d/%Y"))+";"+ "\t"+ "Start time: "+ str(time.strftime("%H:%M:%S"))+"\n"+
                    "##########################################################"+"\n"+
                    "Max_primer_pair_Tm_diff    =   "+str(pair_Tm_diff)+"\n"+
                    "Min_primer_misprimeTm_diff =   "+str(pair_misprimeTm_diff)+"\n"+
                    "Min_Amplicon_size          =   "+str(amplicon_size_min)+"\n"+
                    "Max_Amplicon_size          =   "+str(amplicon_size_max)+"\n"+
                    "amplicon_gap_filter        =   "+str(amplicon_gap_filter)+ "(1:Yes; 0:No)" +"\n"+
                    "min_misprime_dg            =   "+str(min_misprime_dg)+"\n"+
                    "multiplex_primers          =   "+str(multiplex_primers)+"(1:Yes; 0:No)" +"\n"+
                    "multiplex_Tm_min           =   "+str(multiplex_Tm_min)+"\n"+
                    "multiplex_Tm_max           =   "+str(multiplex_Tm_max)+"\n"+
                    "Multiplx_calcscores        =   "+str(multiplx_calcscores)+"\n"+
                    "Grouping_stringency        =   "+str(grouping_stringency)+"\n"+
                    "##########################################################"+"\n"
                )
                
    primer_picking_output = pick_primer_pairs(input_file)
    mplex_ampl_coords_list, Gfor_subsetting, amplicon_len_list, amplicon_coords_list, primer_pair_coords, no_unique_primers_picked = primer_picking_output['mplex_ampl_coords_list'], primer_picking_output['Gfor_subsetting'], primer_picking_output['amplicon_len_list'], primer_picking_output['amplicon_coords_list'], primer_picking_output['primer_pair_coords'], primer_picking_output['no_unique_primers_picked']
    
    ############################################################
    # Graph analysis to identify minimum tiling path
    ############################################################
    sorted_amplicon_coords_list =   sorted(list(set(mplex_ampl_coords_list)))
    ### build network of all possible amplicons (connections are drawn based on overlap)
    for i in xrange(len(sorted_amplicon_coords_list) -1):
        query_node  = sorted_amplicon_coords_list[i]
        query_amplicon_start, query_amplicon_stop =  int(query_node[0]), int(query_node[1])
        for j in xrange(len(sorted_amplicon_coords_list)-(i+1)):
            target_node  = sorted_amplicon_coords_list[i+j+1]
            target_amplicon_start, target_amplicon_stop =  int(target_node[0]), int(target_node[1])
            if (query_amplicon_stop >= target_amplicon_start) and (target_amplicon_stop > query_amplicon_stop) and (target_amplicon_start > query_amplicon_start):
                Gfor_subsetting.add_edge(query_node,target_node)

    ### pick all sub networks
    graphs = list(nx.connected_component_subgraphs(Gfor_subsetting))
    bestshortest_paths_list  =   []
    all_short_paths = [] 
    f_primer_info_all   =  nx.get_node_attributes(Gfor_subsetting,'f_primer_info')                                              # to get info. about all f primers from the orginal network 
    r_primer_info_all   =  nx.get_node_attributes(Gfor_subsetting,'r_primer_info')                                              # to get info. about all r primers from the orginal network 
    amplicon_info_all   =  nx.get_node_attributes(Gfor_subsetting,'amplicon_info')                                              # to get info. about all amplicons from the orginal network 
    all_nodes_connectedcomponents = [] 
    graph_no    =   0
    subnetwork_stats = {}
    counter = 1
    subnetwork_stats_included = []
    for graph in graphs:
        amplicon_coords = sorted(nx.nodes(graph))
        amplicon_coords_sorted_by_second = sorted(amplicon_coords, key=lambda tup: tup[1])
        all_nodes_connectedcomponents.append(amplicon_coords)
        start_node  =  amplicon_coords[0] 
        stop_node   = amplicon_coords_sorted_by_second[-1]  

        ### record the start stop pos for each subnetwork
        if len(amplicon_coords) > 1:
            graph_no    += 1
            subnetwork_stats['subnetwork# '+ str(graph_no)]  = ['start = ' + str(start_node[0]), 'stop = ' + str(stop_node[1]), 'length = '+str(int(stop_node[1]) - int(start_node[0]) + 1)+ ' bp', '# primer pairs = ' + str(len(amplicon_coords))]
            subnetwork_stats_included.append((start_node[0], stop_node[1]))
        
        
        max_coverage_len = 0 
        DG=nx.DiGraph()
        DG.add_nodes_from(amplicon_coords)
        
        # first iteration to get max edge weight (penalized cumulative coverage)
        for i in xrange(len(amplicon_coords)-1):
                start_pos_n0,stop_pos_n0 = amplicon_coords[i][0], amplicon_coords[i][1]
                for j in xrange(len(amplicon_coords[i+1:])):
                    start_pos_n1,stop_pos_n1 = amplicon_coords[i+j+1][0], amplicon_coords[i+j+1][1]
                    if (stop_pos_n0 >= start_pos_n1) and (stop_pos_n1 > stop_pos_n0) and (start_pos_n1 > start_pos_n0):
                        coverage    = ((stop_pos_n1 - start_pos_n0)+1) - ((stop_pos_n0 - start_pos_n1) +1)
                        if coverage > max_coverage_len:
                            max_coverage_len = coverage
                            
        # second iteration with adjusted edge weight (penalized cumulative coverage) and draw edged with the adjusted edge weight                
        for i in xrange(len(amplicon_coords)-1):
                start_pos_n0,stop_pos_n0 = amplicon_coords[i][0], amplicon_coords[i][1]
                for j in xrange(len(amplicon_coords[i+1:])):
                    start_pos_n1,stop_pos_n1 = amplicon_coords[i+j+1][0], amplicon_coords[i+j+1][1]
                    if (stop_pos_n0 >= start_pos_n1) and (stop_pos_n1 > stop_pos_n0) and (start_pos_n1 > start_pos_n0):    
                        coverage    = max_coverage_len - (((stop_pos_n1 - start_pos_n0)+1) - ((stop_pos_n0 - start_pos_n1) +1)) 
                        DG.add_weighted_edges_from([(amplicon_coords[i],amplicon_coords[i+j+1],coverage)])
   
        best_shortest_path    =    nx.dijkstra_path(DG,start_node, stop_node)
        bestshortest_paths_list.append(best_shortest_path)
    bestshortest_paths_list = dijkstra_max_cover_amplicon_sets(bestshortest_paths_list)
    
    for item in bestshortest_paths_list:
        if len(item) == 1:
            graph_no    += 1
            subnetwork_stats['subnetwork# '+ str(graph_no)]  = ['start = ' + str(item[0][0]), 'stop = ' + str(item[0][1]), 'length = '+str(int(item[0][1]) - int(item[0][0]) + 1)+ ' bp', '# primer pairs = ' + str(len(item))]
    
    all_nodes_connectedcomponents   = [item for sublist in all_nodes_connectedcomponents for item in sublist]
    bestshortest_paths_list         = [item for sublist in bestshortest_paths_list for item in sublist]
       
    ### separate odd and even elements of bestshortest_paths_list so that overlapping amplicons are not multiplexed together
    bestshortest_paths_list = sorted(bestshortest_paths_list)
    set1    = bestshortest_paths_list[::2]    # odd nodes
    set2    = bestshortest_paths_list[1::2]   # even nodes

    ### test multiplx_input_write() function
    if set1:
        multiplx_input_create(set1, multiplx_input_set1_1)
    if set2: 
        multiplx_input_create(set2, multiplx_input_set2_1)
    # end graph analysis

    ### summary stat for amplicon lengths
    if len(amplicon_len_list) >= 1:
            min, max, mean, median  =   int(min(amplicon_len_list)), int(max(amplicon_len_list)), int(round(sum(amplicon_len_list)/float(len(amplicon_len_list)))), int(round(np.median(amplicon_len_list)))
    else:
            min, max, mean, median  =   0, 0, 0, 0

    ### cumulative amplicon coverage calculation for primer pair picking step
    cumulative_amplicon_ranges      = spf.range_overlap_adjust(primer_pair_coords)                  # adjusts for the overlap between amplicons
    cumulative_amplicon_coverage    = spf.cumulative_ranges_coverage(cumulative_amplicon_ranges)    # caclulate cumulative amplicon coverage

    ### cumulative amplicon coverage calculation for multiplex step; calculate the cumulative coverage by selected amplicons in set1 & set2
    whole_list_mplex_amplicons = bestshortest_paths_list
    whole_list_mplex_amplicons = sorted(whole_list_mplex_amplicons)
    cumulative_mplex_amplicon_ranges      = spf.range_overlap_adjust(whole_list_mplex_amplicons)                # adjusts for the overlap between amplicons
    cumulative_mplex_amplicon_coverage    = spf.cumulative_ranges_coverage(cumulative_mplex_amplicon_ranges)    # caclulate cumulative amplicon coverage

    ### primer3 thermodynamic parameters based calculating score table and saving it to file
    if set1:
        ### real input file set1
        p1 = sp.Popen(["%scmultiplx" %multiplx_path, "-primers", "%s%s_multiplx_input_set1_1.txt" %(path,Time_stamp), "-thermodynamics", "%sthermodynamics.txt.primer3" %multiplx_path, "-csalt", "%s" %monovalent_salts, "-cmg", "%s" %Mg,"-cdna", "%s" %primer_conc,"-calcscores", "%s" %multiplx_calcscores, "-savescores", "%s%s_score_set1_1" %(path,Time_stamp)])
        output, error = p1.communicate()
    if set2:
        ### real input file set2
        p1 = sp.Popen(["%scmultiplx" %multiplx_path, "-primers", "%s%s_multiplx_input_set2_1.txt" %(path,Time_stamp), "-thermodynamics", "%sthermodynamics.txt.primer3" %multiplx_path, "-csalt", "%s" %monovalent_salts, "-cmg", "%s" %Mg,"-cdna", "%s" %primer_conc,"-calcscores", "%s" %multiplx_calcscores, "-savescores", "%s%s_score_set2_1" %(path,Time_stamp)])
        output, error = p1.communicate()

    ### Calculating multiplex groups from score table
    if set1:
        ### set1
        f0 = open(os.devnull, 'w')
        p2 = sp.Popen(["%scmultiplx" %multiplx_path, "-primers", "%s%s_multiplx_input_set1_1.txt" %(path,Time_stamp), "-scores", "%s%s_score_set1_1" %(path,Time_stamp), "-thermodynamics", "%sthermodynamics.txt.primer3" %multiplx_path, "-csalt", "%s" %monovalent_salts, "-cmg", "%s" %Mg,"-cdna", "%s" %primer_conc, "-calcgroups", "%s" %maxgroups, "%s" %maxitemsingroup, "-stringency", "%s" %grouping_stringency, "-savegroups", "%s%s_multiplex_groups_set1_1.txt" %(path,Time_stamp)],stderr=f0)
        output, error = p2.communicate()
    if set2:
        ### set2
        p2 = sp.Popen(["%scmultiplx" %multiplx_path, "-primers", "%s%s_multiplx_input_set2_1.txt" %(path,Time_stamp), "-scores", "%s%s_score_set2_1" %(path,Time_stamp), "-thermodynamics", "%sthermodynamics.txt.primer3" %multiplx_path, "-csalt", "%s" %monovalent_salts, "-cmg", "%s" %Mg,"-cdna", "%s" %primer_conc, "-calcgroups", "%s" %maxgroups, "%s" %maxitemsingroup, "-stringency", "%s" %grouping_stringency, "-savegroups", "%s%s_multiplex_groups_set2_1.txt" %(path,Time_stamp)],stderr=f0)
        output, error = p2.communicate()

    ### read multiplex_groups file to create a single out put file with primer info
    if set1:
        multiplex_groups_set1_1  = path + Time_stamp + "_" + 'multiplex_groups_set1_1.txt'
    if set2: 
        multiplex_groups_set2_1  = path + Time_stamp + "_" + 'multiplex_groups_set2_1.txt'
    if set1: 
        ### print output files for multiplex_group_output for first set and record last primer_pair_counter value
        primer_pair_counter = multiplex_group_output(multiplex_groups_set1_1, 1, 0)
    if set2:
        ### print output files for multiplex_group_output for second set and continue with primer_pair_counter
        primer_pair_counter = multiplex_group_output(multiplex_groups_set2_1, 2, primer_pair_counter)

    ############################################################
    # Time to run the code: end timer
    ############################################################
    t1 = time.time() 
    total = t1-t0
    total = ("{0:.2f}".format(round(total,2)))
    parameters_used.write(      
                                "no. of primers picked for all possible pairs : " + str(no_unique_primers_picked)+'\n'
                                "target region length                         : " + str(locus_len) + ' bp' + '\n'
                                "cumulative amplicon coverage (primer pairs)  : " + str(cumulative_amplicon_coverage)+' bp'+'\n'
                                "amplicon size stats                          : " + "min length = " + str(min) +"; max length = " + str(max) +"; mean length = " + str(mean) + "; median length = " + str(median) +' (bp)'+'\n'
                                "% of target region covered (primer pairs)    : " + str(("{0:.2f}".format(round((float(cumulative_amplicon_coverage)/locus_len)*100,2))))+' %'+'\n'
                                "subnetwork coverage stats                    : " + str(subnetwork_stats) + '\n'
                                "# minimum tiling primers                     : " + str(primer_pair_counter) + '\n'
                                "cumulative amplicon coverage (mplex)         : " + str(cumulative_mplex_amplicon_coverage)+' bp'+'\n'
                                "% of target region covered (mplex)           : " + str(("{0:.2f}".format(round((float(cumulative_mplex_amplicon_coverage)/locus_len)*100,2))))+' %'+'\n'
                                "### PPS run duration           : " + str(total) + " seconds"+'\n'
                                "##########################################################"+"\n"+
                                "\n"+"\n"
                        )
    parameters_used.close()

### Change log
### v1.02 -> v1.03
    # uses updated dependency thal.c code to predict heterodimer_dg
### v1.03 -> v1.04
    # Minimum tiling path identification using shortest path in graph approach


    
    
    
