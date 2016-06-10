### Repository for sequence processing functions (spf) for UOD, PSE and PPS algorithms
### Version-1.05: 06/10/2016
### Author: Felix Francis (felixfrancier@gmail.com); Under the guidance and financial support of Randall J Wisser (rjw@udel.edu) 


### function to get reverse complement of an input sequence
def rev_complement(seq):
    seq = seq.upper()
    basecomplement = {'A':'T', 
                      'C':'G', 
                      'G':'C', 
                      'T':'A', 
                      '-':'-', 
                      'N':'N'}
    letters = list(seq)
    letters = [basecomplement[base] for base in letters]
    complement = (''.join(letters))
    return complement[::-1]
    
    
### Get sequence from a fasta file (use for FASTA file with a single sequence)
def fasta_to_seq(input_file):
    with open (input_file) as fasta_data:
        line = fasta_data.read()
        lines = line.split("\n")
        sequence = lines[1]
    return sequence

### function to calculate gc content
def gc_content(sequence):
    sequence = sequence.upper()
    gc_count =  sequence.count('G')+sequence.count('C')
    sequence_len = len(sequence)
    gc_cont    = (float(gc_count)/sequence_len)*100
    return round(gc_cont, 2)


### get sequence based on given corrdinates, from an input fasta file
def seq_extraction_loci(loci_fasta_seq, actual_locus_start_pos, genomic_start_coord, genomic_stop_coord):  
    loci_start_coord    =   (genomic_start_coord - actual_locus_start_pos)+1
    loci_stop_coord     =   (genomic_stop_coord - actual_locus_start_pos)+1
    loci_seq            =   fasta_to_seq(loci_fasta_seq)
    select_seq          =   loci_seq [(loci_start_coord-1):(loci_stop_coord)]
    return select_seq

### range overlap function
def range_overlap_adjust(list_ranges):
    overlap_corrected   =   []
    for start, stop in sorted(list_ranges):
        if  overlap_corrected and (start-1 <= overlap_corrected [-1][1] and stop >= overlap_corrected [-1][1]):
            overlap_corrected [-1] = min(overlap_corrected [-1][0], start), stop
        elif overlap_corrected and (start <= overlap_corrected [-1][1] and stop <= overlap_corrected [-1][1]):
            pass
        else:
            overlap_corrected.append((start,stop))
    return overlap_corrected

### cumulative list of ranges coverage
def cumulative_ranges_coverage(list_ranges):
    sum = 0
    for start, stop in (list_ranges):
        sum += stop-start+1
    return sum




