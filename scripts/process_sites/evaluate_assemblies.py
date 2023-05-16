import os
import pandas as pd
import numpy as np
from itertools import groupby
import gzip

"""
Functions fasta_iter(), read_genome() and calculate_stats() are from:
https://github.com/MikeTrizna/assembly_stats/blob/master/assembly_stats/assembly_stats.py
fasta_iter() has been modified by me to handle gzipped files.
"""

def fasta_iter(fasta_file):
    with gzip.open(fasta_file, mode='rt') as fh:
        fa_iter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
        for header in fa_iter:
            header = next(header)[1:].strip()
            seq = "".join(s.upper().strip() for s in next(fa_iter))
            yield header, seq

def read_genome(fasta_file):
    """Takes a FASTA file, and produces 2 lists of sequence lengths. It also
    calculates the GC Content, since this is the only statistic that is not
    calculated based on sequence lengths.

    Parameters
    ----------
    fasta_file : str
        The file location of the FASTA file

    Returns
    -------
    contig_lens: list
        A list of lengths of all contigs in the genome.
    scaffold_lens: list
        A list of lengths of all scaffolds in the genome.
    gc_cont: float
        The percentage of total basepairs in the genome that are either G or C.
    """

    gc = 0
    total_len = 0
    contig_lens = []
    scaffold_lens = []
    for _, seq in fasta_iter(fasta_file):
        scaffold_lens.append(len(seq))
        if "NN" in seq:
            contig_list = seq.split("NN")
        else:
            contig_list = [seq]
        for contig in contig_list:
            if len(contig):
                gc += contig.count('G') + contig.count('C')
                total_len += len(contig)
                contig_lens.append(len(contig))
    gc_cont = (gc / total_len) * 100
    return contig_lens, scaffold_lens, gc_cont


def calculate_stats(seq_lens, gc_cont):
    stats = {}
    seq_array = np.array(seq_lens)
    stats['sequence_count'] = seq_array.size
    stats['gc_content'] = gc_cont
    sorted_lens = seq_array[np.argsort(-seq_array)]
    stats['longest'] = int(sorted_lens[0])
    stats['shortest'] = int(sorted_lens[-1])
    stats['median'] = np.median(sorted_lens)
    stats['mean'] = np.mean(sorted_lens)
    stats['total_bps'] = int(np.sum(sorted_lens))
    csum = np.cumsum(sorted_lens)
    for level in [10, 20, 30, 40, 50]:
        nx = int(stats['total_bps'] * (level / 100))
        csumn = min(csum[csum >= nx])
        l_level = int(np.where(csum == csumn)[0])
        n_level = int(sorted_lens[l_level])
        stats['L' + str(level)] = l_level + 1
        stats['N' + str(level)] = n_level
    return stats


asm_dirs_path='/home/jay/master_thesis/results/process_sites/03_megahit'
#create pandas dict with headers assembly, lenght, GC, nb_contigs other things?
asm_stats = pd.DataFrame(columns=['assembly', 'sequence_count', 'gc_content', 
                                  'longest_contig', 'shortest_contig', 'median', 'mean', 
                                  'total_bps', 'L10', 'N10', 'L20', 'N20', 'L30',
                                  'N30', 'L40', 'N40', 'L50', 'N50'])

for dir in os.listdir(asm_dirs_path):
    for asm in os.listdir(asm_dirs_path+'/'+dir):     
        contig_lens, scaffold_lens, gc_cont = read_genome(asm_dirs_path+'/'+dir+'/'+asm)
        contig_stats = calculate_stats(contig_lens, gc_cont)
        new_row = {'assembly':dir,'sequence_count':contig_stats['sequence_count'],
                   'gc_content':contig_stats['gc_content'],
                   'longest_contig':contig_stats['longest'],
                   'shortest_contig':contig_stats['shortest'],
                   'median':contig_stats['median'],
                   'mean':contig_stats['mean'],
                   'total_bps':contig_stats['total_bps'],
                   'L10':contig_stats['L10'],
                   'N10':contig_stats['N10'],
                   'L20':contig_stats['L20'],
                   'N20':contig_stats['N20'],
                   'L30':contig_stats['L30'],
                   'N30':contig_stats['N30'],
                   'L40':contig_stats['L40'],
                   'N40':contig_stats['N40'],
                   'L50':contig_stats['L50'],
                   'N50':contig_stats['N50']}
        asm_stats = pd.concat([asm_stats, pd.DataFrame([new_row])], ignore_index=True)

asm_stats.to_csv("/home/jay/master_thesis/results/process_sites/sites_asm_stats.tsv", sep='\t', index=False)

