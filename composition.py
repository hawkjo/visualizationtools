import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
from misctools import gzip_friendly_open
from Bio import SeqIO

def get_composition_consensus_quality(fpath=None, records=None):
    if fpath is not None:
        with gzip_friendly_open(fpath) as f:
            rec = next(SeqIO.parse(f, 'fastq'))
            seq_len = len(rec)
        records_iterator = SeqIO.parse(gzip_friendly_open(fpath), 'fastq')
    else:
        assert records is not None
        seq_len = len(records[0])
        records_iterator = records

    bases = 'ACGT'
    composition_dict = {b: np.zeros((seq_len,)) for b in bases}
    quality_counts = np.zeros((41, seq_len))
    num_seqs = 0
    for rec in records_iterator:
        num_seqs += 1
        for b in bases:
            composition_dict[b] += np.array([int(c == b) for c in str(rec.seq)])
        for i, q in enumerate(rec.letter_annotations['phred_quality']):
            quality_counts[q, i] += 1
    for b in bases:
        composition_dict[b] = composition_dict[b] / float(num_seqs)
        
    consensus = []
    for i in range(len(composition_dict['A'])):
        cons_base = max(bases, key=lambda b: composition_dict[b][i])
        if composition_dict[cons_base][i] < 0.5:
            cons_base = 'N'
        consensus.append(cons_base)
        
    weighted_quality = np.array([quality_counts[i]*i for i in range(41)])
    avg_quality = weighted_quality.sum(axis=0)/num_seqs
    
    return composition_dict, ''.join(consensus), avg_quality


def composition_plot(intended=None, xlim=(0, 140), run_name='', **kwargs):
    """
    Draw composition plot. Must either contain named argument

        fpath:  Path to R1 read file.

    or all three named arguments

        composition_dict, consensus, avg_quality:   from get_composition_consensus_quality
    """
    try:
        composition_dict = kwargs['composition_dict']
        consensus = kwargs['consensus']
        avg_quality = kwargs['avg_quality']
    except:
        fpath = kwargs['fpath']
        composition_dict, consensus, avg_quality = get_composition_consensus_quality(fpath)

    fig, ax = plt.subplots(figsize=(20, 7))
    bases = 'ACGT'
    for b in bases:
        ax.plot(composition_dict[b], label=b)
    ax.set_xticks(range(len(consensus)))
    ax.set_xticklabels(consensus)
    if intended is None:
        xlabel = 'Consensus Sequence'
    elif consensus.startswith(intended):
        xlabel = 'Consensus Sequence (Matches intended sequence)'
    else:
        xlabel = 'Consensus Sequence (Does NOT match intended sequence)'
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Fraction')

    ax2 = ax.twinx()
    ax2.plot(avg_quality, '--', color='grey', alpha=0.5, label='Phred Score')
    ax2.set_ylabel('Phred Scores')

    ax.set_xlim(xlim)
    ax.grid(False)
    ax.legend()
    ax.set_title('%s Base Composition and Quality' % run_name)
    return ax
