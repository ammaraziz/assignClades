#!/usr/bin/env python3

# TO DO:
#   1. add checks to input options to help user
#   2. split the main function into parsing and running

import argparse, sys, os
from Bio import SeqIO, Seq, SeqRecord
from codon_align import align_pairwise, get_cds, codon_align
from utils import  load_features, read_in_clade_definitions, is_node_in_clade, safe_translate

def msg(name = None):
    return ''' Assign clades to HA influenza sequences. 
    
    Usage: assignClades.py -s seq.fasta -l h3n2 -b outputName
    -s, --sequence          Path to input sequence [YourSequences.fasta]
    -l, --lineage           Lineage of input strains [h1n1, h3n2, vic, yam]
    -b, --batchName         The batch name
    
    Example usage: 
        assign_clades.py --sequences Batch999_01Jan20.fasta --lineage h1n1 --batchName Batch999_results.txt
    
    Output: 
        batchName_clades.txt        complete clade provenance
        batchName_results.txt       current clade and vaccine result 
    
    Clade defintions are stored in /config/{lineage}.tsv
    '''
# get clade internal clade and likeness
clades_relatives = {}
with open(f"config/clades_relative.tsv") as f:
    from csv import DictReader
    data = DictReader(f, delimiter = "\t", )
    for entry in data:
        clades_relatives[entry['clade']] = entry['relative']

internalDict = {
    "A1b/131K":"3C.2a1b+131K",
    "A1b/135K":"3C.2a1b+135K",
    "3c3.A":"3C.3a",
    "3c2.A1":"3C.2a1",
    "3c2.A1A":"3C.2a1a",
    "A3":"A1b/137F"}

class tmpNode(object):
    def __init__(self):
        self.sequences = {}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = "Assign clades to sequences",
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        usage = msg()
    )
    parser.add_argument("-s", "--sequences", required=True, help = "FASTA file of HA sequences")
    parser.add_argument("-l", "--lineage", required=True, help = "Lineage of the sequences supplied - h1n1, h3n2, vic, yam")
    parser.add_argument("-b", "--batchName", required=True, help = "The name of the batch")
    args = parser.parse_args()

    refname = f"config/reference_{args.lineage}_ha.gb"
    seqs = SeqIO.parse(args.sequences, 'fasta')
    ref = SeqIO.read(refname, 'genbank')
    features = load_features(refname)
    clade_designations = read_in_clade_definitions(f"config/clades_{args.lineage}_ha.tsv")
    results_out = args.batchName + "_results.txt"
    clades_out = args.batchName + "_clades.txt"

    # get sequence as string, CDS seq, amino acid sequence, and start/end pos
    refstr, refCDS, refAA, cds_start, cds_end = get_cds(ref)

    alignment = []
    for seq in seqs:
        seq_container = tmpNode()
        seq_aln = codon_align(seq,  refstr, refAA, cds_start, cds_end)
        if seq_aln is None:
            print(f"{seq.id}\tnot translatable", file=sys.stdout)
            continue

        seq_container.sequences['nuc'] = {i:c for i,c in enumerate(seq_aln)}
        for fname, feat in features.items():
            if feat.type != 'source':
                seq_container.sequences[fname] = {i:c for i,c in enumerate(safe_translate(feat.extract(seq_aln)))}

        matches = []
        for clade_name, clade_alleles in clade_designations.items():
            if is_node_in_clade(clade_alleles, seq_container, ref):
                matches.append(clade_name)

        #write out results
        with open(clades_out, 'a') as f:
            print(f"{seq.description}\t{', '.join(matches)}", file = f)

        # compress clade and find likeness
        if matches:
            clade_final = matches.pop(-1) # get the last clade
            if clade_final in clades_relatives.keys():
                like = clades_relatives[clade_final]
                out = f"{seq.description}\t{clade_final}\t{like}"
            elif clade_final in internalDict:
                tmp_like = internalDict[clade_final]
                like = clades_relatives[tmp_like]
                out = f"{seq.description}\t{clade_final}\t{like}"
            else:
                out = "unknown clade! check clade definitions"

            with open(results_out, 'a') as f:
                print(out, file = f)
            print(out, file=sys.stdout)
