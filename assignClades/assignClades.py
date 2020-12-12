#!/usr/bin/env python3

# TO DO:
#   1. add checks to input options to help user

import argparse, sys, os
from collections import defaultdict
from Bio import SeqIO
from codon_align import get_cds, codon_align
from utils import load_features, load_relatives, read_in_clade_definitions, is_node_in_clade, safe_translate

cwd = os.path.dirname(os.path.realpath(__file__))

def msg(name=None):
    return ''' Assign clades to HA influenza sequences.

    Usage: assignClades.py -s seq.fasta -l h3n2 -b outputName
    -s, --sequence          Path to input sequence [YourSequences.fasta]
    -l, --lineage           Lineage of input strains [h1n1, h3n2, vic, yam]
    -b, --batchName         The prefix (or batch name) which will be appended to the output files

    Example usage:
        assign_clades.py --sequences Batch999_01Jan20.fasta --lineage h1n1 --batchName Batch999_results.txt

    Output:
        batchName_cladeResults.txt          complete clade provenance
        batchName_provenance.txt            current clade and vaccine result
        batchName_error.txt                 sequences that did not align or translate

    Clade defintions are stored in /config/{lineage}.tsv
    '''

def register_arguments():
    parser = argparse.ArgumentParser(
        description="Assign clades to sequences",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        usage=msg()
    )

    parser.add_argument("-s", "--sequences", required=True, help="FASTA file of HA sequences")
    parser.add_argument("-l", "--lineage", required=True,
                        help="Lineage of the sequences supplied - h1n1, h3n2, vic, yam")
    parser.add_argument("-b", "--batchName", required=True, help="The name/path to write results")
    args = parser.parse_args()

    # lineage
    if args.lineage not in ['h1n1', 'h3n2', 'vic', 'yam']:
        raise ValueError(
            f"Incorrect lineage entered: {args.lineage}. Choose one: h1n1, h3n2, vic or yam")
    return(args)

class tmpNode(object):
    def __init__(self):
        self.sequences = {}

def get_provanence(seq_aln, features, clade_designations, reference):
    '''
    get_provanence 

    Parameters
    ----------
    seq_aln : str
        sequence aligned to reference
    features : 
        genbank annotations
    clade_designations 
        clade designations as defined by tsv file (nextstrain)
    
    Returns
    ----------
    clade_matches : list
        list of matched clades
    '''

    seq_container = tmpNode()

    seq_container.sequences['nuc'] = {i: c for i, c in enumerate(seq_aln)}
    for fname, feat in features.items():
        if feat.type != 'source':
            seq_container.sequences[fname] = {
                i: c for i, c in enumerate(safe_translate(feat.extract(seq_aln)))}
    clade_matches = []
    for clade_name, clade_alleles in clade_designations.items():
        if is_node_in_clade(clade_alleles, seq_container, reference):
            clade_matches.append(clade_name)
    
    return(clade_matches)

def get_likeness(seq, provanence, clades_relatives, internal_clades):
    '''
    get_likeness - finds closest virus relative (-like)

    Parameters
    ----------
    seq : SeqRecord
        raw seq record (not aligned) 
    provanence : list of str
        clade provanence in a list
    clades_relatives : dict
        dictionary of relatives (-like) viruses
    internal_clades : dict
        For matching specific clades that are only used in the centre
    
    Returns
    -------
    out : str
        string formatted for output
    '''    
    clade_final = provanence.pop(-1)  # get the last clade

    if clade_final in clades_relatives.keys():
        like = clades_relatives[clade_final]
        out = f"{seq.description}\t{clade_final}\t{like}"

    elif clade_final in internal_clades:
        tmp_like = internal_clades[clade_final]
        like = clades_relatives[tmp_like]
        out = f"{seq.description}\t{clade_final}\t{like}"

    else:
        out = f"{seq.description}\tUnknown clade found. Check clade defintions in config/clade_relatives.tsv"
    
    print(out, file = sys.stdout)
    return(out)

def write_results(provancence, results_relative, output_file):
    
    with open(output_file, 'a') as f:
        print(provancence, file = f)
        print(f"{', '.join(provancence)}\t{results_relative}", file = f)

def main():
    # register arguments
    args = register_arguments()

    input_sequences = SeqIO.parse(args.sequences, 'fasta')
    clade_designations = read_in_clade_definitions(f"config/clades_{args.lineage}_ha.tsv")

    refname = (f"config/reference_{args.lineage}_ha.gb")
    ref = SeqIO.read(refname, 'genbank')
    features = load_features(refname)
    refstr, refCDS, refAA, cds_start, cds_end = get_cds(ref)

    # get clade internal clade and likeness
    clades_relatives, internal_clades = load_relatives()

    results_out = args.batchName + "_cladeResults.txt"
    errors_out = args.batchName + "_error.txt"

    for seq in input_sequences:

        seq_container = tmpNode()
        seq_aln = codon_align(seq, refstr, refAA, cds_start, cds_end)
        # error checking
        if seq_aln is None:
            print(f"{seq.id}\tError translating, check lineage and correct", file=sys.stdout)
            with open(errors_out, 'a') as ef:
                print(f"{seq.id}\tError! Sequence was not translatable. Did you select the correct lineage?", file=ef)
            continue

        seq_container.sequences['nuc'] = {i: c for i, c in enumerate(seq_aln)}
        for fname, feat in features.items():
            if feat.type != 'source':
                seq_container.sequences[fname] = {
                    i: c for i, c in enumerate(safe_translate(feat.extract(seq_aln)))}

        matches = []
        for clade_name, clade_alleles in clade_designations.items():
            if is_node_in_clade(clade_alleles, seq_container, ref):
                matches.append(clade_name)

        # write out results
        with open(clades_out, 'a') as cf:
            print(f"{seq.description}\t{', '.join(matches)}", file=cf)

        # Get most recent clade  and find likeness ("Result")
        if matches:
            clade_final = matches.pop(-1)  # get the last clade

            if clade_final in clades_relatives.keys():
                like = clades_relatives[clade_final]
                out = f"{seq.description}\t{clade_final}\t{like}"

            elif clade_final in internal_clades:
                tmp_like = internal_clades[clade_final]
                like = clades_relatives[tmp_like]
                out = f"{seq.description}\t{clade_final}\t{like}"

            else:
                out = f"{seq.description}\tUnknown clade found. Check clade defintions in config/clade_relatives.tsv"

            with open(results_out, 'a') as rf:
                print(out, file=rf)
            print(out, file=sys.stdout)

if __name__ == '__main__':
    main()
