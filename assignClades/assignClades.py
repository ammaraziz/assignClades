#!/usr/bin/env python3

# TO DO:
#   1. add checks to input options to help user

import argparse
import sys
import os
from Bio import SeqIO
from codon_align import get_cds, codon_align
from utils import load_features, load_relatives, read_in_clade_definitions, is_node_in_clade, safe_translate

cwd = os.path.dirname(os.path.realpath(__file__))
print(cwd)

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

    reference : seq.record
        reference (same as alignment)
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

    print(out, file=sys.stdout)
    return(clade_final, like)


class ResultsBucket:

    def __init__(self):
        import pandas as pd
        results_bucket = pd.DataFrame({'Seq No': [],
                                       'HA Clade': [],
                                       'Result (relative)': [],
                                       'Ha Clade Provanence': [],
                                       'Amino Acid Mutation': [],
                                       'H275Y': [],
                                       'S31N': [],
                                       'I38X': []
                                       })

        results_bucket.set_index('Seq No', inplace=True)
        self.df = results_bucket

    def __str__(self):
        return(self.df.to_string())

    __repr__ = __str__

    def add_result(self, seqno, ha_clade, result, prov=None, aa_mut=None, h275y=None, s31n=None, i38x=None):
        '''
        add row of results, everything is a str
        '''
        if self.seq_in_data(seqno) is True:
            raise ValueError("Entry exists. Check code")
        self.df.loc[seqno] = [ha_clade, result, prov, aa_mut, h275y, s31n, i38x]

    def mod_result(self, seqno, **values):
        '''
        updates row where *column are the columns to update and *values are the entries
        values is a dictionary
        '''
        if self.seq_in_data(seqno) is False:
            raise ValueError("Entry does not exist. Check input.")
        else:
            for col, val in values.items():
                self.df[col].loc[seqno] = val

    def seq_in_data(self, item):
        return(item in self.df.index)

    def write_results(provancence, results_relative, output_file):
        with open(output_file, 'a') as f:
            print(provancence, file=f)
            print(f"{', '.join(provancence)}\t{results_relative}", file=f)


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

    # output files
    prov_out = args.batchName + "_provanence.txt"
    results_out = args.batchName + "_cladeResults.txt"
    errors_out = args.batchName + "_error.txt"

    # results
    results_bucket = ResultsBucket()

    for seq in input_sequences:

        seq_container = tmpNode()
        seq_aln = codon_align(seq, refstr, refAA, cds_start, cds_end)
        # error checking
        if seq_aln is None:
            print(f"{seq.id}\tError translating, check lineage and correct", file=sys.stdout)
            with open(errors_out, 'a') as ef:
                print(f"{seq.id}\tError translating, check lineage and correct", file=ef)
            continue

        clade_provanence = get_provanence(seq_aln, features, clade_designations, ref)
        # write out results
        with open(prov_out, 'a') as cf:
            print(f"{seq.description}\t{', '.join(clade_provanence)}", file=cf)

        clade_final = get_likeness(seq, clade_provanence, clades_relatives, internal_clades)
        ha_clade, fuz_result = get_likeness(seq, clade_provanence, clades_relatives, internal_clades)
        with open(results_out, 'a') as rf:
            print(clade_final, file=rf)
        print(clade_final, file=sys.stdout)

        results_bucket.add_result(seqno=seq.description, ha_clade=ha_clade, result=fuz_result,
                                  prov=', '.join(clade_provanence))
        print(results_bucket)


if __name__ == '__main__':
    main()
