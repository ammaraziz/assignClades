import os
from Bio import SeqIO, SeqFeature, Seq
from collections import defaultdict
from pandas import read_csv

def load_features(reference, feature_names=None):
    #read in appropriately whether GFF or Genbank
    #checks explicitly for GFF otherwise assumes Genbank
    if not os.path.isfile(reference):
        print("ERROR: reference sequence not found. looking for", reference)
        return None

    features = {}
    if '.gff' in reference.lower():
        #looks for 'gene' and 'gene' as best for TB
        try:
            from BCBio import GFF #Package name is confusing - tell user exactly what they need!
        except ImportError:
            print("ERROR: Package BCBio.GFF not found! Please install using \'pip install bcbio-gff\' before re-running.")
            return None
        limit_info = dict( gff_type = ['gene'] )

        with open(reference, encoding='utf-8') as in_handle:
            for rec in GFF.parse(in_handle, limit_info=limit_info):
                for feat in rec.features:
                    if feature_names is not None: #check both tags; user may have used either
                        if "gene" in feat.qualifiers and feat.qualifiers["gene"][0] in feature_names:
                            fname = feat.qualifiers["gene"][0]
                        elif "locus_tag" in feat.qualifiers and feat.qualifiers["locus_tag"][0] in feature_names:
                            fname = feat.qualifiers["locus_tag"][0]
                        else:
                            fname = None
                    else:
                        if "gene" in feat.qualifiers:
                            fname = feat.qualifiers["gene"][0]
                        else:
                            fname = feat.qualifiers["locus_tag"][0]
                    if fname:
                        features[fname] = feat

            if feature_names is not None:
                for fe in feature_names:
                    if fe not in features:
                        print("Couldn't find gene {} in GFF or GenBank file".format(fe))

    else:
        for feat in SeqIO.read(reference, 'genbank').features:
            if feat.type=='CDS':
                if "locus_tag" in feat.qualifiers:
                    fname = feat.qualifiers["locus_tag"][0]
                    if feature_names is None or fname in feature_names:
                        features[fname] = feat
                elif "gene" in feat.qualifiers:
                    fname = feat.qualifiers["gene"][0]
                    if feature_names is None or fname in feature_names:
                        features[fname] = feat
            elif feat.type=='source': #read 'nuc' as well for annotations - need start/end of whole!
                features['nuc'] = feat

    return features
    
def safe_translate(sequence, report_exceptions=False):
    """Returns an amino acid translation of the given nucleotide sequence accounting
    for gaps in the given sequence.

    Optionally, returns a tuple of the translated sequence and whether an
    exception was raised during initial translation.

    >>> safe_translate("ATG")
    'M'
    >>> safe_translate("ATGGT-")
    'MX'
    >>> safe_translate("ATG---")
    'M-'
    >>> safe_translate("ATGTAG")
    'M*'
    >>> safe_translate("")
    ''
    >>> safe_translate("ATGT")
    'MX'
    >>> safe_translate("ATG", report_exceptions=True)
    ('M', False)
    >>> safe_translate("ATGA-G", report_exceptions=True)
    ('MX', True)
    """
    from Bio.Data.CodonTable import TranslationError
    from Bio.Seq import CodonTable
    translation_exception = False

    #sequences not mod 3 give messy BiopythonWarning, so avoid by padding.
    if len(sequence)%3:
        sequence_padded = sequence + "N"*(3-len(sequence)%3)
    else:
        sequence_padded = sequence
    try:
        # Attempt translation by extracting the sequence according to the
        # BioPhython SeqFeature in frame gaps of three will translate as '-'
        translated_sequence = str(Seq.Seq(sequence_padded).translate(gap='-'))
    except TranslationError:
        translation_exception = True
        # Any other codon like '-AA' or 'NNT' etc will fail. Translate codons
        # one by one.
        codon_table  = CodonTable.ambiguous_dna_by_name['Standard'].forward_table
        str_seq = str(sequence_padded)
        codons = np.frombuffer(str_seq[:len(str_seq) - len(str_seq) % 3].encode(), dtype='S3').astype("U")
        assert len(codons) > 0
        aas = []

        for c in codons:
            # Parse result of single codon translation, add amino acids as
            # appropriate.
            try:
                aa = codon_table.get(c)
                if aa is None:
                    if c == '---':
                        aas.append('-')
                    else:
                        aas.append('X')
                else:
                    aas.append(aa)
            except (TranslationError, ValueError):
                aas.append('X')

        translated_sequence = "".join(aas)

    if report_exceptions:
        return translated_sequence, translation_exception
    else:
        return translated_sequence
        

def read_in_clade_definitions(clade_file):
    '''
    Reads in tab-seperated file that defines clades by amino acid or nucleotide mutations

    Format
    ------
    clade    gene    site alt
    Clade_1    ctpE    81  D
    Clade_2    nuc 30642   T
    Clade_3    nuc 444296  A
    Clade_4    pks8    634 T

    Parameters
    ----------
    clade_file : str
        meta data file

    Returns
    -------
    dict
        clade definitions as :code:`{clade_name:[(gene, site, allele),...]}`
    '''

    clades = defaultdict(list)
    df = read_csv(clade_file, sep='\t' if clade_file.endswith('.tsv') else ',')
    for index, row in df.iterrows():
        allele = (row.gene, row.site-1, row.alt)
        clades[row.clade].append(allele)
    clades.default_factory = None

    return clades


def is_node_in_clade(clade_alleles, node, ref):
    '''
    Determines whether a node matches the clade definition based on sequence
    For any condition, will first look in mutations stored in node.sequences,
    then check whether a reference sequence is available, and other reports 'non-match'

    Parameters
    ----------
    clade_alleles : list
        list of clade defining alleles
    node : Phylo.Node
        node to check, assuming sequences (as mutations) are attached to node
    ref : str/list
        positions

    Returns
    -------
    bool
        True if in clade

    '''
    conditions = []
    for gene, pos, clade_state in clade_alleles:
        if gene in node.sequences and pos in node.sequences[gene]:
            state = node.sequences[gene][pos]
        elif ref and gene in ref:
            state = ref[gene][pos]
        else:
            state = ''

        conditions.append(state==clade_state)

    return all(conditions)