"""Extract sequences from a given FASTA file that match the given list of sample names.
"""
from numpy import array, inf
from utils import safe_translate

scoring_params = {"score_match": 3, "score_mismatch": -
                  1, "score_gapext": -1, "score_gapopen": -10}


def align_pairwise(seq1, seq2):
    from Bio import pairwise2
    aln = pairwise2.align.globalms(seq1, seq2,
                                   scoring_params['score_match'], scoring_params['score_mismatch'],
                                   scoring_params['score_gapopen'], scoring_params['score_gapext'],
                                   penalize_end_gaps=False, one_alignment_only=True)[0]
    return aln[2], aln[0], aln[1]


def get_cds(ref):
    '''
    assuming there is one contiguous coding region which might be
    split into multiple sub-proteins like HA1 and HA2.
    loop over all features, pull out min and max of their union
    '''
    cds_start, cds_end = inf, 0
    for feature in ref.features:
        if feature.type == 'CDS':
            if feature.location.start < cds_start:
                cds_start = feature.location.start
            if feature.location.end > cds_end:
                cds_end = feature.location.end

    refstr = str(ref.seq).upper()
    refCDS = refstr[cds_start:cds_end]
    refAA = safe_translate(refstr[cds_start:cds_end])
    return refstr, refCDS, refAA, cds_start, cds_end


def codon_align(seq, refstr, refAA, cds_start, cds_end):
    seqstr = str(seq.seq).upper()
    score, refaln, seqaln = align_pairwise(refstr, seqstr)
    if score < 0:  # did not align
        return None
    ref_aln_array = array(list(refaln))
    seq_aln_array = array(list(seqaln))

    # stip gaps
    ungapped = ref_aln_array != '-'
    #ref_aln_array_ungapped = ref_aln_array[ungapped]
    seq_aln_array_ungapped = seq_aln_array[ungapped]

    seq5pUTR = "".join(seq_aln_array_ungapped[:cds_start])
    seq3pUTR = "".join(seq_aln_array_ungapped[cds_end:])
    seqCDS = "".join(seq_aln_array_ungapped[cds_start:cds_end])
    seqCDS_ungapped = seqCDS.replace('-', '')
    seqAA = safe_translate(seqCDS_ungapped)

    scoreAA, refalnAA, seqalnAA = align_pairwise(refAA, seqAA)
    if scoreAA < 0 or sum(seqAA.count(x) for x in ['*', 'X']) > 5 or refalnAA.count('-') > 5:
        #print(seq.id, "didn't translate properly", file=sys.stderr)
        return None

    seqCDS_aln = seq5pUTR
    pos = 0
    for aa_ref, aa_seq in zip(refalnAA, seqalnAA):
        if aa_seq == '-':
            seqCDS_aln += '---'
            # if the nucleotide sequence is gapped
            # (i.e. because of missing data at the 5p and 3p end, advance pos)
            if seqCDS_ungapped[pos:pos+3] == '---':
                pos += 3
        else:
            if len(seqCDS_ungapped) >= pos+3:
                seqCDS_aln += seqCDS_ungapped[pos:pos+3]
            else:
                seqCDS_aln += '---'
            pos += 3

    return ''.join(seqCDS_aln)+seq3pUTR
