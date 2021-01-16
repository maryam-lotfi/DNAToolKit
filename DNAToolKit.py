# DNA toolkit file
import collections

from structures import *


# check the sequence to make sure it is a DNA string
def validateseq(dna_seq):
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in nucleotides:
            return False
    return tmpseq


# DNA Frequency counter
def countnucfreq(seq):
    tmpfreqdict = {"A": 0, "T": 0, "G": 0, "C": 0}
    for nuc in seq:
        tmpfreqdict[nuc] += 1
    return tmpfreqdict


# DNA -> RNA Transcription
def transcription(seq):
    return seq.replace("T", "U")


# DNA complement
def complement(seq):
    return ''.join([Nucleotides[nuc] for nuc in seq])


# DNA reverse complement
def reverse_complement(seq):
    return ''.join([Nucleotides[nuc] for nuc in seq])[::-1]


# faster approach
# def reverse(seq):
#  mapping = str.maketrans('ATCG', 'TAGC')
#  return seq.translate(mapping)

# gc content calculation

def gc_content(seq):
    return round(seq.count('C') + seq.count('G') / len(seq) * 100)


# calculate GC content in particular part of the genome subseq
def gc_content_subsec(seq, k=20):
    res = []
    for i in range(0, len(seq) - k + 1, k):
        subseq = seq[i:i + k]
        res.append(gc_content(subseq))
    return res


# translate a DVA sequence into an amino acid sequence
def translate_seq(seq, init_pos=0):
    return [DNA_Codons[seq[pos:pos + 3]] for pos in range(init_pos, len(seq) - 2, 3)]


# provide the frequency of each a particular codon in a DNA sequence
# provide the frequency of each codon encoding a given amino acid in a DNA sequence
def codon_usage(seq, aminoacid):
    tmplist = []
    for i in range(0, len(seq) - 2, 3):
        if DNA_Codons[seq[i:i + 3]] == aminoacid:
            tmplist.append(seq[i:i + 3])

    freq_dict = dict(collections.Counter(tmplist))
    total = sum(freq_dict.values())
    for seq in freq_dict:
        freq_dict[seq] = round(freq_dict[seq] / total, 2)
    return freq_dict


# Generate the six reading frames od a DNA sequence including the reverse complement
def gen_reading_frames(seq):
    frames = []
    frames.append(translate_seq(seq, 0))
    frames.append(translate_seq(seq, 1))
    frames.append(translate_seq(seq, 2))
    frames.append(translate_seq(reverse_complement(seq), 0))
    frames.append(translate_seq(reverse_complement(seq), 1))
    frames.append(translate_seq(reverse_complement(seq), 2))

    return frames

# Compute all possible proteins in an aminoacid sequence and return a list of possible proteins
def proteins_from_rf(this_seq):
    current_protein = []
    proteins = []
    for a in this_seq:
        if a == "_":

            # STOP accumulating amino acids if _ was found
            if current_protein:
                for p in current_protein:
                    proteins.append(p)
                current_protein = []
        else:
            # START accumulating aminoacids if M - start was found
            if a == "M":
                current_protein.append("")
            for i in range(len(current_protein)):
                current_protein[i] += a
    return proteins


# Compute all possible proteins for all open reading frames
def all_proteins_rf(seq, startpos=0, endpos=0, ordered=False):
    if endpos > startpos:
        rfs = gen_reading_frames(seq[startpos: endpos])
    else:
        rfs = gen_reading_frames(seq)

    result = []
    for rf in rfs:
        proteins = proteins_from_rf(rf)
        for p in proteins:
            result.append(p)

    if ordered:
        return sorted(result, key=len, reverse=True)
    return result
