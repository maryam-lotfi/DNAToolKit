# DNA Toolset/Code testing file
from DNAToolKit import *
# from utilities import colored
from structures import *
import random
# rnd_dna_str = "acctttacgaatcccttcacgggttatataccctcatctcttacggtctg"
# creating a random DNA sequence for testing
rnd_dna_str = ''.join([random.choice(nucleotides)
                       for nuc in range(30)])

DNAstr = validateseq(rnd_dna_str)

print(f'\nSequence: {DNAstr}\n')
print(f'[1] + Sequence Length: {len(DNAstr)}\n')
print(f'[2] + Nucleotides Frequency: {countnucfreq(DNAstr)}\n')
print(f'[3] + DNA -> RNA: {transcription(DNAstr)}\n')
print(f"[4] + DNA String + complement + Reverse Complement:\n5' {DNAstr} 3'")
print(f"   {''.join(['|' for c in range(len(DNAstr))])}")
print(f"3' {complement(DNAstr)} 5'")
print(f"5' {reverse_complement(DNAstr)} 3'\n")
print(f'[5] + GC Content: {gc_content(DNAstr)}%\n')
print(f'[6] + GC content subsec k=5: {gc_content_subsec(DNAstr, k=5)}\n')
print(f'[7] + Aminoacids Sequence from DNA: {translate_seq(DNAstr)}\n')
print(f'[8] + Codon usage (L): {codon_usage(DNAstr, "L")}\n')
print('[9] + Reading frames:')
for frame in gen_reading_frames(DNAstr):
    print(frame)
print('\n[10] + All proteins in six open reading frames:')
for protein in all_proteins_rf(DNAstr, 0, 0, True):
    print(f'{protein}')
