# Defining codon table in dictionary
codontable = {
                'UUU': 'F',      'CUU': 'L',      'AUU': 'I',      'GUU': 'V',
                'UUC': 'F',      'CUC': 'L',      'AUC': 'I',      'GUC': 'V',
                'UUA': 'L',      'CUA': 'L',      'AUA': 'I',      'GUA': 'V',
                'UUG': 'L',      'CUG': 'L',      'AUG': 'M',      'GUG': 'V',
                'UCU': 'S',      'CCU': 'P',      'ACU': 'T',      'GCU': 'A',
                'UCC': 'S',      'CCC': 'P',      'ACC': 'T',      'GCC': 'A',
                'UCA': 'S',      'CCA': 'P',      'ACA': 'T',      'GCA': 'A',
                'UCG': 'S',      'CCG': 'P',      'ACG': 'T',      'GCG': 'A',
                'UAU': 'Y',      'CAU': 'H',      'AAU': 'N',      'GAU': 'D',
                'UAC': 'Y',      'CAC': 'H',      'AAC': 'N',      'GAC': 'D',
                'UAA': 'Stop',   'CAA': 'Q',      'AAA': 'K',      'GAA': 'E',
                'UAG': 'Stop',   'CAG': 'Q',      'AAG': 'K',      'GAG': 'E',
                'UGU': 'C',      'CGU': 'R',      'AGU': 'S',      'GGU': 'G',
                'UGC': 'C',      'CGC': 'R',      'AGC': 'S',      'GGC': 'G',
                'UGA': 'Stop',   'CGA': 'R',      'AGA': 'R',      'GGA': 'G',
                'UGG': 'W',      'CGG': 'R',      'AGG': 'R',      'GGG': 'G'
}
# Open file and setting parameters
Seqfile = open('rosalind_prot.txt', 'r')
RNAseq = Seqfile.read()
protseq = ''
# Running codon check for every 3 base
for i in range(0, len(RNAseq), 3):
    codon = RNAseq[i:i+3]
# return AA from codontable
    AA = codontable[codon]
# stop translation upon stop codon
    if AA == 'Stop':
        break
# otherwise add AA sequence to protseq
    else:
        protseq += AA
# output result
print protseq
