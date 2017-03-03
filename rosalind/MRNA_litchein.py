# Read file for sequence and add Stop as *
input_file = open('rosalind_mrna.txt').read().splitlines()
prot_seq = ''.join(input_file)
prot_seq += '*'

# Populate amino codon frequency table from codon table
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
amino_codon_f = {}
for codon in codontable:
    if codontable[codon] == 'Stop':
        codontable[codon] = '*'
    if codontable[codon] in amino_codon_f.keys():
        amino_codon_f[codontable[codon]] += 1
    else:
        amino_codon_f[codontable[codon]] = 1

# Calculate # of RNA strings, modulo 1,000,000
rna_string = 1
mod = 1000000
for amino in prot_seq:
    rna_string *= amino_codon_f[amino]
    rna_string %= mod

# Output results
print rna_string


