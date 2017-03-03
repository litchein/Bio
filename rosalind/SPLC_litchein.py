# Read file and parse sequences
input_file = open('rosalind_splc.txt').read().split('>')[1:]
seq_list = []
for i in input_file:
    gene = i.splitlines()[0]
    seq = ''.join(i.splitlines()[1:])
    seq_list.append([gene, seq])
template = seq_list[0]
introns = seq_list[1:]

# Removing all introns from template
for i in introns:
    i_seq = i[1]
    template[1] = template[1].replace(i_seq, '')


# Translate into protein
def trans_prot(template_seq):
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
    template_rna = template_seq.replace('T', 'U')
    prot = ''
    for n in xrange(0, len(template_rna), 3):
        codon = template_rna[n:n+3]
        if codontable[codon] == 'Stop':
            break
        prot += codontable[codon]
    return prot

# Output protein sequence
print trans_prot(template[1])
