from string import maketrans
# Setting variables
codontable = {
'UUU' :'F',      'CUU': 'L',      'AUU': 'I',      'GUU': 'V',
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

# Reading and arranging sequence and reverse complement
input_file = open('rosalind_orf.txt').read().split('>')[1:]
gene = input_file[0].splitlines()[0]
seq = ''.join(input_file[0].splitlines()[1:])
seq_revc = seq.translate(maketrans('ATCG','TAGC'))[::-1].replace('T', 'U')
seq = seq.replace('T', 'U')
seq = [seq, seq_revc]
# Search and store list of possible proteins
prot_list =[]
for i in seq:
    for n in xrange(0, len(i)-3):
        prot_seq = ''
        if i[n:n+3] == 'AUG':
            for m in xrange(n, len(i)-3, 3):
                codon = i[m:m+3]
                if codon > 3:
                    if codontable[codon] == 'Stop':
                        prot_list.append(prot_seq)
                        break
                    else:
                        prot_seq += codontable[codon]
                else:
                    prot_list.append(prot_seq)
                    break

# Prints unique list of proteins
for i in set(prot_list):
    print i