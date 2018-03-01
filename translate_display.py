from string import maketrans

def reverse_complement(dna):
    com = dna.translate(maketrans('ATCG', 'TAGC'))
    return com[::-1]

def gc(dna):
    gc_count = 0.0
    for base in dna:
        if base == 'G' or base == 'C':
            gc_count += 1
    return gc_count / len(dna) * 100

def tm(dna):
    #      Tm= 64.9 +41*(yG+zC-16.4)/(wA+xT+yG+zC)
    count_a = dna.count('A')
    count_t = dna.count('T')
    count_c = dna.count('C')
    count_g = dna.count('G')
    tm = 64.9 + (41*(count_c + count_g - 16.4) / (len(dna)))
    return tm


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
ste20_dna = open('ste20.txt').read()
ste20_rna = ste20_dna.replace('T', 'U')
protseq = ''
aa_seq = ''
for i in range(0, len(ste20_rna), 3):
    codon = ste20_rna[i:i+3]
# return AA from codontable
    AA = codontable[codon]
    if AA == 'Stop':
        break
    else:
        protseq += ' %s ' % AA
        aa_seq += AA

# output result
print protseq
print '\n'

def getdomain(seq,start,end,name='Domain'):
    print '%s is %s' % (name, seq[start-1:end])
    print '%s residues' % len(seq[start-1:end])
    loc = (start -1) * 3
    dna_len = len(seq[start-1:end]) * 3
    print 'DNA sequence: %s \n' % ste20_dna[loc:loc+dna_len]


def display(DNA,PROT,char=33,out_file='output_disp.txt'):
    output = open(out_file, 'w')
    for x in range(0,len(DNA),char):
        output.write('%s\n' % DNA[x:x+char])
        output.write('%s\n' % PROT[x:x+char])
    output.close()


def make_del_primer(DNA, aa_start, aa_end,length=25):
    start = (aa_start * 3) - 2
    end = aa_end * 3
    fwd_primer = DNA[start-8:start-1]+DNA[end:end+(length-7)]
    rev_primer = reverse_complement(DNA[start-(length-7):start-1]+DNA[end:end+8])
    print (length-16)*' ','%s' % fwd_primer
    print '%s' % rev_primer[::-1]
    print '\n'
    print 'forward: \t %s' % fwd_primer
    print '\tGC: %.2f \tTm: %.2f\t%5s mer' % (gc(fwd_primer), tm(fwd_primer), len(fwd_primer))

    print 'reverse: \t %s' % rev_primer
    print '\tGC: %.2f \tTm: %.2f\t%5s mer' % (gc(rev_primer), tm(rev_primer), len(rev_primer))
    print '\n\n'


# Middle domain is 818-916
# Domain 4 is 926-1025
getdomain(aa_seq,818,916,'Middle domain')
getdomain(aa_seq,926,1025,'Domain 4')
# display(DNAseq,protseq,66)
print 'Generating delta_M primers..\n'
make_del_primer(ste20_dna, 818, 916, 30)
print 'Generating delta_D4 primers..\n'
make_del_primer(ste20_dna,926, 1025, 30)
print 'Generating delta_M_D4 primers..\n'
make_del_primer(ste20_dna, 818, 1025, 30)
