from string import maketrans
#setting parameters and determine reverse complement seqence
seqfile = open('rosalind_revp.txt', 'r').read().splitlines()
dna = ''.join(seqfile[1:])
com = dna.translate(maketrans('ATCG','TAGC'))
revcom = com[::-1]

#search for palindrome in 4, 6, 8, 10 ,12 bases
for i in xrange(4,12,2):
    for y in xrange( len(dna) - i +1):
        pos = y-len(dna)
        if dna[pos:pos+i] == com[pos+i-1:pos-1:-1]:
            print str(y+1)+' '+str(i)
