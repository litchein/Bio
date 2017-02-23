from string import maketrans
#setting parameters and determine reverse complement seqence
seqfile = open('rosalind_revp.txt', 'r').read().splitlines()
dna = ''.join(seqfile[1:])
com = maketrans('ATCG','TAGC')

for i in xrange(4,13,2):
    for j in range(len(dna) - i +1):
        #generate slice and reverse complement
        forward = dna[j:j+i]
        complement = forward.translate(com)
        reverse = complement[::-1]


        if forward == reverse:
            print str(j+1) + ' ' + str(i)
