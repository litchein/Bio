import urllib
#open file and retrieve seq into list
accessID = open('rosalind_mprt.txt','r').read().splitlines()
protseq = []
for ID in accessID:
    web = 'http://www.uniprot.org/uniprot/%s.fasta' % ID
    seqfile = urllib.urlopen(web)
#skips first line in fasta file
    seqfile.readline()
#to remove \n, read into list with splitlines and join back
    seq = seqfile.read().splitlines()
    seq = ''.join(seq)
#build in to protseq list
    protseq.append([ID, seq])

#setting parameters
output = []
i = 0

#search motif N{P}[ST]{P} and record matching positions
for entry in protseq:
    match = []
    for j in xrange(len(protseq[i][1])-3):
        if protseq[i][1][j] == 'N':
            if protseq[i][1][j+1] != 'P':
                if protseq[i][1][j+2] == 'S' or protseq[i][1][j+2] == 'T':
                    if protseq[i][1][j+3] != 'P':
                        match.append(str(j+1))

#only add to output if there is a match
    if match != []:
        output.append([ protseq[i][0], match ])
    i += 1

#print output
for k in xrange(len(output)):
#prints ID
    print output[k][0]
#prints positions without brackets
    print ' '.join(output[k][1])
