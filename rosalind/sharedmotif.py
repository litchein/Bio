#open file and build list for FASTA
Seqfile = open('test.txt', 'r')
Seqlist = Seqfile.read().split('>')
#skip first empty >
Seqlist = Seqlist[1:]
genelist = []
#sorting genes into genelist by [ gene , seq ]
for gene in Seqlist:
    gene = gene.splitlines()
    gene= [ gene[0] , ''.join( gene[1:] ) ]
    genelist.append(gene)
#define the shortest sequence
short = min( genelist , key=lambda k:k[1])
maxlen = len( max(genelist , key=lambda k:k[1]) )
minlen = len( short[1] )
#remove shortest sequence from query genelist
for x in xrange( len(genelist) ):
    if genelist[x] == short:
        del genelist[x]

print short
print minlen
print genelist
print short[1][0:0+minlen]
print '\n'
match = 0
pos = 0
#Keep searching until match found
while match != len(genelist):
    match = 0
#Stop operation if not found
    if minlen < 2:
        print 'Not found!'
        break
#Going through the gene
    print 'Position:' + str(pos)
    for y in xrange(maxlen-pos):
        query = short[1][pos-1:pos+minlen-1]
        runs = len(gene[1]) - minlen +1
        for i in xrange(0,runs):
                for gene in genelist:
                    print 'Gene length:' + str(len(gene[1]))
                    print 'Minimum length:'+str(minlen)
                    print 'No. of runs:' + str(runs) +'\n'
                    print 'Template:' + gene[0]
                    print 'Searching:' +gene[1][i:i+minlen]
                    print 'Query:'+query
                if gene[1][i:i+minlen] == query:
                    match += 1
                print 'Match:' + str(match) +'\n'
        if match == len(genelist):
            print 'Match found! ::' + query
        else:
            minlen -= 1
            pos += 1
        print 'End of search, match:' + str(match) +'\n'
        match = 0
