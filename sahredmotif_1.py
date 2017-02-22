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
maxlen = len( max(genelist , key=lambda k:k[1])[1] )
minlen = len( short[1] )
#remove shortest sequence from query genelist
for x in xrange( len(genelist) ):
    if genelist[x] == short:
        del genelist[x]
match = 0
print 'Minimum length: %s' % minlen

def searchmotif(query,template,match):
    for i in xrange(1+len(template[1])-len(query)):
        search = template[1][i:i+len(query)]
        print 'Query: %s' % query
        print 'Search: %s ID: %s' % (search, template[0])
        if search == query:
            print 'Match \n'
            match += 1
        else:
            print 'No Match! \n'
            
while minlen > 2:
    for combi in xrange( 1 + len(short[1]) - minlen ):
        query = short[1][combi:combi+minlen+1]
        print 'New query: %s \n' % query
        for gene in genelist:
            searchmotif(query,gene,0)
            if match >= 1:
                continue
            else:
                minlen -= 1
