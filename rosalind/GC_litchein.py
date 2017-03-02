#open file and build list, split by > for FASTA
Seqfile = open('rosalind_gc.txt', 'r')
Seqlist = Seqfile.read().split('>')
#skip first empty >
Seqlist = Seqlist[1:]
genelist = []
#sorting genes into genelist by [ gene , seq ]
for gene in Seqlist:
    gene = gene.splitlines()
    gene= [ gene[0] , ''.join( gene[1:] ) ]
    genelist.append(gene)

#counting GC for each gene
for gene in xrange( len(genelist) ):
    #reset counter and count GC
    GC = 0
    for base in genelist[gene][1]:
        if base == 'G' or base == 'C':
                GC += 1
#calculate GC% and append to genelist, [ gene , seq , GCcontent ]
    GCcontent = float(GC) / len( genelist[gene][1] ) * 100
    genelist[gene].append(GCcontent)

#return max GC value from genelist
maxGC = max( genelist, key=lambda k:k[2])
print maxGC[0]
print maxGC[2]
