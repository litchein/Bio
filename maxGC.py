from Bio import SeqIO
#build empty dictionary and parse FASTA format
dictGC={}
for gene in SeqIO.parse('rosalind_gc.txt','fasta'):
#reset GC values
    GC=0
#Calculate each GC values
    for base in gene.seq:
        if base == 'G':
            GC+=1
        if base == 'C':
            GC+=1
    GCcontent = float(GC) / len(gene) * 100
#store the GC values into dictionary
    dictGC[gene.id] = GCcontent
#Search for the gene with highest GC values
maxGC = max(dictGC.items(), key=lambda k: k[1])
#print out gene and GC value
print maxGC[0]
print maxGC[1]
