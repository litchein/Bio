def parsegff(gff_file):
    """
    Parse the specified gff file and assign the values of each gene(line) to a dictionaries,
    return as a list of dictionaries
    """
    data = []
    gff_f = open( gff_file , 'r')
    gff = gff_f.read().splitlines()
    gff_f.close()
    print 'Parsing GFF..'
    for line in gff:
        if line.startswith('##'):
            pass
        else:
            field = line.split()
            entry = {}
            entry['seqid'] = field[0]
            entry['source'] = field[1]
            entry['type'] = field[2]
            entry['start']  = int(field[3])
            entry['end'] = int(field[4])
            entry['score'] = field[5]
            entry['strand'] = field[6]
            entry['phase'] = field[7]
            entry['attrib'] = field[8].split(';')
            data.append(entry)

    print '  Done!'
    return data

def sortbychrom(genome):
    """
    Initialize a dictionary for all ID sorted by chromosomes and return as a dictionary of chromosomes,
    and each chromosome is a dictionary of genes with start and end position values

    cdslist = { chromosome: {'XXXXXXX' : [ start , end ] } }
    """
    cdslist = {}
    print 'Initializing dictionaries..'
    print '  Total entries: %s' % len(genome)
    for entry in genome:
        IDloc = {}
        chromidloc = {}

        chrom = entry['seqid']
        ID = entry['attrib'][0][18:]
        start = int(entry['start'])
        end = int(entry['end'])

        IDloc[ID] =[ [start , end] ]
        chromidloc[chrom] = IDloc
    #check if chromosome or ID exist in dictionary
        if chrom in cdslist:
            if ID in cdslist[chrom]:
                cdslist[chrom][ID].append( [start , end] )
            else:
                cdslist[chrom].update(IDloc)
        else:
            cdslist.update(chromidloc)
    print '  Done!'
    return cdslist

def genestart(cdslist,chrom,gene):
    """
    Compare and return the start position of target gene
    """
    return min( cdslist[chrom][gene] , key=lambda k:k[0])[0]

def geneend(cdslist,chrom,gene):
    """
    Compare and return the end position of target gene
    """
    return max( cdslist[chrom][gene] , key=lambda k:k[1])[1]

def searchnested(cdslist):
    """
    Search for nested genes (with exons) by each chromosome and return in a list of ID
    nestedgene = [ [nested] , [nest] ]
    """
    nestedgene = []
    nest =[]
    nested = []
    for chrom in cdslist:
        count = 0
        search1 = cdslist[chrom].keys()
        search2 = cdslist[chrom].keys()
        print 'Searching chromosome %s' % chrom
        print 'No. of CDS: %s' % len(search1)
        for ID1 in search1:
            ID1start = genestart(cdslist,chrom,ID1)
            ID1end = geneend(cdslist,chrom,ID1)
            ID1exon = len(cdslist[chrom][ID1])
            search2.remove(ID1)
            count += 1
            print 'Processing %s..          \r' % count,
            for ID2 in search2:
                ID2start = genestart(cdslist,chrom,ID2)
                ID2end = geneend(cdslist,chrom,ID2)
                ID2exon = len(cdslist[chrom][ID2])
                if ID1start > ID2start and ID1end < ID2end:
                    if ID1 not in nested and ID1exon > 1:
                        nested.append(ID1)
                        nest.append(ID2)

                elif ID2start > ID1start and ID2end < ID1end:
                    if ID2 not in nested and ID2exon > 1:
                        nested.append(ID2)
                        nest.append(ID1)
    print '  Done!'
    nestedgene = [ nested , nest ]
    return nestedgene

genome = parsegff('exon_transcript.gff3')
cdslist = sortbychrom(genome)
nestedgene = searchnested(cdslist)

print ''
print 'Total number of nested genes: %s  ' % len(nestedgene[0])
print 'List of nested genes:        '
#output ID of nestedgene
for i in xrange(len(nestedgene[0])):
    print '%s is nested in %s' % (nestedgene[0][i],nestedgene[1][i])
print ''
