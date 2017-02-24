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
        ID = entry['attrib'][0][7:]
        start = int(entry['start'])
        end = int(entry['end'])

        IDloc[ID] = [ start , end ]
        chromidloc[chrom] = IDloc
    #check if chromosome or ID exist in dictionary
        if chrom in cdslist:
            if ID in cdslist[chrom]:
                if start < cdslist[chrom][ID][0]:
                    cdslist[chrom][ID][0] = start
                if end > cdslist[chrom][ID][1]:
                    cdslist[chrom][ID][1] = end
            else:
                cdslist[chrom].update(IDloc)
        else:
            cdslist.update(chromidloc)
    print '  Done!'
    return cdslist

def searchnested(cdslist):
    """
    Search for nested genes by each chromosome and return in a list of ID
    """
    nest =[]
    nested = []
    for chrom in cdslist:
        search1 = cdslist[chrom].items()
        search2 = cdslist[chrom].items()
        print 'Searching chromosome %s' % chrom
        print 'No. of CDS: %s' % len(search1)
        for ID1 in search1:
            search2.remove(ID1)
            for ID2 in search2:
                if ID1[1][0] > ID2[1][0] and ID1[1][1] < ID2[1][1]:
                    if ID1[0] not in nested:
                        nested.append(ID1[0])
                        nest.append(ID2[0])

                if ID2[1][0] > ID1[1][0] and ID2[1][1] < ID1[1][1]:
                    if ID2[0] not in nested:
                        nested.append(ID2[0])
                        nest.append(ID1[0])

    print 'Done!'
    nestedgene = [ nested , nest ]
    return nestedgene

genome = parsegff('Coding_transcript_CDS.gff3')
cdslist = sortbychrom(genome)
nestedgene = searchnested(cdslist)

print ''
print 'Total number of nested genes: %s  ' % len(nestedgene[0])
print 'List of nested genes:        '
#output ID of nestedgene
for i in xrange(len(nestedgene[0])):
    print '%s is nested in %s' % (nestedgene[0][i],nestedgene[1][i])
print ''
