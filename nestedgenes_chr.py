#setting parameters
cdslist = {}
nestedgene = []
#genome =[]
#open file and parse info into list of dictionaries (genome)
def parsegff(gff_file):
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

    return data

genome = parsegff('Coding_transcript_CDS.gff3')

#initialize a dictionary for all ID
#  cdslist = { chromosome: {'ID=CDS:XXXXXXX' : [ start , end ] } }
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

for chrom in cdslist:
    for ID1 in cdslist[chrom].iteritems():
        for ID2 in cdslist[chrom].iteritems():

            if ID1 != ID2:
                if ID1[1][0] > ID2[1][0] and ID1[1][1] < ID2[1][1]:
                    if ID1[0] in nestedgene:
                        pass
                    else:
                        nestedgene.append(ID1[0])


print 'Total number of nested genes: %s  ' % len(nestedgene)
print 'List of nested genes:        '
#output unique ID of nestedgene
for i in nestedgene:
    print i
print ''
