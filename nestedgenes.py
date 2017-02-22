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

genome = parsegff('sampleCDS.gff3')

#initialize a dictionary for all ID
#  { 'ID=CDS:XXXXXXX' : [ chrom , start , end ] }
print 'Initializing dictionaries..'
print '  Total entries: %s' % len(genome)
for entry in genome:
    cdslist[entry['attrib'][0][7:]] = [ entry['seqid'] , int(entry['start']) , int(entry['end'])]

#build a dictionary of start stop for each gene in genome
#compare start, end for each duplicate ID and compile
count = 0
print 'Building dictionaries..'
for entry in genome:
    count += 1
    print count,' processed..                                       \r',
    attrib = entry['attrib'][0][7:]
    cds={}
    cds[attrib] = [ entry['seqid'] , int(entry['start']) , int(entry['end']) ]

    for ID in cdslist.keys():
        if ID == attrib:
            if cds[attrib][1] < cdslist[ID][1]:
                cdslist[ID] = [ cdslist[ID][0] , cds[attrib][1] , cdslist[ID][2] ]
            if cds[attrib][2] > cdslist[ID][2]:
                cdslist[ID] = [ cdslist[ID][0] , cdslist[ID][1] , cds[attrib][2] ]
        print ID , cdslist[ID]
# search in cdslist for nested genes
del genome
print'  Total CDS: %s          ' % len(cdslist)
print 'Searching for nested genes..'
count = 0

for key,value in cdslist.iteritems():
    count += 1
    print count,' processed..                                       \r',
#Compares only from the same chromosome
    for x in cdslist.iteritems():
        if value[0] == x[1][0] and value[1] > x[1][1] and value[2] < x[1][2]:
            nestedgene.append(key)

print 'Completed!          \n\n\n'
print 'Total number of nested genes: %s  ' % len(set(nestedgene))
print 'List of nested genes:        '
#output unique ID of nestedgene
for i in set(nestedgene):
    print i
print ''
