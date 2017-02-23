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
            entry['chrom'] = field[0]
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

def sortintron(parsedgff):
    """
    Sorting intron information into dictionaries
    :param parsedgff: parsed gff file using parsegff()
    :return: sorted dictionary of introns

    """

intronfile = parsegff('sampleintron.gff3')
