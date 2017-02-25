def parse(gff_file):
    """
    Parse the specified gff file and assign the values of each gene(line) to a dictionaries,
    return as a list of dictionaries
    """
    data = []
    gff_f = open(gff_file, 'r')
    gff = gff_f.read().splitlines()
    gff_f.close()
    print 'Parsing GFF..'
    for line in gff:
        gene_entry = {}
        if line.startswith('##'):
            pass
        else:
            field = line.split()
            gene_entry['chrom'] = field[0]
            gene_entry['source'] = field[1]
            gene_entry['id_type'] = field[2]
            gene_entry['start'] = int(field[3])
            gene_entry['end'] = int(field[4])
            gene_entry['score'] = field[5]
            gene_entry['strand'] = field[6]
            gene_entry['phase'] = field[7]
            gene_entry['attrib'] = field[8].split(';')
            if gene_entry['id_type'] == 'intron' or gene_entry['id_type'] == 'exon':
                gene_entry['gene_id'] = gene_entry['attrib'][0][18:]
            elif gene_entry['id_type'] == 'mRNA':
                gene_entry['gene_id'] = gene_entry['attrib'][0][14:]
        data.append(gene_entry)
    print '  Done!'
    return data


def sort(parsedgff):
    """
    Sorting intron information into dictionaries
    :param parsedgff: parsed gff file using parse()
    :return: sorted dictionary

    """
    gene_dict = {}
    for gene_entry in parsedgff:
        chrom = gene_entry['chrom']
        start = gene_entry['start']
        end = gene_entry['end']
        gene_id = gene_entry['gene_id']
        id_type = gene_entry['id_type']
        if id_type == 'mRNA':
            geneinfo = {id_type: [start, end], 'strand': gene_entry['strand']}
        else:
            geneinfo = {id_type: [[start, end]], 'strand': gene_entry['strand']}
        chrominfo = {gene_id: geneinfo}
        chrom_entry = {chrom: chrominfo}
        if chrom in gene_dict:
            if gene_id in gene_dict[chrom]:
                if id_type in gene_dict[chrom][gene_id]:
                    gene_dict[chrom][gene_id][id_type].append([start, end])
                else:
                    gene_dict[chrom][gene_id].update(geneinfo)
            else:
                gene_dict[chrom].update(chrominfo)
        else:
            gene_dict.update(chrom_entry)
    return gene_dict


def getdict(gff_file):
    output = parse(gff_file)
    return sort(output)


def filterintron(gene_dict):
    """
    remove entries without introns
    :param gene_dict: sorted gene_dict from gff
    :return: gene_dict of genes with introns
    """
    output = gene_dict.copy()
    for chrom in gene_dict:
        gene_id = gene_dict[chrom].keys()
        for gene_id in gene_id:
            if 'intron' not in gene_dict[chrom][gene_id]:
                del output[chrom][gene_id]
    return output


def getintrongene(gff_file):
    output = getdict(gff_file)
    output = filterintron(output)
    return output


def filterlong(gene_dict):
    """
    Filter and remove short isoforms from dictionary
    :param gene_dict:
    :return: filtered dictionary
    """
    alpha = 'abcdefghijklmnopqrstuvwxyz'
    intron_2 = gene_dict.copy()
    for chrom in intron_2:
        gene_id = intron_2[chrom].keys()
        filter_isoform = []
        for id in gene_id:
            #build list of id with isoforms
            if id.endswith(tuple(alpha)) or id.count('.') == 2:
                filter_isoform.append(id)
        filter_isoform.sort()
        isoforms = []
        long_isoform = ''
        for num, id in enumerate(filter_isoform):
            id_num = num + 1
            if id in isoforms or id is long_isoform or id_num == len(filter_isoform):
                pass
            else:
                isoforms = [id]
                parent = id
                if id.count('.') == 2:
                    parent = parent[:-2]
                if parent.endswith(tuple(alpha)):
                    parent = parent[:-1]
                for i in xrange(id_num, len(filter_isoform)):
                    query = filter_isoform[i]
                    if query.count('.') == 2:
                        query = query[:-2]
                    if query.endswith(tuple(alpha)):
                        query = query[:-1]
                    if parent == query:
                        isoforms.append(filter_isoform[i])
                    else:
                        break
                #process only when there is more than one isoform with introns
                if len(isoforms) > 1:
                    compare_iso = []
                    for iso1 in isoforms:
                        iso_start = intron_2[chrom][iso1]['mRNA'][0]
                        iso_end = intron_2[chrom][iso1]['mRNA'][1]
                        iso_length = iso_end - iso_start
                        compare_iso.append([iso1, iso_length])
                    #omit the longest isoform and remove the rest from the dictionary
                    long_isoform = max(compare_iso, key=lambda k: k[1])[0]
                    isoforms.remove(long_isoform)
                    for iso2 in isoforms:
                        del intron_2[chrom][iso2]
    return intron_2


#reorganize and compute intron information
def get_intron_table(intron_dict):
    intron_table = []
    for i in intron_dict:
        listed = intron_dict[i].keys()
        listed.sort()
        for j in listed:
            introns = intron_dict[i][j]['intron']
            mrna = intron_dict[i][j]['mRNA']
            strand = intron_dict[i][j]['strand']
            introns.sort()
            mrna.sort()
            if strand == '-':
                introns = introns[::-1]
                for k in xrange(len(introns)):
                    introns[k]= introns[k][::-1]
                mrna = mrna[::-1]
            for x in xrange(len(introns)):
                num = x+1
                intron_id = 'intron_%s_%s' % (num, j)
                if strand =='+':
                    intron_dist = introns[x][0] - mrna[0]
                    intron_size = introns[x][1] - introns[x][0]
                else:
                    intron_dist = mrna[0] - introns[x][0]
                    intron_size = introns[x][0] - introns[x][1]
                gene_id = j
                intron_table.append({'IntronID': intron_id, 'Distance': intron_dist, 'GeneID': gene_id, 'Size': intron_size})
    return intron_table



f = getintrongene('wormbase_all.gff3')
print 'Filtering..'
f = filterlong(f)
print '  Done!'
print 'Getting intron table..'
f = get_intron_table(f)
for n in range(10):
    print f[n]
