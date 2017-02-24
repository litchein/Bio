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
wormbase = getdict('wormbase_all.gff3')
intron_genes = filterintron(wormbase)