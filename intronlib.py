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
            if gene_entry['id_type'] == 'mRNA':
                gene_entry['gene_id'] = gene_entry['attrib'][0][14:]
            if gene_entry['id_type'] == 'intron' or gene_entry['id_type'] == 'exon':
                gene_entry['gene_id'] = gene_entry['attrib'][0][18:]
            if gene_entry['id_type'] == 'CDS':
                gene_entry['gene_id'] = gene_entry['attrib'][1][18:]
                if gene_entry['gene_id'].count('.') > 1:
                    gene_names = gene_entry['attrib'][1].split('=')[1].split(',')
                    if len(gene_names) == 1:
                        gene_entry['gene_id'] = gene_names[0][11:]
                    else:
                        id_list = []
                        for n in gene_names:
                            id_list.append(n[11:])
                        gene_entry['gene_id'] = id_list
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
        strand = gene_entry['strand']
        if id_type == 'mRNA':
            geneinfo = {id_type: [start, end], 'strand': strand}
        else:
            geneinfo = {id_type: [[start, end]], 'strand': strand}

        if type(gene_id) is list:
            for i in gene_id:
                chrominfo = {i: geneinfo}
                chrom_entry = {chrom: chrominfo}
                if chrom in gene_dict:
                    if i in gene_dict[chrom]:
                        if id_type in gene_dict[chrom][i]:
                            gene_dict[chrom][i][id_type].append([start, end])
                        else:
                            gene_dict[chrom][i].update(geneinfo)
                    else:
                        gene_dict[chrom].update(chrominfo)
                else:
                    gene_dict.update(chrom_entry)
        else:
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
    output_dict = parse(gff_file)
    return sort(output_dict)


def filterintron(gene_dict):
    """
    remove entries without introns
    :param gene_dict: sorted gene_dict from gff
    :return: gene_dict of genes with introns
    """
    output_dict = gene_dict.copy()
    for chrom in gene_dict:
        gene_id = gene_dict[chrom].keys()
        for gene in gene_id:
            if 'intron' not in gene_dict[chrom][gene]:
                del output_dict[chrom][gene]
    return output_dict


def getintrongene(gff_file):
    output_dict = getdict(gff_file)
    output_dict = filterintron(output_dict)
    return output_dict


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
        for transcript in gene_id:
            # build list of id with isoforms
            if transcript.endswith(tuple(alpha)) or transcript.count('.') == 2:
                filter_isoform.append(transcript)
        filter_isoform.sort()
        isoforms = []
        long_isoform = ''
        for num, transcript in enumerate(filter_isoform):
            id_num = num + 1
            if transcript in isoforms or transcript is long_isoform or id_num == len(filter_isoform):
                pass
            else:
                isoforms = [transcript]
                parent = transcript
                if transcript.count('.') == 2:
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
                # process only when there is more than one isoform with introns
                if len(isoforms) > 1:
                    compare_iso = []
                    for iso1 in isoforms:
                        iso_start = intron_2[chrom][iso1]['mRNA'][0]
                        iso_end = intron_2[chrom][iso1]['mRNA'][1]
                        iso_length = iso_end - iso_start
                        compare_iso.append([iso1, iso_length])
                    # omit the longest isoform and remove the rest from the dictionary
                    long_isoform = max(compare_iso, key=lambda k: k[1])[0]
                    isoforms.remove(long_isoform)
                    for iso2 in isoforms:
                        del intron_2[chrom][iso2]
    return intron_2


# reorganize and compute intron information
def get_intron_table(intron_dict):

    output_intron_table = []
    for chrom in intron_dict:
        listed = intron_dict[chrom].keys()
        for gene_id in listed:
            introns = intron_dict[chrom][gene_id]['intron']
            mrna = intron_dict[chrom][gene_id]['mRNA']
            cds = intron_dict[chrom][gene_id]['CDS']
            cds.sort()
            strand = intron_dict[chrom][gene_id]['strand']
            introns.sort()
            mrna.sort()
            atg_start = min(min(cds))
            if strand == '-':
                introns = introns[::-1]
                for k in xrange(len(introns)):
                    introns[k] = introns[k][::-1]
                mrna = mrna[::-1]
                atg_start = max(max(cds))
            for i in xrange(len(introns)):
                num = i+1
                intron_id = '%s.i%s' % (gene_id, num)
                intron_start = introns[i][0]
                intron_end = introns[i][1]
                mrna_start = mrna[0]
                if strand == '+':
                    intron_dist_mrna = intron_start - mrna_start
                    intron_size = intron_end - intron_start + 1
                    intron_dist_atg = intron_start - atg_start
                else:
                    intron_dist_mrna = mrna_start - intron_start
                    intron_size = intron_start - intron_end + 1
                    intron_dist_atg = atg_start - intron_start
                output_intron_table.append({'IntronID': intron_id, 'Dist_mRNA': intron_dist_mrna,
                                            'Dist_ATG': intron_dist_atg, 'GeneID': gene_id,
                                            'Size': intron_size, 'IntronStart': intron_start,
                                            'IntronEnd': intron_end, 'Strand': strand, 'Chrom': chrom})
    return output_intron_table

# change gff3 file if required
all_genes = getdict('wormbase_all.gff3')
all_genes_with_introns = filterintron(all_genes)
print 'Getting intron table..'
intron_table = get_intron_table(all_genes_with_introns)

# writing to file
output = open('intron_table.txt', 'w')
output.write('## IntronID\tChrom\tIntronStart\tIntronEnd\tStrand\tIntronSize\tDist_mRNA\tDist_ATG\tGeneID\n')
for x in range(len(intron_table)):
    output.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (intron_table[x]['IntronID'], intron_table[x]['Chrom'],
                                                           intron_table[x]['IntronStart'], intron_table[x]['IntronEnd'],
                                                           intron_table[x]['Strand'], intron_table[x]['Size'],
                                                           intron_table[x]['Dist_mRNA'], intron_table[x]['Dist_ATG'],
                                                           intron_table[x]['GeneID']))

output.close()
print '  Done!'
