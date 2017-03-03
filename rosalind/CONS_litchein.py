# Parse file into seq_list
input_file = open('rosalind_cons.txt').read().split('>')[1:]
seq_list = []
for i in input_file:
    seq = ''.join(i.splitlines()[1:])
    seq_list.append(seq)
seq_len = len(seq_list[0])

# Populate lists for counting
array = []
alist = ['A:']
clist = ['C:']
glist = ['G:']
tlist = ['T:']
for i in range(seq_len):
    alist.append(0)
    clist.append(0)
    glist.append(0)
    tlist.append(0)
array.append(alist)
array.append(clist)
array.append(glist)
array.append(tlist)

# Count bases in each sequence
for seq in seq_list:
    for n in xrange(seq_len):
        base = seq[n]
        pos = n + 1
        if base == 'A':
            array[0][pos] += 1
        elif base == 'C':
            array[1][pos] += 1
        elif base == 'G':
            array[2][pos] += 1
        elif base == 'T':
            array[3][pos] += 1

# Determine consensus sequence and print
cons = ''
for n in range(seq_len):
    pos = n + 1
    cons += max(array, key=lambda k: k[pos])[0][0]
print cons

# Display matrix
for each in array:
    print ' '.join(str(s) for s in each)
