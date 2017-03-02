# setting parameters for input
seqfile = open('rosalind.txt', 'r').read().splitlines()
s = str(seqfile[0])
t = str(seqfile[1])
pos = []
# compare seqs
for i in xrange(0, len(s)):
    if s[i:i+len(t)] == t:
        # append pos to list as strings
        pos.append(str(i+1))
# output positions and remove brackets
print ' '.join(pos)
