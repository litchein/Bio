# Input for permutation calculation
n = 5 # input
list_n = range(1,n+1)
# Compute permutations
perm_n = 1
for i in list_n:
    perm_n *= i
print perm_n

# Generate permutation lists
perm_list = []
counter = 0
def perm(n, i=0):
    if i == len(n) - 1:
        print ' '.join(str(s) for s in n)
    else:
        for j in range(i, len(n)):
            n[i], n[j] = n[j], n[i]
            perm(n, i + 1)
            n[i], n[j] = n[j], n[i]
perm(list_n)