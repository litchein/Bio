# Input for calculation
k = 23   # AA
m = 25   # Aa
n = 22   # aa

# Calculate probabilities
total = float(k + m + n)
# Probability of having 1 dominant allele = 1 - probability of having 2 recessive
# 2 recessive = aa*aa + Aa*aa + Aa*Aa

aa_aa = n / total * (n-1) / (total-1)
Aa_aa = (m / total * n / (total-1)) * (0.5) * 2
Aa_Aa = m / total * (m-1) / (total-1) * (0.25)

print 1-(aa_aa + Aa_Aa + Aa_aa)

