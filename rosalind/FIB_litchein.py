n = 36   # months
k = 3   # litter size
# Only one pair if less than 3 months
if  n < 3:
    print 1
else:
    mating_adult = 1
    adult = 0
    baby = 0
    for i in xrange(n-2):
        if adult != 0:
            mating_adult += adult
            adult = 0

        if baby != 0:
            adult += baby
            baby = 0

        baby += mating_adult * k


print (mating_adult+adult+baby)

