import RNA

seq = 'GGGGGGGGGGAAAGGGGGGGGGG&CCCCCCCCCCUUUCCCCCCCCCC'

# create fold_compound data structure
fc = RNA.fold_compound(seq)

e = RNA.cofold(seq)
print(e)