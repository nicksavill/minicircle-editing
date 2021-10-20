
"""
Test if thermodynamic co-folding of a gRNA and all possible mRNA edits predicts
the correct gRNA based on minimum free energy.
The test is the initiator gRNA of COX3. All possible edits can be constructed from 
the anchor region onwards (389,168 unique possibilities).
The edits and gRNA are cofolded to obtain MFEs and secondary structures.
Edits can be discarded if there are too many mismatches to reduce the number of cofoldings
which are time consuming.
These correlate well with the base-pairing weighting as used by Josie and Paul in 
their dissertations.   
Conclusion:
Rather than use Josie's method of constructing editing blocks using arbitrary weightings
we instead construct all possible edits and select the most probable (eg those with 
probability ratios greater than 0.1 relative to the edit with the lowest MFE)
"""
import RNA
import numpy as np
import multiprocessing as mp
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

weights = {
    'GC': 30,
    'CG': 30,
    'AU': 20,
    'UA': 20,
    'GU': 15,
    'UG': 15
}
pairings = {'GC':'|', 'CG':'|', 'AU':'|', 'UA':'|', 'GU':':', 'UG':':'}

def get_pairing(mrna, grna):
    return ''.join([pairings[mb+gb] if mb+gb in pairings else '.' for mb, gb in zip(mrna, grna)])
    
def get_weight(mrna, grna):
    """
    Weights as used in Josie's and Paul's disertations
    """
    return -sum([weights.get(mb+gb, 0) for mb, gb in zip(mrna, grna)])

def fold(e_mrna, grna, anchor):
    """
    Cofold mRNA and gRNA.
    Anchor is added so that basepairing occurs correctly just after anchor
    Return edited mRNA, secondary structure and MFE
    """
    seq = e_mrna + anchor[0] + '&' + anchor[1][::-1] + grna[::-1]
    md = RNA.md()
    fc = RNA.fold_compound(seq, md)
    structure, MFE = fc.mfe_dimer()
    return e_mrna, structure, MFE

def next_grna_base(u_mrna, grna, e_mrna, mrna_idx, grna_idx, length, mismatches, e_mrna_set):
    """
    Add next edit to a mRNA-gRNA alignment
    mRNAs and gRNAs are aligned in direction 3' to 5' for convenience
    Disregard edits with too many mismatches
    """
    if get_pairing(e_mrna.upper(), grna).count('.') > mismatches:
        return
    if grna_idx == length:
        # return edited mRNA with 5' to 3' direction
        e_mrna_set.add(e_mrna[::-1])
        return
        
    mbase = u_mrna[mrna_idx]
    gbase = grna[grna_idx]
  
    # no edit
    next_grna_base(u_mrna, grna, e_mrna+mbase, mrna_idx+1, grna_idx+1, length, mismatches, e_mrna_set)
    
    # U-insertion
    if gbase in 'AG':
        next_grna_base(u_mrna, grna, e_mrna+'u', mrna_idx, grna_idx+1, length, mismatches, e_mrna_set)
        
    # U-deletion
    if mbase == 'U':
        next_grna_base(u_mrna, grna, e_mrna+u_mrna[mrna_idx+1], mrna_idx+2, grna_idx+1, length, mismatches, e_mrna_set)

# the alignment 5' to 3'
anchor = ['UAGUUUGUAG', 'AUCAAACAUC']
mrna = 'uuAuGuGuuAuGuAuuuGuGuGuGuAAuuuuAuuGGuGUUUuuuuAGUUGuuGAuuAGuuAAGuuuuAuuuGG'
grna = 'AUUAGUGAUUGAUCAAUUCGAGAUAGACU'

mismatches = 5
# length = 20
length = len(grna)
grna = grna[-length:]

e_mrna_set = set() # unique edited mRNAs
u_mrna = mrna.replace('u', '')

# create all edits with fewer than 6 mismatches
next_grna_base(u_mrna[::-1], grna[::-1], '', 0, 0, length, mismatches, e_mrna_set)
print(len(e_mrna_set))

edits = pd.DataFrame(list(e_mrna_set), columns=['e_mrna'])
edits['pairing'] = edits['e_mrna'].apply(lambda x: get_pairing(x.upper(), grna))
edits['mismatches'] = edits['pairing'].str.count('\.')
edits['weight'] = edits['e_mrna'].apply(lambda x: get_weight(x.upper(), grna))
edits = edits.set_index('e_mrna')

# parallel cofold of gRNA with mRNA to get MFE
with mp.Pool(mp.cpu_count()) as p:
    results = p.starmap(fold, [(e_mrna, grna, anchor) for e_mrna in e_mrna_set])

edits = edits.join(pd.Series([i[1] for i in results], index=edits.index, name='structure'))
edits = edits.join(pd.Series([i[2] for i in results], index=edits.index, name='mfe'))

# probability of edit relative to edit with minimum MFE
minmfe = edits['mfe'].min()
k = 1.986e-3  # boltzmann constant, kcal/(mol.K)
T = 310       # mammalian body temperature, K
edits['ratio'] = np.exp((minmfe-edits['mfe'])/(k*T))

edits = edits.sort_values('mfe')
edits.head(100).to_csv('edits.csv')

sns.scatterplot(x='mfe', y='weight', hue='mismatches', data=edits)
plt.show()
