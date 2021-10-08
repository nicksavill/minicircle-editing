"""
This file contains functions to test whether a gRNA and a mRNA match at the anchor region.
At the moment, there is just one methods, but different ones may be added in the future.
For example it could be interesting to allow a maximum number of GU pairs, or to use a minimum thermodynamic stability.
"""


# TODO: remove testing lines (currently commented out)
# TODO: anchor range assertions
# TODO: may rename position vars

# The most simple approach: Check for perfect matches (CW-base pairing) only, given a specific anchor size.
# gRNA bases must be upper and can be U or T.
#
# @author
# Josie Everatt (significantly adapted from original function by Paul Frohn)
#
# @param
# anchor_size: size of anchor
# mRNA: str of mRNA in 3' to 5' order
# gRNA: str of gRNA in 5' to 3' order
# position: binding position on mRNA (where is the first gRNA base?)
#
# @return
# True if it is a match, False else.
#
def perfect(anchor_min, anchor_max, mRNA, gRNA, position):
    assert isinstance(anchor_min,
                      int) and anchor_min > 1, "Please provide a valid minimum anchor size. Must be int and >1"
    assert isinstance(anchor_max,
                      int) and anchor_max > 1, "Please provide a valid maximum anchor size. Must be int and >1"
    assert isinstance(mRNA, str), "Please provide a valid mRNA String."
    assert isinstance(gRNA, str), "Please provide a valid gRNA String."
    assert isinstance(position, int), "Please provide a valid int binding-position."

    anchor_size = 0
    # anchor = ""
    # mRNAs = ""

    for current_base in range(anchor_max):
        #print("testing anchor size", anchor_size)
        if gRNA[current_base] == "A" and not (mRNA[position + current_base] == "U" or mRNA[position + current_base] == "u"):
            break
        if (gRNA[current_base] == "U" or gRNA[current_base] == "T") and mRNA[position + current_base] != "A":
            break
        if gRNA[current_base] == "G" and mRNA[position + current_base] != "C":
            break
        if gRNA[current_base] == "C" and mRNA[position + current_base] != "G":
            break
        else:
            # anchor = anchor + gRNA[current_base]
            # mRNAs = mRNAs + mRNA[position+current_base]
            anchor_size += 1

    if anchor_size >= anchor_min:
        return anchor_size
    else:
        return False
