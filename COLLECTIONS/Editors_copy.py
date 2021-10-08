'''
This file contains definitions of the 4 different editing approaches:
1. U-insertion: Always inserting a u if a mismatch can not be resolved.
2. Mismatch: Ignoring mismatches that can not be resolved.
3. Falloff: Stop editing at mismatches that can not be resolved.
4. SlidingThermo: Create the thermodynamically most stable mRNA(s)
'''

from math import exp
#from scipy.constants import k


# Edits the mRNA.
# Creates the all possible edit options with a sliding window approach,
# and returns the final top n edits:
# A window size is provided. The editor goes through the duplex base per base,
# constructing all possible edits within the window, with
# the window starting at the current position. E.g.:
# assume the alignment is this:
#    ____
# 3' ACGUGUACG  mRNA
# 5' UGCAAAA    gRNA
#
# where the anchor is the area with the line above.
# Then, the first window if windowsize == 3 will look at:
# 3' GUA    mRNA
# 5' AAA    gRNA
# And build a tree of possible edits. This tree will be:
#
# initial seq:                                GUA
# 1st position options:           GUA                     uGUA
# 2nd position options:     GUA         GAC         uGU         uuGUA
# 3rd position options:   GUA GUu     GAC GAu     uGU uGC     uuG   uuu
#
# From this tree, the highest scoring leaves are selected and used as initial sequence for the next iteration,
# then starting at the next position (2nd position in the tree)
#
# Editing stops when the optimal solution for the window includes more mismatches than allowed.
#
# @param
# mRNA: The mRNA to edit
# gRNA: The gRNA to edit with
# position: The position of the mRNA where to align the gRNA
# anchorsize: Size of the anchor region
# windowsize: Size of the window where to find the thermodynamical best duplex
#
# @return
# a dict of the possible edited mRNAs. Structure:
# {<mRNA-seq>: <position-after-editing>}
#
def new_slidingthermo_editing(mRNA, gRNA, position, anchorsize, windowsize, mismatches_currently, mismatches_allowed,
                              value_GC, value_AU, value_GU, value_mismatch):
    # TODO assertions
    # Necessary functions for slidingthermo_editing().

    # A function that tests which kind of match occurs.
    # returns 0 for mismatch, 1 for GC match, 2 for AU match, 3 for GU match
    def which_match(base1, base2):
        if (base1 == "G" and base2 == "C") or (base1 == "C" and base2 == "G"):
            return 1
        elif (base1 == "A" and (base2 == "T" or base2 == "u" or base2 == "U")) or \
                ((base1 == "T" or base1 == "U" or base1 == "u") and base2 == "A"):
            return 2
        elif ((base1 == "u" or base1 == "U" or base1 == "T") and base2 == "G") or \
                (base1 == "G" and (base2 == "u" or base2 == "U" or base2 == "T")):
            return 3
        else:
            return 0

    # Counts the number of mismatches between two sequences. seq1 and seq2 have to be of the same length
    # Returns int number of mismatches between two sequences
    def count_mismatches(seq1, seq2):
        assert len(seq1) == len(seq2), "Please provide sequences of the same length for comparison."

        mismatches = 0

        for current_base in range(len(seq1)):
            if which_match(seq1[current_base], seq2[current_base]) == 0:
                mismatches += 1

        return mismatches

    # Before filtering and returning edits, normalise their probabilities.
    # dict_of_all_edits: new_mRNAs
    # Structure: {mRNA-seq.upper():
    #             [actual-mRNA-seq, nr-of-mismatches-currently, position-after-editing, thermodynamic-probability]}
    # Returns list of highest Ps in descending order # TODO: check has normalised dict - should do so in place w/o needing to return it
    # Would potentially be cleaner as two separate functions but performance is an issue so wanted to avoid another for
    def normalise_probabilities(new_mRNAs):
        sum_probabilities = 0
        for mRNA_info in new_mRNAs.values():
            sum_probabilities += mRNA_info[3]

        top_probabilities = []

        # normalise and replace the old probability
        for mRNA_info in new_mRNAs.values():
            mRNA_info[3] = mRNA_info[3] / sum_probabilities
            top_probabilities.append(mRNA_info[3])

        # get descending order of Ps
        top_probabilities = sorted(top_probabilities, reverse=True)

        return top_probabilities

    # The function builds a tree to calculate all possible best edits.
    # Node class is needed to build the tree.
    class Node():
        def __init__(self, position, mRNA):

            # Information about this node that is given
            self.mRNA = mRNA
            self.position = position

            # Information about this node that is calculated
            self.value = None
            self.mismatches = 0
            self.created_by_deletion = False

            # The three different kinds of children
            self.stay = None
            self.insert = None
            self.delete = None
            # All children
            self.children = []

            self.leaves = []
            self.probability = None

            # TODO that is pretty ugly. I should either do it without the children list or only with the children list.

        # Add a child without editing at current position
        def add_stay(self):
            self.stay = Node(self.position + 1, self.mRNA)
            # If this position was created by deletion, pass this on to stay child to prevent re-insertion
            # TODO: should created_by_deletion be passed on? Surely stay is moving to the next position so it
            # doesn't matter?
            if self.created_by_deletion: #and (self.mRNA[self.position + 1] == "U" or self.mRNA[self.position + 1] == "u"):
                self.stay.created_by_deletion = True
            self.children.append(self.stay)

        # Add a child by inserting u at current position
        def add_insert(self):
            self.insert = Node(self.position + 1, self.mRNA[:self.position] + "u" + self.mRNA[self.position:])
            self.children.append(self.insert)

        # Add a child by deleting u at current position
        def add_delete(self):
            self.delete = Node(self.position, self.mRNA[:self.position] + self.mRNA[self.position + 1:])
            self.delete.created_by_deletion = True
            self.children.append(self.delete)

        # Building the tree with self as root
        # Recursive
        # Performs every edit regardless of complementarity to the gRNA, then scores the results based on the grNA
        def build_tree(self, depth):

            # Check if editing is wanted and possible
            if depth < windowsize and len(self.mRNA) > self.position:

                # Staying is always an option.
                self.add_stay()
                self.stay.build_tree(depth + 1)

                # If current position is not an U and if a U was not deleted from this position before,
                # inserting an u is an option.
                # At first I was wondering about this but I guess it makes sense - if a U is needed (i.e. the gRNA
                # is G or A), there is a U there already. If not, it can be deleted. Inserting one would be a waste.
                if not (self.mRNA[self.position] == "U" or self.mRNA[self.position] == "u") and not \
                        self.created_by_deletion:
                    self.add_insert()
                    self.insert.build_tree(depth + 1)

                # If current position is an U, deletion is an option.
                if self.mRNA[self.position] == "U" or self.mRNA[self.position] == "u":
                    self.add_delete()
                    # as deletion brings a new base to the current position, don't move to next position.
                    self.delete.build_tree(depth)

        # append all leaves to a given list.
        def get_leaves(self, list_of_leaves):
            # could we check by doing if children? i.e. list not empty?
            at_least_one_child = False

            if self.delete is not None:
                self.delete.get_leaves(list_of_leaves)
                at_least_one_child = True

            if self.stay is not None:
                self.stay.get_leaves(list_of_leaves)
                at_least_one_child = True

            if self.insert is not None:
                self.insert.get_leaves(list_of_leaves)
                at_least_one_child = True

            # It's a leaf!
            # Recursion base case
            if not at_least_one_child:
                list_of_leaves.append(self)

        # Normalises thermodynamic score using the Boltzmann constant
        def boltzmann_normalise_weights(self, root_probability):

            # 37C in Kelvin
            T = 310.15
            R = 8.314
            # Partition function - used to normalise for Boltzmann distribution
            Z = 0

            # first need to calculate value of each leaf
            for leaf in self.leaves:
                leaf.calculate_value(gRNA, value_GC, value_AU, value_GU, value_mismatch)

            # Z = sum over all pairs(e^(-G/kBT))
            for leaf in self.leaves:
                Z += exp(-((-leaf.value/100) / R * T))

            # Calculate probability using free energy and the Boltzmann disn
            for leaf in self.leaves:
                leaf.probability = (exp(-((-leaf.value/100) / R * T)) / Z) * root_probability

        # Normalises thermodynamic score by dividing the leaf's value by the total of all leaf values
        # Multiplies probability by that of the parent edit to allow tracking probability along the window
        def normalise_weights(self, root_probability):

            # first need to calculate value of each leaf
            for leaf in self.leaves:
                leaf.calculate_value(gRNA, value_GC, value_AU, value_GU, value_mismatch)

            total_thermo_score = 0
            for leaf in self.leaves:
                total_thermo_score += leaf.value
            if total_thermo_score == 0:
                total_thermo_score = 0.1
            for leaf in self.leaves:
                leaf.probability = (leaf.value / total_thermo_score) * root_probability
                # print("probability calculated", leaf.value, leaf.probability, root_probability)

        # calculate the value of the node and update the number of mismatches.
        # Only valid for leaves!
        def calculate_value(self, gRNA, value_GC, value_AU, value_GU, value_mismatch):

            value = 0

            # go through all bases and sum up their scores
            for current_position in range(windowsize):

                # what kind of base-pairing is present?
                # 0 for mismatch, 1 for GC, 2 for AU, 3 for GU
                basepair = which_match(self.mRNA[self.position - windowsize + current_position],
                                       gRNA[self.position - windowsize + current_position])

                if basepair == 1:
                    value += value_GC
                elif basepair == 2:
                    value += value_AU
                elif basepair == 3:
                    value += value_GU
                else:
                    value += value_mismatch
                    self.mismatches += 1

            # set node value
            self.value = value

    # Constructing the thermodynamic best sequences in the given window.
    #
    # @param
    # mRNA: string, the mRNA beginning at the binding position of the gRNA
    # gRNA: string, full gRNA
    # current_position: int, start position of editing window
    #
    def get_all_edits(mRNA_editingregion, gRNA, current_position, current_mismatches, root_probability=1):

        # Build the tree of possible edits for the current window
        # print(mRNA_editingregion)
        root = Node(current_position, mRNA_editingregion)
        root.probability = root_probability
        root.mismatches = current_mismatches
        root.build_tree(0)

        # Get the leaves = possible results
        leaves = []
        root.get_leaves(leaves)
        root.leaves = leaves
        root.boltzmann_normalise_weights(root.probability)

        # Calculate the number of mismatches in the given mRNA before the window.
        mismatches_before = count_mismatches(gRNA[anchorsize: current_position],
                                             mRNA_editingregion[anchorsize: current_position])

        # A dict to store the results
        edited_mRNAs = {}

        # Get a list of all values where the number of mismatches is within the limit
        for leaf in leaves:
            if mismatches_before + leaf.mismatches <= mismatches_allowed:
                edited_mRNAs[leaf.mRNA] = [leaf.mismatches, leaf.probability]

        # return dict of best mRNAs with number of mismatches and thermodynamic [probability]
        # or emtpy dict if not possible with given maximum mismatches.
        return edited_mRNAs

    # == End of definitions for slidingthermo_editing(). Now the actual program starts.

    # Current position to investigate ON THE gRNA.
    current_position = anchorsize

    # The generated mRNAs with the number of mismatches.
    # Keys are the uppercase version of the edit - this allows rapid searching to see if the edit has been produced
    # before or if it's new.
    # structure: {mRNA-seq.upper():
    #               [actual-mRNA-seq, nr-of-mismatches-currently, position-after-editing, thermodynamic-probability]}
    new_mRNAs = {mRNA.upper(): [mRNA, mismatches_currently, position + anchorsize, 1]}  # Initial P = 1 until assigned
    # print("starting with", new_mRNAs)
    # 1st loop: Going through all positions until either gRNA or mRNA ends or too many mismatches. Here: gRNAs
    while current_position + windowsize - 1 < len(gRNA):

        if len(new_mRNAs) == 0:
            break
        # Store the newly-generated mRNAs. Same structure as new_mRNAs.
        next_mRNAs = {}

        # print("new round. Investigating", len(new_mRNAs))
        # count = 1

        # Each round, iterate through all generated sequences to produce all new possible edits.
        # Reminder: values structure: [mRNA-seq, current-mismatches, position-after-editing, thermodynamic-probability]
        for uppercase_mRNA, values in new_mRNAs.items():
            # print(count)
            current_mRNA = values[0]

            # Go through all positions until gRNA or mRNA ends, or too many mismatches.
            # If too many/wrong position: moves to next mRNA in for loop
            if position + current_position + windowsize - 1 < len(current_mRNA) and \
                    new_mRNAs[uppercase_mRNA][1] <= mismatches_allowed:

                # Get all possible edits
                # Structure {edited-mRNA-sequence : [mismatches, thermodynamic-probability]}
                # get_all_edits params: mRNA_editingregion, gRNA, current_position, current_mismatches,
                # root_probability
                edits = get_all_edits(current_mRNA[position:], gRNA, current_position, new_mRNAs[uppercase_mRNA][1],
                                      new_mRNAs[uppercase_mRNA][3])
                # print("Returned", len(edits), "edits at position", current_position)

                # If the edit exists, increase its probability.
                if edits != {}:
                    for new_mRNA, details in edits.items():
                        # next_mRNAs keys are upper case edits because U and u are biologically equivalent
                        # If the edit has already been generated, increase its probability
                        if (uppercase_mRNA[:position] + new_mRNA.upper()) in next_mRNAs:
                            # print("already there")
                            # print(details)
                            # increase probability
                            # print("probability before", next_mRNAs[uppercase_mRNA[:position] + new_mRNA.upper()][3])
                            # print("probability to add", details[1])
                            next_mRNAs[uppercase_mRNA[:position] + new_mRNA.upper()][3] += details[1]
                            # print("probability after", next_mRNAs[uppercase_mRNA[:position] + new_mRNA.upper()][3])

                        else:
                            # A new mRNA was generated. Hooray! Now store it with the correct structure:
                            # {mRNA-seq.upper():
                            # [mRNA-seq, nr-mismatches-currently, position-after-editing, thermodynamic-probability]}
                            next_mRNAs[uppercase_mRNA[:position] + new_mRNA.upper()] = \
                                [current_mRNA[:position] + new_mRNA, details[0],
                                 position + current_position + windowsize, details[1]]
                            # print("adding new one. Next now length", len(next_mRNAs))
               # else:
                    # print("no edits made")
                    # May still be edited further downstream so add to list to be checked in future
                #    next_mRNAs[uppercase_mRNA] = new_mRNAs[uppercase_mRNA]
                    # Make sure not passed along in duplicate by exceeding mismatches?
                 #   new_mRNAs[uppercase_mRNA][1] = mismatches_allowed + 1  # TODO: WHY?
            # count += 1

        # Set the just generated mRNAs as initial mRNAs for the next iteration and go to next position.
        # Reset new mRNAs. So new mrnas are the newly produced best edits from each round.
        # They are initially stored in next mrnas when made. Once duplicates are accounted for, next becomes new
        # i.e. newly created and ready to be edited again.
        # print("updating for next round. Will investigate", len(next_mRNAs))
        # print("new", new_mRNAs)
        # print("next", next_mRNAs)
        new_mRNAs = next_mRNAs
        # print(next_mRNAs)
        current_position += 1

    # Once out of while loop, have edited and generated all possible sequences. Normalise their probabilities and
    # proceed with the <user-selected> most likely sequences.
    top_probabilities = normalise_probabilities(new_mRNAs)

    # Final sequences to return - don't need uppercase anymore, only interested in actual edit, position, probability
    # Structure: {mRNA-seq: [position-after-editing, thermodynamic-probability]}
    # print("ended editing with", len(new_mRNAs), "final edits")
    final_mRNAs = {}

    # get top 3 probabilities
    top_probabilities = top_probabilities[:3]

    for mRNA_info in new_mRNAs.values():
        if mRNA_info[3] in top_probabilities:
            final_mRNAs[mRNA_info[0]] = mRNA_info[2:]

    # Return all generated sequences.
    return final_mRNAs


def middle_slidingthermo_editing(mRNA, gRNA, position, anchorsize, windowsize, mismatches_currently, mismatches_allowed,
                                 value_GC, value_AU, value_GU, value_mismatch):
    # TODO assertions

    # Necessary functions for slidingthermo_editing().

    # Necessary functions for get_best_edits().
    # A function that tests which kind of match occurs.
    # returns 0 for mismatch, 1 for GC match, 2 for AU match, 3 for GU match
    def which_match(base1, base2):
        if (base1 == "G" and base2 == "C") or (base1 == "C" and base2 == "G"):
            return 1
        elif (base1 == "A" and (base2 == "T" or (base2 == "u" or base2 == "U"))) or (
                (base1 == "T" or (base1 == "U" or base1 == "u")) and base2 == "A"):
            return 2
        elif ((base1 == "u" or base1 == "U" or base1 == "T") and base2 == "G") or (
                base1 == "G" and (base2 == "u" or base2 == "U" or base2 == "T")):
            return 3
        else:
            return 0

    # Counting the number of mismatches between two sequences. nuc1 and nuc2 have to be of the same length
    def count_mismatches(nuc1, nuc2):
        assert len(nuc1) == len(nuc2), "Please provide sequences of the same length for comparison."

        mismatches = 0

        for current_base in range(len(nuc1)):
            if which_match(nuc1[current_base], nuc2[current_base]) == 0:
                mismatches += 1

        return mismatches

    # The function builds a tree to calculate all possible best edits.
    # A node-class is needed to build the tree.
    # Inherits from Node class in BF_Networks.py
    class Node():
        def __init__(self, position, mRNA):

            # Information about this node that is given
            self.mRNA = mRNA
            self.position = position

            # Information about this node that is calculated
            self.value = None
            self.mismatches = 0
            self.created_by_deletion = False

            # The three different kinds of children
            self.stay = None
            self.insert = None
            self.delete = None
            # All children
            self.children = []

            self.leaves = []
            self.probability = None

            # TODO that is pretty ugly. I should either do it without the children list or only with the children list.

        # Add a child without editing at current position
        def add_stay(self):
            self.stay = Node(self.position + 1, self.mRNA)
            # If this position was created by deletion, pass this on to stay child to prevent re-insertion
            # TODO: should created_by_deletion be passed on? Surely stay is moving to the next position so it
            # doesn't matter?
            if self.created_by_deletion and (self.mRNA[self.position+1] == "U" or self.mRNA[self.position+1] == "u"):
                self.stay.created_by_deletion = True
            self.children.append(self.stay)

        # Add a child by inserting u at current position
        def add_insert(self):
            self.insert = Node(self.position + 1, self.mRNA[:self.position] + "u" + self.mRNA[self.position:])
            self.children.append(self.insert)

        # Add a child by deleting u at current position
        def add_delete(self):
            self.delete = Node(self.position, self.mRNA[:self.position] + self.mRNA[self.position + 1:])
            self.delete.created_by_deletion = True
            self.children.append(self.delete)

        # Building the tree with self as root
        # Recursive
        # Performs every edit regardless of complementarity to the gRNA, then scores the results based on the grNA
        def build_tree(self, depth):

            # Check if editing is wanted and possible
            if depth < windowsize and len(self.mRNA) > self.position:

                # Staying is always an option.
                self.add_stay()
                self.stay.build_tree(depth + 1)

                # If current position is not an U and if a U was not deleted from this position before,
                # inserting an u is an option.
                # At first I was wondering about this but I guess it makes sense - if a U is needed (i.e. the gRNA
                # is G or A), there is a U there already. If not, it can be deleted. Inserting one would be a waste.
                if not (self.mRNA[self.position] == "U" or self.mRNA[self.position] == "u") and not \
                        self.created_by_deletion:
                    self.add_insert()
                    self.insert.build_tree(depth + 1)

                # If current position is an U, deletion is an option.
                if self.mRNA[self.position] == "U" or self.mRNA[self.position] == "u":
                    self.add_delete()
                    # as deletion brings a new base to the current position, don't move to next position.
                    self.delete.build_tree(depth)

        # append all leaves to a given list.
        def get_leaves(self, list_of_leaves):
            # could we check by doing if children? i.e. list not empty?
            at_least_one_child = False

            if self.delete is not None:
                self.delete.get_leaves(list_of_leaves)
                at_least_one_child = True

            if self.stay is not None:
                self.stay.get_leaves(list_of_leaves)
                at_least_one_child = True

            if self.insert is not None:
                self.insert.get_leaves(list_of_leaves)
                at_least_one_child = True

            # It's a leaf!
            # Recursion base case
            if not at_least_one_child:
                list_of_leaves.append(self)

        # Normalises thermodynamic score using the Boltzmann constant
        def boltzmann_normalise_weights(self, root_probability):

            # 37C in Kelvin
            T = 310.15
            R = 8.314
            # Partition function - used to normalise for Boltzmann distribution
            Z = 0

            # first need to calculate value of each leaf
            for leaf in self.leaves:
                leaf.calculate_value(gRNA, value_GC, value_AU, value_GU, value_mismatch)

            # Z = sum over all pairs(e^(-G/kBT))
            for leaf in self.leaves:
                Z += exp(-((-leaf.value / 100) / R * T))

            # Calculate probability using free energy and the Boltzmann disn
            for leaf in self.leaves:
                leaf.probability = (exp(-((-leaf.value / 100) / R * T)) / Z) * root_probability

        # Normalises thermodynamic score by dividing the leaf's value by the total of all leaf values
        # Multiplies probability by that of the parent edit to allow tracking probability along the window
        def normalise_weights(self, root_probability):

            # first need to calculate value of each leaf
            for leaf in self.leaves:
                leaf.calculate_value(gRNA, value_GC, value_AU, value_GU, value_mismatch)

            total_thermo_score = 0
            for leaf in self.leaves:
                total_thermo_score += leaf.value
            if total_thermo_score == 0:
                total_thermo_score = 0.1
            for leaf in self.leaves:
                leaf.probability = (leaf.value / total_thermo_score) * root_probability

        # calculate the value of the node and update the number of mismatches.
        # Only valid for leaves!
        def calculate_value(self, gRNA, value_GC, value_AU, value_GU, value_mismatch):

            value = 0

            # go through all bases and sum up their scores
            for current_position in range(windowsize):

                # what kind of base-pairing is present?
                # 0 for mismatch, 1 for GC, 2 for AU, 3 for GU
                basepair = which_match(self.mRNA[self.position - windowsize + current_position],
                                       gRNA[self.position - windowsize + current_position])

                if basepair == 1:
                    value += value_GC
                elif basepair == 2:
                    value += value_AU
                elif basepair == 3:
                    value += value_GU
                else:
                    value += value_mismatch
                    self.mismatches += 1

            # set node value
            self.value = value

    # Constructing the thermodynamical best sequences in the given window.
    #
    # @param
    # mRNA: the mRNA beginning at the binding position of the gRNA
    # gRNA: full gRNA
    #
    def get_best_edits(mRNA_editingregion, gRNA, current_position, current_mismatches, root_probability=1):

        # === End of definitions for get_best_edits(). Now the actual get_best_edits function starts.

        # Build the tree of possible edits
        root = Node(current_position, mRNA_editingregion)
        root.probability = root_probability
        root.mismatches = current_mismatches
        root.build_tree(0)

        # Get the leaves = possible results
        leaves = []
        root.get_leaves(leaves)
        root.leaves = leaves
        root.boltzmann_normalise_weights(root.probability)

        # Calculate the number of mismatches in the given mRNA before the window.
        mismatches_before = count_mismatches(gRNA[anchorsize: current_position],
                                             mRNA_editingregion[anchorsize: current_position])

        # Get a list of all values where the number of mismatches is within the limit
        leaf_probabilities = [leaf.probability for leaf in leaves
                              if mismatches_before + leaf.mismatches <= mismatches_allowed]

        # A dict to store the results
        best_mRNAs = {}

        # Check if editing can occur with the given limit of mismatches.
        # If yes, fill dict with the best mRNAs: {<new_mRNA-subsequence>: [<nr-of-mismatches>, <thermodynamic score>]}
        # Do we need to check the inner if statement? If it got this far surely the mms etc. have been checked already?
        # ALT when check mms above, use that to at the same time remove leaves with too many?
        # because this is looping thru leaves, not leaf_values. So while leaf_values are those with < max mms, leaves
        # contain all leaves still, even if > max mms
        if leaf_probabilities:
            max_probability = max(leaf_probabilities)
            for leaf in leaves:
                if leaf.probability == max_probability and leaf.mismatches <= mismatches_allowed:
                    best_mRNAs[leaf.mRNA] = [leaf.mismatches, leaf.probability]

        # return dict of best mRNAs with number of mismatches and thermodynamic score
        # or emtpy dict if not possible with given maximum mismatches.
        return best_mRNAs

    # == End of definitions for slidingthermo_editing(). Now the actual program starts.

    # Current position to investigate ON THE gRNA.
    current_position = anchorsize
    # The generated mRNAs with the number of mismatches.
    # structure: {uppercase-mRNA-seq: [actual-mRNA-seq, nr-of-mismatches-currently, position-after-editing-if-stopped-here, thermodynamic score]}
    new_mRNAs = {
        mRNA.upper(): [mRNA, mismatches_currently, position + anchorsize, 1]}  # Initial probability 0 until assigned

    # 1st loop: Going through all positions until either gRNA or mRNA ends or too many mismatches. Here: gRNAs
    while current_position + windowsize - 1 < len(gRNA):

        # if new_mRNAs == {}:
        #    break
        # A dict to store the new generated mRNAs. Same structure as new_mRNAs.
        next_mRNAs = {}

        # If there are multiple sequences with the same score, I have to iterate through them.
        for uppercase_mRNA, values in new_mRNAs.items():
            current_mRNA = values[0]

            # 2nd loop: Going through all positions until either gRNA or mRNA ends or too many mismatches. Here: mRNA
            if position + current_position + windowsize - 1 < len(current_mRNA) and \
                    new_mRNAs[uppercase_mRNA][1] <= mismatches_allowed:

                # Getting the best possible edits at this position (most stable duplex within the window)
                # Structure {mRNA sequence : [mismatches, thermodynamic_probability]}
                # Returned mRNA sequence is the mRNA from the editing block onwards - not the entire mRNA strand.
                best_edits = get_best_edits(current_mRNA[position:], gRNA, current_position,
                                            new_mRNAs[uppercase_mRNA][1],
                                            new_mRNAs[uppercase_mRNA][3])

                # If editing is possible, test if the result already exists and store it if not.
                if best_edits != {}:
                    for new_mRNA, details in best_edits.items():
                        if uppercase_mRNA[:position] + new_mRNA.upper() in next_mRNAs:
                            next_mRNAs[uppercase_mRNA[:position] + new_mRNA.upper()][3] += details[1]
                        else:
                            # A new mRNA was generated. Hooray! Now store it with the correct structure:
                            # {uppercase-mRNA-seq: [actual-mRNA-seq, nr-mismatches-currently, position-after-editing, thermodynamic-probability]}
                            next_mRNAs[uppercase_mRNA[:position] + new_mRNA.upper()] = \
                                [current_mRNA[:position] + new_mRNA, details[0],
                                 position + current_position + windowsize, details[1]]
                # else:
                #    next_mRNAs[current_mRNA] = new_mRNAs[current_mRNA]
                #    new_mRNAs[current_mRNA][0] = mismatches_allowed + 1  # TODO: WHY?

        # Set the just generated mRNAs as initial mRNAs for the next iteration and go to next position.
        # Reset new mRNAs. So new mrnas are the newly produced best edits from each round.
        # They are initially stored in next mrnas when made. Once duplicates are accounted for, next becomes new
        # i.e. newly created and ready to be edited again.
        # print("updating for next round. Will investigate", len(next_mRNAs))
        new_mRNAs = next_mRNAs
        current_position += 1

    # Delete double mRNAs
    all_keys = list(new_mRNAs.keys())
    for key1 in all_keys:
        occurrences = 0
        for key2 in new_mRNAs.keys():
            if key1.upper() == key2.upper():
                occurrences += 1
        if occurrences > 1:
            new_mRNAs.pop(key1)

    final_mRNAs = {}

    for mRNA_info in new_mRNAs.values():
        final_mRNAs[mRNA_info[0]] = mRNA_info[2:]

    # Return all generated sequences.
    return final_mRNAs


# TODO: decide if want to remove separate in/del function. Faster without them (see fiddle.py)

# If mismatch (every mismatch!) and no U on mRNA, a U is inserted.
# Allowes GU pairing when editing
# Stops when either mRNA or gRNA ends.
# Anchor is not edited.
# Accepts only upper letters except u in mRNA (for speed reasons)
#
# @author
# Created by Paul, modified by Josie
#
# @param
# mRNA: str of mRNA in 3' to 5' order!
# gRNA: str of gRNA in 5' to 3' order!
# position: the position on the mRNA where the gRNA did bind (1st base).
# anchorsize: anchorsize
#
# @return
# (sequence of the edited mRNA, position on mRNA after editing)
def overhang_uinsert(mRNA, gRNA, position, anchorsize):
    assert isinstance(anchorsize, int) and anchorsize > 1, "Please provide a valid int anchorsize >1."
    assert isinstance(mRNA, str), "Please provide a valid mRNA String."
    assert isinstance(gRNA, str), "Please provide a valid gRNA String."
    assert isinstance(position, int), "Please provide a valid int binding-position."

    # Anchor is not edited, so skip it.
    current_position = anchorsize
    edit_position = current_position + position

    # Iterating through the gRNA until either gRNA or mRNA is finished.
    while current_position < len(gRNA) and edit_position < len(mRNA):

        # case A at gRNA mismatch: Insert u
        if (gRNA[current_position] == "A") and (
                mRNA[edit_position] != "U" or mRNA[edit_position] != "u"):
            mRNA = insertion(mRNA, edit_position)
            current_position += 1

        # case G at gRNA mismatch: Insert u
        elif (gRNA[current_position] == "G") and not (
                mRNA[edit_position] == "U" or mRNA[edit_position] == "u" or mRNA[
            edit_position] == "C"):
            mRNA = insertion(mRNA, edit_position)
            current_position += 1

        # case U at mRNA mismatch. Delete U
        elif (mRNA[edit_position] == "U" or mRNA[edit_position] == "u") and not (
                gRNA[current_position] == "A" or gRNA[current_position] == "G"):
            mRNA = deletion(mRNA, edit_position)

        # case unresolvable mismatch with C at gRNA: Insert u
        elif (gRNA[current_position] == "C" and not (
                mRNA[edit_position] == "U" or mRNA[edit_position] == "u" or mRNA[
            edit_position] == "G")):
            mRNA = insertion(mRNA, edit_position)
            current_position += 1

        # case unresolvable mismatch with U/T at gRNA: Insert u
        elif ((gRNA[current_position] == "U" or gRNA[current_position] == "T") and not (
                mRNA[edit_position] == "A" or mRNA[edit_position] == "G")):
            mRNA = insertion(mRNA, edit_position)
            current_position += 1

        # case match
        else:
            current_position += 1

    # Returning the edited mRNA and position after editing when done.
    # (usually position + len(gRNA), but not if mRNA-end is reached first)
    return (mRNA, edit_position)


def insertion(mRNA, edit_position):
    mRNA = mRNA[:edit_position] + "u" + mRNA[edit_position:]
    return mRNA


def deletion(mRNA, edit_position):
    mRNA = mRNA[:edit_position] + mRNA[edit_position + 1:]
    return mRNA


# Edits the mRNA. If resolvable mismatch and no U on mRNA, a U is inserted.
# Mismatches that can't be resolved via U indel editing are ignored.
# Allowes Gu pairing when editing.
# Stops when either mRNA or gRNA ends.
# Anchor is not edited.
# Accepts only upper letters except u in mRNA (for speed reasons)
#
# @param
# mRNA: str of mRNA in 3' to 5' order!
# gRNA: str of gRNA in 5' to 3' order!
# position: the position on the mRNA where the gRNA did bind (1st base).
# anchorsize: anchorsize
#
# @return
# (sequence of the edited mRNA, position on mRNA after editing)


def overhang_ignore_mismatch(mRNA, gRNA, position, anchorsize):
    assert isinstance(anchorsize, int) and anchorsize > 1, "Please provide a valid int anchorsize >1."
    assert isinstance(mRNA, str), "Please provide a valid mRNA String."
    assert isinstance(gRNA, str), "Please provide a valid gRNA String."
    assert isinstance(position, int), "Please provide a valid int binding-position."

    # Anchor is not edited, so skip it.
    current_position = anchorsize
    edit_position = current_position + position

    # Iterating through the gRNA until either gRNA or mRNA is finished.
    while current_position < len(gRNA) and edit_position < len(mRNA):

        # case A at gRNA mismatch: Insert u
        if (gRNA[current_position] == "A") and (mRNA[edit_position] != "U" or mRNA[edit_position] != "u"):
            mRNA = insertion(mRNA, edit_position)
            current_position += 1

        # case G at gRNA mismatch: Insert u
        elif (gRNA[current_position] == "G") and not (mRNA[edit_position] == "U" or mRNA[edit_position] == "u" or
                                                      mRNA[edit_position] == "C"):
            mRNA = insertion(mRNA, edit_position)
            current_position += 1

        # case U at mRNA mismatch. Delete U
        elif (mRNA[edit_position] == "U" or mRNA[edit_position] == "u") and not (gRNA[current_position] == "A" or
                                                                                 gRNA[current_position] == "G"):
            mRNA = deletion(mRNA, edit_position)

        # case unresolvable mismatch with C at gRNA: Ignore it and go to next position
        elif (gRNA[current_position] == "C" and not (mRNA[edit_position] == "U" or mRNA[edit_position] == "u" or
                                                     mRNA[edit_position] == "G")):
            current_position += 1

        # case unresolvable mismatch with U/T at gRNA: Ignore it and go to next position
        elif ((gRNA[current_position] == "U" or gRNA[current_position] == "T") and not (
                mRNA[edit_position] == "A" or mRNA[edit_position] == "G")):
            current_position += 1

        # case match
        else:
            current_position += 1

    # Returning the edited mRNA and position after editing when done.
    # (usually position + len(gRNA), but not if mRNA-end is reached first)
    return (mRNA, edit_position)


# Edits the mRNA.
# If there is an unresolvable mismatch, editing stops and gRNA "falls off" (-> next editing block is started).
# Allowes Gu pairing when editing.
# Stops when either mRNA or gRNA ends or unresoilvable mismatch occurs.
# Anchor is not edited.
# Accepts only upper letters except u in mRNA (for speed reasons)
#
# @param
# mRNA: str of mRNA in 3' to 5' order!
# gRNA: str of gRNA in 5' to 3' order!
# position: the position on the mRNA where the gRNA did bind (1st base).
# anchorsize: anchorsize
#
# @return
# (sequence of the edited mRNA, position on mRNA after editing)
#
def overhang_falloff(mRNA, gRNA, position, anchorsize):
    assert isinstance(anchorsize, int) and anchorsize > 1, "Please provide a valid int anchorsize >1."
    assert isinstance(mRNA, str), "Please provide a valid mRNA String."
    assert isinstance(gRNA, str), "Please provide a valid gRNA String."
    assert isinstance(position, int), "Please provide a valid int binding-position."

    # Anchor is not edited, so skip it.
    current_position = anchorsize
    edit_position = current_position + position

    # Iterating through the gRNA until either gRNA or mRNA is finished.
    while current_position < len(gRNA) and edit_position < len(mRNA):

        # case A at gRNA mismatch: Insert u
        if (gRNA[current_position] == "A") and (
                mRNA[edit_position] != "U" or mRNA[edit_position] != "u"):
            mRNA = insertion(mRNA, edit_position)
            current_position += 1

        # case G at gRNA mismatch: Insert u
        elif (gRNA[current_position] == "G") and not (
                mRNA[edit_position] == "U" or mRNA[edit_position] == "u" or mRNA[edit_position] == "C"):
            mRNA = insertion(mRNA, edit_position)
            current_position += 1

        # case U at mRNA mismatch. Delete U
        elif (mRNA[edit_position] == "U" or mRNA[edit_position] == "u") and not (
                gRNA[current_position] == "A" or gRNA[current_position] == "G"):
            mRNA = deletion(mRNA, edit_position)

        # case unresolvable mismatch with C at gRNA: Fall off
        elif (gRNA[current_position] == "C" and not (mRNA[edit_position] == "u" or
                                                     mRNA[edit_position] == "U" or mRNA[edit_position] == "G")):
            break

        # case unresolvable mismatch with U/T at gRNA: Fall off
        elif ((gRNA[current_position] == "U" or gRNA[current_position] == "T") and not (
                mRNA[edit_position] == "A" or mRNA[edit_position] == "G")):
            break

        # case match
        else:
            current_position += 1

    # Returning the edited mRNA and position after editing when done.
    return (mRNA, edit_position)


# Edits the mRNA.
# Creates the thermodynamically best option with a sliding window approach:
# A window size is provided. The editor goes through the duplex base per base,
# always constructing the thermodynamically best option within the window, with
# the window starting at the current position. E.g.:
# assume the alignment is this:
#    ____
# 3' ACGUGUACG  mRNA
# 5' UGCAAAA    gRNA
#
# where the anchor is the area with the line above.
# Then, the first window if windowsize == 3 will look at:
# 3' GUA    mRNA
# 5' AAA    gRNA
# And build a tree of possible edits. This tree will be:
#
# initial seq:                                GUA
# 1st position options:           GUA                     uGUA
# 2nd position options:     GUA         GAC         uGU         uuGUA
# 3rd position options:   GUA GUu     GAC GAu     uGU uGC     uuG   uuu
#
# From this tree, the highest scoring leaves are selected and used as initial sequence for the next iteration,
# then starting at the next position (2nd position in the tree)
#
# Editing stops when the optimal solution for the window includes more mismatches than allowed.
#
# @param
# mRNA: The mRNA to edit
# gRNA: The gRNA to edit with
# position: The position of the mRNA where to align the gRNA
# anchorsize: Size of the anchor region
# windowsize: Size of the window where to find the thermodynamical best duplex
#
# @return
# a dict of the possible edited mRNAs. Structure:
# {<mRNA-seq>: <position-after-editing>}
#
def slidingthermo_editing(mRNA, gRNA, position, anchorsize, windowsize, mismatches_currently, mismatches_allowed,
                          value_GC, value_AU, value_GU, value_mismatch):
    # TODO assertions

    # Necessary functions for slidingthermo_editing().

    # Constructing the thermodynamical best sequences in the given window.
    #
    # @param
    # mRNA: the mRNA beginning at the binding position of the gRNA
    # gRNA: full gRNA
    #
    def get_best_edits(mRNA_editingregion, gRNA, current_position, windowsize, current_mismatches, value_GC, value_AU,
                       value_GU, value_mismatch):

        # Necessary functions for get_best_edits().

        # A function that tests which kind of match occurs.
        # returns 0 for mismatch, 1 for GC match, 2 for AU match, 3 for GU match
        def which_match(base1, base2):
            if (base1 == "G" and base2 == "C") or (base1 == "C" and base2 == "G"):
                return 1
            elif (base1 == "A" and (base2 == "T" or (base2 == "u" or base2 == "U"))) or (
                    (base1 == "T" or (base1 == "U" or base1 == "u")) and base2 == "A"):
                return 2
            elif ((base1 == "u" or base1 == "U" or base1 == "T") and base2 == "G") or (
                    base1 == "G" and (base2 == "u" or base2 == "U" or base2 == "T")):
                return 3
            else:
                return 0

        # Counting the number of mismatches between two sequences. nuc1 and nuc2 have to be of the same length
        def count_mismatches(nuc1, nuc2):
            assert len(nuc1) == len(nuc2), "Please provide sequences of the same length for comparison."

            mismatches = 0

            for current_base in range(len(nuc1)):
                if which_match(nuc1[current_base], nuc2[current_base]) == 0:
                    mismatches += 1

            return mismatches

        # The function builds a tree to calculate all possible best edits.
        # A node-class is needed to build the tree.
        # Possible improvement: change this to Node from Networks.py
        class Node():
            def __init__(self, position, mRNA):

                # Information about this node that is given
                self.position = position
                self.mRNA = mRNA

                # Information about this node that is calculated
                self.value = None
                self.mismatches = 0
                self.created_by_deletion = False

                # The three different kinds of children
                self.stay = None
                self.insert = None
                self.delete = None
                # All children
                self.children = []

                # TODO that is pretty ugly. I should either do it without the children list or only with the children list.

            # Add a child without editing at current position
            def add_stay(self):
                self.stay = Node(self.position + 1, self.mRNA)
                if self.created_by_deletion and (self.mRNA[self.position] == "U" or self.mRNA[self.position] == "u"):
                    self.stay.created_by_deletion = True
                self.children.append(self.stay)

            # Add a child by inserting u at current position
            def add_insert(self):
                self.insert = Node(self.position + 1, self.mRNA[:self.position] + "u" + self.mRNA[self.position:])
                self.children.append(self.insert)

            # Add a child by deleting u at current position
            def add_delete(self):
                self.delete = Node(self.position, self.mRNA[:self.position] + self.mRNA[self.position + 1:])
                self.delete.created_by_deletion = True
                self.children.append(self.delete)

            # Building the tree with self as root
            def build_tree(self, depth):

                # Check if editing is wanted and possible
                if depth < windowsize and len(self.mRNA) > self.position:

                    # Staying is always an option.
                    self.add_stay()
                    self.stay.build_tree(depth + 1)

                    # If current position is not an U and if a U was not deleted from this position before,
                    # inserting an u is an option.
                    if not (self.mRNA[self.position] == "U" or self.mRNA[self.position] == "u") and not \
                            self.created_by_deletion:
                        self.add_insert()
                        self.insert.build_tree(depth + 1)

                    # If current position is an U, deletion is an option.
                    if self.mRNA[self.position] == "U" or self.mRNA[self.position] == "u":
                        self.add_delete()
                        # as deletion brings a new base to the current position, don't move to next position.
                        self.delete.build_tree(depth)

            # append all leaves to a given list.
            def get_leaves(self, list):
                # could we check by doing if children? i.e. list not empty?
                at_least_one_child = False

                if self.delete is not None:
                    self.delete.get_leaves(list)
                    at_least_one_child = True

                if self.stay is not None:
                    self.stay.get_leaves(list)
                    at_least_one_child = True

                if self.insert is not None:
                    self.insert.get_leaves(list)
                    at_least_one_child = True

                # It's a leaf!
                if not at_least_one_child:
                    list.append(self)

            # calculate the value of the node and update the number of mismatches.
            # Only valid for leaves!
            def calculate_value(self, gRNA, value_GC, value_AU, value_GU, value_mismatch):

                value = 0

                # go through all bases and sum up their scores
                for current_position in range(windowsize):

                    # what kind of basepairing do we have here?
                    # 0 for mismatch, 1 for GC match, 2 for AU match, 3 for GU match
                    basepair = which_match(self.mRNA[self.position - windowsize + current_position],
                                           gRNA[self.position - windowsize + current_position])

                    if basepair == 1:
                        value += value_GC
                    elif basepair == 2:
                        value += value_AU
                    elif basepair == 3:
                        value += value_GU
                    else:
                        value += value_mismatch
                        self.mismatches += 1

                # set node value
                self.value = value

        # === End of definitions for get_best_edits(). Now the actual Program starts.

        # Build the tree of possible edits
        root = Node(current_position, mRNA_editingregion)
        root.mismatches = current_mismatches
        root.build_tree(0)

        # Get the leaves = possible results
        leaves = []
        root.get_leaves(leaves)

        # Calculate the values for the leaves
        [leaf.calculate_value(gRNA, value_GC, value_AU, value_GU, value_mismatch) for leaf in leaves]

        # Calculate the number of mismatches in the given mRNA before the window.
        mismatches_before = count_mismatches(gRNA[anchorsize: current_position],
                                             mRNA_editingregion[anchorsize: current_position])

        # Get a list of all values where the number of mismatches is within the limit
        leaf_values = [leaf.value for leaf in leaves if mismatches_before + leaf.mismatches <= mismatches_allowed]

        # A dict to store the results
        best_mRNAs = {}

        # Check if editing can occur with the given limit of mismatches.
        # If yes, fill dict with the best mRNAs: {<new_mRNA-subsequence>: [<nr-of-mismatches>, <thermodynamic score>]}
        # Do we need to check the inner if statement? If it got this far surely the mms etc. have been checked already?
        # ALT when check mms above, use that to at the same time remove leaves with too many?
        # because this is looping thru leaves, not leaf_values. So while leaf_values are those with < max mms, leaves
        # contain all leaves still, even if > max mms
        if leaf_values:
            max_value = max(leaf_values)
            for leaf in leaves:
                if leaf.value == max_value and leaf.mismatches <= mismatches_allowed:
                    best_mRNAs[leaf.mRNA] = [leaf.mismatches, leaf.value]

        # return dict of best mRNAs with number of mismatches and thermodynamic score
        # or emtpy dict if not possible with given maximum mismatches.
        return best_mRNAs

    # == End of definitions for slidingthermo_editing(). Now the actual program starts.

    # Current position to investigate ON THE gRNA.
    current_position = anchorsize
    # The generated mRNAs with the number of mismatches.
    # structure: {mRAN-seq: [nr-of-mismatches-currently, position-after-editing-if-stopped-here, thermodynamic score]}
    # Initial score of 0 until actual score assigned
    new_mRNAs = {mRNA: [mismatches_currently, position + anchorsize, 0]}

    # 1st loop: Going through all positions until either gRNA or mRNA ends or too many mismatches. Here: gRNAs
    while current_position + windowsize - 1 < len(gRNA):

        # A dict to store the new generated mRNAs. Same structure as new_mRNAs.
        next_mRNAs = {}

        # If there are multiple sequences with the same score, I have to iterate through them.
        for current_mRNA in new_mRNAs.keys():

            # 2nd loop: Going through all positions until either gRNA or mRNA ends or too many mismatches. Here: mRNA
            if position + current_position + windowsize - 1 < len(current_mRNA) and \
                    new_mRNAs[current_mRNA][0] <= mismatches_allowed:

                # Getting the best possible edits at this position (most stable duplex within the window)
                best_edits = get_best_edits(current_mRNA[position:], gRNA, current_position, windowsize,
                                            new_mRNAs[current_mRNA][0], value_GC, value_AU, value_GU, value_mismatch)

                # If editing is possible, test if the result already exists and store it if not.
                if best_edits != {}:
                    for new_mRNA, details in best_edits.items():
                        already_there = False
                        for key in next_mRNAs.keys():
                            if key.upper() == (current_mRNA[:position].upper() + new_mRNA.upper()):
                                already_there = True
                        if not already_there:
                            # A new mRNA was generated. Hooray! Now store it with the correct structure:
                            # {mRNA-seq: [nr-of-mismatches-currently, position-after-editing-if-stopped-here]}
                            next_mRNAs[current_mRNA[:position] + new_mRNA] = \
                                [details[0], position + current_position + windowsize, details[1]]
                else:
                    next_mRNAs[current_mRNA] = new_mRNAs[current_mRNA]
                    new_mRNAs[current_mRNA][0] = mismatches_allowed + 1

        # Set the just generated mRNAs as initial mRNAs for the next iteration and go to next position.
        # Reset new mRNAs. So new mrnas are the newly produced best edits from each round.
        # They are initially ka next mrnas when made. Once duplicates are accounted for, next becomes new
        # i.e. newly created and ready to be edited again.
        new_mRNAs = next_mRNAs
        current_position += 1

    # Delete double mRNAs
    all_keys = list(new_mRNAs.keys())
    for key1 in all_keys:
        occurrences = 0
        for key2 in new_mRNAs.keys():
            if key1.upper() == key2.upper():
                occurrences += 1
        if occurrences > 1:
            new_mRNAs.pop(key1)

    # As the number of mismatches is not interesting as output, delete it from dict.
    # Only position after editing stays: {mRAN-seq: position-after-editing-if-stopped-here}
    for key in new_mRNAs.keys():
        new_mRNAs[key] = new_mRNAs[key][1:]

    # Return all generated sequences.
    return new_mRNAs
