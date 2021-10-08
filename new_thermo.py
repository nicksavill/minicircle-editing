import timeit


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
def new_slidingthermo_editing(mRNA, gRNA, position, anchorsize, windowsize, mismatches_currently, mismatches_allowed,
                              value_GC, value_AU, value_GU, value_mismatch):
    # TODO assertions

    # Necessary functions for slidingthermo_editing().

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

        # Normalises thermodynamic score by dividing the leaf's value by the total of all leaf values
        # Multiplies probability by that of the parent edit to allow tracking probability along the window
        def normalise_weights(self, root_probability):

            # first need to calculate value of each leaf
            for leaf in self.leaves:
                leaf.calculate_value(gRNA, value_GC, value_AU, value_GU, value_mismatch)

            total_thermo_score = 0
            for leaf in self.leaves:
                total_thermo_score += leaf.value
            for leaf in self.leaves:
                leaf.probability = (leaf.value / total_thermo_score) * root_probability
                #print("probability calculated", leaf.value, leaf.probability, root_probability)

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
    def get_all_edits(mRNA_editingregion, gRNA, current_position, windowsize, current_mismatches, value_GC, value_AU,
                      value_GU, value_mismatch, root_probability=1):

        # === End of definitions for get_all_edits(). Now the actual get_all_edits function starts.

        # Build the tree of possible edits
        root = Node(current_position, mRNA_editingregion)
        root.probability = root_probability
        root.mismatches = current_mismatches
        root.build_tree(0)

        # Get the leaves = possible results
        leaves = []
        root.get_leaves(leaves)
        root.leaves = leaves
        root.normalise_weights(root.probability)

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
    #               [actual-mRNA-seq, nr-of-mismatches-currently, position-after-editing, thermodynamic probability]}
    new_mRNAs = {mRNA.upper(): [mRNA, mismatches_currently, position + anchorsize, 1]}  # Initial P = 1 until assigned
    #print("starting with", new_mRNAs)
    # 1st loop: Going through all positions until either gRNA or mRNA ends or too many mismatches. Here: gRNAs
    while current_position + windowsize - 1 < len(gRNA):

        # Store the newly-generated mRNAs. Same structure as new_mRNAs.
        next_mRNAs = {}

        print("new round. Investigating", len(new_mRNAs))
        #count = 1

        # Each round, iterate through all generated sequences to produce all new possible edits.
        # Reminder: values structure: [mRNA-seq, current-mismatches, position-after-editing, thermodynamic-probability]
        for uppercase_mRNA, values in new_mRNAs.items():
            #print(count)
            current_mRNA = values[0]

            # Go through all positions until gRNA or mRNA ends, or too many mismatches.
            # If too many/wrong position: moves to next mRNA in for loop
            if position + current_position + windowsize - 1 < len(current_mRNA) and \
                    new_mRNAs[uppercase_mRNA][1] <= mismatches_allowed:

                # Get all possible edits
                # Structure {edited-mRNA-sequence : [mismatches, thermodynamic-probability]}
                # get_all_edits params: mRNA_editingregion, gRNA, current_position, windowsize, current_mismatches,
                # value_GC, value_AU, value_GU, value_mismatch, root_probability
                edits = get_all_edits(current_mRNA[position:], gRNA, current_position, windowsize,
                                      new_mRNAs[uppercase_mRNA][1], value_GC, value_AU, value_GU, value_mismatch,
                                      new_mRNAs[uppercase_mRNA][3])
                #print("Returned", len(edits), "edits at position", current_position)

                # If the edit exists, increase its probability.
                if edits != {}:
                    for new_mRNA, details in edits.items():
                        # next_mRNAs keys are upper case edits because U and u are biologically equivalent
                        if (uppercase_mRNA[:position] + new_mRNA.upper()) in next_mRNAs:
                            #print("already there")
                            #print(details)
                            # increase probability
                            #print("probability before", next_mRNAs[uppercase_mRNA[:position] + new_mRNA.upper()][3])
                            #print("probability to add", details[1])
                            next_mRNAs[uppercase_mRNA[:position] + new_mRNA.upper()][3] += details[1]
                            #print("probability after", next_mRNAs[uppercase_mRNA[:position] + new_mRNA.upper()][3])

                        else:
                            # A new mRNA was generated. Hooray! Now store it with the correct structure:
                            # {mRNA-seq.upper():
                            # [mRNA-seq, nr-mismatches-currently, position-after-editing, thermodynamic-probability]}
                            next_mRNAs[uppercase_mRNA[:position] + new_mRNA.upper()] = \
                                [current_mRNA[:position] + new_mRNA, details[0], position + current_position + windowsize, details[1]]
                            #print("adding new one. Next now length", len(next_mRNAs))
                else:
                    #print("no edits made")
                    # May still be edited further downstream so add to list to be checked in future
                    next_mRNAs[uppercase_mRNA] = new_mRNAs[uppercase_mRNA]
                    # Make sure not passed along in duplicate by exceeding mismatches?
                    new_mRNAs[uppercase_mRNA][1] = mismatches_allowed + 1  # TODO: WHY?
            #count += 1

        # Set the just generated mRNAs as initial mRNAs for the next iteration and go to next position.
        # Reset new mRNAs. So new mrnas are the newly produced best edits from each round.
        # They are initially stored in next mrnas when made. Once duplicates are accounted for, next becomes new
        # i.e. newly created and ready to be edited again.
        print("updating for next round. Will investigate", len(next_mRNAs))
        #print("new", new_mRNAs)
        #print("next", next_mRNAs)
        new_mRNAs = next_mRNAs
        # print(next_mRNAs)
        current_position += 1

    # Final sequences to return - don't need uppercase anymore, only interested in actual edit, position, probability
    # Structure: {mRNA-seq: [position-after-editing, thermodynamic-probability]}
    final_mRNAs = {}
    top_probabilities = []
    for mRNA_info in new_mRNAs.values():
        top_probabilities.append(mRNA_info[3])

    # get 3 most probable edits
    top_probabilities = sorted(top_probabilities, reverse=True)[:5]

    for mRNA_info in new_mRNAs.values():
        if mRNA_info[3] in top_probabilities:
            final_mRNAs[mRNA_info[0]] = mRNA_info[1:]

    # Return all generated sequences.
    return final_mRNAs


def old_slidingthermo_editing(mRNA, gRNA, position, anchorsize, windowsize, mismatches_currently, mismatches_allowed,
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
                       value_GU, value_mismatch, root_probability=1):

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
                # TODO: took out self.mRNA bc I assume that'll be accessible from the parent class
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

            # Normalises thermodynamic score by dividing the leaf's value by the total of all leaf values
            # Multiplies probability by that of the parent edit to allow tracking probability along the window
            def normalise_weights(self, root_probability):

                # first need to calculate value of each leaf
                for leaf in self.leaves:
                    leaf.calculate_value(gRNA, value_GC, value_AU, value_GU, value_mismatch)

                total_thermo_score = 0
                for leaf in self.leaves:
                    total_thermo_score += leaf.value
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
        root.normalise_weights(root.probability)

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
    # structure: {mRNA-seq: [nr-of-mismatches-currently, position-after-editing-if-stopped-here, thermodynamic score]}
    new_mRNAs = {mRNA: [mismatches_currently, position + anchorsize, 1]}  #  Initial probability 0 until assigned

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
                # Structure {mRNA sequence : [mismatches, thermodynamic_probability]}
                best_edits = get_best_edits(current_mRNA[position:], gRNA, current_position, windowsize,
                                            new_mRNAs[current_mRNA][0], value_GC, value_AU, value_GU, value_mismatch,
                                            new_mRNAs[current_mRNA][2])

                # If editing is possible, test if the result already exists and store it if not.
                if best_edits != {}:
                    for new_mRNA, details in best_edits.items():
                        already_there = False
                        # Can't just test if new_mrNA in next_mRNAs because CuUG and CUuG are equivalent but will be
                        # seen as different. So need to compare uppercase
                        for key in next_mRNAs.keys():
                            if key.upper() == (current_mRNA[:position].upper() + new_mRNA.upper()):
                                already_there = True
                        if not already_there:
                            # A new mRNA was generated. Hooray! Now store it with the correct structure:
                            # {mRNA-seq: [nr-mismatches-currently, position-after-editing, thermodynamic-probability]}
                            next_mRNAs[current_mRNA[:position] + new_mRNA] = \
                                [details[0], position + current_position + windowsize, details[1]]
                else:
                    next_mRNAs[current_mRNA] = new_mRNAs[current_mRNA]
                    new_mRNAs[current_mRNA][0] = mismatches_allowed + 1  # TODO: WHY?

        # Set the just generated mRNAs as initial mRNAs for the next iteration and go to next position.
        # Reset new mRNAs. So new mrnas are the newly produced best edits from each round.
        # They are initially stored in next mrnas when made. Once duplicates are accounted for, next becomes new
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
    # Final structure: {mRNA-seq: [position-after-editing-if-stopped-here, thermodynamic-score]}
    #for key in new_mRNAs.keys():
    #    new_mRNAs[key] = new_mRNAs[key][1:]

    # Return all generated sequences.
    return new_mRNAs


unedited = "AGAAAUAUGAAUAUGUGUAUGAUAUAUAAAAACAAGGAUUUUUUGGGGGUUUAGGGACAGAGGGUUUAUUUUUGAGGAUUUUAGGAGGAGAAAAGGGAUGGGAAACAGAAGGACAUAAGAAAAGUUUCGUUAUUAGAUUAAAAAAGUAUGCAAAUAAUUUUUGU"
unedited = unedited[::-1]
gRNA = "UUAAUCUAAUAACGGAGACUUGUGUAUGAUAGUAAUGGUGAUAUUUGUAUA"
gRNA = gRNA

start = timeit.default_timer()
new_new_mRNAs = new_slidingthermo_editing(unedited, gRNA, 23, 6, 3, 0, 3, 30, 20, 15, 0)
end = timeit.default_timer()

print("new time:", end - start, "seconds")
print("Returned", len(new_new_mRNAs), "edits")


start = timeit.default_timer()
old_new_mRNAs = old_slidingthermo_editing(unedited, gRNA, 23, 6, 3, 0, 3, 30, 20, 15, 0)
end = timeit.default_timer()

print("old time:", end - start, "seconds")
print("Returned", len(old_new_mRNAs), "edits")

print("new")
for k,v in new_new_mRNAs.items():
    print(k,v)

print("old")
for k,v in old_new_mRNAs.items():
    print(k,v)

for k,v in old_new_mRNAs.items():
    for k2, v2 in new_new_mRNAs.items():
        if k.upper() == k2.upper():
            print(v[2], v2[2])



# start = timeit.default_timer()
# new_mRNAs = old_slidingthermo_editing(unedited, gRNA, 0, 6, 3, 0, 3, 30, 20, 15, 0)
# end = timeit.default_timer()

# print(len(new_mRNAs), "old time:", end - start, "seconds")
