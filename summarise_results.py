#from ThermoOnly.COLLECTIONS import Stats
from COLLECTIONS import BF_Stats

"""
This file will generate a statistical summary of the results from the simulation.
You must provide:
    1. The path to the output .DATA.txt file
    2. The path to a directory containing the fully-edited sequences.
    3. The names of all genes of interest.
"""

# The names of genes of interest.
genes = ["RPS12", "ND9", "ND8_v2", "ND8_v1", "ND7", "ND3", "MURF2", "CYB", "CR4", "CR3", "COX3", "A6_v1", "A6_v2"]
genes2 = ["ND9", "ND8_v2", "ND8_v1", "ND7"]

for gene in genes:
    # The path to the directory containing genes of interest.
    # It is assumed that all edited mRNAs are stored with the same naming convention.
    filepath_mRNA = "./FILES/new_edited_mRNAs"

    # Where are they stored?
    # The .DATA.txt file for the genes of interest.
    # It is assumed that all simulations were run under the same parameters
    # i.e. the gene name is the only thing that should change in the filepath.
    output_directory = "./ThermoOnly/EditingArea/ColourOnly/All"

    # What anchor sizes do you have?
    anchormin = 6
    anchormax = 21

    # Editing type. Allowed: Uinsert, Falloff, Mismatch, SlidingThermo
    editing_type = "SlidingThermo"

    # Graph type. Allowed: All, ColourOnly, Pruned
    graph_type = "ColourOnly"

    # Simulation type. Allowed: All, Window, EditingArea
    simulation_type = "EditingArea"

    # The type of gRNA data used. E.g. "Real", "Artificial", "RealSubset"...
    used_gRNAs = "RealAll_midthermo"

    # Maximum depth of recursion.
    # If a highly ramified Graph is expected, this can limit calculation time to a responsible limit.
    maximum_depth = None

    # which position should be the first to investigate?
    starting_position = 0

    threshold = 1e-50

    # Overlap. How far to 3' should the just edited sequence be investigated?
    overlap = 40
    # Overhang. How far to the 5' end should be investigated?
    overhang = 10

    # Values for the Thermodynamics
    value_GC = 30
    value_AU = 20
    value_GU = 15
    windowsize = 3
    mismatches_allowed = 2

    if editing_type == "SlidingThermo":
        if used_gRNAs.find("Real") > -1:
            if gene == 'ND7':
                simulation_type = "All"
            file = output_directory + "/Simulation=" + simulation_type + "_Editing=" + editing_type + "_GC=" + str(
                value_GC) + "_AU=" + str(value_AU) + "_GU=" + str(value_GU) + "_windowsize=" + str(
                windowsize) + "mismatches=" + str(mismatches_allowed) + "_Graph=" + graph_type + "_Anchor=" + str(
                anchormin) + "-" + str(anchormax) + "_Depth=" + str(maximum_depth) + "_threshold=" + str(threshold) +\
                   "_" + gene + "_" + used_gRNAs \
                + ".DATA.txt"
        else:
            file = output_directory + "/Simulation=" + simulation_type + "_Editing=" + editing_type + "_GC=" + str(
                value_GC) + "_AU=" + str(value_AU) + "_GU=" + str(value_GU) + "_windowsize=" + str(
                windowsize) + "mismatches=" + str(mismatches_allowed) + "_Graph=" + graph_type + "_Anchor=" + \
                str(anchormin) + "-" + str(anchormax) + "_Depth=" + str(maximum_depth) + "_" + gene + "_" + used_gRNAs + \
                ".DATA.txt"
    else:
        if used_gRNAs.find("Real") > -1:
            file = output_directory + "/Simulation=" + simulation_type + "_Editing=" + editing_type + "_Graph=" + \
                   graph_type + "_Anchor=" + str(anchormin) + "-" + str(anchormax) + "_Depth=" + str(maximum_depth) + "_" \
                   + gene + "_" + used_gRNAs + gene + ".DATA.txt"
        else:
            file = output_directory + "/Simulation=" + simulation_type + "_Editing=" + editing_type + "_Graph=" + \
                   graph_type + "_Anchor=" + str(anchormin) + "-" + str(anchormax) + "_Depth=" + str(maximum_depth) + \
                   "_" + gene + "_" + used_gRNAs + ".DATA.txt"

    # Read in expected mRNA sequence.
    with open(filepath_mRNA + "/edited_" + gene + ".fa") as f:
        data = f.read()
        data = data.strip().split("\n")
        full_edited_mRNA = data[1]
    #full_edited_mRNA = open(filepath_mRNA + "/edited_" + gene + ".fa").read().strip().split("\n")[1]

    # Create stats generator
    my_stats_generator = BF_Stats.Stats_Generator()

    # Initialize
    my_stats_generator.read_graph(file)

    # Generate desired statistical output.
    positions = len(full_edited_mRNA) - my_stats_generator.get_correct_coverage_position(full_edited_mRNA)
    coverage = my_stats_generator.get_correct_coverage_percentage(full_edited_mRNA)
    correct_node = my_stats_generator.get_most_correct_node(full_edited_mRNA)
    worst_node, best_node, leaf_Ps = my_stats_generator.get_probability_nodes(full_edited_mRNA)
    print("2m, 3w, midthermo, boltzmann, all")
    print("Gene", gene, ":", str(coverage) + "% coverage", "\nProbability of most correct node:", correct_node.probability,
          "\nsequence:", correct_node.mRNA, "\nname:", correct_node.name, "\nHighest probability node:",
          best_node.probability, "\nsequence:", best_node.mRNA, "\nname:", best_node.name, "\nLowest probability node:",
          worst_node.probability, "\nsequence:", worst_node.mRNA, "\nname:", worst_node.name,)


# correct COX3     my_stats_generator.read_graph("./OUTPUT/new_seqs/thermo_colour_4mm/Simulation=All_Editing=SlidingThermo_GC=30_AU=20_GU=15_windowsize=3mismatches=4_Graph=ColourOnly_Anchor=6-21_Depth=None_COX3_NewReal_0threshold.DATA.txt")
