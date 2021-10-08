from COLLECTIONS import Programs_round2
'''
This is the file to execute when you want to run a simulation. 
To do so, you must 
 1. provide the paths to the input files
 2. specify the type of simulation you want
 3. specify the parameters of your simulation
 4. specify the type of output you want and where to store it. 
Then, simply run  this program. 

I am sure there are still several bugs in this project. If you have any questions, 
feel free to contact me at info@belafrohn.de. 
'''

# Leaving COX3 out for now bc takes so long
genes = ["RPS12", "ND9", "ND8_v2", "ND8_v1", "ND7", "ND3", "MURF2", "CR4", "CR3", "CYB", "COX3", "A6_v1", "A6_v2"]
genes2 = [ "RPS12","ND9","ND8_v2","ND8_v1", "ND7"]
genes3 = ["CR4"]
round2 = ["RPS12_340_o", "RPS12_281_o", "RPS12_282_o", "ND8_v2_181_o", "ND8_v2_184_o", "ND8_v2_3575_o", "ND8_v1_208_o", "ND8_v1_211_o", "ND8_v1_3724_o", "ND7_1493_o", "ND7_1495_o", "ND7_2532_o", "ND3_898_o", "ND3_899_o", "ND3_900_o", "CYB_4052_o", "CYB_4053_o", "CYB_4054_o", "CR4_1028_o", "CR4_1029_o", "CR4_3107_o", "CR3_88_o", "CR3_115_o", "CR3_116_o", "COX3_790_o", "COX3_1489_o", "COX3_1491_o", "A6_v1_636_o", "A6_v1_1547_o", "A6_v1_1848_o", "A6_v2_673_o", "A6_v2_1575_o", "A6_v2_1876_o"]

simulationtypes = ["EditingArea"]

for gene in round2:

    for type in simulationtypes:

        # Step 1: Provide the paths to the input files

        # The unedited sequence. Must be a fasta file. Sequence must be from 5' to 3'.
        filepath_gene = "./FILES/SecondRound/ColourOnly/"+gene+".fa"
        #if gene.find("ND8") > -1 or gene == 'ND3':
        #    filepath_gene = "../STEP/FILES/new_unedited_mRNAs/unedited_"+gene.split("_")[0]+"_further.fa"

        # The fully edited sequence. Must be a fasta file. Sequence must be from 5' to 3'.
        if gene.find("_v") > -1:
            filepath_control = "../STEP/FILES/new_edited_mRNAs/edited_"+gene.split("_")[0]+"_"+gene.split("_")[1]+".fa"
        else:
            filepath_control = "../STEP/FILES/new_edited_mRNAs/edited_"+gene.split("_")[0]+".fa"

        # The set of gRNAs. Must be a fasta file, each gRNA an entry. Sequence must be from 5' to 3' (and the sequence that forms the duplex, so after all it resembles the reverse complement of the fully edited mRNA).
        # The names of the gRNAs must follow this convention:
        # >mO_005(II)_gA6(691-741)
        # mO_005(II)_gA6 - name of gRNA
        # (691-741) - position of gRNA on fully edited sequence
        #if (gene.find("ND") > -1 and gene != "ND3") or (gene == "RPS12") or (gene =="MURF2"):
        #    filepath_gRNAs = "./FILES/new_grnas/ideal_"+gene+"_grnas.fa"
        #else:
        #    filepath_gRNAs = "./FILES/new_grnas/"+gene+"_grnas.fa"
        filepath_gRNAs = "./FILES/new_grnas/ideal_all_grnas.fa"


        # The legend for the graph that provides sequences for node names. Must follow this convention:
        # <nodename>,<seq>\n
        # E.g.:
        # 1,AUGCAG
        # 2,UUCAACC
        if gene.find("_v") > -1:
            filepath_nodes = "./FILES/legends/" + gene.split("_")[0] + "_" + gene.split("_")[1] + ".txt"
        else:
            filepath_nodes = "./FILES/legends/"+gene.split("_")[0]+".txt"


        # Step 2: Specify the type of simulation you want to run.

        # Editing type. Allowed: Uinsert, Falloff, Mismatch, SlidingThermo
        editing_type = "SlidingThermo"

        # Graph type. Allowed: All, ColourOnly, Pruned
        graph_type = "ColourOnly"

        # Simulation type. Allowed: All, Window, EditingArea
        simulation_type = type


        # Step 3: Specify the parameters of the simulation you want. Which ones you have to specify depends on your Step 2.

        # A: Specifications always to make.

        # Anchorsize. Multiple values possible, resulting in multiple runs and outputs.
        anchormin = 6
        anchormax = 21

        # Maximum depth of recursion. If a highly ramified Graph is expected, this can limit calculation time to a responsible limit.
        maximum_depth = None

        # Probability below which editing will not occur
        threshold = 1e-150

        # which position should be the first to investigate?
        starting_position = 0

        # B: Specifications depending on your choice of simulation or editing type type.

        # B.1: If you chose Window as simulation type

        # Overlap. How far to 3' should the just edited sequence be investigated?
        overlap = 20
        # Overhang. How far to the 5' end should be investigated?
        overhang = 0

        # B.2: If you chose SlidingThermo as editing type

        # Values for the Thermodynamics
        value_GC = 30
        value_AU = 20
        value_GU = 15
        value_mismatch = 0
        windowsize = 3
        mismatches_allowed = 2

        second_round = True


        # Step 4: Specify your output.

        # Path to directory where to store the generated data.
        directorypath_output = "./OUTPUT/SecondRound"

        # The type of gRNA data used. E.g. "Real", "Artificial", "RealSubsetA6"...
        used_gRNAs = "RealAll_midthermo_round2"

        # The used gene.
        used_gene = gene




        # ======== RUNNING THE SIMULATION =========

        myprogram = Programs_round2.Program()
        myprogram.set_filepaths(filepath_gene, filepath_control, filepath_gRNAs, filepath_nodes, directorypath_output)

        if gene == "ND7":
            simulation_type = "All"
        myprogram.set_parameters(editing_type, simulation_type, graph_type, anchormin, anchormax, starting_position, overlap, overhang, value_GC, value_AU, value_GU, value_mismatch, windowsize, mismatches_allowed, threshold)

        # run the simulation. OLD: at the moment returns name of simulation. Do I want that? No.
        myprogram.run_simulation(used_gRNAs, used_gene, maximum_depth, second_round)

            # OLD: clear program? Not necessary as new Simulation is created for each run_simulation()

        # save the generated data
        myprogram.save_data(directorypath_output)

        # display and safe the generated graph in graphviz / png format.
        myprogram.save_and_show_graphviz_windows(directorypath_output)

