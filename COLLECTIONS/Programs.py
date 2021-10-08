from COLLECTIONS import BF_Simulations_anchor
import os

# TODO: anchor assertions

'''
This class wraps the whole program. It stores the parameters and outputs, 
so you can run the simulation or store results. 
Different methods allow a range of things to do with the simulation. 
Mainly, this class exists to keep files to generate data (e.g. main_explained.py) non-messy and consistent 
while allowing a broad range of different actions. 
'''

class Program():
    def __init__(self):

        # Paths to input files
        self.filepath_gene = None
        self.filepath_control = None
        self.filepath_gRNAs = None
        self.filepath_nodes = None

        # Boolean to check if input was read already
        self.paths_set = False

        # Specification of simulation type
        self.editing_type = None
        self.simulation_type = None
        self.graph_type = None

        # Parameters of the simulation
        self.anchor_min = None
        self.anchor_max = None
        self.starting_position = None
        self.overlap = None
        self.overhang = None
        self.maximum_depth = None
        self.value_GC = None
        self.value_AU = None
        self.value_GU = None
        self.windowsize = None
        self.mismatches_allowed = None
        self.threshold = None

        # Boolean to check if parameters are set
        self.params_set = False

        # Names of data used; for naming output
        self.used_gRNAs = None
        self.used_gene = None
        self.directorypath_output = None

        # Boolean to check if output was specified
        self.names_set = False

        # The simulation of this program.
        self.simulation = None

        # Information if simulation did run
        self.simulation_finished = False

    # set the paths to the necessary files.
    def set_filepaths(self, gene, control, gRNAs, nodes, output):
        assert isinstance(gene, str) and os.path.exists(gene), "Please provide a valid path to the gene file. Wrong path: "+gene
        assert isinstance(control, str) and os.path.exists(control), "Please provide a valid path to the control file. Wrong path: "+control
        assert isinstance(gRNAs, str) and os.path.exists(gRNAs), "Please provide a valid path to the gRNAs file. Wrong path: "+gRNAs
        assert isinstance(nodes, str) and os.path.exists(nodes), "Please provide a valid path to the nodes file. Wrong path: "+nodes
        assert isinstance(output, str) and os.path.exists(output), "Please provide a valid path to the output directory. Wrong path: "+output

        self.filepath_gene = gene
        self.filepath_control = control
        self.filepath_gRNAs = gRNAs
        self.filepath_nodes = nodes
        self.directorypath_output = output

        self.paths_set = True

    # specify parameters.
    def set_parameters(self, editing_type, simulation_type, graph_type, anchor_min, anchor_max, starting_position,
                       overlap=None, overhang=None, value_GC=None, value_AU=None, value_GU=None, windowsize=None,
                       mismatches_allowed=None, threshold=None):
        assert editing_type in ("Uinsert", "Falloff", "Mismatch", "SlidingThermo"), "Please provide a valid editing " \
            "type. Allowed types are Uinsert, Falloff, Mismatch, SlidingThermo. "
        assert simulation_type in ("All", "Window", "EditingArea"), "Please provide a valid simulation type. Allowed " \
            "types are All, Window, EditingArea."
        assert graph_type in ("All", "ColourOnly", "Pruned"), "Please provide a valid graph type. Allowed types are " \
            "All, ColourOnly, Pruned. "
        assert anchor_min > 1, "Please provide an anchor_min > 1."
        assert anchor_max > 1, "Please provide an anchor_max > 1"
        assert starting_position >= 0, "Please provide a starting position >= 0."
        assert (overlap is None and overhang is None) or overlap + overhang > 0, "Please provide a valid window or no window."
        assert (value_GU == value_AU == value_GC == windowsize == mismatches_allowed is None) or \
               (isinstance(value_GC, int) and isinstance(value_AU, int) and isinstance(value_GU, int) and
                isinstance(windowsize, int) and isinstance(mismatches_allowed, int)), \
            "Please provide valid Thermodynamics parameters"

        self.editing_type = editing_type
        self.simulation_type = simulation_type
        self.graph_type = graph_type
        self.anchor_min = anchor_min
        self.anchor_max = anchor_max
        self.starting_position = starting_position
        self.overlap = overlap
        self.overhang = overhang
        self.value_GC = value_GC
        self.value_AU = value_AU
        self.value_GU = value_GU
        self.windowsize = windowsize
        self.mismatches_allowed = mismatches_allowed
        self.threshold = threshold

        self.params_set = True

    # run the simulation if all parameters are set.
    def run_simulation(self, used_gRNAs, used_gene, maximum_depth=None):
        assert self.paths_set, "You must provide paths to the necessary files first."
        assert self.params_set, "You must provide parameters first."
        assert isinstance(used_gRNAs, str), "Please provide a valid name for the gRNAs you used."
        assert isinstance(used_gene, str), "Please provide a valid name for the gene you are investigating."
        assert maximum_depth is None or maximum_depth > 0, "Please provide a maximum depth >0."

        self.used_gRNAs = used_gRNAs
        self.used_gene = used_gene
        self.maximum_depth = maximum_depth

        self.names_set = True

        # specify the name of the simulation. This will be the filename if the generated data will be stored or printed.
        if self.simulation_type == "All" or self.simulation_type == "EditingArea":
            if self.editing_type == "SlidingThermo":
                simulation_name = "Simulation=" + self.simulation_type + "_Editing=" + self.editing_type + \
                                  "_GC="+str(self.value_GC)+"_AU="+str(self.value_AU)+"_GU="+str(self.value_GU) + \
                                  "_windowsize="+str(self.windowsize)+"mismatches="+str(self.mismatches_allowed) + \
                                  "_Graph=" + self.graph_type + "_Anchor=" + str(self.anchor_min) + "-" + \
                                  str(self.anchor_max) + "_Depth="+str(self.maximum_depth)+ "_threshold=" + \
                                  str(self.threshold) + "_" + self.used_gene + "_" + self.used_gRNAs
            else:
                simulation_name = "Simulation=" + self.simulation_type + "_Editing=" + self.editing_type + "_Graph=" + \
                                  self.graph_type + "_Anchor=" + str(self.anchor_min) + "-" + str(self.anchor_max) + \
                                  "_Depth="+str(self.maximum_depth) + "_threshold=" + str(self.threshold) + "_" + \
                                  self.used_gene + "_" + self.used_gRNAs
        elif self.simulation_type == "Window":
            if self.editing_type == "SlidingThermo":
                simulation_name = "Simulation=" + self.simulation_type + "_Editing=" + self.editing_type + "_GC=" + \
                                  str(self.value_GC)+"_AU="+str(self.value_AU)+"_GU="+str(self.value_GU)+"_windowsize="\
                                  +str(self.windowsize)+"mismatches="+str(self.mismatches_allowed)+ "_Graph="+\
                                  self.graph_type + "_Anchor=" + str(self.anchor_min) + "-" + str(self.anchor_max) + \
                                  "_Depth="+str(self.maximum_depth)+ "_lap,hang=" + str(self.overlap) + "," + \
                                  str(self.overhang) + "_threshold=" + str(self.threshold) + "_" + self.used_gene + \
                                  "_" + self.used_gRNAs
            else:
                simulation_name = "Simulation=" + self.simulation_type + "_Editing=" + self.editing_type + "_Graph=" + \
                    self.graph_type + "_Anchor=" + str(self.anchor_min) + "-" + str(self.anchor_max) + "_Depth=" + \
                    str(self.maximum_depth) + "_lap,hang=" + str(self.overlap) + "," + str(self.overhang) + \
                    "_threshold=" + str(self.threshold) + "_" + self.used_gene + "_" + self.used_gRNAs

        # Create simulation
        self.simulation = BF_Simulations_anchor.Simulation(name=simulation_name)

        # initialize simulation
        self.simulation.initialize(self.filepath_gene, self.filepath_control, self.filepath_gRNAs, self.filepath_nodes)

        # run simulation
        self.simulation.run(editing_type=self.editing_type, simulation_type=self.simulation_type,
                            graph_type=self.graph_type, starting_position=self.starting_position,
                            anchor_min=self.anchor_min, anchor_max=self.anchor_max, overlap=self.overlap,
                            overhang=self.overhang, maximum_depth=self.maximum_depth, threshold=self.threshold,
                            windowsize=self.windowsize, mismatches_allowed=self.mismatches_allowed,
                            value_GC=self.value_GC, value_AU=self.value_AU, value_GU=self.value_GU, value_mismatch=0,
                            end_of_editingarea=len(self.simulation.mRNA_unedited))

        # Save information that simulation was executed to allow printing and saving.
        self.simulation_finished = True

    # TODO: save_data()

    def save_and_show_graphviz_windows(self, directorypath_output):
        assert self.simulation_finished, "You must run a simulation before you can display any results."
        assert self.simulation is not None, "Simulation not found."
        assert isinstance(directorypath_output, str) and os.path.exists(directorypath_output), "Please provide a valid output directory."

        self.simulation.graphviz_windows(directorypath_output)

    def save_data(self, directorypath_output):
        assert self.simulation_finished, "You must run a simulation before you can display any results."
        assert self.simulation is not None, "Simulation not found."
        assert isinstance(directorypath_output, str) and os.path.exists(directorypath_output), "Please provide a valid output directory."

        self.simulation.save_data(directorypath_output)




