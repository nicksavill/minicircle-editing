from COLLECTIONS import AnchorAligners, Editors, BF_Network_white
import sys, time, os

'''
This file contains the Simulation class, which stores all the generated data. 
At the moment, it is just one class, but it might be interesting to add different simulations. 
E.g. a simulation with a completely different data structure. 
'''


# TODO: track used gRNAs so can stop them re-binding; track position

class Simulation():
    def __init__(self, name=None):

        self.name = name
        self.gene_name = None

        # A dict to store the legend: {<seq>:<nodename> ... }
        self.legend = None
        # Strings of the initial / correct mRNA-sequences
        self.mRNA_unedited = None
        self.mRNA_edited = None
        # A list to store gRNAs: [ [<gRNA1_name>,<gRNA1_seq>], [<gRNA2_name>,<gRNA2_seq>], ... ]
        self.gRNAs = None
        # Start position if second round of editing
        self.second_round_start = None

        # parameters for Thermodynamic approach
        self.windowsize = None
        self.mismatches_allowed = None
        self.value_GC = None
        self.value_AU = None
        self.value_GU = None
        self.value_mismatch = None

        # Objects to run the simulation
        self.network = None
        self.root = None

        # Parameters to control graph creation
        self.to_visit = []
        self.visited = []
        self.best_whites = []

        # path to the legend file to store new nodes there
        self.nodes_filepath = None

        # booleans to check what was done already
        self.initialized = False
        self.finished = False

    # Add a network to the simulation
    def add_network(self, network):
        assert isinstance(network, BF_Network_white.Network), "Please provide a valid Network object. Current type is " + str(
            type(network)) + "."
        self.network = network

    # Add a unedited mRNA to the simulation
    def add_mRNA_unedited(self, mRNA_unedited):
        assert isinstance(mRNA_unedited, str), "Please provide a valid mRNA-unedited sequence."
        self.mRNA_unedited = mRNA_unedited

    # Add a fully edited mRNA to the simulation
    def add_mRNA_edited(self, mRNA_edited):
        assert isinstance(mRNA_edited, str), "Please provide a valid mRNA-edited sequence."
        self.mRNA_edited = mRNA_edited

    # Add a set of gRNAs in the required structure (see above) to the simulation
    def add_gRNAs(self, gRNAs):
        assert isinstance(gRNAs, list), "Please provide a valid list of gRNAs. " \
                                        "Must be of structure [ [<gRNA1_name>,<gRNA1_seq>], [<gRNA2_name>,<gRNA2_seq>], ... ]."
        for gRNA in gRNAs:
            assert isinstance(gRNA, list), "Please provide a valid list of gRNAs. " \
                                           "Must be of structure [ [<gRNA1_name>,<gRNA1_seq>], [<gRNA2_name>,<gRNA2_seq>], ... ]."
        self.gRNAs = gRNAs

    # Add a path to the global nodes file to the simulation
    def add_nodes_filepath(self, nodes_filepath):
        assert isinstance(nodes_filepath, str) and os.path.exists(
            nodes_filepath), "Please provide a valid path to legend file."
        self.nodes_filepath = nodes_filepath

    # Add a legend of existing nodes to the simulation
    def add_legend(self, legend):
        assert isinstance(legend, dict), "Please provide a valid dict for the node legend. "
        self.legend = legend

    # Add the root to the simulation
    def add_root(self, root):
        assert isinstance(root, BF_Network_white.Node), "Only Node-objects can be root."  # maybe specific graphnode object?
        self.root = root
        self.to_visit.append(root)
        root.probability = 1

    # Initialize a simulation
    def initialize(self, filepath_gene, filepath_control, filepath_RNAs, filepath_nodes, second_round=False):
        assert isinstance(filepath_gene, str) and os.path.exists(
            filepath_gene), "Provide a valid filepath for the gene."
        assert isinstance(filepath_control, str) and os.path.exists(
            filepath_control), "Provide a valid filepath for the control."
        assert isinstance(filepath_RNAs, str) and os.path.exists(
            filepath_RNAs), "Provide a valid filepath for the gRNAs."
        assert isinstance(filepath_nodes, str) and os.path.exists(
            filepath_nodes), "Provide a valid filepath for the nodes."

        # add mRNAs
        if second_round:
            with open(filepath_gene) as f:
                mRNA_info = f.read()
                self.gene_name = mRNA_info.split("\n")[0].split(">")[1]
                mRNA_unedited = mRNA_info.split("\n")[1][::-1]  # Reverse because editing is 3' to 5'.
                self.second_round_start = int(mRNA_info.split("\n")[0].split("last-block=")[1])
            self.add_mRNA_unedited(mRNA_unedited)
        else:
            with open(filepath_gene) as f:
                mRNA_unedited = f.read()
                self.gene_name = mRNA_unedited.split("\n")[0].split(">")[1]
                mRNA_unedited = mRNA_unedited.split("\n")[1][::-1]  # Reverse because editing is 3' to 5'.
            self.add_mRNA_unedited(mRNA_unedited)
        # print(self.mRNA_unedited[87])

        with open(filepath_control) as f:
            mRNA_edited = f.read()
            mRNA_edited = mRNA_edited.split("\n")[1][::-1]  # Easier to compare to generated sequences.
        self.add_mRNA_edited(mRNA_edited)
        # print(self.mRNA_edited)

        # add gRNAs
        with open(filepath_RNAs) as f:
            gRNAs = f.read()
            gRNAs = gRNAs.strip().split(">")
            gRNAs.pop(0)
            gRNAs = [gRNA.strip().split("\n") for gRNA in gRNAs]
        self.add_gRNAs(gRNAs)

        # create legend
        # Format {'<node mRNA sequence>' : 'node name'}
        with open(filepath_nodes) as f:
            nodes = f.read()
            nodes = nodes.strip().split("\n")
            nodes = [node.strip().split(",") for node in nodes]
            legend = {}
            for node in nodes:
                # sequences are stored in uppercase to avoid redundancy.
                legend[node[1].upper()] = node[0]
        self.add_legend(legend)
        self.add_nodes_filepath(filepath_nodes)

        # create network, its root (= node with unedited seq) and add root to the network.
        network = BF_Network_white.Network(name=self.name)
        # Allow first edit to occur anywhere along its length
        # Used length instead of None to avoid if statement in simulation_editing_area_only
        root = BF_Network_white.Node(name='"unedited mRNA"', mRNA=self.mRNA_unedited, edited_until=len(self.mRNA_unedited))
        network.add_node(root)
        self.add_network(network)
        self.add_root(root)

        # confirm initialization
        self.initialized = True

    # Save graph as png and gv file and display it.
    def graphviz_windows(self, directorypath_output):
        assert isinstance(directorypath_output, str) and os.path.exists(
            directorypath_output), "Please provide a valid output directory."
        assert self.finished, "The simulation didn't run yet."

        # add again! You can out-comment this to make things faster by skipping the generation of pictures by graphviz.
        # self.network.graphviz_windows(directorypath_output)
        # print("Graph picture saved to two files:", self.name + ".gv(.png) .")

    # Saving the generated data to provided directory.
    def save_data(self, directorypath_output):
        assert self.finished, "The simulation didn't run yet."
        assert isinstance(directorypath_output, str) and os.path.exists(
            directorypath_output), "Please provide a valid output directory."

        self.network.save_data(directorypath_output)
        print("Graph data saved to file:", self.name + ".DATA.txt .")

    # Run the simulation
    def run(self, editing_type, simulation_type, graph_type, starting_position, anchor_min, anchor_max, overlap,
            overhang, maximum_depth, threshold, windowsize=None, mismatches_allowed=None, value_GC=None, value_AU=None,
            value_GU=None, value_mismatch=None, end_of_editingarea=None):
        assert self.initialized, "Simulation must be initialized first."
        assert (editing_type == "Uinsert" or editing_type == "Falloff" or editing_type == "Mismatch" or
                editing_type == "SlidingThermo"), \
            "Provide a correct editing type. Allowed types are: Uinsert, Falloff, Mismatch, SlidingThermo. "
        assert (simulation_type == "All" or simulation_type == "EditingArea" or simulation_type == "Window"), \
            "Provide a correct simulationi type. Allowed types are: All, Window, EditingArea. "
        assert (graph_type == "All" or graph_type == "ColourOnly" or graph_type == "Pruned"), \
            "Provide a correct graph type. Allowed types are: All, ColourOnly, Pruned. "
        assert maximum_depth is None or (isinstance(maximum_depth, int) and maximum_depth > 0), \
            "Provide a valid maximum depth."
        assert isinstance(starting_position, int) and starting_position >= 0, "Provide a valid starting position."
        assert isinstance(anchor_min, int) and anchor_min > 1, \
            "Please provide a valid minimum anchor size. Must be int and >1"
        assert isinstance(anchor_max, int) and anchor_max > 1, \
            "Please provide a valid maximum anchor size. Must be int and >1"
        # assert isinstance(threshold, int) and 0 <= threshold <= 1, \
        #    "Please provide a valid threshold. Must be int, 0 <= threshold <= 1"
        assert windowsize is None or (isinstance(windowsize, int) and windowsize > 0), \
            "Please provide a valid windowsize for the SlidingThermo editing approach. Must be int and >0."
        assert end_of_editingarea is None or (isinstance(end_of_editingarea, int) and end_of_editingarea > 0), \
            "Please provide a valid end_of_editing_area position. Must be int and >0. "
        assert (value_GU == value_AU == value_GC == windowsize == mismatches_allowed == None) or \
               (isinstance(value_GC, int) and isinstance(value_AU, int) and isinstance(value_GU, int) and
                isinstance(windowsize, int) and isinstance(mismatches_allowed, int) and mismatches_allowed >= 0), \
            "Please provide valid Thermodynamics parameters"

        # set parameters of Thermodynamic approach
        if editing_type == "SlidingThermo":
            self.windowsize = windowsize
            self.mismatches_allowed = mismatches_allowed
            self.value_GC = value_GC
            self.value_AU = value_AU
            self.value_GU = value_GU
            self.value_mismatch = value_mismatch

        # if second round of editing: start a little ways upstream of previous editing block.
        if self.second_round_start:
            starting_position = self.second_round_start - 20

        # Tell what's going on and get starting time to compare it to end time later.
        print("Starting simulation...")
        start_time = time.time()

        # run the simulation
        if simulation_type == "All":
            self.simulation_all(editing_type, graph_type, starting_position, anchor_min, anchor_max, 1,
                                maximum_depth, threshold)

        #        if simulation_type == "Window":
        #            self.simulation_window(editing_type, graph_type, starting_position, anchor_min, anchor_max, 1,
        #                                   maximum_depth, threshold, overlap, overhang)

        if simulation_type == "EditingArea":
            self.simulation_editing_area_only(editing_type, graph_type, starting_position, anchor_min,
                                              anchor_max, 1, maximum_depth, threshold, end_of_editingarea)

        # get end time and print feedback.
        end_time = time.time()
        print("Simulation finished.")
        print("Time consumed:", end_time - start_time, "seconds")

        # store that simulation was executed.
        self.finished = True

    # Simulate all possible edits: Searching every position of the provided mRNA for matching gRNAs.
    # If a match is found, the mRNA is edited with the given editing type and this new mRNA is introduced as a new node
    # to the network.
    # The new node is searched again, resulting in a recursive building of the network.
    # But the starting position of the new node, so where to start the investigation, is position+1.
    # Thus, leaves are not re-introduced (the process does not start all over again)
    # E.g.:
    # If this is your current mRNA:
    # 3' UGCAACACGUGUACGACUAGCUGACG
    # and your gRNA binds at the first G, say like this:
    # 3' UGCAACACGUGUACGACUAGCUGACG
    # 5'  cguugugccaa
    # which then results in this edited mRNA
    # 3' UGCAACACGGUuACGACUAGCUGACG
    # 5'  cguugugccaa
    # Then the next node will have this mRNA as sequence, but the first base that is allowed for gRNAs to bind will be
    # the C (position + 1).
    #
    # @ param
    # editing type: Which editing scheme should be used. Choices: Uinsert, Falloff, Mismatch, SlidingThermo
    # graph_type: Which type of graph should be generated? Choices: All, ColourOnly, Pruned
    # node: the current node to investigate in the network
    # position: The position to start investigation on the mRNA of the node
    # anchor_min: The allowed minimum anchor size
    # anchor_max: The allowed maximum anchor size
    # depth: How many editing blocks have there been including this one?
    # maximum_depth: Sometimes the depth has to be limited to make computation possible
    #
    # @ return
    # nothing, because the given node is changed.
    #
    def simulation_all(self, editing_type, graph_type, position, anchor_min, anchor_max, depth, maximum_depth,
                       threshold):

        # If a maximum depth is given, stop if current depth is deeper.
        if maximum_depth is not None and depth > maximum_depth:
            return

        # Continue editing while there are still nodes to be visited
        while self.to_visit:
            # print(self.to_visit)

            # Work from end of list to keep visiting depth-first
            current_node = self.to_visit.pop()
            # Track which nodes have been visited to prevent re-visiting if not needed
            self.visited.append(current_node)

            # node.gRNAs: dict of gRNAs used to create this node : position the gRNA bound
            # node.gRNAs.values(): [[bound_posn_1], [bound_posn_2, bound_posn_3]]
            # Get last editing position to prevent starting from the beginning (i.e. simulate 3'-5')
            # if loop accounts for root, where no gRNAs list yet
            if current_node.latest_start is not None:
                position = current_node.latest_start + 1

            # Find all gRNAs that can bind and edit this mRNA
            # Stops when no base is left after the anchor alignment on the mRNA.
            for gRNA in self.gRNAs:

                # Set current position to initial position for every gRNA
                current_position = position
                # 2nd loop for: Go through all gRNAs and all valid positions (thus checking for all possible matches).
                # Stops when no base is left after the anchor alignment on the mRNA.
                # TODO: decide if len - anchor_min is correct way to go about this
                while current_position < len(current_node.mRNA) - anchor_max:

                    # Check if current gRNA and mRNA match at current position
                    anchor_size = AnchorAligners.perfect(anchor_min, anchor_max, current_node.mRNA, gRNA[1],
                                                         current_position)
                    if anchor_size:

                        # Do the correct type of editing
                        if editing_type == "Uinsert":
                            # edit
                            (new_mRNA, position_after_editing) = Editors.overhang_uinsert(current_node.mRNA, gRNA[1],
                                                                                          current_position, anchor_size)
                            # add node and go to next editing block
                            self._add_child(new_mRNA, current_node, gRNA, current_position=current_position,
                                            position_after_editing=position_after_editing, thermo_score=None,
                                            graph_type=graph_type)

                        elif editing_type == "Mismatch":
                            # edit
                            (new_mRNA, position_after_editing) = Editors.overhang_ignore_mismatch(current_node.mRNA,
                                                                                                  gRNA[1],
                                                                                                  current_position,
                                                                                                  anchor_size)
                            # add node and go to next editing block
                            self._add_child(new_mRNA, current_node, gRNA, current_position=current_position,
                                            position_after_editing=position_after_editing, thermo_score=None,
                                            graph_type=graph_type)

                        elif editing_type == "Falloff":
                            # edit
                            (new_mRNA, position_after_editing) = Editors.overhang_falloff(current_node.mRNA, gRNA[1],
                                                                                          current_position, anchor_size)
                            # add node and go to next editing block
                            self._add_child(new_mRNA, current_node, gRNA, current_position=current_position,
                                            position_after_editing=position_after_editing, thermo_score=None,
                                            graph_type=graph_type)

                        elif editing_type == "SlidingThermo":
                            # edit
                            # new_mRNAs: {<mRNA sequence> : [position_after_editing, thermodynamic_probability]}
                            new_mRNAs = Editors.middle_slidingthermo_editing(current_node.mRNA, gRNA[1],
                                                                             current_position,
                                                                             anchor_size, self.windowsize, 0,
                                                                             self.mismatches_allowed, self.value_GC,
                                                                             self.value_AU, self.value_GU,
                                                                             self.value_mismatch)

                            # Thermodynamic approach can generate multiple mRNAs with equivalent scores.
                            # Thus, all of them are used.
                            # details: [position_after_editing, thermo_score]
                            for new_mRNA, details in new_mRNAs.items():
                                # u or U is no difference. Therefore, legend only stores U sequences to avoid redundancy.
                                # IDEA: Might be interesting: do such sequences show a higher abundance in real data?
                                # new_mRNA = new_mRNA.upper()
                                # add node and go to next editing block
                                self._add_child(new_mRNA, current_node, gRNA, details[1], current_position,
                                                details[0], graph_type)

                    # go to next position (while loop, all approaches)
                    current_position += 1
            position += 1

            current_node.normalise_weights()

            # only add the highest probability white node
            max_P = 0
            best_white_node = None
            for edge in current_node.outgoing_edges.values():
                if not edge.target.green and not edge.target.orange:
                    if edge.target.probability > max_P:
                        max_P = edge.target.probability
                        best_white_node = edge.target

            if best_white_node:
                #print("white", best_white_node)
                if best_white_node not in self.to_visit and best_white_node not in self.visited and best_white_node.probability > threshold:
                    self.to_visit.append(best_white_node)
                    self.best_whites.append(best_white_node)
                elif best_white_node not in self.to_visit and best_white_node not in self.visited:
                    self.best_whites.append(best_white_node)

            to_remove = []
            for name, edge in current_node.outgoing_edges.items():
                if edge.target.green or edge.target.orange:
                #print("before", edge.target.green, edge.target.orange)
                    if edge.target.probability > threshold:
                        if edge.target not in self.to_visit and edge.target not in self.visited:
                            #print("colour test", edge.target.green, edge.target.orange)
                            self.to_visit.append(edge.target)
                            assert edge.target.green or edge.target.orange
                # Remove excess white nodes so graph isn't too cluttered
                else:
                    #print(edge.target.name, edge.target in self.network.nodes, type(edge.target.name))
                    if edge.target not in self.best_whites:
                        if edge.target.name in self.network.nodes:
                            self.network.remove_node(edge.target)
                        to_remove.append(name)
                        #print("remove", to_remove)

            for white_edge in to_remove:
                #print(len(current_node.outgoing_edges))
                del current_node.outgoing_edges[white_edge]



    # Simulate possible edits with only allowing new gRNAs to bind in areas that have just been edited or were anchor.
    def simulation_editing_area_only(self, editing_type, graph_type, position, anchor_min, anchor_max, depth,
                                     maximum_depth, threshold, end_of_editingarea=None):

        # If a maximum depth is given, stop if current depth is deeper.
        if maximum_depth is not None and maximum_depth < depth:
            return

        # Continue editing while there are still nodes to be visited
        while self.to_visit:
            # print(self.to_visit)
            # Work from end of list to keep visiting depth-first
            current_node = self.to_visit.pop()
            end_of_editingarea = current_node.edited_until
            # Track which nodes have been visited to prevent re-visiting if not needed
            self.visited.append(current_node)

            # TODO: is this the error?
            # node.gRNAs: dict of gRNAs used to create this node : position the gRNA bound
            # node.gRNAs.values(): [[bound_posn_1], [bound_posn_2, bound_posn_3]]
            # Get last editing start position to prevent starting from the beginning (i.e. simulate 3'-5')
            #if current_node.latest_start is not None:
            #    position = current_node.latest_start + 1 
            if current_node.gRNAs.values():
                position = min(min(current_node.gRNAs.values())) + 1

            print("next_round")

            # 1st loop for: Go through all gRNAs and all valid positions (thus checking for all possible matches).
            # Stops when no base is left after the anchor alignment on the mRNA.
            for gRNA in self.gRNAs:
                # Set current position to initial position for every gRNA
                current_position = position

                # 2nd loop for: Go through all gRNAs and all valid positions(thus checking for all possible matches).
                # valid positions are all positions that were just edited or anchor.
                # Stops when no base is left after the anchor alignment on the mRNA.
                # TODO: is anchor_min correct?
                while current_position < len(current_node.mRNA) - anchor_max and current_position < end_of_editingarea:

                    # Check if current gRNA and mRNA match at current position
                    anchor_size = AnchorAligners.perfect(anchor_min, anchor_max, current_node.mRNA, gRNA[1],
                                                         current_position)
                    if anchor_size:

                        # Do the correct type of editing
                        if editing_type == "Uinsert":
                            # edit
                            (new_mRNA, position_after_editing) = Editors.overhang_uinsert(current_node.mRNA, gRNA[1],
                                                                                          current_position, anchor_size)
                            # add node and go to next editing block
                            self._add_child(new_mRNA, current_node, gRNA, current_position=current_position,
                                            position_after_editing=position_after_editing, thermo_score=None,
                                            graph_type=graph_type)

                        elif editing_type == "Mismatch":
                            # edit
                            (new_mRNA, position_after_editing) = Editors.overhang_ignore_mismatch(current_node.mRNA,
                                                                                                  gRNA[1],
                                                                                                  current_position,
                                                                                                  anchor_size)
                            # add node and go to next editing block
                            self._add_child(new_mRNA, current_node, gRNA, current_position=current_position,
                                            position_after_editing=position_after_editing, thermo_score=None,
                                            graph_type=graph_type)

                        elif editing_type == "Falloff":
                            # edit
                            (new_mRNA, position_after_editing) = Editors.overhang_falloff(current_node.mRNA, gRNA[1],
                                                                                          current_position, anchor_size)
                            # add node and go to next editing block
                            self._add_child(new_mRNA, current_node, gRNA, current_position=current_position,
                                            position_after_editing=position_after_editing, thermo_score=None,
                                            graph_type=graph_type)

                        elif editing_type == "SlidingThermo":
                            # edit
                            # new_mRNAs: {<mRNA sequence> : [position_after_editing, thermodynamic_score]}
                            new_mRNAs = Editors.middle_slidingthermo_editing(current_node.mRNA, gRNA[1],
                                                                             current_position,
                                                                             anchor_size, self.windowsize, 0,
                                                                             self.mismatches_allowed, self.value_GC,
                                                                             self.value_AU, self.value_GU,
                                                                             self.value_mismatch)
                            # print(new_mRNAs)

                            # Thermodynamic approach can generate multiple mRNAs with equivalent scores.
                            # Thus, all of them are used.
                            #print(len(new_mRNAs), "edits returned")
                            for new_mRNA, details in new_mRNAs.items():
                                # u or U is no difference. Therefore, legend stores U sequences to avoid redundancy.
                                # IDEA: Might be interesting: do such sequences show a higher abundance in real data?
                                # new_mRNA = new_mRNA.upper()
                                # add node and go to next editing block
                                self._add_child(new_mRNA, current_node, gRNA, details[1], current_position,
                                                details[0], graph_type)

                    # go to next position (while loop, all approaches)
                    current_position += 1
            position += 1

            current_node.normalise_weights()

            # only add the highest probability white node
            max_P = 0
            best_white_node = None
            print(len(current_node.outgoing_edges))
            for edge in current_node.outgoing_edges.values():
                if not edge.target.green and not edge.target.orange:
                    if edge.target.probability > max_P:
                        max_P = edge.target.probability
                        best_white_node = edge.target

            if best_white_node:
                # Allow editing to continue if white node exceeds threshold and has not been seen before.
                if best_white_node not in self.to_visit and best_white_node not in self.visited and best_white_node.probability > threshold:
                    self.to_visit.append(best_white_node)
                    self.best_whites.append(best_white_node)
                # Allow white node to be included but not edited further if below threshold.
                elif best_white_node not in self.to_visit and best_white_node not in self.visited:
                    self.best_whites.append(best_white_node)

            # All other white nodes and their incoming edges need to be removed.
            edges_to_remove = []
            for name, edge in current_node.outgoing_edges.items():
                # Green and orange always allowed, provided they exceed the threshold.
                if edge.target.green or edge.target.orange:
                    if edge.target.probability > threshold:
                        if edge.target not in self.to_visit and edge.target not in self.visited:
                            self.to_visit.append(edge.target)
                            assert edge.target.green or edge.target.orange
                else:
                    if edge.target not in self.best_whites:
                        # some edges point to same node so may have already been removed - prevents key not found error
                        if edge.target.name in self.network.nodes:
                            self.network.remove_node(edge.target)
                        edges_to_remove.append(name)

            for white_edge in edges_to_remove:
                #print(len(current_node.outgoing_edges))
                del current_node.outgoing_edges[white_edge]

            for edge_name in list(current_node.outgoing_edges.keys()):
                assert current_node.outgoing_edges[edge_name].target.name in self.network.nodes
                if current_node.outgoing_edges[edge_name].target.name not in self.network.nodes:
                    del current_node.outgoing_edges[edge_name]
                    print('deletion')

    # TODO: description and assertions
    def _add_child(self, new_mRNA, node, gRNA, thermo_score, current_position, position_after_editing, graph_type):

        # Check if the just generated mRNA already exists in the network
        existing_node = self.network.check_mRNA_upper(new_mRNA)

        if existing_node:
            # Add edge to current Node that points to existing Node.
            node.add_edge(existing_node, gRNA[0], thermo_score)

            # TODO: increment weight, check if was terminal and now not (i.e. A)
            # TODO: edge case A
            # Tell the existing Node that it was created using the current gRNA
            existing_node.add_gRNA(gRNA[0], current_position)
            existing_node.update_edited_position(current_position, position_after_editing)

            # TODO: what if it's orange from one path and green from another?
            if not existing_node.orange:

                # Assign green if editing results in correct sequence over full len(gRNA)
                if new_mRNA[:position_after_editing].upper() == self.mRNA_edited[:position_after_editing].upper():
                    existing_node.green = True

                # If this gRNA is used at the correct position, tell the node that it is orange.
                # If correctly edited until this block (i.e. block itself contains mis-edit/s) but gRNA at correct posn
                # Nodes can be green and/or orange - if both, green will be used when visualising
                if new_mRNA[:current_position].upper() == self.mRNA_edited[:current_position].upper() and \
                        not existing_node.green and \
                        current_position == len(self.mRNA_edited) - int(gRNA[0].split("-")[1].split(")")[0]):
                    if gRNA[0].find(self.gene_name) > -1:
                        existing_node.orange = True

            assert existing_node.name in self.network.nodes

        # If it does not exist yet, introduce a new node (depending on graph_type)
        else:
            # If mRNA is known in the legend, use name specified there.
            # If not known, add it to the legend and the legend file.
            if new_mRNA.upper() in self.legend.keys():
                new_node_name = self.legend[new_mRNA.upper()]
            else:
                new_node_name = str(len(self.legend) - 1)
                self.legend[new_mRNA.upper()] = new_node_name
                # stored with small us!
                # TODO:
                with open(self.nodes_filepath, "a") as f:
                    f.write(new_node_name + "," + new_mRNA + "\n")

            # Check colour of new node
            green = new_mRNA[:position_after_editing].upper() == self.mRNA_edited[:position_after_editing].upper()
            # NB: Colour checking relies on gRNAs following the correct naming convention!
            # Prevent non-cognate gRNAs creating orange nodes because they happen to bind at the position in their name
            if gRNA[0].find(self.gene_name) > -1:
                #print(gRNA[0], self.gene_name)
                orange = new_mRNA[:current_position].upper() == self.mRNA_edited[:current_position].upper() and current_position == len(self.mRNA_edited) - int(gRNA[0].split("-")[1].split(")")[0])
            else:
                orange = False

            # Create Node objects (always allowed for coloured nodes).
            if green:
                new_node = BF_Network_white.Node(name=new_node_name, mRNA=new_mRNA, gRNA=gRNA[0], latest_start=current_position,
                                           edited_until=position_after_editing, gRNA_position=current_position,
                                           green=True)
                # Printing feedback to command line that something was found
                print("green node found:", new_node_name)
            elif orange:
                new_node = BF_Network_white.Node(name=new_node_name, mRNA=new_mRNA, gRNA=gRNA[0], latest_start=current_position,
                                           edited_until=position_after_editing, gRNA_position=current_position,
                                           orange=True)
                # Printing feedback to command line that something was found
                print("orange node found:", new_node_name)
                # Create Node object if white nodes are allowed by graph_type
            else:
                if graph_type != "ColourOnly":
                    # print("white node found:", new_node_name)
                    new_node = BF_Network_white.Node(name=new_node_name, mRNA=new_mRNA, gRNA=gRNA[0], latest_start=current_position,
                                               edited_until=position_after_editing, gRNA_position=current_position)

            # Add new node to network (if allowed by graph_type)
            if green or orange or graph_type != "ColourOnly":
               # print("before", new_node.name in self.network.nodes, new_node.name)
                self.network.add_node(new_node)
                #print("after", new_node.name in self.network.nodes, new_node.name)
                # Add edge from current node to this new node
                node.add_edge(new_node, gRNA[0], thermo_score)
                assert new_node.name in self.network.nodes
                assert new_node in list(self.network.nodes.values())
