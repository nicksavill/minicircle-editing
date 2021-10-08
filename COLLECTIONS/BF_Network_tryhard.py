#import graphviz
import os, sys

"""
This file contains the classes that are needed to construct the graph. 
A Network represents the whole editing process, so it consists of Nodes and Edges
A Node represents a (intermediate or leaf) mRNA
An Edge represents an editing block
"""


# Nodes representing mRNAs.
class Node():

    def __init__(self, name=None, mRNA=None, gRNA=None, gRNA_position=None, edited_until=None,
                 orange=False, green=False):

        # Just strings
        self.name = name
        self.mRNA = mRNA

        # Edges stored in a dict: {<edge-name>: <Edge-object>}
        self.outgoing_edges = {}
        # For white node simulations - to control which edges get permanently added.
        self.temp_edges = {}

        # gRNAs that lead to this Node in a dict: {<gRNA-name>: [<gRNA-start-position>]}
        self.gRNAs = {}

        # Booleans for the colour.
        # Green = Correctly edited until last base of causing gRNA
        # Orange = Correctly edited until gRNA binding position and correct gRNA, but mis-edited in editing block.
        self.green = green
        self.orange = orange

        # Pass cumulative weights between edges
        self.probability = None

        # Track end of last editing block
        self.edited_until = edited_until

        if gRNA is not None:
            self.gRNAs[gRNA] = [gRNA_position]

    # Add an (outgoing) edge.
    def add_edge(self, target_node, new_gRNA_name, thermo_probability=None, normalised_probability=None):
        assert isinstance(target_node, Node), \
            "Please provide a valid Node object as target for the new edge. Currently: " + str(type(target_node))
        assert isinstance(new_gRNA_name, str), \
            "Please provide a valid gRNA name String for the new edge. Currently: " + str(type(new_gRNA_name))

        new_edges_name = self.name + "^to^" + target_node.name + "^with^" + new_gRNA_name + "^thermo_probability" + \
                         "_" + str(thermo_probability)

        # If edge already exists increase number of occurrences. If not, add edge.
        # TODO: update weight of edge?
        if new_edges_name in list(self.outgoing_edges.keys()):
            new_edge = self.outgoing_edges[new_edges_name]
            new_edge.increment_occurrence()
        else:
            self.outgoing_edges[new_edges_name] = Edge(start=self, target=target_node, gRNA=new_gRNA_name,
                                                       thermo_probability=thermo_probability, normalised_probability=normalised_probability)

    # Add a temporary edge (for white node simulations, so can control how many edges added)
    def add_temp_edge(self, target_node, new_gRNA_name, thermo_probability=None):
        assert isinstance(target_node, Node), "Please provide a valid Node object as target for new edge. Currently:" + str(type(target_node))
        assert isinstance(new_gRNA_name, str), "Please provide a valid String name for the gRNA. Currently:" + str(type(new_gRNA_name))

        temp_edge_name = self.name + "^to^" + target_node.name + "^with^" + new_gRNA_name

        if temp_edge_name in list(self.temp_edges.keys()):
            temp_edge = self.temp_edges[temp_edge_name]
        else:
            self.temp_edges[temp_edge_name] = Edge(start=self, target=target_node, gRNA=new_gRNA_name, thermo_probability=thermo_probability)

    # TODO: add position after editing into gRNAs.
    # add position into add_gRNA
    # modify all code to reflect new structure of node.gRNAs
    # modify editingarea to use this to restrict editing

    # Add a new gRNA that can lead to this Node
    def add_gRNA(self, gRNA_name, gRNA_position):
        assert isinstance(gRNA_name, str), "Please provide a valid gRNA name String to add. Currently: " + str(
            type(gRNA_name))
        assert isinstance(gRNA_name, str), "Please provide a valid gRNA position String to add. Currently: " + str(
            type(gRNA_position))

        # If gRNA already exists, add the new position for this gRNA. If not, add gRNA.
        if gRNA_name in self.gRNAs.keys():
            self.gRNAs[gRNA_name].append(gRNA_position)
        else:
            self.gRNAs[gRNA_name] = [gRNA_position]

    #
    def update_edited_position(self, new_edited_until):
        self.edited_until = new_edited_until

    # Use weights of all daughter edges to normalise and find probabilities
    def normalise_weights(self):
        thermo_total = 0
        for edge in self.temp_edges.values():
            thermo_total += edge.thermo_probability
        for edge in self.temp_edges.values():
            edge.normalised_probability = edge.thermo_probability / thermo_total
            edge.name = edge.name + "_normalised_p" + str(edge.normalised_probability)
            # print(edge.probability)
            #if edge.target.green or edge.target.orange:
            assigned_probability = edge.normalised_probability * self.probability
            #else:
            #    assigned_probability = edge.normalised_probability * self.probability * 0.1  # penalise white nodes
            # update probability of target node
            if edge.target.probability is None:
                edge.target.probability = assigned_probability
            else:
                edge.target.probability += assigned_probability
            # print(edge.target.name, "probability assigned", edge.probability * self.probability)


# Edges between Nodes representing editing blocks.
class Edge():
    def __init__(self, start=None, target=None, gRNA=None, thermo_probability=None, normalised_probability=None):
        assert start is None or isinstance(start, Node), "Please provide a valid start Node for the new edge."
        assert target is None or isinstance(target, Node), "Please provide a valid target Node for the new edge."
        assert gRNA is None or isinstance(gRNA, str), "Please provide a valid gRNA String for the new edge."

        # Node object where edge comes from
        self.start = start
        # Node object where edge leads to
        self.target = target
        # Used gRNA name
        self.gRNA = gRNA
        # How often is this gRNA used?
        self.occurrences = 1
        # Name of edge
        self.name = self.start.name + "^to^" + self.target.name + "^with^" + self.gRNA + "^thermo_probability" + \
                    "_" + str(thermo_probability)
        # Weights and probability
        self.thermo_probability = thermo_probability
        self.normalised_probability = normalised_probability

    # Increasing the number of occurrences
    def increment_occurrence(self):
        self.occurrences += 1


# The network that represents the editing process in the cell.
class Network():
    def __init__(self, name=None, nodes=None):
        assert nodes is None or isinstance(nodes, list), \
            "Only a list of Node object is allowed to be added. Currently: " + str(type(nodes))
        if nodes is not None:
            for node in nodes:
                assert isinstance(node, Node), \
                    "Only a list of Node object is allowed to be added. At least one entry is: " + str(type(node))

        # Name of the Network (for saving files)
        self.name = name

        # All nodes of the Network in a dict: {<node-name>: <Node-object>}
        self.nodes = {}

        # If a list of nodes is given, add it.
        if nodes is not None:
            for node in nodes:
                self.add_node(node)

    # Add a new node
    def add_node(self, node):
        assert isinstance(node, Node)

        self.nodes[node.name] = node

    # Remove an exising node - specifically for white node simulations
    def remove_node(self, node):
        assert isinstance(node, Node)

        del self.nodes[node.name]

    # Check if a specific mRNA exists in the network already.
    # If yes, return the Node, if not, return None.
    def check_mRNA(self, mRNA):
        for node in self.get_nodes():
            if node.mRNA.upper() == mRNA.upper():
                return node
        return None

    # TODO: Check this wasn't needed for one of the stats functions. If not, fix return.
    # Check if a given mRNA exists in upper() already in the Network.
    # Returns True/False
    def check_mRNA_upper(self, mRNA):
        for node in self.get_nodes():
            if node.mRNA.upper() == mRNA.upper():
                return node
        return False

    # return a list of all nodes (Node objects)
    def get_nodes(self):
        return list(self.nodes.values())

    # return a list of all node names (str objects)
    # TODO: decide if want to just names = [name for name in self.nodes.keys()] - surely they're the same names?
    def get_node_names(self):
        names = [node.name for node in self.nodes.values()]
        return names

    # save the graph as .gv and .png and display it.
    def graphviz_windows(self, directorypath_output):
        assert isinstance(directorypath_output, str) and os.path.exists(
            directorypath_output), "Please provide a valid output directory."

        # Neccessary for windows
        os.environ["PATH"] += os.pathsep + 'C:/Program Files (x86)/Graphviz2.38/bin/'

        # Create graph in PNG format
        graph = graphviz.Digraph(name=self.name, format="png")

        # Go through all nodes and add them to the graph
        for node in self.get_nodes():
            if node.green:
                graph.node(name=node.name, label=node.name + "_p" + str(node.probability), style="filled",
                           color="green")
            elif node.orange:
                graph.node(name=node.name, label=node.name + "_p" + str(node.probability), style="filled",
                           color="orange")
            else:
                graph.node(name=node.name, label=node.name + "_p" + str(node.probability))

            # add edges of the node to the graph
            for edge in list(node.outgoing_edges.values()):
                if edge.occurrences > 1:
                    graph.edge(node.name, edge.target.name, label=edge.gRNA.split("_")[1] + "x" + str(edge.occurrences)
                               + "_p" + str(edge.probability))
                else:
                    graph.edge(node.name, edge.target.name, label=edge.gRNA.split("_")[1] +
                               "_p" + str(edge.probability))

        # Saving graph to given directory. Then displaying it with standard png-program.
        graph.render(directory=directorypath_output)
        # graph.view(directory=directorypath_output)

    # Saving the data in this network to a file within the given directory.
    def save_data(self, directorypath_output):
        assert isinstance(directorypath_output, str) and os.path.exists(
            directorypath_output), "Please provide a valid output directory."

        # open file and change output from command line to it.
        original_stdout = sys.stdout
        outfile = open(directorypath_output + "/" + self.name + ".DATA.txt", "w")
        sys.stdout = outfile

        # All information is within nodes (and edges, which are stored in nodes).
        for node in self.get_nodes():

            # Head of each entry: Name of node
            print(">" + node.name + "_p" + str(node.probability))

            # The sequence of the node ***in 5' to 3' direction***. Thus, has to be reversed!
            print(node.mRNA[::-1])

            # The gRNAs that lead to this node. Format: {<gRNA-name>, <number-of-occurrences>}
            print(node.gRNAs)

            # The outgoing edges of this node. Each edges has four characteristics.
            # Four of them are stored in the edge's name: <startname^to^targetname^with^gRNAname^weight^weight>
            # This can easily be split by ^. Thus:
            # {<name-of-edge>: occurrences}
            edges = {}
            for edge in list(node.outgoing_edges.values()):
                edges[edge.name] = edge.occurrences
            print(edges)

            # colour status of this node
            print("green:", node.green)
            print("orange:", node.orange)

        # close file and change output to command line again.
        outfile.close()
        sys.stdout = original_stdout
