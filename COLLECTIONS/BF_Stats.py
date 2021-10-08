from STEP.COLLECTIONS import BF_Network
import os
import numpy as np

'''
This file stores the Stats_Generator class, which can read graphs from .DATA.txt files and work with them. 
Other classes might be added if they are necessary to generate more complex stats. 
'''


# A class to read a .DATA.txt graph file and generate some meaningful statistics.
class Stats_Generator():
    def __init__(self):

        # The graph to investigate.
        self.graph = None

    # Reading the graph from a .DATA.txt file.
    def read_graph(self, filepath_graph):
        assert os.path.exists(filepath_graph), "Please provide a path to an existing graph DATA file"

        # Create graph for generating stats
        self.graph = BF_Network.Network(name="Graph to investigate")

        # Reading the provided file into the graph

        # splitting into node-entries
        with open(filepath_graph) as f:
            nodes = f.read()
            nodes = nodes.split(">")
            # getting rid of first empty entry
            nodes.pop(0)
            # separate the information about a node
            nodes = [node.strip().split("\n") for node in nodes]

        # index-legend:
        # 0: node-name
        # 1: mRNA-sequence
        # 2: gRNAs that generated this node. {<gRNA-name>: <gRNA-position>}
        # 3: edges: {<edge-name>: <occurrences>}
        # 4: green
        # 5: orange

        # Make strings to dicts / lists / booleans (as eval is used, this is a security risk).
        # List of lists. Each node is contained its own list, within the outer nodes list.
        nodes = [[node[0], node[1], eval(node[2]), eval(node[3]), eval(node[4].split(": ")[1]),
                  eval(node[5].split(": ")[1])] for node in nodes]
        # print(nodes)

        # Add nodes and their gRNAs to graph. Before edges can be added, all nodes must be stored.
        # Node[0] is the node's name, in convention number_in_legend>_p<probability_of_node>
        for node in nodes:
            self.graph.add_node(BF_Network.Node(name=node[0].split("_p")[0], mRNA=node[1],
                                                green=node[4], orange=node[5]))
            self.graph.nodes[node[0].split("_p")[0]].gRNAs = node[2]
            self.graph.nodes[node[0].split("_p")[0]].probability = node[0].split("_p")[1]

        # Add edges to nodes:
        for node in nodes:

            # Go through all edges
            for edge_name, edge_occurrences in node[3].items():

                # Information about edge is in its name
                edge_content = edge_name.split("^")
                edge_target = self.graph.nodes[edge_content[2]]
                edge_gRNA = edge_content[4]
                #edge_anchor_weight = edge_content[6].split("_")[-2]
                edge_thermo_weight = edge_content[5].split("_")[-1]
                #edge_abundance_weight = edge_content[4].split("_")[-1]
                edge_probability = edge_content[6].split("_p")[1]

                # Add edge to node as often as necessary
                i = 0
                while i < edge_occurrences:
                    self.graph.nodes[node[0].split("_p")[0]].add_edge(edge_target, edge_gRNA, edge_thermo_weight)
                    i += 1

    # Return the last correctly edited position (in python indexing) of best generated mRNA and the node associated
    # with this position
    #
    # @ param
    # fully_edited_mRNA: 3' to 5' !
    def get_correct_coverage_position(self, fully_edited_mRNA):
        assert self.graph is not None, "You must read a graph first."

        # The best position
        correct_position_maximum = 0

        for node in self.graph.get_nodes():

            # If node has been edited and is orange or green,
            # get the "worst" position of all gRNAs generating this node.
            if len(node.gRNAs) != 0 and (node.orange or node.green):
                correct_until = min(list(node.gRNAs.values())[0])
            else:
                correct_until = 0

            # Compare positions after correct_until
            while correct_until < len(fully_edited_mRNA):
                # have to compare from 3' to 5', so [::-1]
                if node.mRNA[::-1][correct_until].upper() == fully_edited_mRNA[::-1][correct_until].upper():
                    correct_until += 1
                else:
                    break

            # Compare found value to current correct_position_maximum
            if correct_until > correct_position_maximum:
                correct_position_maximum = correct_until

        # return maximum
        return correct_position_maximum

    # Return up to which percentage the correct mRNA was generated.
    #
    # @ param
    # fully_edited_mRNA: 3' to 5' !
    def get_correct_coverage_percentage(self, fully_edited_mRNA):
        assert self.graph is not None, "You must read a graph first."

        # Get the last correctly edited position and divide it by the length of the correct edited mRNA.
        return self.get_correct_coverage_position(fully_edited_mRNA) * 100 / len(fully_edited_mRNA)

    def get_most_correct_node(self, fully_edited_mRNA):
        assert self.graph is not None, "You must read in a graph first."

        # The best position
        correct_position_maximum = 0
        most_correct_node = None

        for node in self.graph.get_nodes():

            # If node has been edited and is orange or green,
            # get the "worst" position of all gRNAs generating this node.
            if len(node.gRNAs) != 0 and (node.orange or node.green):
                correct_until = min(list(node.gRNAs.values())[0])
            else:
                correct_until = 0

            # Compare positions after correct_until
            while correct_until < len(fully_edited_mRNA):
                # have to compare from 3' to 5', so [::-1]
                if node.mRNA[::-1][correct_until].upper() == fully_edited_mRNA[::-1][correct_until].upper():
                    correct_until += 1
                else:
                    break

            # Compare found value to current correct_position_maximum
            if correct_until > correct_position_maximum:
                correct_position_maximum = correct_until
                most_correct_node = node

        # return maximum
        return most_correct_node

    def get_probability_nodes(self, fully_edited_mRNA):
        assert self.graph is not None, "You must read in a graph first."

        # The best position
        lowest_probability = 1
        lowest_probability_node = None

        highest_probability = 0
        highest_probability_node = None

        leaf_probabilities = []

        for node in self.graph.get_nodes():

            # collect Ps of terminal nodes
            if node.outgoing_edges == {}:
                leaf_probabilities.append(node.probability)

            # Compare found value to current correct_position_maximum
            if float(node.probability) < lowest_probability:
                lowest_probability = float(node.probability)
                lowest_probability_node = node

            if float(node.probability) > highest_probability and node.name != '"unedited mRNA"':
                highest_probability = float(node.probability)
                highest_probability_node = node

        #leaf_probabilities = sorted(leaf_probabilities, reverse=True)

        # return maximum, minimum, all terminal Ps in descending order
        return lowest_probability_node, highest_probability_node, leaf_probabilities

    # Return the names of all used gRNAs in a list
    def get_used_gRNAs(self):
        assert self.graph is not None, "You must read a graph first."

        # Create a set to avoid redundancy
        all_gRNAs = {""}

        # Go through all gRNAs of all nodes.
        for node in self.graph.nodes.values():
            for gRNA in node.gRNAs:
                all_gRNAs.add(gRNA)

        # Remove the initial empty entry
        all_gRNAs.remove("")

        # Change set to list and return
        return list(all_gRNAs)

    # Return the ratio of green vs orange nodes in the simulation (root is not counted)
    # Returns -1 if no oranges are found (avoiding dividing by 0)
    def get_green_to_orange_ratio(self):
        assert self.graph is not None, "You must read a graph first."

        greens = 0
        oranges = 0

        # Go through all nodes
        for node in self.graph.get_nodes():
            if node.green:
                greens += 1
            elif node.orange:
                oranges += 1

        # Avoid dividing by 0
        if oranges == 0:
            return -1

        return (greens / oranges)

    # Get the average incoming edges per node
    def get_average_incoming_edges(self):
        assert self.graph is not None, "You must read a graph first."

        total_edges = 0

        # Count all incoming edges
        for node in self.graph.get_nodes():
            if len(node.gRNAs) != 0:

                # Some gRNAs are used multiple times. As their position is stored, count these positions.
                for times_used in node.gRNAs.values():
                    total_edges += len(times_used)

        return total_edges / len(self.graph.nodes)

    # Get the average outgoing edges per node. Leaves are not included!
    def get_average_outgoing_edges(self):
        assert self.graph is not None, "You must read a graph first."

        total_edges = 0

        # count all incoming edges
        for node in self.graph.get_nodes():
            if len(node.outgoing_edges) != 0:

                # check whether an edge occurs multiple times
                for edge in list(node.outgoing_edges.values()):
                    total_edges += edge.occurrences

        return total_edges / (len(self.graph.nodes) - self.count_leaves())

    # Get the average number of different gRNAs used for incoming edges per node
    def get_average_incoming_edges_different_gRNAs(self):
        assert self.graph is not None, "You must read a graph first."

        total_edges = 0

        # Count all incoming edges
        for node in self.graph.get_nodes():
            if len(node.gRNAs) != 0:
                # Some gRNAs are used multiple times. As their position is stored, count these positions.
                total_edges += len(node.gRNAs)

        return total_edges / len(self.graph.nodes)

    # Return the average number of incoming edges for white nodes only
    # Returns -1 if no white nodes.
    def get_average_incoming_edges_white_nodes(self):
        assert self.graph is not None, "You must read a graph first."

        total_edges = 0
        white_nodes = 0

        # Count all incoming edges
        for node in self.graph.get_nodes():

            # Count white nodes
            if not (node.orange or node.green):
                if len(node.gRNAs) != 0:
                    white_nodes += 1
                    # Some gRNAs are used multiple times from different previous mRNAs.
                    # As their position is stored, count these positions.
                    for times_used in node.gRNAs.values():
                        total_edges += len(times_used)

        if white_nodes == 0:
            return -1

        return total_edges / white_nodes

    # Get the number of occurrences for different numbers of incoming edges
    # Returns a dict: {<# incoming edges>: <# occurrences>}
    def get_occurrences_incoming_edges_white_nodes(self):
        assert self.graph is not None, "You must read a graph first."

        # A dict storing the number of occurrences for different numbers of incoming edges
        occurrences = {}

        # Count all incoming edges
        for node in self.graph.get_nodes():
            # Count white nodes
            if not (node.orange or node.green) and len(node.gRNAs) != 0:
                number_of_edges = 0
                for times_used in node.gRNAs.values():
                    number_of_edges += len(times_used)

                # Add to dict
                if not number_of_edges in occurrences.keys():
                    occurrences[number_of_edges] = 1
                else:
                    occurrences[number_of_edges] += 1

        return occurrences

    # Get the number of occurrences for different numbers of incoming edges
    # Returns a dict: {<# incoming edges>: <# occurrences>}
    def get_occurrences_incoming_edges_coloured_nodes(self):
        assert self.graph is not None, "You must read a graph first."

        # A dict storing the number of occurrences for different numbers of incoming edges
        occurrences = {}

        # Count all incoming edges
        for node in self.graph.get_nodes():
            # Count white nodes
            if (node.orange or node.green) and len(node.gRNAs) != 0:
                number_of_edges = 0
                for times_used in node.gRNAs.values():
                    number_of_edges += len(times_used)

                # Add to dict
                if not number_of_edges in occurrences.keys():
                    occurrences[number_of_edges] = 1
                else:
                    occurrences[number_of_edges] += 1

        return occurrences

    # Get the average number of different gRNAs used for outgoing edges per node. Leaves not included!
    def get_average_outgoing_edges_different_gRNAs(self):
        assert self.graph is not None, "You must read a graph first."

        total_edges = 0

        # count all incoming edges
        for node in self.graph.get_nodes():
            if len(node.outgoing_edges) != 0:

                used_gRNAs = []

                # check whether an edge occurs multiple times
                for edge in list(node.outgoing_edges.values()):

                    if edge.gRNA not in used_gRNAs:
                        total_edges += 1
                        used_gRNAs.append(edge.gRNA)

        return total_edges / (len(self.graph.nodes) - self.count_leaves())

    # Return the ratio of green or oragne vs white nodes in the simulation (root is not counted)
    # Returns -1 if no whites are found (avoiding dividing by 0)
    def get_colour_to_white_ratio(self):
        assert self.graph is not None, "You must read a graph first."

        # Whites = -1 because root is always white
        whites = -1
        coloured = 0
        for node in self.graph.get_nodes():
            if node.green or node.orange:
                coloured += 1
            else:
                whites += 1
        if whites == 0:
            return -1
        return (coloured / whites)

    # Simply return # of white nodes. Root is ignored.
    # If no white nodes except root returns -1
    def get_number_of_white_nodes(self):
        assert self.graph is not None, "You must read a graph first."

        # whites = -1 because root is always white
        whites = -1
        for node in self.graph.get_nodes():
            if not (node.green or node.orange):
                whites += 1
        if whites == 0:
            return -1
        return whites

    # Simply return # of coloured nodes. Root is ignored.
    # If no coloured nodes nodes except root returns -1
    def get_number_of_coloured_nodes(self):
        assert self.graph is not None, "You must read a graph first."

        # whites = -1 because root is always white
        coloured_nodes = 0
        for node in self.graph.get_nodes():
            if (node.green or node.orange):
                coloured_nodes += 1
        if coloured_nodes == 0:
            return -1
        return coloured_nodes

    # Return the percentage of leaves among orange nodes.
    # Returns -1 if no orange nodes.
    def get_leaves_among_oranges(self):
        assert self.graph is not None, "You must read a graph first."

        total_oranges = 0
        leaf_oranges = 0

        # Go throug all nodes and increment fitting variables
        for node in self.graph.get_nodes():
            if node.orange:
                total_oranges += 1
                if len(node.outgoing_edges) == 0:
                    leaf_oranges += 1

        if total_oranges == 0:
            return -1

        return (leaf_oranges / total_oranges) * 100

    def create_node_coverage(self, fully_edited_mRNAuence, node_coverage_vector):
        assert self.graph is not None, "You must read a graph first."
        assert isinstance(fully_edited_mRNAuence, str), "Please provide a valid control mRNA."
        assert isinstance(node_coverage_vector, np.ndarray), "Please provide a valid array to store the coverage."
        assert len(fully_edited_mRNAuence) <= len(
            node_coverage_vector), "Array must be same <= control mRNA. Currently " + str(len(node_coverage_vector))

        for node in self.graph.get_nodes():

            # only investigate coloured nodes
            if node.green or node.orange:

                position = min(node.gRNAs.values())[0]

                while position < len(fully_edited_mRNAuence) and node.mRNA[-1 - position].upper() == \
                        fully_edited_mRNAuence[-1 - position].upper():
                    node_coverage_vector[position] += 1
                    position += 1

    # Return the number of nodes in the simulation
    def count_nodes(self):
        assert self.graph is not None, "You must read a graph first."

        return len(self.graph.get_nodes())

    # return the number of leaves in the simulation
    def count_leaves(self):
        assert self.graph is not None, "You must read a graph first."

        leaves = 0
        for node in self.graph.get_nodes():
            if len(node.outgoing_edges) == 0:
                leaves += 1

        return leaves

    # Return dict of gRNA : number-of-times-used
    # Total number of times gRNA was used, so if same gRNA
    def count_gRNAs_used_for_white_nodes(self):
        assert self.graph is not None, "You must read a graph first."

        gRNAs_used = {}

        for node in self.graph.get_nodes():
            if not (node.orange or node.green) and node.name != '"unedited mRNA"':
                for gRNA, times_used in node.gRNAs.items():
                    if gRNA in gRNAs_used:
                        gRNAs_used[gRNA] += len(times_used)
                    else:
                        gRNAs_used[gRNA] = len(times_used)

        return gRNAs_used

    def characterise_white_node_children(self):

        total_nodes = self.count_nodes()
        white_nodes_list = self.graph.get_white_nodes()

        # Features of interest
        terminal_nodes = 0
        green_edges = 0
        orange_edges = 0
        white_edges = 0
        total_outgoing_edges = 0
        total_children = 0
        # Prevent same node from being counted twice
        children_seen_before = []
        # Children of each colour
        greens = []
        oranges = []
        whites = []
        # white nodes leading to coloured nodes
        leads_green = []
        leads_orange = []
        for node in white_nodes_list:
            # Count terminal white nodes
            if node.outgoing_edges == {}:
                terminal_nodes += 1
            # If not terminal, collect info on the types of nodes the white nodes give rise to
            else:
                for edge in node.outgoing_edges.values():
                    total_outgoing_edges += 1
                    # Make sure children nodes are not counted twice (bc edges can lead to the same nodes)
                    if edge.target not in children_seen_before:
                        total_children += 1
                        children_seen_before.append(edge.target)
                    # for each colour, store the nodes and count the edges that can lead to these nodes
                    # nodes should only be counted once, but may have multiple incoming edges
                    if edge.target.green:
                        green_edges += 1
                        if edge.target not in greens:
                            greens.append(edge.target)
                        if edge.start not in leads_green:
                            leads_green.append(edge.start)
                    elif edge.target.orange:
                        orange_edges += 1
                        if edge.target not in oranges:
                            oranges.append(edge.target)
                        if edge.start not in leads_orange:
                            leads_orange.append(edge.start)
                    else:
                        white_edges += 1
                        if edge.target not in whites:
                            whites.append(edge.target)

        terminal_pc = terminal_nodes / len(white_nodes_list) * 100

        print("Total nodes:", total_nodes)
        print("Total white nodes:", len(white_nodes_list))
        print("Total children of white nodes:", total_children)
        print("Total outgoing edges of white nodes:", total_outgoing_edges)
        print("Percentage of white nodes nodes that are terminal:", terminal_pc, "total:", terminal_nodes)
        print("Percentage of white node children that are white nodes:", len(whites) / total_children * 100, "total:",
              len(whites))
        print("Percentage of edges leading to white node children:", white_edges / total_outgoing_edges * 100, "total:",
              white_edges)
        print("Percentage of white node children that are orange nodes:", len(oranges) / total_children * 100, "total:",
              len(oranges))
        print("Percentage of edges leading to orange node children:", orange_edges / total_outgoing_edges * 100,
              "total:",
              orange_edges)
        print("Percentage of white node children that are green nodes:", len(greens) / total_children * 100, "total:",
              len(greens))
        print("Percentage of edges leading to green node children:", green_edges / total_outgoing_edges * 100, "total:",
              green_edges)

        print("Orange nodes arising from white nodes:")
        if oranges:
            for node in oranges:
                print(node.name)
        else:
            print('None')

        print("White nodes leading to orange nodes:")
        if leads_orange:
            for node in leads_orange:
                print(node.name)
        else:
            print('None')

        print("Green nodes arising from white nodes:")
        if greens:
            for node in greens:
                print(node.name)
        else:
            print('None')

        print("White nodes leading to green nodes:")
        if leads_green:
            for node in leads_green:
                print(node.name)
        else:
            print('None')

        return terminal_nodes, total_nodes, whites, oranges, greens, white_nodes_list
