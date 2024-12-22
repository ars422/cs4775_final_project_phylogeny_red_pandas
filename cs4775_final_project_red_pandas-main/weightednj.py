#!/usr/bin/env python3

"""
Generates a Newick formatted tree using the Weighted Neighbor-Joining algorithm, with an outgroup species specified.

Arguments:
    -f: The path to the distance matrix file (a symmetric matrix with 0 on the diagonal).
        (Default: dist10.txt)
    -nwk: The path to the output file to store the nwk format 
        (Default: tree.nwk)
    -o: Name of the outgroup species to root the tree.

Outputs:
    A Newick formatted tree.

Example usage:
    python weighted_nj.py -f dist10.txt -o Green_monkey -nwk tree.nwk
"""

import argparse
import numpy as np
from scipy.special import erfc


def read_data(distances_file):
    with open(distances_file, "r") as f:
        lines = [l.strip().split() for l in f.readlines()]
        mapping = {i: s for i, s in enumerate(lines[0])}
        lines = [l[1:] for l in lines[1:]]
        D = {i: {} for i in range(len(lines))}
        for i, l in enumerate(lines):
            for j, sval in enumerate(l):
                D[i][j] = float(sval)
    return D, mapping

"""
    Computes the variance based on a given evolutionary distance
    Args:
        distance (float): Evolutionary distance between two species
    Returns:
        float: Calculated variance
"""
def compute_variance(distance):
    D = 0.75 * (1 - np.exp(-4 * distance / 3))
    variance = np.exp(8 * distance / 3) * D * (1 - D) / 438
    return variance

"""
    Computes the positivity score using the variance and a z-score
    Args:
        dik, djk, diP, djP, dPQ (float): Distances between nodes
        variance (float): Variance of the distance
    Returns:
        float: Positivity score
"""
def compute_positivity(dik, djk, diP, djP, dPQ, variance):
    z_score = dPQ / (2 ** 0.5 * variance ** 0.5)
    positivity = -np.log(erfc(z_score))
    return positivity

"""
    Computes the additivity score for a pair of nodes (i, j)
    Args:
        i, j (int): Indices of the nodes being evaluated
        corrected_distance (dict): Distance matrix with corrections
    Returns:
        float: Additivity score for the pair.
"""
def compute_additivity(i, j, corrected_distance):
    scores = []
    for k in corrected_distance.keys():
        if k != i and k != j:
            dik = corrected_distance.get(i, {}).get(k, None)
            djk = corrected_distance.get(j, {}).get(k, None)
            if dik is None or djk is None:
                continue 
            dij = corrected_distance.get(i, {}).get(j, None)
            if dij is None:
                continue
            diP = (dik + dij - djk) / 2
            djP = (djk + dij - dik) / 2
            dPk = (dik - djk) - (diP - djP)
            add_score = dPk ** 2
            scores.append(add_score)
    return np.sum(scores)

"""
    Combines the additivity and positivity scores to compute total cost
    Args:
        i, j (int): Indices of the nodes being evaluated
        corrected_distance (dict): Distance matrix with corrections
        variances (dict): Variance for each pair of nodes
    Returns:
        float: Total cost for the pair.
"""
def compute_combined_cost(i, j, corrected_distance, variances):
    add_score = compute_additivity(i, j, corrected_distance)
    pos_score = 0
    for k in corrected_distance.keys():
        if k != i and k != j:
            dik = corrected_distance.get(i, {}).get(k, None)
            djk = corrected_distance.get(j, {}).get(k, None)
            if dik is None or djk is None:
                continue #for missing distances
            diP = (dik + corrected_distance.get(i, {}).get(j, 0) - djk) / 2
            djP = (djk + corrected_distance.get(i, {}).get(j, 0) - dik) / 2
            dPQ = corrected_distance.get(i, {}).get(j, 0) - diP - djP
            variance = variances.get(i, {}).get(k, None)
            if variance is not None:
                pos_score += compute_positivity(dik, djk, diP, djP, dPQ, variance)
    return add_score + pos_score

"""
    Computes the total distances for each node by summing the pairwise distances
    to all other nodes in the current distance matrix
    Args:
        nodes (list): List of node indices currently in the tree
        updated_distances (dict): Distance matrix with current distances
    Returns:
        dict: A dictionary where each key is a node, and the value is the total distance
"""
def compute_total_distances(nodes, updated_distances):
    total_distances = {}
    for i in nodes:
        total_distances[i] = sum(updated_distances[i][j] for j in nodes if j != i)
    return total_distances

"""
    Selects the pair of nodes that minimizes the combined cost based on additivity
    and positivity scores
    Args:
        nodes (list): List of node indices currently in the tree
        updated_distances (dict): Distance matrix with current distances
        variances (dict): Variance matrix for all pairs of nodes
    Returns:
        tuple: The pair of nodes (i, j) with the lowest combined cost
"""
def select_best_pair_with_cost(nodes, updated_distances, variances):
    best_pair = None
    min_cost = float('inf')
    for i in nodes:
        for j in nodes:
            if i != j:
                cost = compute_combined_cost(i, j, updated_distances, variances)
                if cost < min_cost:
                    min_cost = cost
                    best_pair = (i, j)
    return best_pair

"""
    Updates the distance matrix after clustering two nodes into a new node
    Args:
        nodes (list): List of node indices currently in the tree
        updated_distances (dict): Distance matrix with current distances
        cluster_sizes (dict): Dictionary tracking the sizes of clusters
        x, y (int): Indices of the nodes being clustered
        new_node (int): Index of the new node formed from clustering
"""
def update_distances(nodes, updated_distances, cluster_sizes, x, y, new_node):
    updated_distances[new_node] = {}
    for node in nodes:
        if node not in (x, y):
            updated_distances[new_node][node] = updated_distances[node][new_node] = (
                (cluster_sizes[x] * updated_distances[x][node] + cluster_sizes[y] * updated_distances[y][node])
                / (cluster_sizes[x] + cluster_sizes[y])
            )

"""
    Implements the Weighted Neighbor-Joining algorithm to construct a phylogenetic tree
    Args:
        D (dict): Initial distance matrix
        og (int): Index of the outgroup node
    Returns:
        tuple: Edges of the tree, branch lengths, and the fake root node index
"""
def neighbor_join(D, og):
    edges = []
    updated_distances = {key: value.copy() for key, value in D.items()}
    nodes = list(updated_distances.keys())
    cluster_sizes = {node: 1 for node in nodes}
    node_count = max(nodes) + 1
    branch_lengths = {}

    while len(nodes) > 2:
        total_distances = compute_total_distances(nodes, updated_distances)
        variances = {i: {j: compute_variance(updated_distances[i][j]) for j in updated_distances[i]} for i in nodes}
        x, y = select_best_pair_with_cost(nodes, updated_distances, variances)
        delta_x = 0.5 * updated_distances[x][y] + (total_distances[x] - total_distances[y]) / (2 * (len(nodes) - 2))
        delta_y = updated_distances[x][y] - delta_x
        new_node = node_count
        node_count += 1
        cluster_sizes[new_node] = cluster_sizes[x] + cluster_sizes[y]
        branch_lengths[(new_node, x)] = delta_x
        branch_lengths[(new_node, y)] = delta_y
        update_distances(nodes, updated_distances, cluster_sizes, x, y, new_node)
        nodes.remove(x)
        nodes.remove(y)
        nodes.append(new_node)
        updated_distances.pop(x, None)
        updated_distances.pop(y, None)
        for d in updated_distances.values():
            d.pop(x, None)
            d.pop(y, None)
        edges.extend([(new_node, x), (new_node, y)])

    a, b = nodes
    branch_length = updated_distances[a][b] / 2
    branch_lengths[(a, b)] = branch_length
    branch_lengths[(b, a)] = branch_length
    edges.append((a, b))
    fake_root = node_count

    if og in nodes:
        other_node = nodes[0] if nodes[1] == og else nodes[1]
        edges.remove((a, b))
        edges.extend([(fake_root, og), (fake_root, other_node)])
        branch_lengths[(fake_root, og)] = branch_lengths[(og, fake_root)] = branch_length
        branch_lengths[(fake_root, other_node)] = branch_lengths[(other_node, fake_root)] = branch_length
    else:
        fake_root = nodes[0]

    return edges, branch_lengths, fake_root

"""
    Constructs a tree representation from the edges of the phylogenetic tree
    Args:
        fake_root (int): Index of the root node
        E (list): List of edges in the tree
    Returns:
        dict: A nested dictionary representing the tree structure
"""
def assemble_tree(fake_root, E):
    tree_structure = {}

    def build_tree(current_node, parent_node=None):
        if current_node not in tree_structure:
            tree_structure[current_node] = []
        for edge in E:
            node1, node2 = edge
            if node1 == current_node and node2 != parent_node:
                tree_structure[current_node].append(node2)
                build_tree(node2, current_node)
            elif node2 == current_node and node1 != parent_node:
                tree_structure[current_node].append(node1)
                build_tree(node1, current_node)

    build_tree(fake_root)
    return tree_structure

"""
    Generates a Newick formatted string for the tree
    Args:
        fake_root (int): Index of the root node
        tree_map (dict): Tree structure as a nested dictionary
        branch_lengths (dict): Dictionary of branch lengths
        mapping (dict, optional): Mapping of node indices to species names
    Returns:
        str: Newick formatted string representing the tree
"""
def generate_newick(fake_root, tree_map, branch_lengths, mapping=None):
    def construct_newick(node, parent=None):
        name = mapping.get(node, '')
        children = [child for child in tree_map.get(node, []) if child != parent]
        if not children:
            return name
        else:
            newick_children = []
            for child in children:
                branch_length = branch_lengths.get((node, child), 0.0)
                child_str = construct_newick(child, node)
                newick_children.append('%s:%.6f' % (child_str, branch_length))
            return '(%s)%s' % (','.join(newick_children), name if name else "")

    return construct_newick(fake_root) + ';'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Weighted Neighbor-Joining on a Set of Sequences')
    parser.add_argument('-f', action="store", dest="f", type=str, default='distancematrix.txt')
    parser.add_argument('-nwk', action="store", dest="nwk", type=str, default='tree.nwk')
    parser.add_argument('-o', action="store", dest="o", type=str, default='HQ992981.1')
    args = parser.parse_args()
    distances_file = args.f
    og_ = args.o
    nwk_ = args.nwk
    D, mapping = read_data(distances_file)
    og = dict(map(reversed, mapping.items()))[og_]
    E, branch_lengths, fake_root = neighbor_join(D, og)
    mapping[fake_root] = ""
    tree_map = assemble_tree(fake_root, E)
    nwk_str = generate_newick(fake_root, tree_map, branch_lengths, mapping)
    print(nwk_str)
    with open(nwk_, "w") as nwk_file:
        print(nwk_str, file=nwk_file)