import argparse
import numpy as np
from newick import loads
import time

''' Reads in the sequences from the fasta files. From A4

Arguments:
    filename: which filename to read sequences from
Returns:
    output: list of sequences in fasta file
'''
def read_fasta(filename):
    with open(filename, "r") as f:
        output = []
        ids = []
        s = ""
        for l in f.readlines():
            temp_lst = l.strip()
            if temp_lst[0] == ">":
                # skip the line that begins with ">"
                ids.append(temp_lst[1:11])
                if s == "": continue
                output.append(s)
                s = ""
            # keep appending to the current sequence when the line doesn't begin
            # with ">"
            else: s += temp_lst
        output.append(s)
        return output,ids
    
''' Reads in haplotype and ID mapping from txt file

Arguments: 
    filename: filename which has the mapping between haplotypes and IDs
Returns:
    output: map_hap: a dictionary which maps ID -> haplotype
'''
def read_map_txt(filename):
    map_hap = {}
    with open(filename, "r") as f:
      for l in f.readlines():
        temp_lst = l.split()
        map_hap[temp_lst[1]] = temp_lst[0]
      return map_hap
    
''' Maps haplotypes to their sequences

Arguments:
    filename: filename which has the sequences in it
    mapper: maps IDs to haplotypes
Returns:
    output: seq: dictionary haplotype -> sequence
'''
def read_txt(filename, mapper):
    with open(filename, "r") as f:
      seq = {}
      for l in f.readlines():
        temp_lst = l.split()
        if len(temp_lst) > 0:
          temp_key = mapper[temp_lst[0]]
          if mapper[temp_lst[0]] in seq.keys():
            seq[temp_key] += temp_lst[2]
          else:
            seq[temp_key] = temp_lst[2]
      return seq
    
''' Returns the IDs of a nodes left and right children

Arguments:
    key: the name of the parent node (string)
    t_dict: the dictionary mapping nodes to their children
    map_hap: the dictionary mapping haplotypes to their ID
Returns:
    kid1: the ID of the left child
    kid2: the ID of the right child
'''
def get_kids(key, t_dict, map_hap):
    kid1 = t_dict[key][0]
    kid2 = t_dict[key][1]

    if kid1[0] != 'I':
        kid1 = map_hap[kid1]
    if kid2[0] != 'I':
        kid2 = map_hap[kid2]

    return kid1, kid2

''' Returns the node for the given ID

Arguments:
    kid: the ID of a node (string)
    I: a list of internal nodes
    H: a list of haplotype leaf nodes
Returns:
    kid: the node corresponding to that ID
'''
def kid_node(kid, I, H):
    id = int(kid[1:])

    if kid[0] =='I':
        kid = I[id]
    else: kid = H[id]

    return kid

class Tree():
    ''' Initializes a tree with given parameters.

    Arguments:
        count: starting count for internal nodes (int)
        root: root of tree (Node)
        format_dict: each key is a label for a node mapping to a list of 
            the names of its children nodes
            (dict with keys as strings and items as lists of strings)
    '''
    def __init__(self, count, root):
        self.count = count
        self.root = root
        self.format_dict = {}

    ''' Recursively creates a format_dict from the given root node

    Arguments:  
        parent: the root node of the tree (Node from the newick library)
    '''
    def create_tree(self, parent):
        if parent != None:
            name = parent.name

            if name == None:
                name = "I" + str(self.count)

            self.format_dict[name] = []

            for n in parent.descendants:
                if n.name == None:
                    self.count = self.count + 1
                self.format_dict[name].append(self.create_tree(n))

            return name
            
        else:
            return []

class Node():
    ''' Initializes a node with given parameters. From A3 2a.py

    Arguments:
        name: name of node (only relevant for leaves)
        left: left child (Node)
        right: right child (Node)
        branch_length: length of branch that leads to this node (float)
        branch_id: id of branch that leads to this node (int)
        probs: probability of observed bases beneath this node
                [list of 4 probs for 'ACGT'] (initialized to None)
    '''
    def __init__(self, name, left, right): 
        self.name = name
        self.left = left
        self.right = right

''' Traverse a tree in post order

Arguments:
    root: the root of the tree to traverse (Node)
Returns:
    a list of nodes in post order
'''
def post_order(root):
    if root != None:
        return post_order(root.left) + post_order(root.right) + [root]
    else:
        return []

''' Maps sequences to their haplotype labels

Arguments:
    seqs: list of sequences
    ids: corresponding list of IDs
    mapper: maps ID -> haplotype
Returns:
    D: dictionary of sequences with haplotype labels as the keys (H1, H2 etc.)
'''
def sequence_dict(seqs, ids, mapper):
    D = {}

    for i in range(len(seqs)):
        D[mapper[ids[i]]] = seqs[i]
    return D

''' Calculates parsimony scores for each state for every node using Sankoff's algorithm

Arguments:
    D: dictionary of sequences
    lst: Post order traversal
    index: index of sequence being examined
Returns:
    min_score_matrix: dictionary of dictionaries: minimum score matrix (4 states per internal node)
'''
def sankoff(D, lst, index):
    char_list = ['A', 'C', 'G', 'T', '-']
    min_score_matrix = {}

    # Traverse through nodes in post order
    for node in lst:
        min_score_matrix[node.name] = {'A' : np.inf, 'C' : np.inf, 'G' : np.inf, 'T' : np.inf, '-' : np.inf}

        # Check if the node is a leaf
        if node.left is None and node.right is None:
            # Set the value for its state to be 0
            char = D[node.name][index]
            min_score_matrix[node.name][char] = 0
        else: 
            # recurrence relation: s(i) = min_j[c_ij + s_l(j)] + min_k[c_ik + s_r(k)]
            for char in char_list:
                # Get min left sum if there is a left child
                min_left = 0
                if node.left is not None:
                    min_left = np.inf
                    for left_char in char_list:
                        left_sum = min_score_matrix[node.left.name][left_char]
                        if char != left_char: left_sum += 1
                        min_left = min(min_left, left_sum)

                # Get min right sum if there is a right child
                min_right = 0
                if node.right is not None:
                    min_right = np.inf
                    for right_char in char_list:
                        right_sum = min_score_matrix[node.right.name][right_char]
                        if char != right_char: right_sum += 1
                        min_right = min(min_right, right_sum)

                min_score_matrix[node.name][char] = min_left + min_right

    return min_score_matrix

''' Backtraces from Sankoff's algorithm to assign states to each node

Arguments:
    min_score_matrix: a dictionary of node name -> bases -> scores
    root: the node that refers to the top of the tree
Returns: 
    state_matrix: node name -> min parsimony state
    matrix_only_min: dictionary node name -> min score
'''
def backtrace(min_score_matrix, root):
    state_matrix = {}
    matrix_only_min = {}
    matrix_only_min['I1'] = min(min_score_matrix['I1'].values())
    char_list = ['A', 'C', 'G', 'T', '-']

    state = None
    min_score = np.inf
    score_dict = min_score_matrix[root.name]

    # Get state with minimum parismony score for the root node
    for char in score_dict:
        curr_score = score_dict[char]
        if curr_score < min_score:
            min_score = curr_score
            state = char
        
    state_matrix[root.name] = state
    # Backtrace on the rest of the tree
    def recursive_call(node, score, char):
        # Loop through pairs of characters for children
        for left_char in char_list:
            for right_char in char_list:
                left_score = 0
                left_cost = 0
                right_score = 0
                right_cost = 0

                # Get the score for the children 
                if node.left is not None:
                    left_score = min_score_matrix[node.left.name][left_char]
                    if left_char != char: left_cost = 1
                if node.right is not None:
                    right_score = min_score_matrix[node.right.name][right_char]
                    if right_char != char: right_cost = 1

                # Check that the sum of the children's score equals this nodes min score
                if left_score + left_cost + right_score + right_cost == score:
                    # Set scores for children and recurse
                    if node.left is not None:
                        state_matrix[node.left.name] = left_char
                        matrix_only_min[node.left.name] = left_score
                        recursive_call(node.left, left_score, left_char)

                    if node.right is not None:
                        state_matrix[node.right.name] = right_char
                        matrix_only_min[node.right.name] = right_score
                        recursive_call(node.right, right_score, right_char)
         

    recursive_call(root, min_score, state)

    return state_matrix, matrix_only_min

''' Gets total parsimony score of a given tree

Arguments:
    seq_dict: dictionary of sequences with haplotype labels as the keys (H1, H2 etc.)
    root: the root Node of a phylogenetic tree
Returns:
    total_score: the total parsimony score of a given tree

'''
def get_total_score(seq_dict, root):
    # Get root node and traverse tree in post order
    post_order_lst = post_order(root)

    min_score_matrix = {}
    state_matrix = {}
    matrix_only_min = {}

    for index in range(len(seq_dict['H1'])):
        # Run sankoff's algorithm to get the state for each node
        min_score_matrix = sankoff(seq_dict, post_order_lst, index)

        state_matrix[index],matrix_only_min[index] = backtrace(min_score_matrix, root)

    total_score = 0
    for i in range(len(matrix_only_min)):

        total_score += matrix_only_min[i][root.name]
    
    return total_score

''' Uses Nearest Neighbors Interchange to find neighboring Trees (no siblings of internal and leaves)

Argument:
    root: the root Node of a phylogenetic tree
Returns:
    neighbors: a list of Nodes which contain neighboring trees to the given tree

'''
def get_neighbors(root):
    neighbors = []
    
    # makes a copy of the tree (so we don't have to worry about static vs dynamic typing)
    def make_a_copy(node):
        if node == None:
            return None
        return Node(node.name, make_a_copy(node.left), make_a_copy(node.right))
    
    copied_root = make_a_copy(root)

    def find_sibling(root, node):
        nodes = post_order(root)
        for n in nodes:
            if n.right == node:
                return n.left
            elif n.left == node:
                return n.right
        return None
            
    def recurse_tree(node):
        if node == None or node.left == None or node.right == None:
            return
        sib = find_sibling(copied_root, node)
        
        # note: it does not account for siblings that are leaves
        if sib != None and sib.left != None and sib.right != None:
            # swap 1: switching left subtrees between siblings
            temp_sib_left = sib.left
            temp_node_left = node.left
            sib.left = temp_node_left
            node.left = temp_sib_left
            neighbors.append(make_a_copy(copied_root))
            
            #undo swap
            sib.left = temp_sib_left
            node.left = temp_node_left
            # swap 2: switching right subtrees between siblings
            temp_sib_right = sib.right
            temp_node_right = node.right
            sib.right = temp_node_right
            node.right = temp_sib_right
            neighbors.append(make_a_copy(copied_root))
        recurse_tree(node.left)
        recurse_tree(node.right)

    recurse_tree(copied_root)
    return neighbors

''' Generates the Newick of a tree without distances

Arguments:
    root: the root Node of a phylogenetic tree
    mapping: an optional mapping of Node to a different output
Returns:
    newick: a string which is the Newick format of the tree inputted
'''
# modified to not have distances
def generate_newick(root, mapping = None):
    ''' Complete this function. '''
    def display(root): #this is from discussion
        if root.right is None and root.left is None:
            if mapping == None: 
                return str(root.name)
            else: 
                return mapping[root]
        elif root.right is None or root.left is None:
            return '%s' % (generate_newick(display[root][0]))
        [left_child, right_child] = [root.left, root.right]
        return '(%s, %s)' % (display(left_child), display(right_child))
    return display(root) + ';'

''' Uses Beta hill climbing by exploring neighbors of tree. If a neighboring tree has a lower parsimony
score, then that tree is jumped to and we explore the neighbors of that one.

Arguments:
    seq_dict: dictionary of sequences with haplotype labels as the keys (H1, H2 etc.)
    tree: the root Node of a phylogenetic tree
Returns:
    tree: returns the neighbor found with the lowest parsimony score
    curr_score: returns the parsimony score found with the tree found
'''
# uses beta hill climbing by exploring neighbors
def hill_climbing(seq_dict, tree):
    curr_score = get_total_score(seq_dict, tree)
    
    # keeps finding neighbors and evaluating them
    while True:
        N = get_neighbors(tree)
        best_score = curr_score
        for neighbor in N:
            print(generate_newick(neighbor))
            temp_score = get_total_score(seq_dict, neighbor)
            print(temp_score)
            if best_score > temp_score:
                best_score = temp_score
                best_n = neighbor
                print(best_score)
        if best_score < curr_score:
            tree = best_n
            curr_score = best_score
        else:
            break

    return tree,curr_score


def main():
    # parsing in information
    parser = argparse.ArgumentParser(description = 'Sequence Stage to Fasta')
    
    parser.add_argument('-m', action="store", dest = "m", type=str, default='mapping_haplotypes.txt')
    parser.add_argument('-f', action="store", dest="f", type=str, default='clustalw_slow.fasta')
    parser.add_argument('-nwk', action="store", dest="nwk", type=str, default=None)
    args = parser.parse_args()
    map_file = args.m
    
    fasta_file = args.f
    newick_file = args.nwk
    seq,ids = read_fasta(fasta_file)
    mapper = read_map_txt(map_file)
    seq_dict = sequence_dict(seq,ids,mapper)

    # creates tree from Joshi et al.
    H = [] # indexed at 0, so it has a dummy H0 that will not be used
    I = []
    for i in range(54):
        H.append(Node("H" + str(i), None, None))
    i15 = Node("I15", H[4], H[6])
    i14 = Node("I14", i15, H[8])

    i10 = Node("I10", H[1],H[5])
    i9 = Node("I9", i10, H[2])
    i8 = Node("I8", i9, H[19])
    i5 = Node("I5", i8, H[3])

    i12 = Node("I12", H[16], H[7])
    i13 = Node("I13", H[11], H[12])
    i11 = Node("I11", i13, H[13])
    i7 = Node("I7", i12, i11)
    i6 = Node("I6", i7, H[14])
    i4 = Node("I4", H[17], i6)

    i3 = Node("I3", i5, i4)
    i2 = Node("I2", i14, i3)
    i39 = Node("I39", H[33], H[9])
    i40 = Node("I40", H[34], H[43])
    i43 = Node("I43", H[27], H[35])
    i44 = Node("I44", H[41], H[32])
    i46 = Node("I46", H[36], H[42])
    i27 = Node("I27", H[23], H[25])
    i47 = Node("I47", H[45], H[51])
    i48 = Node("I48", H[52], H[53])
    i52 = Node("I52", H[24], H[20])
    i51 = Node("I51", H[18], H[21])
    i49 = Node("I49", H[47], H[39])
    i35 = Node("I35", H[49], H[44])
    i36 = Node("I36", H[40], H[50])
    i37 = Node("I37", H[48], H[30])
    i38 = Node("I38", H[15], H[28])

    i42 = Node("I42", H[22], i43)
    i41 = Node("I41", i42, i44)
    i26 = Node("I26", i40, i41)
    i21 = Node("I21", i39, H[46])
    i45 = Node("I45", H[31], i46)
    i23 = Node("I23", H[38], i27)
    i28 = Node("I28", i47, i48)
    i50 = Node("I50", i52, H[26])
    i29 = Node("I29", i50, i51)
    i24 = Node("I24", i29, i49)
    i20 = Node("I20", i28, i24)
    i19 = Node("I19", i23, i20)

    i22 = Node("I22", i26, i45)
    i18 = Node("I18", i21, i22)
    i17 = Node("I17", i18, i19)

    i34 = Node("I34", i35, H[37])
    i33 = Node("I33", i34, i36)
    i32 = Node("I32", H[29], i33)
    i30 = Node("I30", H[10], i32)
    i31 = Node("I31", i37, i38)
    i25 = Node("I25", i30, i31)

    i16 = Node("I16", i17, i25)
    i1 = Node("I1", i2, i16)
    root = i1

    if newick_file != None:
        with open(newick_file, "r") as file: nwk = file.readline().strip()
    
        t = Tree(1, loads(nwk)[0])
        t.create_tree(t.root)
        t_dict = t.format_dict

        # H stores haplotype nodes
        H = [None] * 54
        # I stores internal nodes
        I = [None] * (len(t_dict.keys()) - 51)

        # Traverse tree to create and store nodes with IDs as names
        for key in list(reversed(t_dict.keys())):
            if key[0] == 'I':
                kid1, kid2 = get_kids(key, t_dict, mapper)
                kid1 = kid_node(kid1, I, H)
                kid2 = kid_node(kid2, I, H)
                I[int(key[1:])] = Node(key, kid1, kid2)
            else:
                hap = mapper[key]
                H[int(hap[1:])] = Node(hap, None, None)
        
        root = I[1]
    
    print(get_total_score(seq_dict, root))
    out = hill_climbing(seq_dict, root)
    print(generate_newick(out[0]))

if __name__ == "__main__":
    main()
    