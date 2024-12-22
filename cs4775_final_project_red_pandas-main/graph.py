import matplotlib.pyplot as plt
from Bio import Phylo
from io import StringIO

from weightednj import read_data, neighbor_join, assemble_tree, generate_newick

def parse_distance_matrix(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    species = lines[0].strip().split('\t')[1:]
    matrix = [list(map(float, line.strip().split('\t')[1:])) for line in lines[1:]]
    return species, matrix

def create_phylo_tree(newick_str):
    tree = Phylo.read(StringIO(newick_str), "newick")
    tree.root_at_midpoint()
    
    return tree

def visualize_phylo_tree(tree):
    fig, ax = plt.subplots(figsize=(12, 9))
    def format_branch_length(clade):
        return f'{clade.branch_length:.6f}' if clade.branch_length else ''
    Phylo.draw(tree, axes=ax, branch_labels=format_branch_length)
    plt.title('Midpoint Rooted Weighted Neighbor-Joining Phylogenetic Tree', fontsize=16)
    plt.xlabel('Evolutionary Distance', fontsize=12)
    plt.ylabel('Species', fontsize=12)
    plt.show()

def main():
    distances_file = 'distancematrix.txt'
    D, mapping = read_data(distances_file)
    
    edges, branch_lengths, fake_root_node = neighbor_join(D, None) 
    mapping[fake_root_node] = ""
    tree_map = assemble_tree(fake_root_node, edges)
    nwk_str = generate_newick(fake_root_node, tree_map, branch_lengths, mapping)
    tree = create_phylo_tree(nwk_str)
    visualize_phylo_tree(tree)

if __name__ == "__main__":
    main()
