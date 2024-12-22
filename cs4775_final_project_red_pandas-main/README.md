# cs4775_final_project_red_pandas
CS 4775 Final Project on Phylogenetic Trees of Red Pandas: Alexa Sheldon (ars422), Bella Falkenberg (alf253), Megan Adams (maa367), Sandee Basco (myb22)

Find implementations:
- Weighted Neighbor Joining: weightednj.py
- Parsimony: parsimony/parsimony.py

Parsimony - Joshi et al. (2020) tree
- enter the parsimony folder
- type: python parsimony.py
- the defaults set this up

Parsimony - Other Trees
- To specify a mapping file (ie genbank accession numbers to Haplotypes) use the '-m' argument
- To specify the fasta file you want to work with (the sequences) use the '-f' argument
- To specify the Tree you want to test, create a file with the newick and then pass file as the '-nwk' argument
- We already have newicks for NJ and Weighbor titled 'NJ.nwk' and 'WNJ.nwk' respectively
- type in the parsimony folder: python parsimony.py -m <some_mapping_file> -f <some_fasta_file_with_sequences> -nwk <some_phylogenetic_tree_newick>

Parsimony - Output
- The first output will be the parsimony score of the tree you inputted.
- Then, it will output the neighbors' parsimony scores and their Newick formats
- The last parsimony score and Newick format is the tree found by hill-climbing and NNI


Weighted Neighbor-Joining

To run the algorithm: 
python weightednj.py -f <distance_matrix_file> -o <outgroup_species> -nwk <output_newick_file>

Arguments (default is to use sequences listed in paper):
- -f: Path to the distance file matrix file (default: distancematrix.txt)
- -o: Name of the outgroup species to root the tree (default: HQ992981.1)
- -nwk: Path to the output file to store the Newick formatted tree (default: tree.nwk)

Output:
A Newick formatted tree saved in the specified output file.
