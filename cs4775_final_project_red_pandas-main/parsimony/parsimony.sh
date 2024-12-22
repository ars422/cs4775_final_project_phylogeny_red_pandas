#!/usr/bin/env bash
# In the terminal, execute `bash parsimony/parsimony.sh to run all 3 trees.

# Joshi et al.
python3 parsimony.py -m mapping_haplotypes.txt -f clustalw_slow.fasta
echo -e "----------------------------------------------\n"

# NJ
python3 parsimony.py -m mapping_haplotypes.txt -f clustalw_slow.fasta -nwk NJ.nwk
echo -e "----------------------------------------------\n"

# Weighbor
python3 parsimony.py -m mapping_haplotypes.txt -f clustalw_slow.fasta -nwk WNJ.nwk
echo -e "----------------------------------------------\n"




