import argparse

def read_map_txt(filename):
    map_hap = {}
    with open(filename, "r") as f:
      for l in f.readlines():
        temp_lst = l.split()
        map_hap[temp_lst[1]] = temp_lst[0]
      return map_hap
'''
  input:  filename: filename which has the sequences in it
          mapper: maps IDs to haplotypes
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

def main():
    parser = argparse.ArgumentParser(description = 'Sequence Stage to Fasta')
    parser.add_argument('-t', action="store", dest = "t", type=str, default='staging_sequences_to_fasta.txt')
    parser.add_argument('-m', action="store", dest = "m", type=str, default='mapping_haplotypes.txt')
    args = parser.parse_args()
    txt_file = args.t
    map_file = args.m
    mapper = read_map_txt(map_file)
    sequences = read_txt(txt_file,mapper)

if __name__ == "__main__":
    main()
