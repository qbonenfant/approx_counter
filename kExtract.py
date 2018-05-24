import sys

filename = sys.argv[1]

with open(filename, 'r') as f:
    for i,line in enumerate(f):
        seq,nb = line.rstrip("\n").split("\t")
        print(">Kmer_" +str(i)+"_"+nb)
        print(seq)            