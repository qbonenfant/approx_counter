import sys

filename = sys.argv[1]

with open(filename, 'r') as f:
    seq = ""
    ident = ""
    count = ""
    for i,line in enumerate(f):
        if(i%2==0):
            print("\t".join(line[1:].rstrip("\n").split("_")))
