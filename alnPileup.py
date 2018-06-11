#coding=utf-8
# Generate a pileup from an aln alignement file


import sys

fileName = sys.argv[1]
blocs = [[]]
SEP = "\t"
alpha = ["A","T","C","G"]

with open(fileName,'r') as f:
    #skipping first line
    next(f)
    loop = True
    currentBloc = 0
    while(loop):
        try:
            line = next(f).rstrip("\n")
        except StopIteration:
            loop = False
        else:
            # if we are in bloc
            if(line and line[0] != " "):
                idt,seq = line.split()
                blocs[currentBloc].append(seq)
                
            # if it is a new bloc
            elif blocs[-1] != []:
                blocs.append([])
                currentBloc +=1

# Once all sequences are loaded.
# Start counting
for nb,bloc in enumerate( [b for b in  blocs if b !=[] ]):
    
    print("BLOC NUMBER: " +str(nb))
    print(SEP.join(alpha))
    # for each bloc of sequences
    for i in range(len(bloc[0])):
        #count per nucleotides: ATCG
        count = [0,0,0,0]
        for seq in bloc:
            if(seq[i]!="-"):
                count[alpha.index(seq[i].upper())] += 1
        print(SEP.join(str(elem) for elem in count))
    print("\n")

            
            