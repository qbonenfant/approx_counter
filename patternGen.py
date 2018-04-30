# coding=utf-8

import random

# Pattern we will try to find
pattern = "ATGCCTCTCTATCTAGCAGCTACGACATTCTATCGCTAGCATGCTAGCA"

global alpha
alpha = set(["A", "T", "C", "G"])


def mutateSeq(seq, treshold):
    newSeq = ""
    for c in seq:
        if(random.random() <= treshold):
            mutType = random.randint(0, 2)
            if(mutType == 1):
                newSeq += list(alpha - set(c))[random.randint(0, 2)]
            elif(mutType == 2):
                newSeq += c
                newSeq += list(alpha)[random.randint(0, 3)]
        else:
            newSeq += c
    return(newSeq)


def paddSides(seq, size):
    ln = len(seq)
    cut = random.randint(0, size - ln)
    insert = ""
    for i in range(size - ln):
        insert += list(alpha)[random.randint(0, 3)]
        # insert+= "x"
    return(insert[0:cut] + seq + insert[cut:])


out = open("generatedRandomStarts10k.txt", "w")
for i in range(10000):
    out.write(">Generated_" + str(i) + "_random")
    out.write("\n")
    out.write(paddSides(mutateSeq(pattern, 0.15), 100))
    out.write("\n")
out.close()
