# coding=utf-8


import sys

fileName = sys.argv[1] if(len(sys.argv) >= 2) else "nope"
minExact = int(sys.argv[2]) if(len(sys.argv) >= 3) else 10
kmerCount = {}
kmerList = []


def haveOverlap(seq1, seq2):
    minOverlap = min(len(seq1), len(seq2)) - 1
    return(seq1[-minOverlap:] == seq2[:minOverlap])


with open(fileName, 'r') as f:
    loop = True
    ident = ""
    seq = ""
    it = 0
    next(f)  # skipping first line
    while(loop):
        try:
            line = next(f)
        except StopIteration:
            loop = False
        else:

            noErr, oneErr, twoErr, km = line.rstrip("\n").split("\t")
            kmerCount[km] = (noErr, oneErr, twoErr)
            print(">kmer_"+ str(it) + "_" + "_".join(str(el)
                                      for el in [noErr, oneErr, twoErr]))
            print(km)
            it+=1