# coding=utf-8


import sys

fileName = sys.argv[1] if(len(sys.argv) >= 2) else "nope"
minExact = int(sys.argv[2]) if(len(sys.argv) >= 3) else 10
kmerCount = {}
kmerList = []


def haveOverlap(seq1, seq2):
    minOverlap = min(len(seq1), len(seq2)) - 1
    if(seq1[-minOverlap:] == seq2[:minOverlap]):
        return(seq1 + seq2[minOverlap:])
    else:
        return("")


with open(fileName, 'r') as f:
    loop = True
    ident = ""
    seq = ""
    next(f)  # skipping first line
    while(loop):
        try:
            line = next(f)
        except StopIteration:
            loop = False
        else:
            noErr, oneErr, twoErr, km = line.rstrip("\n").split("\t")
            kmerList.append(km)
            kmerCount[km] = (noErr, oneErr, twoErr)

for km in kmerList:
    ov = km
    for km2 in kmerList:

        direct = haveOverlap(ov, km2)
        reverse = haveOverlap(km2, ov)

        if(direct != "" and reverse == ""):
            ov = direct
        elif(reverse != "" and direct == ""):
            ov = reverse

    print(ov)
