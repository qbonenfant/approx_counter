# coding=utf-8


import sys

fileName = sys.argv[1] if(len(sys.argv) >= 2) else "nope"
verbose = True if(len(sys.argv) >= 3 and sys.argv[2] == "v") else False
kmerCount = {}
kmerList = []
adapters = {}


def haveOverlap(seq1, seq2):
    minOverlap = min(len(seq1), len(seq2)) - 1
    if(seq1[-minOverlap:] == seq2[:minOverlap]):
        return(seq1 + seq2[minOverlap:])
    else:
        return("")


def updateCount(kCounter, key):
    if(key in kCounter.keys()):
        kCounter[key] += 1
    else:
        kCounter[key] = 1
    return(kCounter)


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
            km, nb = line.rstrip("\n").split("\t")
            kmerList.append(km)
            kmerCount[km] = nb

for km in kmerList:
    ov = km
    found = True
    used = []
    while(found):
        found = False
        for km2 in kmerList:
            if(km2 not in used):
                direct = haveOverlap(ov, km2)
                reverse = haveOverlap(km2, ov)

                if(direct != "" and reverse == ""):
                    ov = direct
                    found = True
                    used.append(km2)
                    break
                elif(reverse != "" and direct == ""):
                    ov = reverse
                    found = True
                    used.append(km2)
                    break
        if(verbose):
            print(ov)
    if(verbose):
        print("FINAL OVERLAP")
        print(ov)
    adapters = updateCount(adapters, ov)
for value, adapt in sorted([(v, k) for k, v in adapters.items()])[::-1]:
    print(adapt, value)
