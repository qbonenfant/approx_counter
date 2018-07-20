# coding=utf-8


import sys

fileName = sys.argv[1] if(len(sys.argv) >= 2) else "nope"
kmerList = []
adapters = {}
kmerCount = {}


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
            noErr, oneErr, twoErr, km = line.rstrip("\n").split("\t")
            kmerList.append(km)
            kmerCount[km] = (noErr, oneErr, twoErr)

dispArray = []
for km in kmerList[0:1]:
    dispArray.append(
        "FIRST\t" + "\t".join([km] + [elem for elem in kmerCount[km]]))
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
                    dispArray.append("RIGHT\t" +
                                     "\t".join([km2] + [elem for elem in kmerCount[km2]]))
                    break
                elif(reverse != "" and direct == ""):
                    ov = reverse
                    found = True
                    used.append(km2)
                    dispArray.insert(0,"LEFT\t" +
                          "\t".join([km2] + [elem for elem in kmerCount[km2]]))
                    break

    print("\t"+ov+"\t0\t1\t2")
print("\n".join(dispArray))
