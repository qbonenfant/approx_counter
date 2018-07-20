# coding = utf-8
# Compare two kmer count file.
# Displa intersection, difference and unic kmers.

import sys

ref = sys.argv[1]
test = sys.argv[2]


def parseFasta(filename):
    storeSet = set()
    with open(filename) as f:
        seq = ""
        ident = ""
        for line in f:
            if(line[0] == ">"):
                if(seq == ""):
                    ident = line[:-1]
                else:
                    storeSet.add(seq)
                    seq = ""
                    ident = line[:-1]
            else:
                seq += line[:-1]
    storeSet.add(seq)
    return(storeSet)


refKmers = parseFasta(ref)
testKmers = parseFasta(test)
total = len(refKmers | testKmers)
print("Total number of kmers: " + str(total))
print("Number of exact kmers: " + str(len(refKmers)))
print("Number of neighbourgs: " + str(len(testKmers)))
print("New kmers (only in neighbourgs): " +
      str(len(testKmers - (refKmers & testKmers))))
print("Common kmers : " + str(len(refKmers & testKmers)))
