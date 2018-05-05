# coding=utf-8
import sys
global BASES

BASES = set(["A", "T", "C", "G"])


def substitution(sequence):
    allSub = []
    for i, c in enumerate(sequence):
        for n in BASES - set(c):
            allSub.append(sequence[0:i] + n + sequence[i + 1:])

    return(allSub)


def deletion(sequence):
    allDel = []
    for i in range(len(sequence)):
        allDel.append(sequence[0:i] + sequence[i + 1:])
    return(allDel)


def insertion(sequence):
    allIns = []
    for i in range(len(sequence) + 1):
        for n in BASES:
            allIns.append(sequence[0:i] + n + sequence[i:])

    return(allIns)


allKmer = set()

filename = sys.argv[1] if len(sys.argv) >= 2 else "nope.txt"
k = sys.argv[2] if len(sys.argv) >= 3 else 2

with open(filename, "r") as f:
    for seq1 in f:
        seq1 = seq1.rstrip("\n").split(" ")[0]
        oneErr = insertion(seq1) + substitution(seq1) + deletion(seq1)
        for seq in oneErr:
            allKmer |= set(insertion(seq) +
                           substitution(seq) + deletion(seq))

out = open(".".join(filename.split(".")[:-1]) + "_neighbors.fasta", "w")
for i, kmer in enumerate(allKmer):
    out.write(">kmer_" + str(i) + " \n")
    out.write(kmer + " \n")
out.close()
