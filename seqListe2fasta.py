#coding=utf-8
# Convert text file containing one sequence per line to fasta.
import sys

fileName = sys.argv[1]

with open(fileName) as f:
	for i,line in enumerate(f):
		print(">Sequence_" + str(i))
		print(line.split()[0])
