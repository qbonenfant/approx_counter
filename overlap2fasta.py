#coding=utf-8

import sys

fileName = sys.argv[1]

with open(fileName) as f:
	for i,line in enumerate(f):
		print(">Overlap_" + str(i))
		print(line.split()[0])
