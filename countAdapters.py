#coding=utf-8

# Count the number of different adapters  in a file
# and print out the count for each adapter.

import sys

fileName = sys.argv[1]
adapter = {}

with open(fileName,'r') as f:
    for line  in f:
        adapt = line.rstrip("\n")
        if(adapt in adapter):
            adapter[adapt]+=1
        else:
            adapter[adapt] = 1

for elem,value in reversed(sorted([(k,v) for k,v in adapter.items()], key = lambda x: x[1])):
    print(elem,value)

