#!/usr/local/bin/python
import os
import sys
import warnings

filename = str(sys.argv[1])
idx1 = int(sys.argv[2])
idx2 = int(sys.argv[3])

if os.path.exists(filename):
    with open(filename) as f:
        lines = f.readlines()
        counter = 0
        for line in lines:
            counter += 1
            cols=line.split('\t')
            start=int(cols[idx1])
            stop=int(cols[idx2])
            size=stop-start
            if size != 0:
                print line.strip('\n')
            else:
                warnings.warn("Element size 0 for line " + str(counter))
