#!/usr/bin/env python

import sys

"""
Usage: create_testdata.py <interval>
"""

# start is the first binding residue index.
# interval is the number of residues between binding residues.
def generate_pssm(start, interval):
    pssm = []
    for i in xrange(100):
        if i % interval == start:
            pssm.append(map(str, [1]*20))
        elif i % interval == start-1 or i % interval == start+1:
            pssm.append(map(str, [-1]*20))
        else:
            pssm.append(map(str, [0]*20))
    return pssm
 
interval = int(sys.argv[1])
bindres_file = "./bindingData.txt"
pssms_file = "./pssms.txt"
with open(bindres_file, "w") as fp:
    startA = 1
    startB = 2
    startC = 3
    binding_site_indexA = ' '.join(map(str, [i+startA for i in xrange(0, 100, interval)]))
    binding_site_indexB = ' '.join(map(str, [i+startB for i in xrange(0, 100, interval)]))
    binding_site_indexC = ' '.join(map(str, [i+startC for i in xrange(0, 100, interval)]))
    fp.write("http://purl.uniprot.org/uniprot/AAAAAA {}\n".format(binding_site_indexA))
    fp.write("http://purl.uniprot.org/uniprot/BBBBBB {}\n".format(binding_site_indexB))
    fp.write("http://purl.uniprot.org/uniprot/CCCCCC {}\n".format(binding_site_indexC))
with open(pssms_file, "w") as fp:
    fp.write(">http://purl.uniprot.org/uniprot/AAAAAA\n")
    pssm = '\n'.join(map('\t'.join, generate_pssm(startA, interval)))
    fp.write(pssm+"\n")
    fp.write(">http://purl.uniprot.org/uniprot/BBBBBB\n")
    pssm = '\n'.join(map('\t'.join, generate_pssm(startB, interval)))
    fp.write(pssm+"\n")
    fp.write(">http://purl.uniprot.org/uniprot/CCCCCC\n")
    pssm = '\n'.join(map('\t'.join, generate_pssm(startC, interval)))
    fp.write(pssm+"\n")
 
