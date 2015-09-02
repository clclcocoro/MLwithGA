#!/usr/bin/env python

import sys
import random

"""
Usage: create_testdata.py <interval> <sequence_length> [--randomscore]

    <interval>        integer that is bigger than 8
    <sequence_length> integer that is bigger than 9
    --randomscore     score is sampled from discrete random distribution
                      binding residue                             randint(3, 5)
                      the both sides residues of binding residue  randint(-10, -8)
                      non-binding residue                         randint(-2, 0)
"""

# start is the first binding residue index.
# interval is the number of residues between binding residues.
def generate_pssm(start, sequence_length, interval, random_flag=False):
    pssm = []
    for i in xrange(sequence_length):
        if i % interval == start:
            if random_flag:
                pssm.append(map(str, [random.randint(3, 5) for i in xrange(20)]))
            else:
                pssm.append(map(str, [1]*20))
        elif i % interval == start-1 or i % interval == start+1:
            if random_flag:
                pssm.append(map(str, [random.randint(-10, -8) for i in xrange(20)]))
            else:
                pssm.append(map(str, [-1]*20))
        else:
            if random_flag:
                pssm.append(map(str, [random.randint(-2, 0) for i in xrange(20)]))
            else:
                pssm.append(map(str, [0]*20))
    return pssm
 
if sys.argv[1] == "-h" or sys.argv[1] == "-help" or sys.argv[1] == "--help":
    print """
Usage: create_testdata.py <interval> <sequence_length> [--randomscore]

    <interval>        integer that is bigger than 8 
    <sequence_length> integer that is bigger than 9
    --randomscore     score is sampled from discrete random distribution
                      binding residue                             randint(1, 10)
                      the both sides residues of binding residue  randint(-10, -8)
                      non-binding residue                         randint(-7, 0)
"""
    sys.exit(0)
interval = int(sys.argv[1])
if not interval > 8:
    raise ValueError("<interval> must be bigger than 8")
interval += 1 # modify for xrange()
sequence_length = int(sys.argv[2])
if not sequence_length > 9:
    raise ValueError("<sequence_length> must be bigger than 9")
random_flag = False
if len(sys.argv) == 4 and sys.argv[3] == "--randomscore":
    random_flag = True
sequence_length = int(sys.argv[2])
bindres_file = "./bindingData.txt"
if random_flag:
    pssms_file = "./pssms_random_score.txt"
else:
    pssms_file = "./pssms_fixed_score.txt"
with open(bindres_file, "w") as fp:
    startA = 1
    startB = 2
    startC = 3
    binding_site_indexA = ' '.join(map(str, [i+startA for i in xrange(0, sequence_length, interval)]))
    binding_site_indexB = ' '.join(map(str, [i+startB for i in xrange(0, sequence_length, interval)]))
    binding_site_indexC = ' '.join(map(str, [i+startC for i in xrange(0, sequence_length, interval)]))
    fp.write("http://purl.uniprot.org/uniprot/AAAAAA {}\n".format(binding_site_indexA))
    fp.write("http://purl.uniprot.org/uniprot/BBBBBB {}\n".format(binding_site_indexB))
    fp.write("http://purl.uniprot.org/uniprot/CCCCCC {}\n".format(binding_site_indexC))
with open(pssms_file, "w") as fp:
    fp.write(">http://purl.uniprot.org/uniprot/AAAAAA\n")
    pssm = '\n'.join(map('\t'.join, generate_pssm(startA, sequence_length, interval, random_flag)))
    fp.write(pssm+"\n")
    fp.write(">http://purl.uniprot.org/uniprot/BBBBBB\n")
    pssm = '\n'.join(map('\t'.join, generate_pssm(startB, sequence_length, interval, random_flag)))
    fp.write(pssm+"\n")
    fp.write(">http://purl.uniprot.org/uniprot/CCCCCC\n")
    pssm = '\n'.join(map('\t'.join, generate_pssm(startC, sequence_length, interval, random_flag)))
    fp.write(pssm+"\n")
 
