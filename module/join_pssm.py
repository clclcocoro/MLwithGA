import glob
import os
import sys
import base64
import re

def print_pssm(fname):
    with open(fname) as fp:
        for c, line in enumerate(fp):
            if c <= 2:
                continue
            if len(line.rstrip()) == 0:
                break
            rec = line.strip().split()
            if len(rec) == 44:
                print "\t".join(map(lambda x: x.strip(),  rec[2:42]))
            else:
                rec = []
                for i in xrange(9, 67, 3):
                    rec.append(line[i:i+3].strip())
                for i in xrange(70, 150, 4):
                    rec.append(line[i:i+3].strip())
                print "\t".join(rec)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit('Usage: %s pssm directory' % sys.argv[0])
    if not os.path.exists(sys.argv[1]):
        sys.exit('ERROR: Directory %s not found' % sys.argv[1])
    if not os.path.isdir(sys.argv[1]):
        sys.exit('ERROR: %s is not Directory' % sys.argv[1])

    DIR = sys.argv[1]
    if DIR[-1] != "/":
        DIR += "/"
    
    for fname in glob.iglob(DIR + "*.pssm"):
        print ">%s" % base64.b64decode(fname.split("/")[-1].replace(".fasta.pssm",""))
        print_pssm(fname)
