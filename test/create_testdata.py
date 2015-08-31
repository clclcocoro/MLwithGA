#!/usr/bin/env python


def generate_pssm(p): # p = binding site integer between (1, 4)
    pssm = []
    for i in xrange(100):
        if i % 5 == p:
            pssm.append(map(str, [1]*20))
        elif i % 5 == p-1 or i % 5 == p+1:
            pssm.append(map(str, [-1]*20))
        else:
            pssm.append(map(str, [0]*20))
    return pssm
 

bindres_file = "./bindingData.txt"
pssms_file = "./pssms.txt"
with open(bindres_file, "w") as fp:
    binding_site_indexA = ' '.join(map(str, [i+1 for i in xrange(0, 100, 5)]))
    binding_site_indexB = ' '.join(map(str, [i+2 for i in xrange(0, 100, 5)]))
    binding_site_indexC = ' '.join(map(str, [i+3 for i in xrange(0, 100, 5)]))
    fp.write("http://purl.uniprot.org/uniprot/AAAAAA {}\n".format(binding_site_indexA))
    fp.write("http://purl.uniprot.org/uniprot/BBBBBB {}\n".format(binding_site_indexB))
    fp.write("http://purl.uniprot.org/uniprot/CCCCCC {}\n".format(binding_site_indexC))
with open(pssms_file, "w") as fp:
    fp.write(">http://purl.uniprot.org/uniprot/AAAAAA\n")
    pssm = '\n'.join(map('\t'.join, generate_pssm(1)))
    fp.write(pssm+"\n")
    fp.write(">http://purl.uniprot.org/uniprot/BBBBBB\n")
    pssm = '\n'.join(map('\t'.join, generate_pssm(2)))
    fp.write(pssm+"\n")
    fp.write(">http://purl.uniprot.org/uniprot/CCCCCC\n")
    pssm = '\n'.join(map('\t'.join, generate_pssm(3)))
    fp.write(pssm+"\n")
 
