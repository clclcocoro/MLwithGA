#!/usr/bin/env python

import re


class PSSM(object):

    """
    PSSM object created from followng string

    ">http://purl.uniprot.org/uniprot/Q9Y3D8\n
    1\t0\t3...\t-4\n
    1\t0\t3...\t-4\n
    1\t0\t3...\t-4\n
    ...
    1\t0\t3...\t-4\n"
    """

    def __init__(self, raw_pssm):
        self.uniprotURI, self.pssm = self.parse_raw_pssm(raw_pssm)

    def parse_raw_pssm(self, raw_pssm):
        parts = raw_pssm.split('\n')
        uniprotURI = parts[0][1:]
        pssm = []
        for part in parts[1:]:
            if len(part) != 0:
                pssm.append(map(int, part.split('\t')))
        return uniprotURI, pssm

    def get_PSSM(self):
        return self.pssm

    def get_uniprotURI(self):
        return self.uniprotURI


class PSSMData(object):

    """
    PSSMData object is dictionary of PSSM object
    """

    def __init__(self):
        self.pssms = {} 
        self.uniprotURIs = []

    def add_PSSMRecord(self, pssm_obj):
        self.uniprotURIs.append(pssm_obj.uniprotURI)
        self.pssms[pssm_obj.uniprotURI] = pssm_obj

    def get_PSSMRecord(self, uniprotURI):
        return self.pssms.get(uniprotURI)

    def get_uniprotURIs(self):
        return self.uniprotURIs


class BindingResidueData(object):

    """
    BindingResidue object created from followng string

    "http://purl.uniprot.org/uniprot/P00861 384 161 151 303 155 391 154 153 392 152 385 383 389 390 304 156 158\n
     http://purl.uniprot.org/uniprot/P16932 308 80 303 304 306 307 77 78 79 305\n
     http://purl.uniprot.org/uniprot/P30289 124 64 67 68 117\n
     http://purl.uniprot.org/uniprot/P32173 20 56 57 21\n"
    """

    def __init__(self, raw_bindres):
        self.uniprotURIs, self.bindRecords = self.parse_raw_bindres(raw_bindres)

    def parse_raw_bindres(self, raw_bindres):
        parts = raw_bindres.split('\n')
        uniprotURIs = []
        bindRecords = {}
        for part in parts:
            if len(part) == 0:
                continue
            eles = part.split(' ')
            uniprotURIs.append(eles[0])
            buff = []
            for ele in eles[1:]:
                if len(ele) != 0:
                    buff.append(int(ele))
            bindRecords[eles[0]] = set(buff)
        return uniprotURIs, bindRecords

    """
    return set([20, 56, 57, 21])
    """
    def get_bindRecord(self, uniprotURI):
        return self.bindRecords.get(uniprotURI)

    def get_uniprotURIs(self):
        return self.uniprotURIs


def parse_pssms_file(pssms_file):
    pssmData = PSSMData()
    with open(pssms_file) as fp:
        raw_pssm = ""
        for i, line in enumerate(fp):
            if i != 0 and re.match(">", line):
                pssmData.add_PSSMRecord(PSSM(raw_pssm))
                raw_pssm = line
            else:
                raw_pssm += line
        pssmData.add_PSSMRecord(PSSM(raw_pssm))
    return pssmData


def parse_record_files(bindres_file, pssms_file):
    with open(bindres_file) as fp:
        raw_bindres = ''.join(fp.readlines())
        bindingResidueData = BindingResidueData(raw_bindres)
    pssmData = PSSMData()
    with open(pssms_file) as fp:
        raw_pssm = ""
        for i, line in enumerate(fp):
            if i != 0 and re.match(">", line):
                pssmData.add_PSSMRecord(PSSM(raw_pssm))
                raw_pssm = line
            else:
                raw_pssm += line
        pssmData.add_PSSMRecord(PSSM(raw_pssm))
    return bindingResidueData, pssmData


def create_feature_vectors(pssm, window_size):
    """
        terminal spacer, PSSM 
        1 0 0 0 0 0 0 0 0 0 0 ... 0
        1 0 0 0 0 0 0 0 0 0 0 ... 0
        0 -1 -2 -1 3 ...         -3
    """
    feature_vectors = []
    m = pssm.get_PSSM()
    seqlen = len(m)
    for i in xrange(seqlen):
        feature_vector = []
        p = i - window_size 
        if p < 0:
            feature_vector += [1 if j != 0 and (j+1) % 21 == 0 else 0 for j in xrange(21*(window_size-i))]
            p = 0
        if i + window_size <= seqlen - 1:
            while p <= i + window_size:
                feature_vector += m[p]+[0]
                p += 1
        else:
            while p <= seqlen - 1:
                feature_vector += m[p]+[0]
                p += 1
            feature_vector += [1 if j != 0 and (j+1) % 21 == 0 else 0 for j in xrange(21*((i+window_size)-(seqlen-1)))]
        feature_vectors.append(feature_vector)
    return feature_vectors

# Only 5-25th (former or latter) residues from binding residue are used as negative dataset.
# So, 1-4th (former or latter) residues from binding residue must be eliminated from dataset.

def get_negative_data_index_set(bindRecord, sequence_length):
    negative_data_index = set()
    for i in bindRecord:
        for j in xrange(i-25, i-4):
            if j >= 0:
                negative_data_index.add(j)
        for j in xrange(i+5, i+26):
            if j <= sequence_length - 1:
                negative_data_index.add(j)
    for i in bindRecord:
        for j in xrange(i-4, i+5):
            if j in negative_data_index:
                negative_data_index.remove(j)
    return negative_data_index


def create_training_data(bindRecord, feature_vectors, negative_data_index_set):
    positive_data = []
    negative_data = []
    for i, feature_vector in enumerate(feature_vectors):
        if i in bindRecord:
            positive_data.append(feature_vector)
        else:
            if i in negative_data_index_set:
                negative_data.append(feature_vector)
    return positive_data, negative_data


def create_dataset(bindingResidueData, pssmData, window_size):
    positive_dataset = []
    negative_dataset = []
    for uniprotURI in bindingResidueData.get_uniprotURIs():
        pssm = pssmData.get_PSSMRecord(uniprotURI)
        feature_vectors = create_feature_vectors(pssm, window_size)
        bindRecord = bindingResidueData.get_bindRecord(uniprotURI)
        negative_data_index_set = get_negative_data_index_set(bindRecord, len(pssm.get_PSSM()))
        positive_data, negative_data = create_training_data(bindRecord, feature_vectors, negative_data_index_set)
        positive_dataset += positive_data
        negative_dataset += negative_data
    return positive_dataset, negative_dataset


if __name__ == "__main__":
    bindres_file = "/Users/clclcocoro/galaxy/work/data/bindingData.txt"
    pssms_file = "/Users/clclcocoro/galaxy/work/data/pssms.txt"
    bindingResidueData, pssmData = parse_record_files(bindres_file, pssms_file)
    
    """
    print "bindingResidueData"
    print bindingResidueData
    print bindingResidueData.get_uniprotURIs()
    print bindingResidueData.get_bindRecord(bindingResidueData.get_uniprotURIs()[0])
    print "pssmData"
    print pssmData
    print pssmData.get_uniprotURIs()
    print pssmData.get_PSSMRecord(pssmData.get_uniprotURIs()[0]).get_PSSM()[:10]
    """
    
    window_size = 3
    positive_dataset, negative_dataset = create_dataset(bindingResidueData, pssmData, window_size)
    
    print "positive_dataset"
    print positive_dataset
    print "negative_dataset"
    print negative_dataset
