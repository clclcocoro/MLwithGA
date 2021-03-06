#!/usr/bin/env python

import sys
sys.path.append("..")

import unittest
import validate_performance
import feature
import dataset

class TestValidate(unittest.TestCase):

    def test_calculate_TPR_FPR(self):
        TP, FP, TN, FN = 1, 1, 1, 1
        FPR, TPR = validate_performance.calculate_TPR_FPR(TP, FP, TN, FN)
        self.assertEqual(TPR, 0.5)
        self.assertEqual(FPR, 0.5)

        TP, FP, TN, FN = 1, 2, 3, 4
        FPR, TPR = validate_performance.calculate_TPR_FPR(TP, FP, TN, FN)
        self.assertEqual(TPR, 0.2)
        self.assertEqual(FPR, 0.4)

    def test_calculate_AUC(self):
        decision_values = [-1, -0.5, -0.1,  0.1, 0.5, 1]
        correct_labels = [0, 0, 0, 1, 1, 1]
        AUC, mcc = validate_performance.calculate_AUC(decision_values, correct_labels)
        self.assertEqual(AUC, 1.0)

        decision_values = [-1, 0.5, -0.1,  0.1, -0.5, 1]
        correct_labels = [0, 1, 0, 1, 0, 1]
        AUC, mcc = validate_performance.calculate_AUC(decision_values, correct_labels)
        self.assertEqual(AUC, 1.0)

        decision_values = [1.0]*6
        correct_labels = [0, 1, 0, 1, 0, 1]
        AUC, mcc = validate_performance.calculate_AUC(decision_values, correct_labels)
        self.assertEqual(AUC, 0.5)

        decision_values = [-1, -0.5, 0.5, 1]
        correct_labels = [0, 1, 0, 1]
        AUC, mcc = validate_performance.calculate_AUC(decision_values, correct_labels)
        self.assertEqual(AUC, 0.75)

        decision_values = [1.0] * 2 + [-1.0] * 8 + [1.0] * 2 + [0.5] * 8
        correct_labels = [0] * 10 + [1] * 10
        AUC, mcc = validate_performance.calculate_AUC(decision_values, correct_labels)
        self.assertEqual(round(AUC*(10**5))/(10**5), 0.82)

        decision_values = [-1.0] * 8 + [0.0] * 2 + [0.0] * 2 + [0.5] * 8
        correct_labels = [0] * 10 + [1] * 10
        AUC, mcc = validate_performance.calculate_AUC(decision_values, correct_labels)
        self.assertEqual(round(AUC*(10**5))/(10**5), 0.98)

        decision_values = [-1.0] * 8 + [0.2] * 2 + [0.1]*3 +[0.2]*2 + [0.5]*5
        correct_labels = [0] * 10 + [1] * 10
        AUC, mcc = validate_performance.calculate_AUC(decision_values, correct_labels)
        self.assertEqual(round(AUC*(10**5))/(10**5), 0.92)

        decision_values = [-1.0]*6 + [0.0]*2+ [0.2]*2 + [0.0]*2+ [0.1]*1 +[0.2]*2 + [0.5]*5
        correct_labels = [0] * 10 + [1] * 10
        AUC, mcc = validate_performance.calculate_AUC(decision_values, correct_labels)
        self.assertEqual(round(AUC*(10**5))/(10**5), 0.9)


def create_positive_and_negative_dataset(window_size, sequence_length):
        bindres_file = "/tmp/bindingData.txt"
        pssms_file = "/tmp/pssms.txt"
        with open(bindres_file, "w") as fp:
            fp.write("http://purl.uniprot.org/uniprot/AAAAAA 0 1 2\n")
            fp.write("http://purl.uniprot.org/uniprot/BBBBBB 9\n")
            fp.write("http://purl.uniprot.org/uniprot/CCCCCC 7 2\n")
        with open(pssms_file, "w") as fp:
            fp.write(">http://purl.uniprot.org/uniprot/AAAAAA\n")
            pssm = '\n'.join(map('\t'.join, [['1' if i == j else '-1' for i in xrange(20)]+['5' for l in xrange(20)] for j in xrange(sequence_length)]))
            fp.write(pssm+"\n")
            fp.write(">http://purl.uniprot.org/uniprot/BBBBBB\n")
            pssm = '\n'.join(map('\t'.join, [['2' if i == j else '-2' for i in xrange(20)]+['5' for l in xrange(20)] for j in xrange(sequence_length)]))
            fp.write(pssm+"\n")
            fp.write(">http://purl.uniprot.org/uniprot/CCCCCC\n")
            pssm = '\n'.join(map('\t'.join, [['3' if i == j else '-3' for i in xrange(20)]+['5' for l in xrange(20)] for j in xrange(sequence_length)]))
            fp.write(pssm+"\n")
        bindingResidueData, pssmData = feature.parse_record_files(bindres_file, pssms_file)
        positive_dataset, negative_dataset = feature.create_dataset(bindingResidueData, pssmData, window_size)
        return positive_dataset, negative_dataset
 

class TestFeature(unittest.TestCase):

    def test_jensen_shannon_divergence(self):
        JSD = feature.jensen_shennon_divergence(feature.background_amino_acid_probs)
        epsilon = 10**-10
        self.assertTrue(JSD < epsilon)
        bg = [1.0 if i == 1 else 0 for i in xrange(20)]
        JSD = feature.jensen_shennon_divergence([1 if i == 0 else 0 for i in xrange(20)], background_probs=bg)
        epsilon = 10**-5
        self.assertTrue(JSD > 1-epsilon)
        JSD = feature.jensen_shennon_divergence([5 for i in xrange(20)], background_probs=feature.background_amino_acid_probs)
        self.assertEqual(JSD, 0.0375597100129578)

    def test_create_feature_vectors(self):
        create_positive_and_negative_dataset(1, 10)
        pssms_file = "/tmp/pssms.txt"
        conservation = True
        pssmData = feature.parse_pssms_file(pssms_file)
        pssm = pssmData.get_PSSMRecord(pssmData.get_uniprotURIs()[0])
        feature_vectors = feature.create_feature_vectors(pssm, 1, conservation=True)
        self.assertEqual(66, len(feature_vectors[0]))
        self.assertEqual(0, feature_vectors[0][-3])
        self.assertEqual(0, feature_vectors[9][-1])

    def test_create_dataset(self):
        window_size = 1
        positive_dataset, negative_dataset = create_positive_and_negative_dataset(window_size, 10)
        self.assertEqual(len(positive_dataset), 6)
        self.assertEqual(len(negative_dataset), 8)
        correct_positive = [
                    [0]*20              +[0]
                    +[1]+[-1]*19        +[1]
                    +[-1]+[1]+[-1]*18   +[0],

                    [1]+[-1]*19         +[1]
                    +[-1]+[1]+[-1]*18   +[0]
                    +[-1]*2+[1]+[-1]*17 +[0],

                    [-1]+[1]+[-1]*18    +[0]
                    +[-1]*2+[1]+[-1]*17 +[0]
                    +[-1]*3+[1]+[-1]*16 +[0],

                    [-2]*8+[2]+[-2]*11  +[0]
                    +[-2]*9+[2]+[-2]*10 +[1]
                    +[0]*20             +[0],

                    [-3]+[3]+[-3]*18    +[0]
                    +[-3]*2+[3]+[-3]*17 +[0]
                    +[-3]*3+[3]+[-3]*16 +[0],
                                            
                    [-3]*6+[3]+[-3]*13  +[0]
                    +[-3]*7+[3]+[-3]*12 +[0]
                    +[-3]*8+[3]+[-3]*11 +[0]
                ]
        for i in xrange(6):
            for j, ele in enumerate(positive_dataset[i]):
                self.assertEqual(ele, correct_positive[i][j])
        correct_negative = [
                    [-1]*6+[1]+[-1]*13  +[0]
                    +[-1]*7+[1]+[-1]*12 +[0]
                    +[-1]*8+[1]+[-1]*11 +[0], # AAAAAA 7

                    [-1]*7+[1]+[-1]*12  +[0]
                    +[-1]*8+[1]+[-1]*11 +[0]
                    +[-1]*9+[1]+[-1]*10 +[1], # AAAAAA 8
                                            
                    [-1]*8+[1]+[-1]*11  +[0]
                    +[-1]*9+[1]+[-1]*10 +[1]
                    +[0]*20             +[0], # AAAAAA 9
                                            
                    [0]*20              +[0]
                    +[2]+[-2]*19        +[1]
                    +[-2]+[2]+[-2]*18   +[0], # BBBBBB 0
                ]
        index = [0, 1, 2, 3]
        for i in xrange(4):
            for j, ele in enumerate(negative_dataset[index[i]]):
                self.assertEqual(ele, correct_negative[i][j])

        window_size = 2
        positive_dataset, negative_dataset = create_positive_and_negative_dataset(window_size, 10)
        self.assertEqual(len(positive_dataset), 6)
        self.assertEqual(len(negative_dataset), 8)
        correct_positive = [
                    [0]*20              +[0]
                    +[0]*20             +[0]
                    +[1]+[-1]*19        +[1]
                    +[-1]+[1]+[-1]*18   +[0]
                    +[-1]*2+[1]+[-1]*17 +[0],

                    [0]*20              +[0]
                    +[1]+[-1]*19        +[1]
                    +[-1]+[1]+[-1]*18   +[0]
                    +[-1]*2+[1]+[-1]*17 +[0]
                    +[-1]*3+[1]+[-1]*16 +[0],

                    [1]+[-1]*19         +[1]
                    +[-1]+[1]+[-1]*18   +[0]
                    +[-1]*2+[1]+[-1]*17 +[0]
                    +[-1]*3+[1]+[-1]*16 +[0]
                    +[-1]*4+[1]+[-1]*15 +[0],

                    [-2]*7+[2]+[-2]*12  +[0]
                    +[-2]*8+[2]+[-2]*11 +[0]
                    +[-2]*9+[2]+[-2]*10 +[1]
                    +[0]*20             +[0]
                    +[0]*20             +[0],

                    [3]+[-3]*19         +[1]
                    +[-3]+[3]+[-3]*18   +[0]
                    +[-3]*2+[3]+[-3]*17 +[0]
                    +[-3]*3+[3]+[-3]*16 +[0]
                    +[-3]*4+[3]+[-3]*15 +[0],
                                            
                    [-3]*5+[3]+[-3]*14  +[0]
                    +[-3]*6+[3]+[-3]*13 +[0]
                    +[-3]*7+[3]+[-3]*12 +[0]
                    +[-3]*8+[3]+[-3]*11 +[0]
                    +[-3]*9+[3]+[-3]*10 +[1]
                ]
        for i in xrange(6):
            for j, ele in enumerate(positive_dataset[i]):
                self.assertEqual(ele, correct_positive[i][j])
        correct_negative = [
                    [-1]*5+[1]+[-1]*14  +[0]
                    +[-1]*6+[1]+[-1]*13 +[0]
                    +[-1]*7+[1]+[-1]*12 +[0]
                    +[-1]*8+[1]+[-1]*11 +[0]
                    +[-1]*9+[1]+[-1]*10 +[1], # AAAAAA 7

                    [-1]*6+[1]+[-1]*13  +[0]
                    +[-1]*7+[1]+[-1]*12 +[0]
                    +[-1]*8+[1]+[-1]*11 +[0]
                    +[-1]*9+[1]+[-1]*10 +[1]
                    +[0]*20             +[0], # AAAAAA 8
                                            
                    [-1]*7+[1]+[-1]*12  +[0]
                    +[-1]*8+[1]+[-1]*11 +[0]
                    +[-1]*9+[1]+[-1]*10 +[1]
                    +[0]*20             +[0]
                    +[0]*20             +[0], # AAAAAA 9
                                            
                    [0]*20              +[0]
                    +[0]*20             +[0]
                    +[2]+[-2]*19        +[1]
                    +[-2]+[2]+[-2]*18   +[0]
                    +[-2]*2+[2]+[-2]*17 +[0], # BBBBBB 0
                ]
        index = [0, 1, 2, 3]
        for i in xrange(4):
            for j, ele in enumerate(negative_dataset[index[i]]):
                self.assertEqual(ele, correct_negative[i][j])

        window_size = 20
        positive_dataset, negative_dataset = create_positive_and_negative_dataset(window_size, 10)
        self.assertEqual(len(positive_dataset), 6)
        self.assertEqual(len(negative_dataset), 8)
        correct_positive = [0 for i in xrange(21*20)] + [1 if i == j else -1 for i in xrange(10) for j in xrange(21)] + [0 for i in xrange(21*11)]
        for i in xrange(10):
            if i == 0 or i == 9:
                correct_positive[420+21*i+20] = 1
            else:
                correct_positive[420+21*i+20] = 0
        for i, ele in enumerate(positive_dataset[0]):
            self.assertEqual(ele, correct_positive[i])

        window_size = 1
        positive_dataset, negative_dataset = create_positive_and_negative_dataset(window_size, 30)
        self.assertEqual(len(positive_dataset), 6)
        self.assertEqual(len(negative_dataset), 60)
        correct_negative = [-2 for i in xrange(42)] + [0]*21
        for i in xrange(2):
            if i == 1:
                correct_negative[21*i+20] = 1
            else:
                correct_negative[21*i+20] = 0
        for i, ele in enumerate(negative_dataset[41]):
            self.assertEqual(ele, correct_negative[i])


class TestDataset(unittest.TestCase):

    def test_folded_dataset(self):
        window_size = 1
        positive_dataset, negative_dataset = create_positive_and_negative_dataset(window_size, 10)
        foldedDataset = dataset.FoldedDataset(positive_dataset, negative_dataset)
        test_labels, test_dataset, train_labels, train_dataset = foldedDataset.get_test_and_training_dataset(0)
        self.assertEqual(len(test_labels), 4)
        self.assertEqual(len(train_labels), 10)

        test_labels, test_dataset, train_labels, train_dataset = foldedDataset.get_test_and_training_dataset(1)
        self.assertEqual(len(test_labels), 3)
        self.assertEqual(len(train_labels), 11)

        test_labels, test_dataset, train_labels, train_dataset = foldedDataset.get_test_and_training_dataset(2)
        self.assertEqual(len(test_labels), 3)
        self.assertEqual(len(train_labels), 11)

        for i in xrange(3, 5):
            test_labels, test_dataset, train_labels, train_dataset = foldedDataset.get_test_and_training_dataset(i)
            self.assertEqual(len(test_labels), 2)
            self.assertEqual(len(train_labels), 12)

    def test_folding(self):
        window_size = 1
        positive_dataset, negative_dataset = create_positive_and_negative_dataset(window_size, 10)
        foldedDataset = dataset.FoldedDataset(positive_dataset, negative_dataset)
        myDataset = [i for i in xrange(11)]
        size = len(myDataset)
        fold = 5
        folded_dataset = foldedDataset.folding(size, myDataset, fold)
        correct = [[0, 1, 10],
                   [2, 3],
                   [4, 5],
                   [6, 7], 
                   [8, 9]
                ]
        for i in xrange(len(folded_dataset)):
            for j in xrange(len(folded_dataset[i])):
                self.assertEqual(folded_dataset[i][j], correct[i][j])


if __name__ == "__main__":
    unittest.main()
