#!/usr/bin/env python

from pybrain.tools.shortcuts     import buildNetwork
from pybrain.supervised.trainers import BackpropTrainer
from pybrain.datasets            import SupervisedDataSet
from pybrain.structure.modules   import SigmoidLayer
from sklearn.ensemble            import RandomForestClassifier 
from sklearn import svm
import numpy
import feature
import dataset
import validate_performance

class CrossValidation(object):

    """
    Gene Scale in GA has to be (0, n).  n is greater than 0.
    """
    def __init__(self, bindres_file, pssms_file, log_file, method, fold=5, undersampling=True, shuffle=True, maxEpochs_for_trainer=10, geneScale=(0, 10)):
        if geneScale[0] != 0 or geneScale[1] <= geneScale[0]:
            raise ValueError("Gene Scale in GA has to be (0, n).  n is greater than 0.")
        if method != "neuralNetwork" and method != "randomForest" and method != "SVM":
            raise ValueError("method must be neuralNetwork or randomForest or SVM [{}]".format(method))
        self.bindingResidueData, self.pssmData = feature.parse_record_files(bindres_file, pssms_file)
        self.log_file = log_file
        self.method = method
        self.fold = fold
        self.undersampling = undersampling
        self.shuffle = shuffle
        self.maxEpochs_for_trainer=maxEpochs_for_trainer
        self.SVMParamScales = {"cost" : (-10, 10), "gamma" : (-10, 5)} 
        self.NNParamScales = {"node_num" : (5, 50), "learning_rate" : (0.01, 0.1)} 
        self.RFParamScales = {"n_estimators" : (101, 1001), "max_features" : (2, 30)}
        self.windowSizeScales = (1, 19)
        self.geneScale = geneScale
        self.log = {}

    def create_folded_dataset(self, window_size):
        positive_dataset, negative_dataset = feature.create_dataset(self.bindingResidueData, self.pssmData, window_size)
        folded_dataset = dataset.FoldedDataset(positive_dataset, negative_dataset, fold=self.fold,
                                                undersampling=self.undersampling, shuffle=self.shuffle)
        return folded_dataset

    def window_size_scaling(self, window_size_gene):
        return int(round(window_size_gene/float(self.geneScale[1])*(self.windowSizeScales[1]-self.windowSizeScales[0]))+self.windowSizeScales[0])

    def power_two_scaling(self, gene, scale_max, scale_min):
        return 2**int(round(gene/float(self.geneScale[1])*(scale_max - scale_min))+scale_min)

    def int_scaling(self, gene, scale_max, scale_min):
        return int(round(gene/float(self.geneScale[1])*(scale_max - scale_min))+scale_min)

    def float_scaling(self, gene, scale_max, scale_min):
        return gene/float(self.geneScale[1])*(scale_max - scale_min)+scale_min

    # Each gene is between 0 and n. n is greater than 0.
    def decode_chromosome(self, chromosome):
        if self.method == "neuralNetwork":
            node_num = self.int_scaling(chromosome[0], self.NNParamScales["node_num"][1], 
                                        self.NNParamScales["node_num"][0])
            learning_rate = self.float_scaling(chromosome[1], self.NNParamScales["learning_rate"][1],
                                               self.NNParamScales["learning_rate"][0])
            window_size = self.window_size_scaling(chromosome[2])
            return node_num, learning_rate, window_size
        elif self.method == "randomForest":
            n_estimators = self.int_scaling(chromosome[0], self.RFParamScales["n_estimators"][1],
                                            self.RFParamScales["n_estimators"][0])
            max_features = self.int_scaling(chromosome[1], self.RFParamScales["max_features"][1],
                                            self.RFParamScales["max_features"][0])
            window_size = self.window_size_scaling(chromosome[2])
            return n_estimators, max_features, window_size
        elif self.method == "SVM":
            cost = self.power_two_scaling(chromosome[0], self.SVMParamScales["cost"][1], self.SVMParamScales["cost"][0])
            gamma = self.power_two_scaling(chromosome[1], self.SVMParamScales["gamma"][1], self.SVMParamScales["gamma"][0])
            window_size = self.window_size_scaling(chromosome[2])
            return cost, gamma, window_size

    def write_log(self, gene1, gene2, gene3, mean_AUC, mean_decision_value, mean_mcc):
        with open(self.log_file, 'a') as fp:
            #fp.write("{} {} {} {}\n".format(gene1, gene2, gene3, mean_AUC))
            fp.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(gene1, gene2, gene3, mean_AUC, mean_decision_value, mean_mcc))

    def add_log(self, gene1, gene2, gene3, mean_AUC, mean_decision_value, mean_mcc):
        self.log[(gene1, gene2, gene3)] = (mean_AUC, mean_decision_value, mean_mcc)

    def check_log(self, gene1, gene2, gene3):
        if (gene1, gene2, gene3) in self.log:
            return True
        return False

    def get_means_from_log(self, gene1, gene2, gene3):
        return self.log.get((gene1, gene2, gene3))

    def eval_func(self, chromosome):
        if len(chromosome) != 3:
            raise ValueError("len(chromosome) is must be 3 [{}]".format(len(chromosome)))
        if self.method == "neuralNetwork":
            return self.neuralNetwork_eval_func(chromosome)
        elif self.method == "randomForest":
            return self.randomForest_eval_func(chromosome)
        elif self.method == "SVM":
            return self.SVM_eval_func(chromosome)

    def neuralNetwork_eval_func(self, chromosome):
        node_num, learning_rate, window_size = self.decode_chromosome(chromosome)
        if self.check_log(node_num, learning_rate, window_size):
            return self.get_means_from_log(node_num, learning_rate, window_size)[0]
        folded_dataset = self.create_folded_dataset(window_size)
        indim = 21 * (2 * window_size + 1)
        mean_AUC = 0
        mean_decision_value = 0
        mean_mcc = 0
        sample_size_over_thousand_flag = False
        for test_fold in xrange(self.fold):
            test_labels, test_dataset, train_labels, train_dataset = folded_dataset.get_test_and_training_dataset(test_fold)
            if len(test_labels) + len(train_labels) > 1000:
                sample_size_over_thousand_flag = True
            ds = SupervisedDataSet(indim, 1)
            for i in xrange(len(train_labels)):
                ds.appendLinked(train_dataset[i], [train_labels[i]])
            net = buildNetwork(indim, node_num, 1, outclass=SigmoidLayer, bias=True)
            trainer = BackpropTrainer(net, ds, learningrate=learning_rate)
            trainer.trainUntilConvergence(maxEpochs=self.maxEpochs_for_trainer)
            decision_values = [net.activate(test_dataset[i]) for i in xrange(len(test_labels))]
            decision_values = map(lambda x: x[0], decision_values)
            AUC, decision_value_and_max_mcc = validate_performance.calculate_AUC(decision_values, test_labels)
            mean_AUC += AUC
            mean_decision_value += decision_value_and_max_mcc[0]
            mean_mcc += decision_value_and_max_mcc[1]
            if sample_size_over_thousand_flag:
                break
        if not sample_size_over_thousand_flag:
            mean_AUC /= self.fold
            mean_decision_value /= self.fold
            mean_mcc /= self.fold
        self.write_log(node_num, learning_rate, window_size, mean_AUC, mean_decision_value, mean_mcc)
        self.add_log(node_num, learning_rate, window_size, mean_AUC, mean_decision_value, mean_mcc)
        return mean_AUC

    def randomForest_eval_func(self, chromosome):
        n_estimators, max_features, window_size = self.decode_chromosome(chromosome)
        if self.check_log(n_estimators, max_features, window_size):
            return self.get_means_from_log(n_estimators, max_features, window_size)[0]
        folded_dataset = self.create_folded_dataset(window_size)
        indim = 21 * (2 * window_size + 1)
        mean_AUC = 0
        mean_decision_value = 0
        mean_mcc = 0
        sample_size_over_thousand_flag = False
        for test_fold in xrange(self.fold):
            test_labels, test_dataset, train_labels, train_dataset = folded_dataset.get_test_and_training_dataset(test_fold)
            if len(test_labels) + len(train_labels) > 1000:
                sample_size_over_thousand_flag = True
            clf = RandomForestClassifier(n_estimators=n_estimators, max_features=max_features)
            clf.fit(train_dataset, train_labels)
            probas = clf.predict_proba(test_dataset)
            decision_values = map(lambda x: x[1], probas) # Probability of being binding residue
            AUC, decision_value_and_max_mcc = validate_performance.calculate_AUC(decision_values, test_labels)
            mean_AUC += AUC
            mean_decision_value += decision_value_and_max_mcc[0]
            mean_mcc += decision_value_and_max_mcc[1]
            if sample_size_over_thousand_flag:
                break
        if not sample_size_over_thousand_flag:
            mean_AUC /= self.fold
            mean_decision_value /= self.fold
            mean_mcc /= self.fold
        self.write_log(n_estimators, max_features, window_size, mean_AUC, mean_decision_value, mean_mcc)
        self.add_log(n_estimators, max_features, window_size, mean_AUC, mean_decision_value, mean_mcc)
        return mean_AUC

    def SVM_eval_func(self, chromosome):
        cost, gamma, window_size = self.decode_chromosome(chromosome)
        if self.check_log(cost, gamma, window_size):
            return self.get_means_from_log(cost, gamma, window_size)[0]
        folded_dataset = self.create_folded_dataset(window_size)
        indim = 21 * (2 * window_size + 1)
        mean_AUC = 0
        mean_decision_value = 0
        mean_mcc = 0
        sample_size_over_thousand_flag = False
        for test_fold in xrange(self.fold):
            test_labels, test_dataset, train_labels, train_dataset = folded_dataset.get_test_and_training_dataset(test_fold)
            if len(test_labels) + len(train_labels) > 1000:
                sample_size_over_thousand_flag = True
            clf = svm.SVC(C=cost, gamma=gamma, class_weight='auto')
            clf.fit(train_dataset, train_labels)
            decision_values = clf.decision_function(test_dataset)
            if type(decision_values[0]) is list or type(decision_values[0]) is numpy.ndarray:
                decision_values = map(lambda x: x[0], decision_values)
            AUC, decision_value_and_max_mcc = validate_performance.calculate_AUC(decision_values, test_labels)
            mean_AUC += AUC
            mean_decision_value += decision_value_and_max_mcc[0]
            mean_mcc += decision_value_and_max_mcc[1]
            if sample_size_over_thousand_flag:
                break
        if not sample_size_over_thousand_flag:
            mean_AUC /= self.fold
            mean_decision_value /= self.fold
            mean_mcc /= self.fold
        self.write_log(cost, gamma, window_size, mean_AUC, mean_decision_value, mean_mcc)
        self.add_log(cost, gamma, window_size, mean_AUC, mean_decision_value, mean_mcc)
        return mean_AUC
