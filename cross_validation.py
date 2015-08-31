#!/usr/bin/env python

from pybrain.tools.shortcuts     import buildNetwork
from pybrain.supervised.trainers import BackpropTrainer
from pybrain.datasets            import SupervisedDataSet
from sklearn.ensemble            import RandomForestClassifier 
from sklearn import svm
import feature
import dataset
import validate

class CrossValidation(object):

    """
    Gene Scale in GA has to be (0, n).  n is greater than 0.
    """
    def __init__(self, bindres_file, pssms_file, log_file, method, fold=5, undersampling=True, shuffle=True, maxEpochs_for_trainer=30, geneScale=(0, 10)):
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
        self.SVMParamScales = {"cost" : (-10, 10), "gamma" : (-10, 10)} 
        self.NNParamScales = {"node_num" : (5, 50), "learning_rate" : (0.01, 0.1)} 
        self.RFParamScales = {"n_estimators" : (101, 3001), "max_features" : (2, 30)}
        self.windowSizeScales = (1, 19)
        self.geneScale=geneScale

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

    def write_log(self, gene1, gene2, gene3, mean_AUC):
        with open(self.log_file, 'a') as fp:
            #fp.write("{} {} {} {}\n".format(gene1, gene2, gene3, mean_AUC))
            fp.write("{}\t{}\t{}\t{}\n".format(gene1, gene2, gene3, mean_AUC))

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
        folded_dataset = self.create_folded_dataset(window_size)
        indim = 21 * (2 * window_size + 1)
        mean_AUC = 0
        for test_fold in xrange(self.fold):
            test_labels, test_dataset, train_labels, train_dataset = folded_dataset.get_test_and_training_dataset(test_fold)
            ds = SupervisedDataSet(indim, 1)
            for i in xrange(len(train_labels)):
                ds.appendLinked(train_dataset[i], [train_labels[i]])
            net = buildNetwork(indim, node_num, 1)
            trainer = BackpropTrainer(net, ds, learningrate=learning_rate)
            trainer.trainUntilConvergence(maxEpochs=self.maxEpochs_for_trainer)
            decision_values = [net.activate(test_dataset[i]) for i in xrange(len(test_labels))]
            mean_AUC += validate.calculate_AUC(decision_values, test_labels)
        mean_AUC /= self.fold
        self.write_log(node_num, learning_rate, window_size, mean_AUC)
        return mean_AUC

    def randomForest_eval_func(self, chromosome):
        n_estimators, max_features, window_size = self.decode_chromosome(chromosome)
        folded_dataset = self.create_folded_dataset(window_size)
        indim = 21 * (2 * window_size + 1)
        mean_AUC = 0
        for test_fold in xrange(self.fold):
            test_labels, test_dataset, train_labels, train_dataset = folded_dataset.get_test_and_training_dataset(test_fold)
            clf = RandomForestClassifier(n_estimators=n_estimators, max_features=max_features)
            clf.fit(train_dataset, train_labels)
            probas = clf.predict_proba(test_dataset)
            decision_values = map(lambda x: x[1], probas) # Probability of being binding residue
            mean_AUC += validate.calculate_AUC(decision_values, test_labels)
        mean_AUC /= self.fold
        self.write_log(n_estimators, max_features, window_size, mean_AUC)
        return mean_AUC

    def SVM_eval_func(self, chromosome):
        cost, gamma, window_size = self.decode_chromosome(chromosome)
        folded_dataset = self.create_folded_dataset(window_size)
        indim = 21 * (2 * window_size + 1)
        mean_AUC = 0
        for test_fold in xrange(self.fold):
            test_labels, test_dataset, train_labels, train_dataset = folded_dataset.get_test_and_training_dataset(test_fold)
            clf = svm.SVC(C=cost, gamma=gamma, class_weight='auto')
            clf.fit(train_dataset, train_labels)
            decision_values = clf.decision_function(test_dataset)
            mean_AUC += validate.calculate_AUC(decision_values, test_labels)
        mean_AUC /= self.fold
        self.write_log(cost, gamma, window_size, mean_AUC)
        return mean_AUC
