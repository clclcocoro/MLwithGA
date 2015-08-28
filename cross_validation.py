#!/usr/bin/env python

from pybrain.tools.shortcuts     import buildNetwork
from pybrain.supervised.trainers import BackpropTrainer
from pybrain.datasets            import SupervisedDataset
from sklearn.ensemble            import RandomForestClassifier 
import feature
import dataset
import validate

class CrossValidation(object):

    """
    Gene Scale in GA has to be (0, n).  n is greater than 0.
    """
    def __init__(self, bindres_file, pssms_file, method, fold=5, undersampling=True, shuffle=True, maxEpochs_for_trainer=180, geneScale=(0, 10)):
        if geneScale[0] =! 0 or geneScale[1] <= geneScale[0]:
            raise ValueError("Gene Scale in GA has to be (0, n).  n is greater than 0.")
        if method != "neuralNetwork" and method != "randomForest":
            raise ValueError("method must be neuralNetwork or randomForest [{}]".format(method))
        self.bindingResidueData, self.pssmData = feature.parse_record_files(bindres_file, pssms_file)
        self.method = method
        self.fold = fold
        self.undersampling = undersampling
        self.shuffle = shuffle
        self.maxEpochs_for_trainer=maxEpochs_for_trainer
        self.NNParamScales = {"node_num" : (5, 50), "learning_rate" : (-5, 5)} # As for learning_rata, 2**x
        self.RFParamScales = {"n_estimators" : (101, 3001), "max_features" : (2, 30)}
        self.windowSizeScales = (1, 19)
        self.geneScale=geneScale

    def create_folded_dataset(self, window_size):
        positive_dataset, negative_dataset = feature.create_dataset(self.bindingResidueData, self.pssmData, window_size)
        folded_dataset = dataset.FoldedDataset(positive_dataset, negative_dataset, fold=self.fold,
                                                undersampling=self.undersampling, shuffle=self.shuffle)
        return folded_dataset

    # Each gene is between 0 and n. n is greater than 0.
    def decode_chromosome(self, chromosome):
        if self.method == "neuralNetwork":
            node_num = int(round(chromosome[0]*(self.NNParamScales["node_num"][1]/self.geneScale[1])+self.NNParamScales["node_num"][0])
            leaning_rate = 2**int(round(chromosome[1]*(self.NNParamScales["learning_rate"][1]/self.geneScale[1])
                                               + self.NNParamScales["learning_rate"][0]))
            window_size = int(round(chromosome[2]*(self.windowSizeScales[1]/self.geneScale[1])+self.windowSizeScales[0]))
            return node_num, learning_rate, window_size
        elif self.method == "randomForest":
            n_estimators = int(round(chromosome[0]*(self.RFParamScales["n_estimators"][1]/self.geneScale[1])+self.RFParamScales["n_estimators"][0])
            max_features = int(round(chromosome[1]*(self.RFParamScales["max_features"][1]/self.geneScale[1])
                                               + self.RFParamScales["max_features"][0]))
            window_size = int(round(chromosome[2]*(self.windowSizeScales[1]/self.geneScale[1])+self.windowSizeScales[0]))
            return n_estimators, max_features, window_size

    def eval_func(self, chromosome):
        if len(chromosome) != 3:
            raise ValueError("len(chromosome) is must be 3 [{}]".format(len(chromosome)))
        if self.method == "neuralNetwork":
            return self.neuralNetwork_eval_func(chromosome)
        elif self.method == "randomForest":
            return self.randomForest_eval_func(chromosome)
            
    def neuralNetwork_eval_func(self, chromosome):
        node_num, learning_rate, window_size = self.decode_chromosome(chromosome)
        folded_dataset = self.create_folded_dataset(window_size)
        indim = 21 * window_size
        mean_AUC = 0
        for test_fold in xrange(self.fold):
            test_labels, test_dataset, train_labels, train_dataset = folded_dataset.get_test_and_training_dataset(test_fold)
            ds = SupervisedDataset(indim, 1)
            for i in len(train_labels):
                ds.addSample((train_dataset[i]), (train_labels[i],))
            net = buildNetwork(indim, node_num, 1)
            trainer = BackpropTrainer(net, ds, learningrate=learning_rate)
            trainer.trainUntilConvergence(maxEpochs=self.maxEpochs_for_trainer)
            decision_values = [net.activate(test_dataset[i]) for i in len(test_labels)]
            mean_AUC += validate.calculate_AUC(decision_values, test_labels)
        mean_AUC /= self.fold
        return mean_AUC

    def randomForest_eval_func(self, n_estimators, max_features, window_size):
        n_estimators, max_features, window_size = self.decode_chromosome(chromosome)
        folded_dataset = self.create_folded_dataset(window_size)
        indim = 21 * window_size
        mean_AUC = 0
        for test_fold in xrange(self.fold):
            test_labels, test_dataset, train_labels, train_dataset = folded_dataset.get_test_and_training_dataset(test_fold)
            clf = RandomForestClassifier(n_estimators=ntree, max_features=max_features)
            clf.fit(train_dataset, train_labels)
            probas = clf.predict_proba(test_dataset)
            decition_values = map(lambda x: x[1], probas) # Probability of being binding residue
            mean_AUC += validate.calculate_AUC(decision_values, test_labels)
        mean_AUC /= self.fold
        return mean_AUC
