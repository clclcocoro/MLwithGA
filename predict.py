#!/usr/bin/env python

"""Create model with optimized parameters and training data.


Usage:
  create_model.py -i <best_chromosome_file> -p <pssms_file> -m <pickled_model_file> -o <prediction_output_file>
  create_model.py (-h | --help)
  create_model.py --version

Options:
  -i          best chromosome file generated by genetic_algorithm.py. 
  -p          pssms file.
  -m          pickled model file generated by create_model.py
  -o          prediction output.
  -h --help   Show this screen.
  --version   Show version.
"""


import pickle
from docopt import docopt
from pybrain.tools.shortcuts     import buildNetwork
from pybrain.supervised.trainers import BackpropTrainer
from pybrain.datasets            import SupervisedDataSet
from sklearn.ensemble            import RandomForestClassifier 
from sklearn import svm
import feature
import dataset
import common


def predict_with_NN_classifier(net, my_decision_value, test_dataset):
    results = {}
    decision_values = [net.activate(test_dataset[i]) for i in xrange(len(test_dataset))]
    decision_values = map(lambda x: x[0], decision_values)
    predicted_labels = map(lambda x: 1 if x >= my_decision_value else 0, decision_values)
    results['label'] = predicted_labels
    results['decision_value'] = decision_values
    return results


def predict_with_RF_classifier(clf, my_decision_value, test_dataset):
    results = {}
    probas = clf.predict_proba(test_dataset)
    decision_values = map(lambda x: x[1], probas) # Probability of being binding residue
    predicted_labels = [1 if decision_value >= my_decision_value else 0 for decision_value in decision_values]
    results['label'] = predicted_labels
    results['decision_value'] = decision_values
    return results


def predict_with_SVM_classifier(clf, my_decision_value, test_dataset):
    results = {}
    decision_values = clf.decision_function(test_dataset)
    decision_values = map(lambda x: x[0], decision_values)
    predicted_labels = [1 if decision_value >= my_decision_value else 0 for decision_value in decision_values]
    results['label'] = predicted_labels
    results['decision_value'] = decision_values
    return results
 

def predict(method, pickled_model_file, my_decision_value, test_dataset):
    pkl_file = open(pickled_model_file, 'rb')
    clf_or_net = pickle.load(pkl_file)
    pkl_file.close()
    if method == "neuralNetwork":
        return predict_with_NN_classifier(clf_or_net, my_decision_value, test_dataset)
    elif method == "randomForest":
        return predict_with_RF_classifier(clf_or_net, my_decision_value, test_dataset)
    elif method == "SVM":
        return predict_with_SVM_classifier(clf_or_net, my_decision_value, test_dataset)
    else:
        raise ValueError("method must be neuralNetwork or randomForest or SVM [{}]".format(method))


def write_to_output_file(prediction_output_file, user_defined_header, results):
    with open(prediction_output_file, 'a') as fp:
        idx = 0
        for predicted_label, decision_value in zip(results['label'], results['decision_value']):
            fp.write("{} {} {} {}\n".format(user_defined_header, idx, predicted_label, decision_value))
            idx += 1

if __name__ == "__main__":
    arguments = docopt(__doc__)
    best_chromosome_file = arguments['<best_chromosome_file>']
    pssms_file = arguments['<pssms_file>']
    pickled_model_file = arguments['<pickled_model_file>']
    prediction_output_file = arguments['<prediction_output_file>']
    method_genes_decision_value = common.get_method_genes_decision_value(best_chromosome_file)
    method, window_size, my_decision_value = method_genes_decision_value[0], int(method_genes_decision_value[3]), float(method_genes_decision_value[4])
    pssmData = feature.parse_pssms_file(pssms_file)
    for user_defined_header in pssmData.get_uniprotURIs():
        pssm = pssmData.get_PSSMRecord(user_defined_header)
        test_dataset = feature.create_feature_vectors(pssm, window_size)
        results = predict(method, pickled_model_file, my_decision_value, test_dataset)
        write_to_output_file(prediction_output_file, user_defined_header, results)
