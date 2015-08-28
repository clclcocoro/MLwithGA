#!/usr/bin/env python

import random
import copy


class FoldedDataset(object):

    def __init__(self, positive_dataset, negative_dataset, fold=5, undersampling=True, shuffle=True):
        self.original_positive_dataset = copy.deepcopy(positive_dataset)
        self.original_negative_dataset = copy.deepcopy(negative_dataset)
        self.fold = fold
        positive_size = len(positive_dataset)
        negative_size = len(negative_dataset)
        if shuffle:
            my_seed=12345
            random.seed(my_seed)
            random.shuffle(positive_dataset)
            random.seed(my_seed)
            random.shuffle(negative_dataset)
        if undersampling:
            if negative_size > 2 * positive_size:
                negative_dataset = negative_dataset[:2*positive_size]
            if positive_size > 2 * negative_size:
                positive_dataset = positive_dataset[:2*negative_size]
        self.positive_size = len(positive_dataset)
        self.negative_size = len(negative_dataset)
        self.folded_positive_dataset = self.folding(self.positive_size, positive_dataset, fold)
        self.folded_negative_dataset = self.folding(self.negative_size, negative_dataset, fold)

    def folding(self, size, dataset, fold):
        p = size / fold
        r = size % fold
        if r == 0:
            folded_dataset = [dataset[p*i:p*(i+1)] for i in xrange(fold)]
        else:
            folded_dataset = [dataset[p*i:p*(i+1)] for i in xrange(fold)]
            for i in xrange(r):
                folded_dataset[i] += [dataset[fold*p+i]]
        return folded_dataset

    def get_folded_positive_dataset(self):
        return self.folded_positive_dataset

    def get_folded_negative_dataset(self):
        return self.folded_negative_dataset

    def get_original_positive_dataset(self):
        return self.original_positive_dataset

    def get_original_negative_dataset(self):
        return self.original_negative_dataset

    def get_test_and_training_dataset(self, test_fold):
        if test_fold < 0 or self.fold <= test_fold:
            raise ValueError("test_fold [{}] must be between 0 and self.fold-1 [{}]".format(test_fold, self.fold-1))
        test_positive_size = len(self.folded_positive_dataset[test_fold])
        test_negative_size = len(self.folded_negative_dataset[test_fold])
        test_labels = [1] * test_positive_size + [0] * test_negative_size
        test_dataset = self.folded_positive_dataset[test_fold] + self.folded_negative_dataset[test_fold]
        train_labels = [1] * (self.positive_size - test_positive_size) + [0] * (self.negative_size - test_negative_size)
        train_dataset = []
        for i in xrange(self.fold):
            if i != test_fold:
                train_dataset += self.folded_positive_dataset[i]+self.folded_negative_dataset[i]
        return test_labels, test_dataset, train_labels, train_dataset
