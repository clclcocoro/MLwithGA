#!/usr/bin/env python

import random
import copy


class FoldedDataset(object):

    def __init__(self, positive_dataset, negative_dataset, fold=5, undersampling=True, shuffle=True):
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

    def get_test_and_training_dataset(self, test_fold, undersampling_training_dataset=True):
        if test_fold < 0 or self.fold <= test_fold:
            raise ValueError("test_fold [{}] must be between 0 and self.fold-1 [{}]".format(test_fold, self.fold-1))
        test_positive_size = len(self.folded_positive_dataset[test_fold])
        test_negative_size = len(self.folded_negative_dataset[test_fold])
        test_labels = [1] * test_positive_size + [0] * test_negative_size
        test_dataset = self.folded_positive_dataset[test_fold] + self.folded_negative_dataset[test_fold]
        positive_train_labels = [1] * (self.positive_size - test_positive_size)
        negative_train_labels = [0] * (self.negative_size - test_negative_size)
        positive_train_dataset = []
        negative_train_dataset = []
        for i in xrange(self.fold):
            if i != test_fold:
                positive_train_dataset += self.folded_positive_dataset[i]
                negative_train_dataset += self.folded_negative_dataset[i]
        if undersampling_training_dataset:
            if self.positive_size - test_positive_size > 1500:
                positive_train_labels = [1] * 1500
                positive_train_dataset = positive_train_dataset[:1500]
            if self.negative_size - test_negative_size > 1500:
                negative_train_labels = [0] * 1500
                negative_train_dataset = negative_train_dataset[:1500]
        train_labels = positive_train_labels + negative_train_labels
        train_dataset = positive_train_dataset + negative_train_dataset
        return test_labels, test_dataset, train_labels, train_dataset
