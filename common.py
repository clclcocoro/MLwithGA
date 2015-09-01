#!/usr/bin/env python

def get_method_and_genes(best_chromosome_file):
    with open(best_chromosome_file) as fp:
        for i, line in enumerate(fp):
            if i == 1:
                parts = line.rstrip().split('\t')
                return parts[:4] # [method, gene1, gene2, gene3]


def get_method_genes_decision_value(best_chromosome_file):
    with open(best_chromosome_file) as fp:
        for i, line in enumerate(fp):
            if i == 1:
                parts = line.rstrip().split('\t')
                return parts # [method, gene1, gene2, gene3, decision_value]
