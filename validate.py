#!/usr/bin/env python

"""
decision_values = [0.001, -0.3,  ... 0.2, 1.21]
correct_labels = [0, 0, 0, 1, 1, ... 1, 0]
"""

def calculate_TPR_FPR(TP, FP, TN, FN):
    TPR = TP / float(TP + FN)
    FPR = FP / float(FP + TN)
    return TPR, FPR


def update(result, TP, FP, TN, FN):
    if result[1] == 1:
        TP += 1
        FN -= 1
    elif result[1] == 0:
        FP += 1
        TN -= 1
    else:
        raise ValueError("correct_label is not 0 or 1")
    return TP, FP, TN, FN


def calculate_AUC(decision_values, correct_labels):
    results = []
    for i in xrange(len(correct_labels)):
        results.append((decision_values[i], correct_labels[i])))
    results.sort(reverse=True)
    positive_size = correct_labels.count(1)
    negative_size = correct_labels.count(0)
    TP, FP, TN, FN = 0, 0, negative_size, positive_size
    points = []
    prev_decval = -float('inf')
    for i, result in enumerate(results):
        if result[0] != prev_decval or i == len(correct_labels) - 1:
            TP, FP, TN, FN = update(result, TP, FP, TN, FN)
            points.append((calculate_TPR_FPR(TP, FP, TN, FN)))
        else:
            TP, FP, TN, FN = update(result, TP, FP, TN, FN)
    AUC = 0.0
    prev_point = (0, 0)
    for point in points:
        AUC += (point[0]+prev_point[0]) * (point[1]-prev_point[1]) / 2
    return AUC

