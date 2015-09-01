#!/usr/bin/env python

"""
decision_values = [0.001, -0.3,  ... 0.2, 1.21]
correct_labels = [0, 0, 0, 1, 1, ... 1, 0]
"""

import math

def calculate_TPR_FPR(TP, FP, TN, FN):
    FPR = FP / float(FP + TN)
    TPR = TP / float(TP + FN)
    return FPR, TPR


def calculate_MCC(TP, FP, TN, FN):
    if (TP+FN)*(TP+FP)*(TN+FP)*(TN+FN) == 0:
        return 0
    return (TP*TN-FP*FN) / math.sqrt((TP+FN)*(TP+FP)*(TN+FP)*(TN+FN))


def update_decision_value_and_max_mcc(decision_value_and_max_mcc, decision_value, TP, FP, TN, FN):
    current_mcc = calculate_MCC(TP, FP, TN, FN)
    if decision_value_and_max_mcc[1] <= current_mcc:
        decision_value_and_max_mcc = [decision_value, current_mcc]
    return decision_value_and_max_mcc


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
        results.append((decision_values[i], correct_labels[i]))
    results.sort(reverse=True)
    positive_size = correct_labels.count(1)
    negative_size = correct_labels.count(0)
    TP, FP, TN, FN = 0, 0, negative_size, positive_size
    points = []
    prev_decval = float('inf')
    decision_value_and_max_mcc = [0, 0] # [decision_value, max_mcc]
    for i, result in enumerate(results):
        if i == len(correct_labels) - 1: # Final result
            if result[0] != prev_decval:
                points.append((calculate_TPR_FPR(TP, FP, TN, FN)))
                decision_value_and_max_mcc = update_decision_value_and_max_mcc(decision_value_and_max_mcc, result[0], TP, FP, TN, FN)
                TP, FP, TN, FN = update(result, TP, FP, TN, FN)
                points.append((calculate_TPR_FPR(TP, FP, TN, FN)))
                decision_value_and_max_mcc = update_decision_value_and_max_mcc(decision_value_and_max_mcc, result[0], TP, FP, TN, FN)
            else:
                TP, FP, TN, FN = update(result, TP, FP, TN, FN)
                points.append((calculate_TPR_FPR(TP, FP, TN, FN)))
                decision_value_and_max_mcc = update_decision_value_and_max_mcc(decision_value_and_max_mcc, result[0], TP, FP, TN, FN)
            break
        if i != 0 and result[0] != prev_decval:
            points.append((calculate_TPR_FPR(TP, FP, TN, FN)))
            decision_value_and_max_mcc = update_decision_value_and_max_mcc(decision_value_and_max_mcc, result[0], TP, FP, TN, FN)
            TP, FP, TN, FN = update(result, TP, FP, TN, FN)
        else: # the same decision value 
            TP, FP, TN, FN = update(result, TP, FP, TN, FN)
        prev_decval = result[0]
    AUC = 0.0
    prev_point = (0, 0)
    for point in points:
        AUC += (point[1]+prev_point[1]) * (point[0]-prev_point[0]) / 2
        prev_point = point
    return AUC, decision_value_and_max_mcc
