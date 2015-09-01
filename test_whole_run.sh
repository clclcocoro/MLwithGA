#!/bin/bash

python genetic_algorithm.py -m 'neuralNetwork' -b test/bindingData.txt -p test/pssms.txt \
                            -l output/log.txt -o output/neuralNetwork_best_chromosome.tsv
python create_model.py -i output/neuralNetwork_best_chromosome.tsv -b test/bindingData.txt \
                        -p test/pssms.txt -o output/neuralNetwork.pkl
python predict.py -i output/neuralNetwork_best_chromosome.tsv -p test/pssms.txt \
                    -m output/neuralNetwork.pkl -o output/neuralNetwork_prediction_output.txt
python genetic_algorithm.py -m 'SVM' -b test/bindingData.txt -p test/pssms.txt \
                            -l output/log.txt -o output/SVM_best_chromosome.tsv
python create_model.py -i output/SVM_best_chromosome.tsv -b test/bindingData.txt \
                        -p test/pssms.txt -o output/SVM.pkl
python predict.py -i output/SVM_best_chromosome.tsv -p test/pssms.txt \
                    -m output/SVM.pkl -o output/SVM_prediction_output.txt
python genetic_algorithm.py -m 'randomForest' -b test/bindingData.txt -p test/pssms.txt \
                            -l output/log.txt -o output/randomForest_best_chromosome.tsv
python create_model.py -i output/randomForest_best_chromosome.tsv -b test/bindingData.txt \
                        -p test/pssms.txt -o output/randomForest.pkl
python predict.py -i output/randomForest_best_chromosome.tsv -p test/pssms.txt \
                    -m output/randomForest.pkl -o output/randomForest_prediction_output.txt
