#!/bin/bash

function remove_if_exists {
  if [ -f $1 ]; then
    rm $1
  fi
}


if [ $1 ]; then if [ $1 == "-h" ]; then
  echo "Usage: test_whole_run.sh "
  echo "       test_whole_run.sh --randomscore"
  exit
  fi
fi

TESTDATA_DIR="test"
OUTPUT_DIR="output"
cd ".."
remove_if_exists "$OUTPUT_DIR/log.txt"
remove_if_exists "$OUTPUT_DIR/neuralNetwork_prediction_output.txt"
remove_if_exists "$OUTPUT_DIR/SVM_prediction_output.txt"
remove_if_exists "$OUTPUT_DIR/randomForest_prediction_output.txt"

BINDFILE="$TESTDATA_DIR/bindingData.txt"
if [ $1 ]; then if [ $1 == "--randomscore" ]; then
    PSSMFILE="$TESTDATA_DIR/pssms_random_score.txt"
  fi
else
  PSSMFILE="$TESTDATA_DIR/pssms_fixed_score.txt"
fi
echo $PSSMFILE
LOGFILE="$OUTPUT_DIR/log.txt"

SVM_CHROMOSOME_TSV="$OUTPUT_DIR/SVM_best_chromosome.tsv"
SVM_PKL="$OUTPUT_DIR/SVM.pkl"
SVM_PRED_OUTPUT="$OUTPUT_DIR/SVM_prediction_output.txt"
time python genetic_algorithm.py -m 'SVM' -b $BINDFILE -p $PSSMFILE \
                            -l $LOGFILE -o $SVM_CHROMOSOME_TSV
python create_model.py -i $SVM_CHROMOSOME_TSV -b $BINDFILE \
                        -p $PSSMFILE -o $SVM_PKL
python predict.py -i $SVM_CHROMOSOME_TSV -p $PSSMFILE \
                    -m $SVM_PKL -o $SVM_PRED_OUTPUT

NN_CHROMOSOME_TSV="$OUTPUT_DIR/neuralNetwork_best_chromosome.tsv"
NN_PKL="$OUTPUT_DIR/neuralNetwork.pkl"
NN_PRED_OUTPUT="$OUTPUT_DIR/neuralNetwork_prediction_output.txt"
time python genetic_algorithm.py -m 'neuralNetwork' -b $BINDFILE -p $PSSMFILE \
                            -l $LOGFILE -o $NN_CHROMOSOME_TSV
python create_model.py -i $NN_CHROMOSOME_TSV -b $BINDFILE \
                        -p $PSSMFILE -o $NN_PKL
python predict.py -i $NN_CHROMOSOME_TSV -p $PSSMFILE \
                    -m $NN_PKL -o $NN_PRED_OUTPUT

RF_CHROMOSOME_TSV="$OUTPUT_DIR/randomForest_best_chromosome.tsv"
RF_PKL="$OUTPUT_DIR/randomForest.pkl"
RF_PRED_OUTPUT="$OUTPUT_DIR/randomForest_prediction_output.txt"
time python genetic_algorithm.py -m 'randomForest' -b $BINDFILE -p $PSSMFILE \
                            -l $LOGFILE -o $RF_CHROMOSOME_TSV
python create_model.py -i $RF_CHROMOSOME_TSV -b $BINDFILE \
                        -p $PSSMFILE -o $RF_PKL
python predict.py -i $RF_CHROMOSOME_TSV -p $PSSMFILE \
                    -m $RF_PKL -o $RF_PRED_OUTPUT
