#!/usr/bin/env python

"""Run Genetic Algorithm.


Usage:
  genetic_algorithm.py -m <method> -b <binding_residue_file> -p <pssms_file> -o <output_file>
  genetic_algorithm.py (-h | --help)
  calculateAUC_from_blastout.py --version

Options:
  -m          machine learning method ("neuralNetwork" or "randomForest"). 
  -b          binding residue file.
  -p          pssms file.
  -o          output file.
  -h --help   Show this screen.
  --version   Show version.
"""


from pyevolve import GSimpleGA
from pyevolve import G1DList
from pyevolve import Selectors
from pyevolve import Initializators, Mutators
from docopt import docopt
import cross_validation


def run_ga(method, cross_validation, rangemin=0, rangemax=10):
    genome = G1DList.G1DList(3)
    genome.setParams(rangemin=rangemin, rangemax=rangemax)
    if cross_validation.method == "neuralNetwork":
        genome.evaluator.set(cross_validation.neuralNetwork_eval_func)
    elif cross_validation.method == "randomForest":
        genome.evaluator.set(cross_validation.randomForest_eval_func)
    ga = GSimpleGA.GSimpleGA(genome)
    ga.selector.set(Selectors.GRouletteWheel)
    ga.setGenerations(3)
    ga.setPopulationSize(25)
    ga.setMutationRate(0.05)
    ga.setCrossoverRate(0.8)
    ga.evolve(freq_stats=1)

    print ga.bestIndividual()
    best_chromosome = ga.bestIndividual()
    return best_chromosome

if __name__ == "__main__":
    arguments = docopt(__doc__)
    method = arguments['method']
    bindres_file = arguments['binding_residue_file']
    pssms_file = arguments['pssms_file']
    output_file = arguments['output']
    crossValidation = cross_validation.CrossValidation(bindres_file, pssms_file, method)
    best_chromosome = run_ga(crossValidation)
    with open(output_file, "w") as fp:
        if crossValidation.method == "neuralNetwork":
            fp.write("#method\tnode_num\tlearning_rate\twindow_size\n")
        elif crossValidation.method == "randomForest":
            fp.write("#method\tn_estimators\tmax_features\twindow_size\n")
        gene1, gene2, gene3 = cross_validation.decode_chromosome(best_chromosome)
        fp.write("{}\t{}\t{}\t{}".format(method, gene1, gene2, gene3))
