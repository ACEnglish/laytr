"""
Given a som (pkl?) and a kfeat dataframe (joblib), pick the winning neuron. Reports an array saved with joblib
"""
import pickle
import joblib
import argparse

import numpy as np

def parse_args(args):
    """
    Argument parser
    """
    parser = argparse.ArgumentParser(prog="place", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="Input dataframe (.jl)")
    parser.add_argument("-s", "--som", type=str, required=True,
                        help="Self-organizing map (.pkl)")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Output joblib file")
    args = parser.parse_args(args)
    return args


def place_main(args):
    """
    Main
    """
    args = parse_args(args)
    feats = joblib.load(args.input)
    som = pickle.load(open(args.som, 'rb'))
    neurons = np.array([som.winner(_) for _ in feats])
    joblib.dump(neurons, args.output)
