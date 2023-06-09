"""
Given a SOM and a DataFrame of features, pick the winning neuron and saved result with joblib
"""
import joblib
import argparse

import numpy as np

def parse_args(args):
    """
    Argument parser
    """
    parser = argparse.ArgumentParser(prog="map", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="Input dataframe")
    parser.add_argument("-s", "--som", type=str, required=True,
                        help="Self-organizing map")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Output joblib file")
    args = parser.parse_args(args)
    return args

def map_to_som(feats, som):
    """
    Map features to som, return numpy array
    """
    return np.array([som.winner(_) for _ in feats])

def map_main(args):
    """
    Main
    """
    args = parse_args(args)
    data = joblib.load(args.input)
    som = joblib.load(args.som)
    joblib.dump({"index": data["index"], 
                 "map": map_to_som(data["features"], som)},
                args.output)
