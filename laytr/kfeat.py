"""
Given regions and a reference, create a numpy 2D array of the kmer featurization
"""
import re
import argparse
import multiprocessing
from functools import partial
from typing import Iterator, List, Tuple, Type

import pysam
import joblib
import truvari
import numpy as np
from numpy.typing import ArrayLike

NUCS = {'A':0, 'G':1, 'C':2, 'T':3}
NUCFILT = re.compile("[^ATCG]")

def kmer_index(kmer: str) -> int:
    """
    Get the array index of a kmer
    """
    index = 0
    for pos, nuc in enumerate(kmer):
        index += NUCS[nuc] << pos * 2
    return index

def kfeat(seq: str, k: int = 3,
          normalize: bool = True) -> ArrayLike:
    """
    Return the kmer featurization of a sequence
    """
    ret = np.zeros(4 ** k)
    number_of_kmers = len(seq) - k + 1
    for i in range(number_of_kmers):
        ret[kmer_index(seq[i:i+k])] += 1

    if normalize:
        ret /= number_of_kmers
    return ret

def kfeat_seqs(seqs: List[str], k: int = 3,
           normalize: bool = True) -> ArrayLike:
    """
    Return the kmer featurization of sequences
    """
    return np.array([kfeat(_, k, normalize) for _ in seqs])

def get_features(region: Tuple[str, int, int],
                 ref_fn: str, k: int = 3,
                 normalize: bool = True) -> ArrayLike:
    """
    Return the kmer features for a single region
    """
    chrom, start, end = region
    ref = pysam.FastaFile(ref_fn)
    seq = NUCFILT.sub("", ref.fetch(chrom, start, end))
    return kfeat(seq, k, normalize)

class RegionIter:
    """
    Iterate a tab-delimited file with chrom, start, and end columns
    while iterating, will collect and store positions in self.regions
    """
    def __init__(self, fn: str):
        self.file_name = fn
        self.fh = truvari.opt_gz_open(self.file_name)
        self.regions : List[Tuple[str, int, int]] = []

    def __iter__(self) -> Iterator[Tuple[str, int, int]]:
        for line in self.fh:
            data = line.strip().split('\t')
            chrom, start, end = data[:3]
            start = int(start)
            end = int(end)
            self.regions.append((chrom, start, end))
            yield self.regions[-1]

def regions_to_kmers(regions: Type[RegionIter], reference: str,
                     k: int = 3, nproc: int = 1):
    """
    Get kmer featurization for a set of regions over a reference
    """
    m_caller = partial(get_features, ref_fn=reference, k=k)
    with multiprocessing.Pool(nproc) as pool:
        result = [_ for _ in pool.map(m_caller, regions)]
        pool.close()
        pool.join()
    feats = np.vstack([result])
    return feats

def parse_args(args):
    """
    Argument parser
    """
    parser = argparse.ArgumentParser(prog="kfeat", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-r", "--regions", type=str, required=True,
                        help="Bed file of regions to fecth")
    parser.add_argument("-f", "--reference", type=str, required=True,
                        help="Reference genome fasta")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Output joblib file")
    parser.add_argument("-k", "--kmer", type=int, default=3,
                        help="Size of kmer (%(default)s)")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Number of threads to use (%(default)s)")
    args = parser.parse_args(args)
    return args

def kfeat_main(args):
    """
    Main
    """
    args = parse_args(args)
    m_reg = RegionIter(args.regions)
    data = regions_to_kmers(m_reg, args.reference, k=args.kmer, nproc=args.threads)
    joblib.dump({'index': m_reg.regions, 'features':data}, args.output)
