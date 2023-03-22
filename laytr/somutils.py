
import re
import multiprocessing
from typing import List
from functools import partial

import pysam
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

def kfeats(seqs: List[str], k: int = 3,
           normalize: bool = True) -> ArrayLike:
    """
    Return the kmer featurization of sequences
    """
    return np.array([kfeats(_, k, normalize) for _ in seqs])

def get_features(region: tuple[str, int, int],
                 ref_fn: str, k: int = 3,
                 normalize: bool = True) -> ArrayLike:
    """
    Return the kmer features for a single region
    """
    chrom, start, end = region
    ref = pysam.FastaFile(ref_fn)
    seq = NUCFILT.sub("", ref.fetch(chrom, start, end))
    return kfeat(seq, k, normalize)

def iter_regions(fn: str):
    """
    Iterate a tab-delimited file with chrom, start, and end columns
    """
    fh = truvari.opt_gz_open(fn)
    for line in fh:
        data = line.strip().split('\t')
        chrom, start, end = data[:3]
        start = int(start)
        end = int(end)
        yield chrom, start, end

def regions_to_kmers(regions: List[tuple[str, int, int]], reference: str,
                     k: int = 3, ordered: bool = False, nproc: int = 4):
    """
    Get kmer featurization for a set of regions over a reference
    """
    m_caller = partial(get_features, ref_fn=reference, k=k)
    with multiprocessing.Pool(nproc) as pool:
        if ordered:
            result = [_ for _ in pool.map(m_caller, regions)]
        else:
            result = [_ for _ in pool.imap_unordered(m_caller, regions)]
        pool.close()
        pool.join()
    feats = np.vstack([result])
    return feats

