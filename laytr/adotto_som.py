import re
import multiprocessing
from functools import partial

import joblib
import pysam
import truvari
import numpy as np
#from sklearn.preprocessing import normalize

NUCFILT = re.compile("[^ATCG]")

class kmer_featurization:

  def __init__(self, k):
    """
    seqs: a list of DNA sequences
    k: the "k" in k-mer
    """
    self.k = k
    self.letters = ['A', 'T', 'C', 'G']
    self.multiplyBy = 4 ** np.arange(k-1, -1, -1) # the multiplying number for each digit position in the k-number system
    self.n = 4**k # number of possible k-mers

  def obtain_kmer_feature_for_a_list_of_sequences(self, seqs, write_number_of_occurrences=False):
    """
    Given a list of m DNA sequences, return a 2-d array with shape (m, 4**k) for the 1-hot representation of the kmer features.

    Args:
      write_number_of_occurrences:
        a boolean. If False, then in the 1-hot representation, the percentage of the occurrence of a kmer will be recorded; otherwise the number of occurrences will be recorded. Default False.    
    """
    kmer_features = []
    for seq in seqs:
      this_kmer_feature = self.obtain_kmer_feature_for_one_sequence(seq.upper(), write_number_of_occurrences=write_number_of_occurrences)
      kmer_features.append(this_kmer_feature)

    kmer_features = np.array(kmer_features)

    return kmer_features

  def obtain_kmer_feature_for_one_sequence(self, seq, write_number_of_occurrences=False):
    """
    Given a DNA sequence, return the 1-hot representation of its kmer feature.

    Args:
      seq: 
        a string, a DNA sequence
      write_number_of_occurrences:
        a boolean. If False, then in the 1-hot representation, the percentage of the occurrence of a kmer will be recorded; otherwise the number of occurrences will be recorded. Default False.
    """
    number_of_kmers = len(seq) - self.k + 1

    kmer_feature = np.zeros(self.n)

    for i in range(number_of_kmers):
      this_kmer = seq[i:(i+self.k)]
      this_numbering = self.kmer_numbering_for_one_kmer(this_kmer)
      kmer_feature[this_numbering] += 1

    if not write_number_of_occurrences:
      kmer_feature = kmer_feature / number_of_kmers

    return kmer_feature

  def kmer_numbering_for_one_kmer(self, kmer):
    """
    Given a k-mer, return its numbering (the 0-based position in 1-hot representation)
    """
    digits = []
    for letter in kmer:
      digits.append(self.letters.index(letter))

    digits = np.array(digits)

    numbering = (digits * self.multiplyBy).sum()

    return numbering

def get_features(region, kf, ref_fn):
    chrom, start, end = region
    ref = pysam.FastaFile(ref_fn)
    seq = NUCFILT.sub("", ref.fetch(chrom, start, end))
    ret1 = kf.obtain_kmer_feature_for_one_sequence(seq)
    gcpct = (seq.count('G') + seq.count('C')) / len(seq)
    return ret1, gcpct

def iter_regions(fn):
    fh = truvari.opt_gz_open(fn)
    for line in fh:
        data = line.strip().split('\t')
        chrom, start, end = data[:3]
        start = int(start)
        end = int(end)
        yield chrom, start, end


def regions_to_kmers(regions, reference, k=3, ordered=False, nproc=4):
    kf = kmer_featurization(k)
    m_caller = partial(get_features, kf=kf, ref_fn=reference)
    with multiprocessing.Pool(nproc) as pool:
        if ordered:
            result = [_ for _ in pool.map(m_caller, regions)]
        else:
            result = [_ for _ in pool.imap_unordered(m_caller, regions)]
        pool.close()
        pool.join()
    feats = np.vstack([_[0] for _ in result])
    gcs = np.array([_[1] for _ in result])
    return feats, gcs

if __name__ == '__main__':
    ref = "/Users/english/code/references/grch38/GRCh38_1kg_mainchrs.fa"
    regions = "/Users/english/code/adotto/regions/adotto_TRregions_v1.1.bed.gz"
    m_reg = iter_regions(regions)
    data = regions_to_kmers(m_reg, ref, ordered=True, nproc=8)
    joblib.dump(data, "adotto_TRregions_v1.1_3mers.jl")

#3) Try to highlight by GC percent?
