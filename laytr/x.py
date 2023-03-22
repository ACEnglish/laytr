from itertools import product
NUCS = {'A':0, 'G':1, 'C':2, 'T':3}
def kmer_index(kmer: str) -> int:
    """
    Get the array index of a kmer
    """
    index = 0
    for pos, nuc in enumerate(kmer):
        index += NUCS[nuc] << pos * 2
    return index

def generate_kmers(k):
    return [''.join(kmer) for kmer in product(list(NUCS.keys()), repeat=k)]


k = 3
for i in generate_kmers(k):
    print(i, kmer_index(i))

