SOMs

First, you need to create a SOM. You create this off of a set of features. Laytr comes packaged with a way to make
kmer_featurization of genomic regions. See `laytr kfeat`. 


See `notebooks/KmerSOM.ipynb` for an example of how to create a SOM using minisom.
This will create a `.som` file.

We can then map our features to the som to collect the neuron to which each genomic region maps using `laytr map`.

This will create a joblib output with data structure 
{
 "index": 
 "neurons": numpy.array with x/y columns identifying to which neuron each genomic region maps
}

Note that kfeats' index is the genomic regions (chrom, start, end). But any joblib saved data can be put through `laytr map`
as long as it has the structure:
{
 "index":
 "features":
}

From here, go to `notebooks/SOMPlotExample.ipynb` to see how to visualize data.


SurbScore
This is api only.


GIAB


I need to make a tool for truv2ga4gh
Will take the truvari output directory (with/without phab)
And add the needed intermediates such that it can be fed into whatever ga4gh application
links to hap.py and https://github.com/ga4gh/benchmarking-tools/tree/master/doc/ref-impl

