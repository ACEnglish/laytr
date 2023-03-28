Laytr - Library for variant benchmarking stratification

There is no quick start. This is a library for helping bioinformaticians perform analysis. Please read the documentation
and tutorials.

# Quick Start

The core functionality of laytr is the library. See wiki for description of jupyter notebook tutorials.

Additionally, some tools for creating data are available in the command line.
```bash
laytr -h

laytr kfeat regions reference kfeat.jl
# make som > som.pkl
laytr winner kfeat.jl som.pkl kfeat_winners.jl
```
 
Then you can import laytr.plot_hex or laytr.heatmap or umatrix or what?

# Modules

- SeqSom - Self-organizing map based on kmer featurization of DNA.
- TRsom - Self-organizing map based on Adotto TR region features
- surbscore - Exploratory data analysis technique with permutation tests on chi2 scores with multiple test correction

# Use cases -

1. Given regions and a reference, build kmer feats
2. Take those kmer feats and build a som
3. Map regions to the SOM to build columns, Neuron.
4. Describe how to put makers onto neurons to identify neighborhoods
5. Describe umatrix for how to color the neurons
6. Describe how to heatmap the umatrix e.g. average something mnother


# Tutorial

import laytr

laytr.kfeat(seq, k)
laytr.kfeats(seqs, k)

laytr.map(feat, map)
laytr.maps(feats, map)

laytr.hex_plot(

Need:
Loaders for maps (assume pkl files from minisom).

I have the plot code somewhere, just need to clean up that documentation.

And then I want to formalize the custom heatmap stuff.

Once I have all of those, I can try to make instructions for how to do it with
refine.regions.txt

And then I want to make another map of all the TRregion features. See if we got interesting patterns
that will help find good/bad sites?

A data file that holds regions and the neuron they map to.
This will help... so something, I think

# I actually don't want to make a CLI. This is a library for helping not a tool for doing.
seqsom kmerfeats regions reference.bed -k 3 -o data.jl --ordered
turn regions of the reference into kmer features and optionally have them ordered
This data will then need documentation (and an example notebook) to work with



1. layTR
  1. SeqSom
  2. TRregion Som
  3. SurbScore
