Laytr - GIAB TR stratification toolkit


There is no quick start. This is a library for helping bioinformaticians perform analysis. Please read the documentation
and tutorials.

# Quick Start

The core functionality of laytr is the library. See wiki for description of jupyter notebook tutorials.

Additionally, some tools for creating data are available in the command line.
```bash
laytr -h
```

# Modules

- SeqSom - Self-organizing map based on kmer featurization of DNA.
- TRsom - Self-organizing map based on Adotto TR region features
- surbscore - Exploratory data analysis technique with permutation tests on chi2 scores with multiple test correction

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
