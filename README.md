Laytr - Library for variant benchmarking stratification

# Quick Start

The core functionality of laytr is the library. See wiki for description of jupyter notebook tutorials.

There are also some command line tools for creating data.

```bash
laytr kfeat --regions example/example.bed --reference example/reference.fa --output chr22_kfeat.jl
```

You can then create a SOM out of the kmer featurization. See `notebooks/KmerSom.ipynb` for an example.

Next, you can map regions to the SOM.
```bash
laytr map --input chr22_kfeat.jl --som soms/adotto_TRv1.1_3mers.som --output chr22_kfeat_map.jl
```
 
Finally, you can visualize your SOM with `notebooks/SOMPlotExample.ipynb`.

# ToDo
- make notebooks/KmerSom.ipynb
- make SurbScore (no cli)
- 
