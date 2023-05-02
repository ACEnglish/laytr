Laytr - Library for variant benchmarking stratification

# Install

```bash
git clone https://github.com/ACEnglish/laytr.git
cd laytr
python3 -m pip install .
```

# Quick Start

```bash
usage: laytr [-h] CMD ...

laytr v0.0.1 Library for variant benchmarking stratification

Available commands:
    kfeat     Create kmer featuration of genomic regions
    map       Map kfeats to a SOM and report neurons
    tru2ga    Consolidate truvari outputs and annotate with GA4GH intermediates
    giabSV06  GIAB SV v0.6 report on a truvari directory
    giabTR    GIAB TR report on a refine.regions.txt

positional arguments:
  CMD         Command to execute
  OPTIONS     Options to pass to the command
```

# Usage

## SOMs
Create kmer featurization of genomic regions
```bash
laytr kfeat --regions example/example.bed --reference example/reference.fa --output chr22_kfeat.jl
```

You can then create a SOM out of the kmer featurization. See `notebooks/KmerSom.ipynb` for an example.

Next, you can map regions to the SOM.
```bash
laytr map --input chr22_kfeat.jl --som soms/adotto_TRv1.1_3mers.som --output chr22_kfeat_map.jl
```
 
Finally, you can visualize your SOM with `notebooks/SOMPlotExample.ipynb`.

## tru2ga

Consolidate truvari results into GA4GH intermediates.
```bash
laytr tru2ga -i truvari_results/ --with-refine -o result_
# Creates `result_truth.vcf.gz` and `result_query.vcf.gz`
```

## giabTR
Creates an html report from truvari's `refine.regions.txt` on the GIAB TR benchmark.

```bash
# Make a truvari result
truvari bench -b giab_tr.vcf.gz -c tr_caller.vcf.gz --includebed giab_tr.bed -o bench/
truvari refine --reference grch38.fa bench/
laytr giabTR -r bench/refine.regions.txt \
	     -b giab_tr.bed \
	     -t adotto_tr_catalog.bed \
	     -s adotto_TRv1.1_3mers.som # available from laytr repo \
	     -m adotto_TRv1.1_3mers.map # available from laytr repo \
	     -o giabTR_report.html
```
See `examples/giabTR_report.html`
([rendered](https://htmlpreview.github.io/?https://github.com/ACEnglish/laytr/blob/develop/examples/giabTR_report.html))
for a look at the report.

## giabSV06

blank

# ToDo
- make SurbScore (no cli)
