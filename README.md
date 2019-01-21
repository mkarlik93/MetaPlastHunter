# MetaPlastHunter rc
## MetaPlastHunter for read classificaton


The tool for accurate classifing and quantifing reads of plastid origin from metagenomes.

MetaPlstHunter will search the provided sequences (Fastq) or sequence alignment file (SAM) and acurate classify them using MetaPlastHunter algorithm :

MetaPlastHunter provides:

* Searching raw reads against constantly updated and curated DB of protistan plastid genomes, quantifing and visualization using Krona tools
* Analyzing aligned reads aligned DB in SAM format, quantifing and vizualization using Krona tools
* Searching raw reads against DB using heuristic k - mer matching module, quantifing and vizualization using Krona tools

usage: MetaPlastHunter [-h] [--taxonomic_classification]
                       [--rapid_classification] [--settings settings]
                       [--sam_assign] [--in_1 input] [--in_2 [IN_2]]
                       [--output OUTPUT] [--threads [THREADS]] [--check]

Version 1.0.0

MetaPlastHunter -

 The efficient and accurate plastid reads classification pipeline.

Quantitative aproach for eukaryotic metagenomics.

Available workflows:

[--taxonomic_classification/-C] Searching, Classification, Visualization

[--rapid_classification, -Acc] Use it to lunch pipeline with in exact k-mer matching preliminary classification

[--sam_assign, -A]   Sequence alignment file (SAM) classification

[--check]   Check settings and calculate empirical treshold if nessecary



### Manual

A manual is available in the form of the [wiki](https://github.com/mkarlik93/MetaPlastHunter/wiki) here on GitHub.

### License

MetaPlastHunter is licensed under the GNU GPL v3+. See LICENSE.txt for further details. MPH makes use of 18S rRNA sequences sourced from the SILVA database, which employs a dual licensing model. See SILVA.LICENCE.txt for further details.
