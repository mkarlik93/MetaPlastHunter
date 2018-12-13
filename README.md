# MetaPlastHunter rc
## MetaPlastHunter for read classificaton


The tool for accurate classifing and quantifing reads of plastid origin from metagenomes.

MetaPlstHunter will search the provided sequences (Fastq) or sequence alignment file (SAM) and acurate classify them using MetaPlastHunter algorithm :

MetaPlastHunter provides:

* Searching raw reads against constantly updated and curated DB of protistan plastid genomes, quantifing and visualization using Krona tools
* Analyzing aligned reads aligned DB in SAM format, quantifing and vizualization using Krona tools
* Searching raw reads against DB using heuristic k - mer matching module, quantifing and vizualization using Krona tools

### Installation

##### github



##### pip

However, to use all features of MetaPlastHunter extra binary applications are required:

(bbtools)[https://jgi.doe.gov/data-and-tools/bbtools/]

(Krona tools)[https://github.com/marbl/Krona/wiki]

### Manual

A manual is available in the form of the wiki here on GitHub.

### License

MetaPlastHunter is licensed under the GNU GPL v3+. See LICENSE.txt for further details. MPH makes use of 18S rRNA sequences sourced from the SILVA database, which employs a dual licensing model. See SILVA.LICENCE.txt for further details.
