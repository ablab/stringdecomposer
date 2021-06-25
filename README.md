[![Anaconda-Server Badge](https://anaconda.org/bioconda/stringdecomposer/badges/installer/conda.svg)](https://anaconda.org/bioconda/stringdecomposer)

# StringDecomposer

## Version 1.0

As an input StringDecomposer algorithm takes the set of monomers (typically, alpha satellites) and a genomic segment (assembly, Oxford Nanopore or a PacBio HiFi read) that contains a tandem repeat consisting of the given monomers.
StringDecomposer partitions this segment into distinct monomers, providing an accurate translation from the nucleotide alphabet into the monomer alphabet.


## Installation

The recommended way to install StringDecomposer is with conda package manager:
```
conda install -c bioconda stringdecomposer
```


Alternatively, StringDecomposer can be build and installed from source.

Requirements:
- Python3.5+
    - [biopython](https://biopython.org/wiki/Download)
    - [pandas](https://pypi.org/project/pandas/)
    - [python-edlib](https://pypi.org/project/edlib/)
    - [setuptools](https://pypi.org/project/setuptools/)
- g++ (version 5.3.1 or higher)

The required python packages can be installed through conda using 

    conda install --file requirements.txt

Local building without installation:

    git clone https://github.com/ablab/stringdecomposer.git
    cd stringdecomposer
    make

Then, StringDecomposer is available as

    python stringdecomposer.py


Installing from source:

    git clone https://github.com/ablab/stringdecomposer.git
    cd stringdecomposer
    make install

Then, StringDecomposer is available as

    stringdecomposer

Removal of StringDecomposer installed from source:

    make uninstall

## Quick start
The following command assumes that StringDecomposer is either installed through conda or from source.

    stringdecomposer ./test_data/read.fa ./test_data/DXZ1_star_monomers.fa -o ./test_data

The same result can be achieved with `make test_launch` (for local build without installation) and
`make test_launch_install` (for installed from source or via conda).
These `make` rules ensure correctness of StringDecomposer's output on the test dataset.

Testing run results:

    ./test_data/final_decomposition.tsv           final decomposition of sequences to monomer alphabet
    ./test_data/final_decomposition_alt.tsv       final decomposition of sequences to monomer alphabet with alternative monomers for each position
    ./test_data/final_decomposition_raw.tsv       raw decomposition with initial dynamic programming scores instead of identities

Each line in final_decomposition.tsv file has the following form:

    <read-name> <best-monomer> <start-pos> <end-pos> <identity> <second-best-monomer> <second-best-monomer-identity> <homo-best-monomer> <homo-identity> <homo-second-best-monomer> <homo-second-best-monomer-identity> <reliability>

`homo`-related columns represent statistics of the best-scoring (second-best-scoring) monomer after compression of homopolymer runs in both the monomer and the target read.
Reliability is either equal to `?` (signifies unreliable alignment which can be caused by a retrotransposon insertion or a poor quality segment of a read) or `+` (if the alignment is reliable).


## Synopsis

    stringdecomposer [-h] [-t THREADS] [-o OUT_FILE] [-i MIN_IDENTITY] [-s SCORING] [-b BATCH_SIZE] [--fast] sequences monomers

Required arguments:

    sequences                                         fasta-file with long reads or genomic sequences (accepts multiple sequences in one file)
    monomers                                          fasta-file with monomers

Optional arguments:

    -h, --help                                         show this help message and exit

    -t THREADS, --threads THREADS                      number of threads (by default 1)

    -o OUT_FILE, --out-file OUT_FILE                   output tsv-file (by default final_decomposition.tsv)

    -i MIN_IDENTITY, --min-identity MIN_IDENTITY       only monomer alignments with percent identity >= MIN_IDENTITY are printed (by default MIN_IDENTITY=0%)

    -s SCORING, --scoring SCORING                      set scoring scheme for StringDecomposer in the format "insertion,deletion,mismatch,match" (by default "-1,-1,-1,1")

    -b BATCH_SIZE, --batch-size BATCH_SIZE             set size of the batch in parallelization (by default 5000)

    --fast                                             StringDecomposer won't generate <second-best-monomer>, <second-best-monomer-identity>, <reliability> and _homo_-related columns (very useful in case of large number of monomers)

## Latest updates

### StringDecomposer 1.0 release (11 August 2020)

* initial StringDecomposer release
* conda support
* results of StringDecomposer monomer annotation for available centromere assemblies and ONT and Hifi reads of cen6, cen8, and cenX can be found at [Figshare](https://doi.org/10.6084/m9.figshare.12783371)


## Citation

The String Decomposition Problem and its Applications to Centromere Analysis and Assembly. *Tatiana Dvorkina, Andrey V. Bzikadze, Pavel A. Pevzner* Bioinformatics, Volume 36, Issue Supplement_1, July 2020, Pages i93–i101; doi: [https://doi.org/10.1093/bioinformatics/btaa454](https://doi.org/10.1093/bioinformatics/btaa454)

## Contact

In case of any issues please use [issue tracker](https://github.com/ablab/stringdecomposer/issues) or email directly to [t.dvorkina@spbu.ru](mailto:t.dvorkina@spbu.ru)
