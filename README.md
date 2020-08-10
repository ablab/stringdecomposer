[![Anaconda-Server Badge](https://anaconda.org/bioconda/stringdecomposer/badges/installer/conda.svg)](https://anaconda.org/bioconda/stringdecomposer)

# StringDecomposer

## Version 1.0

As an input StringDecomposer (SD) algorithm takes the set of monomers (typically, alpha satellites) and a genomic segment (assembly, Oxford Nanopore or a PacBio HiFi read) that contains a tandem repeat constituted by the given monomers.
StringDecomposer partitions this segment into distinct monomers, providing an accurate translation from the nucleotide alphabet into the monomer alphabet.


## Installation

The recommended way to install StringDecomposer is with conda package manager:
```
conda install -c bioconda stringdecomposer
```


Alternatively, StringDecomposer can be installed from source.

Requirements:
- Linux (OSX is currently not supported)
- Python 3.6 or 3.7
- C++ compliler with C++11 support (g++ version 5.3.1 or higher)
- GNU make

The required python packages can be installed through conda using ```conda install --file requirements.txt```.

Local building without installation:

    git clone https://github.com/ablab/stringdecomposer.git
    cd stringdecomposer
    make

Installing from source:

    git clone https://github.com/ablab/stringdecomposer.git
    cd stringdecomposer
    python setup.py install --record files.txt

Removal of StringDecomposer installed from source:

    xargs rm -rf < files.txt

## Quick start
The following command assumes that StringDecomposer is either installed through conda or from source.

    run_decomposer ./test_data/read.fa ./test_data/DXZ1_star_monomers.fa

Testing run results:

    final_decomposition.tsv           final decomposition of sequences to monomer alphabet
    final_decomposition_alt.tsv       final decomposition of sequences to monomer alphabet with alternative monomers for each position
    final_decomposition_raw.tsv       raw decomposition with initial dynamic programming scores instead of identities

Each line in final_decomposition.tsv file has the following form:

    <read-name> <best-monomer> <start-pos> <end-pos> <identity> <second-best-monomer> <second-best-monomer-identity> <homo-best-monomer> <homo-identity> <homo-second-best-monomer> <homo-second-best-monomer-identity> <reliability>

`homo`-related columns represent statistics of the best-scoring (second-best-scoring) monomer after compression of homopolymer runs in both the monomer and the target read.
Reliability is either equal to `?` (signifies unreliable alignment which can be caused by a retrotransposon insertion or a poor quality segment of a read) or `+` (if the alignment is reliable).


## Synopsis

    run_decomposer [-h] [-t THREADS] [-o OUT_FILE] [-i MIN_IDENTITY] [-s SCORING] [-b BATCH_SIZE] [--fast] sequences monomers

Required arguments:

    sequences                                         fasta-file with long reads or genomic sequences (accepts multiple sequences in one file)
    monomers                                          fasta-file with monomers

Optional arguments:

    -h, --help                                         show this help message and exit

    -t THREADS, --threads THREADS                      number of threads (by default 1)

    -o OUT_FILE, --out-file OUT_FILE                   output tsv-file (by default final_decomposition.tsv)

    -i MIN_IDENTITY, --min-identity MIN_IDENTITY       only monomer alignments with percent identity >= MIN_IDENTITY are printed (by default MIN_IDENTITY=0%)

    -s SCORING, --scoring SCORING                      set scoring scheme for SD in the format "insertion,deletion,mismatch,match" (by default "-1,-1,-1,1")

    -b BATCH_SIZE, --batch-size BATCH_SIZE             set size of the batch in parallelization (by default 5000)

    --fast                                             SD won't generate <second-best-monomer>, <second-best-monomer-identity>, <reliability> and _homo_-related columns (very useful in case of large number of monomers)

## Latest updates

### StringDecomposer 1.0 release (11 August 2020)

* initial StringDecomposer release
* conda support
* results of StringDecomposer monomer annotation for available centromere assemblies and ONT and Hifi reads of cen6, cen8, and cenX can be found at [Figshare](https://doi.org/10.6084/m9.figshare.12783371)


## StringDecomposer+

The set of scripts help to solve decomposition-related problems are placed in `sd/scripts` directory.

Script `extract_hors.py` converts monomer decomposition to HOR decomposition.
Synopsis:

    extract_hors.py <sequences> <monomers>  <decomposition> <output file> [--canonical <txt-file>] [--min-idnt MIN_IDNT] [--min-reliable MIN_RELIABLE] [--min-cnt MIN_CNT] [--min-weight MIN_WEIGHT] [--min-len MIN_LEN] [--max-len MAX_LEN]

File in `--canonical` option represents a list of canonical HORs. Each line represents one HOR, HOR has to be represented as a list of monomer ids from `<monomers>` file, separated by ",".

Script `extract_centromere_related_regions.py` extracts reads (or assembly segment) that are covered by the given set monomers.
Synopsis:

    extract_centromere_related_regions.py  -s <sequences> -m <monomers> -o <output file> [-d <edit-distance>]


Script `convert_identities.py` converts decomposition with arbitrary scores to final_decomposition.tsv (with identites instead of scores).
Synopsis:

    convert_identities.py -s <sequences> -m <monomers> -d <decomposition> [-o <output file>]

## Citation

The String Decomposition Problem and its Applications to Centromere Analysis and Assembly. *Tatiana Dvorkina, Andrey V. Bzikadze, Pavel A. Pevzner* Bioinformatics, Volume 36, Issue Supplement_1, July 2020, Pages i93–i101; doi: [https://doi.org/10.1093/bioinformatics/btaa454](https://doi.org/10.1093/bioinformatics/btaa454)

## Contact

In case of any issues please use [issue tracker](https://github.com/ablab/stringdecomposer/issues) or email directly to [t.dvorkina@spbu.ru](mailto:t.dvorkina@spbu.ru)
