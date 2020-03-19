# String decomposition into monomers

StringDecomposer (SD) algorithm takes the set of monomers and a long error-prone read (or a genomic segment) and partitions this read into distinct monomers, providing an accurate translation of each read from a nucleotide alphabet into a monomer alphabet.


## Installation

    git clone https://github.com/TanyaDvorkina/stringdecomposer.git
    cd src && g++ -o dp main.cpp edlib/src/edlib.cpp -I edlib/include -fopenmp && cd ..

Requirements:
- Python3
    - [biopython](https://biopython.org/wiki/Download)
    - [edlib](https://pypi.org/project/edlib/)
    - [argparse](https://pypi.org/project/argparse/)
- g++ (version 5.3.1 or higher)

Requirements can be installed through Conda as ```conda install --file requirements.txt```.

## Quick start

    run_decomposer.py ./test_data/read.fa ./test_data/DXZ1_star_monomers.fa -r

Testing run results:

    final_decomposition.tsv           final decomposition of sequences to monomer alphabet
    final_decomposition_alt.tsv       final decomposition of sequences to monomer alphabet with alternative monomers for each position
    raw_decomposition.tsv             raw decomposition with initial dynamic programming scores instead of identities

Each line in final_decomposition.tsv file has the following form:

    <read-name> <best-monomer> <start-pos> <end-pos> <identity> <second-best-monomer> <second-best-monomer-identity> <homo-best-monomer> <homo-identity> <homo-second-best-monomer> <homo-second-best-monomer-identity> <reliability>

_homo_-related columns represent statistics of the best-scoring (second-best-scoring) monomer after homopolymer collapsing in both monomer and the target read.


## Synopsis

    run_decomposer.py [-h] [-t THREADS] [-o OUT_FILE] [-i MIN_IDENTITY] [-s SCORING] [-r] sequences monomers

Required arguments:

    sequences                                         fasta-file with long reads or genomic sequences (accepts multiple sequences in one file)
    monomers                                          fasta-file with monomers

Optional arguments:

    -h, --help                                         show this help message and exit

    -t THREADS, --threads THREADS                      number of threads (by default 1)

    -o OUT_FILE, --out-file OUT_FILE                   output tsv-file (by default final_decomposition.tsv)

    -i MIN_IDENTITY, --min-identity MIN_IDENTITY       only monomer alignments with percent identity >= MIN_IDENTITY are printed (by default MIN_IDENTITY=0%)

    -s SCORING, --scoring SCORING                      set scoring scheme for SD in the format "insertion,deletion,mismatch,match" (by default "-1,-1,-1,1")

    -r, --raw                                          save initial monomer decomposition to [OUTPUT_FILE_FOLDER]/raw_decomposition.tsv (by default False)


## Citation

The String Decomposition Problem and its Applications to Centromere Assembly. *Tatiana Dvorkina, Andrey V. Bzikadze, Pavel A. Pevzner* bioRxiv 2019.12.26.888685; doi: [https://doi.org/10.1101/2019.12.26.888685](https://doi.org/10.1101/2019.12.26.888685)

## Contact

In case of any issues please use [issue tracker](https://github.com/ablab/stringdecomposer/issues) or email directly to [t.dvorkina@spbu.ru](mailto:t.dvorkina@spbu.ru)
