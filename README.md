# String decomposition into monomers

StringDecomposer (SD) algorithm takes the set of monomers and a long error-prone read (or a genomic segment) and partitions this read into distinct monomers, providing an accurate translation of each read from a nucleotide alphabet into a monomer alphabet.


## Installation

    git clone https://github.com/TanyaDvorkina/stringdecomposer.git
    cd src && g++ -o dp main.cpp -fopenmp && cd ..

Requirements:
- Python3
    - [biopython](https://biopython.org/wiki/Download)
    - [edlib](https://pypi.org/project/edlib/)
    - [argparse](https://pypi.org/project/argparse/)
- GCC
    - OpenMP

## Running

Synopsis: 

    run_decomposer.py [-h] [-t THREADS] [-o OUT_FILE] [-i MIN_IDENTITY] [-s SCORING] [-r] sequences monomers

Positional arguments:

  `sequences`

fasta-file with long reads or genomic sequences

  `monomers` 

fasta-file with monomers

Optional arguments:

  `-h, --help`            

show this help message and exit

  `-t THREADS, --threads THREADS`

number of threads (by default 1)

  `-o OUT_FILE, --out-file OUT_FILE`

output tsv-file (by default final_decomposition.tsv)

  `-i MIN_IDENTITY, --min-identity MIN_IDENTITY`

only monomer alignments with percent identity >= MIN_IDENTITY are printed (by default MIN_IDENTITY=0%)

  `-s SCORING, --scoring SCORING`

set scoring scheme for SD in the format "insertion,deletion,match,mismatch" (by default "-1,-1,-1,1")

  `-r, --raw`             

save initial monomer decomposition to raw_decomposition.tsv (by default False)


## Output

    final_decomposition.tsv           final sequences to monomer alphabet decomposition
    final_decomposition_alt.tsv       final sequences to monomer alphabet decomposition with alternative monomers for each position
    raw_decomposition.tsv             raw decomposition with initial dynamic programming scores instead of identities


Each line in final_decomposition.tsv file has the following form:

    <read-name> <best-monomer> <start-pos> <end-pos> <identity-score> <second-best-monomer/None> <second-best-monomer-identity/-1>