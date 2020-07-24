# StringDecomposer

StringDecomposer (SD) algorithm takes the set of monomers and a long error-prone read (or a genomic segment) and partitions this read into distinct monomers, providing an accurate translation of each read from a nucleotide alphabet into a monomer alphabet.


## Installation

Requirements:
- Python3.5
    - [biopython](https://biopython.org/wiki/Download)
    - [argparse](https://pypi.org/project/argparse/)
    - [joblib](https://joblib.readthedocs.io/en/latest/installing.html)
    - [networkx](https://pypi.org/project/networkx)
    - [numpy](https://scipy.org/install.html)
    - [matplotlib](https://pypi.org/project/matplotlib/)
    - [pandas](https://pypi.org/project/pandas/)
    - [python-edlib](https://pypi.org/project/edlib/)
    - [seaborn](https://pypi.org/project/seaborn/)
    - [setuptools](https://pypi.org/project/setuptools/)
- g++ (version 5.3.1 or higher)

Requirements can be installed through Conda as ```conda install --file requirements.txt```.

### Local Building (without installation)

    git clone https://github.com/TanyaDvorkina/stringdecomposer.git
    cd stringdecomposer
    make

If you wish to run tests:

    make test

### Installing from source

    git clone https://github.com/TanyaDvorkina/stringdecomposer.git
    cd stringdecomposer
    python setup.py install

## Quick start

    sd/run_decomposer.py ./test_data/read.fa ./test_data/DXZ1_star_monomers.fa

Testing run results:

    final_decomposition.tsv           final decomposition of sequences to monomer alphabet
    final_decomposition_alt.tsv       final decomposition of sequences to monomer alphabet with alternative monomers for each position
    final_decomposition_raw.tsv       raw decomposition with initial dynamic programming scores instead of identities

Each line in final_decomposition.tsv file has the following form:

    <read-name> <best-monomer> <start-pos> <end-pos> <identity> <second-best-monomer> <second-best-monomer-identity> <homo-best-monomer> <homo-identity> <homo-second-best-monomer> <homo-second-best-monomer-identity> <reliability>

_homo_-related columns represent statistics of the best-scoring (second-best-scoring) monomer after homopolymer collapsing in both monomer and the target read. Reliability is either equal to ? (if StringDecomposer suggests that it is a potential gap) or + (if the alignment is reliable).


## Synopsis

    run_decomposer.py [-h] [-t THREADS] [-o OUT_FILE] [-i MIN_IDENTITY] [-s SCORING] [-b BATCH_SIZE] [--fast] sequences monomers

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



## `sd/scripts` folder

The set of scripts help to solve decomposition-related problems.

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
