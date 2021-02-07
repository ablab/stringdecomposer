
# CentromereArchitect

## Version 0.5

CentromereArchitect (CA) is the first tool for annotation of alpha satellite arrays in centromeres of a newly assembled human genome.
CA consists of two modules: 
- Monomer Inference allows extraction of human monomers based on the given alpha-satellite consensus template and centromeric sequence.
- HOR Inference allows extraction of HORs from the centromeric sequences using the inferred monomers.

This branch contains prerelease version of CA that was used to generate results for ISMB-2021 submission. 
We applied CA to the recently generated [complete assembly of a human genome](https://github.com/nanopore-wgs-consortium/CHM13) by the Telomere-to-Telomere consortium (v1.0) and generated the complete set of human monomers and high-order repeats for so-called live centromeres. 

All centromeric sequences, generated monomers and HORs can be found [here](https://figshare.com/s/4e7cd6e7cb3397c6ef6f).

In order to replicate our analysis, please download these data and follow the steps below:

## Install Required Libraries
- Python3.7
    - [biopython](https://biopython.org/wiki/Download)
    - [joblib](https://joblib.readthedocs.io/en/latest/installing.html)
    - [networkx](https://pypi.org/project/networkx)
    - [numpy](https://scipy.org/install.html)
    - [matplotlib](https://pypi.org/project/matplotlib/)
    - [pandas](https://pypi.org/project/pandas/)
    - [python-edlib](https://pypi.org/project/edlib/)
    - [seaborn](https://pypi.org/project/seaborn/)
    - [scikit-learn==0.21.3](https://pypi.org/project/scikit-learn/)
    - [setuptools](https://pypi.org/project/setuptools/)
- g++ (version 5.3.1 or higher)

The required python packages can be installed through conda using ```conda install --file requirements.txt```.

## Build StringDecomposer

    git clone https://github.com/ablab/stringdecomposer.git
    cd stringdecomposer
    make

## Run Monomer Inference

Example of running MonomerInference on cenX (cenXct.fa) with consensus Alpha Satellite (AlphaSat.fa):
```
python3 ca/monomer_inference.py -seq ISMB2021/Centromeres/cenXct.fa -mon ISMB2021/Monomers/AlphaSat.fa -o ISMB2021/Monomers/cenX
```
The result monomers can be found in ```ISMB2021/Monomers/cenX/monomers.fa```.
In the CentromereArchitect paper, MonomerInference step is launched on all centromeric sequences (Centromeres/allCt.fa) and this may be time consuming to replicate.

Example of running MonomerShift on 90 bp left for newly obtained monomers:
```
python3 ca/shift_monomers.py --finalDec ISMB2021/Monomers/cenX/final/final_decomposition.tsv -mn ISMB2021/Monomers/cenX/monomers.fa --shift 90 -o ISMB2021/Monomers/cenX/shift90
```
Shifted monomers can be found in ```ISMB2021/Monomers/cenX/shift90/shifted_mn.fasta```.

## Run HOR Inference

HORDecomposition algorithm applied to cenX using the set of monomers inferred from all centromeres:
```
python3 ca/extract_hors.py ISMB2021/Centromeres/cenXct.fa ISMB2021/Monomers/cenX_mn.fa ISMB2021/MonomerDecomposition/cenX_monodecomposition.tsv ISMB2021/HORDecomposition/cenX_hordecomposition.tsv
```
Resulting HORs can be found in ```ISMB2021/HORDecomposition/cenX_hordecomposition.tsv```.

SuperHORDecomposition algorithm applied to cenX using the set of monomers inferred from all centromeres:
```
python3 ca/extract_hors.py ISMB2021/Centromeres/cenXct.fa ISMB2021/Monomers/cenX_mn.fa ISMB2021/MonomerDecomposition/cenX_monodecomposition.tsv ISMB2021/SuperHORDecomposition/cenX_superhordecomposition.tsv --superhor
```

Resulting superHORs can be found in ```ISMB2021/SuperHORDecomposition/cenX_superhordecomposition.tsv```.


## Contact

In case of any issues please email directly to [t.dvorkina@spbu.ru](mailto:t.dvorkina@spbu.ru)
