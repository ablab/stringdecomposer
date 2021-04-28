#Monomer Graphs
The main script for drawing monomer graph is `BuildDualGraph.py`

Running example
```
python3 BuildDualGraph.py -sdtsv sdres -mon mon -monIA monIA -seq seq  -o graphMS
```

###Input:
* sdres -- folder with results of StringDecomposer for all centromers. The name of SD output for centromere N 
  should have a name "cenN_dec.tsv"(sdres/cen1_dec.tsv, sdres/cen2_dec.tsv...)
    
* mon -- folder with monomers used for StringDecomposer. Files have names cenN_mn.fa(mon/cen1_mn.fa, mon/cen2_mn.fa...)

* monIA -- folder with Ivan's monomoers. Should be present for getting Ivan's names. Files have names 
  cenN_mn.fa(monIA/cen1_mn.fa, monIA/cen2_mn.fa, ...)
  
* seq -- folder with centromere sequence. We only need it if we change the monomer set but resolving Elerov Circulatin. 
  If you don't print monoruns graphs it is not necessary. Files names: "seq/cen1ct.fa", "seq/cen2.ct.fa" etc
  
* graphMS -- output folder for monomers graphs in dot and png formats and also for monoruns graphs.

Files name sre important

###Extra args: 
* --red for showing red edges in graph which represent monomers with small edit distance
* --blue for showing edges from perfect matching
* --norm for normalizing weight on the #the most frequent monomer
* --maxk show De Bruijn Graph for k=1..maxk (default=1). For k=1 is just simple monomer graph
* --monorun also build and show monorun graph

###Other output files
* HORscnt.tsv -- file in directory from which script are running. Contain information about HORs and multiplicity. 
  Each run append to the end of file.
  
* L.csv -- information about vertex in Monoruns graphs. Also root directory and append to the end.