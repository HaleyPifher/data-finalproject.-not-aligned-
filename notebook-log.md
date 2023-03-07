# *Description*
Data on 10 strains of non toxigenic C. difficile with a focus on strain CCUG37785 

-whole genomes 

-reads already assessed using FastQC

-Study relevant to data showed that CCUG37785 spores inoculated orally in mice provided almost full protection against a virulent strain. 

# *Alignment Methods* 
## 1)ClustalW
### Software Description 
ClustalW makes improvements to the progressive multiple alignment method which greatly improve the sensitivity without sacrificing any of the speed and efficiency. ClustalW uses the NJ algorithm by default, but UPGMA can be used for more efficient tree production. The default was used in generating the alignment for this project. 

Method: (i) all pairs of sequences are aligned separately in order to calculate a distance matrix giving the divergence of each pair of sequences; (ii) a guide tree is calculated from the distance matrix; (iii) the sequences are progressively aligned according to the branching order in the guide tree. 
### Strengths & Weaknesses 
In cases where all of the sequences in a data set are very similar (e.g. no pair less than 35% identical), CLUSTAL W will find an alignment which is difficult to improve by eye. In this sense, the alignment is optimal with regard to the alternative of manual
alignment.

Use of different weight matrices as the alignment progresses by-passes the problem of initial choice of weight matrix.

Long processing time  

As more divergent sequences are included, it becomes increasingly difficult to find good alignments and to evaluate them. What we find with CLUSTAL W is that the basic block-like structure of the alignment (corresponding to the major secondary structure elements) is usually recovered, with some of the most divergent sequences misaligned in small regions. This is a very useful starting point for manual refinement, as it helps define the major blocks of similarity. 
### User Choice
In deciding which sequences to align for sequences: CCUG37785, CD37, 5.3 they varied greatly in size(bp), so I found the average of each sequence and then used the following formula in excel to identify the closest size to the average: =INDEX(A2:A28,MATCH(MIN(ABS(D2:D28-D29)),ABS(D2:D28-D29),0)). This was not possible for every sequence because some files were simply too large and not easily transferrable. For sequences: Z31, DS28666, DSM29637, DSM29688, DSM28670, DSM29629, and DSM28669 I randomly selected a sequence to align. *All sequences represent whole genomes for the various strains.* All 10 sequences were then placed in a fasta file and added to the clustalw downloaded software folder. 

Refer to my ClustalW terminal output file for step by step of how the alignment was completed.
### References
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC308517/pdf/nar00046-0131.pdf

http://www.clustal.org/download/clustalw_help.txt

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9045287/#B21

## 2)Muscle
### Software Description 
Following guide tree construction, the fundamental step is pairwise profile alignment, which is used first for progressive alignment and then for refinement. MUSCLE uses two distance measures for a pair of sequences: a k mer distance (for an unaligned pair) and the Kimura distance (for an aligned pair). 
Figure that describes process visually can be found at https://academic.oup.com/view-large/figure/38266783/gkh340f2.jpeg 
### Strengths & Weaknesses
Fast distance estimation using k mers (espec for large sequences)
Compared to other methods (CLUSTALW, T-Coffee, and MAFFT) has average accuracy and high speed for large numbers of sequences 
UPGMA over Neighbor Joining (NJ) is used to build tree because it provides better accuracy at each node by aligning 2 profiles with the fewest differences even if they aren’t evolutionary neighbors.
Lack of statistically significant results in paper (linked in references) ex. Muscle-p (time + space complexity, no refinement included) results are not statistically significant compared to T-Coffee and NWNSI
Stage I of the process emphasizes speed over accuracy 
### User Choice 
Same user choices as described above under alignment method 1 (ClustalW). Additinally, I ran into an error/bug with muscle   "*** ERROR ***  MSA::GetLetter(0/1, 152/24551)=" this was fixed by running the folloing commands: 
*** ERROR ***  MSA::GetLetter(0/1, 2974/115365)='
(base) MacBook-Pro-2:muscle2 haleypifher$ 
(base) MacBook-Pro-2:muscle2 haleypifher$ ls
finalprojectalignmentdata.fasta	muscle3.8.31_i86darwin64	muscle3.8.31_i86darwin64.tar.gz
(base) MacBook-Pro-2:muscle2 haleypifher$ grep "'" finalprojectalignmentdata.fasta 
(base) MacBook-Pro-2:muscle2 haleypifher$ less finalprojectalignmentdata.fasta 
(base) MacBook-Pro-2:muscle2 haleypifher$ ./muscle3.8.31_i86darwin64 
After running Muscle again after fixing error another error appeared: *** ERROR ***  MSA::GetLetter(0/1, 1320/85188)='I'/4294967295 (3/2/23). 
### References 
https://academic.oup.com/nar/article/32/5/1792/2380623?login=false

https://www.drive5.com/muscle/manual/msa_getletter_bug.html

# *Distance and Parsimony Methods*
## 1) R Package - ape for distance-based tree estimation method
### Software Description 
ape is a package written in R for the analysis of phylogenetics and evolution. Functions for reading, writing, plotting, and manipulating phylogenetic trees, analyses of comparative data in a phylogenetic framework, ancestral character analyses, analyses of diversification and macroevolution, computing distances from DNA sequences, reading and writing nucleotide sequences as well as importing from BioConductor, and several tools such as Mantel's test, generalized skyline plots, graphical exploration of phylogenetic data (alex, trex, kronoviz), estimation of absolute evolutionary rates and clock-like trees using mean path lengths and penalized likelihood, dating trees with non-contemporaneous sequences, translating DNA into AA sequences, and assessing sequence alignments. Phylogeny estimation can be done with the NJ, BIONJ, ME, MVR, SDM, and triangle methods, and several methods handling incomplete distance matrices (NJ*, BIONJ*, MVR*, and the corresponding triangle method). 

I created a tree using the classical Neighbor-Joining (NJ) algorithm. 

### Strengths & Weaknesses 
can use the Tamura and Nei 1993 model which allows for different rates of transitions and transversions, heterogeneous base frequencies, and between-site variation of the substitution rate. 
phylogeny estimation can completed with a variety of methods 
R has the ability to match different datasets using labels (names, rownames and colnames), and ape uses this feature to match trees, sequences and other data. 
users can face practical difficulties that prevent efficient analysis 

### User Choice 
In creating tree ran into the following error message: "Error in njs(D) : distance information insufficient to construct a tree, cannot calculate agglomeration criterion." To correct it I refered to page 115 of package 'ape' documentation which has the following note: "If the sequences are very different, most evolutionary distances are undefined and a non-finite value (Inf or NaN) is returned. You may do dist.dna(, model = "raw") to check whether some values are higher than 0.75." None of my distance values were equal to or >.75, but changing the model from TN93 to raw did allow me to create a tree. 

### References 
https://cran.r-project.org/web/packages/ape/ape.pdf

http://ape-package.ird.fr/

https://academic.oup.com/bioinformatics/article/35/3/526/5055127

## 2) R package - phangron 
### Software Description 
phangron is a package for phylogenetic reconstruction and analysis in R. It now offers the possibility of reconstructing phylogenies with distance based methods, maximum parsimony or maximum likelihood (ML) and performing Hadamard conjugation. Extending the general ML framework, this package provides the possibility of estimating mixture and partition models. Furthermore, phangorn offers several functions for comparing trees, phylogenetic models or splits, simulating character data and performing congruence analyses.

I created a tree using the maximum parsimony (MP) method. MP is an optimality criterion for which the preferred tree is the tree that requires the least changes to explain some data. 

### Strengths & Weaknesses 
allows for search and comparison to find the better tree 

parsimony methods have been shown to produce inconsistent trees
### User Choice 
The maximum parsimony method was used to construct the tree.
### Refrences 
https://academic.oup.com/bioinformatics/article/27/4/592/198887
