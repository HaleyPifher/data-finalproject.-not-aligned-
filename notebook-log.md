# *Description*
Data on 10 strains of non toxigenic C. difficile with a focus on strain CCUG37785 
-whole genomes 
-reads already assessed using FastQC
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
# References
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC308517/pdf/nar00046-0131.pdf

http://www.clustal.org/download/clustalw_help.txt
