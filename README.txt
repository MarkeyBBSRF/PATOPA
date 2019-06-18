1.Input Data Sets:

“rectal_table”: All the nonsynonymous mutations in TCGA rectal cancer dataset and its coresponding Polyphen-2 scores. The first column is the sample names, the second column is the mutated genes and the third column is the functional impact score corresponding to each mutations.

“CoreGroup2Gene.txt”: contains the mapping of pathways to genes (genes entrez id). The first column is the name of pathways, the second column is genes included in the pathways.

2.Example files:

The file "example1.R” is an example of R codes which read the data file “rectal_table” and “CoreGroup2Gene.txt” and analyze the 9 pathways in “CoreGroup2Gene.txt” together. The file returns “pairprob.txt” which contains the probabilities of A>B, A<B and A=B for all pairs of pathways.

The file "example2.R” is an example of R codes which read the data file “rectal_table” and “CoreGroup2Gene.txt” and analyze the 9 pathways in “CoreGroup2Gene.txt” pairwisely. The file returns “ijpairprob.txt” which contains the probabilities of P(Pathwayi>Pathwayj), P(Pathwayi>Pathwayj) and P(Pathwayi=Pathwayj) for all pairs of pathways.

The file "example3.R” is an example of R codes which read the data file “rectal_table” and “CoreGroup2Gene.txt” and do 100 simulation of size 50 based on the estimated probabilities  of the mutational order of MAPK signaling pathway and PI3K-Akt signaling pathway.The file returns the average distance of the estimated P(A>B), P(A<B) and P(A=B) from the true values.

3.Function library:

The file "function_library.r" contains R codes used in estimating the order of mutations and simulations. It calls C code "loglik-call2_32.c”. One needs to compile "loglik-call2_32.c” using the following command before calls it:

R CMD SHLIB loglik-call2_32.c


