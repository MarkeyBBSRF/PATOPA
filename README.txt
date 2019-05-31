The file "example.R” is an example of R codes which read the data file “rectal_table” and “CoreGroup2Gene.txt” and returns “pairprob.txt”.

“rectal_table”: All the nonsynonymous mutations in TCGA rectal cancer dataset and its coresponding Polyphen-2 scores.

“CoreGroup2Gene.txt”: contains the mapping of pathways to genes (genes entrez id).


“pairprob.txt" : the estimates of P(A>B), P(A<B),P(A=B) for all pairs of pathway A and pathway B in the CoreGroup2Gene.txt file. 


The file "function_library.r" contains R codes used in estimating the order of mutations and simulations. It calls C code "loglik-call2_32.c”. One needs to compile "loglik-call2_32.c” using the following command before calls it:

R CMD SHLIB loglik-call2_32.c


