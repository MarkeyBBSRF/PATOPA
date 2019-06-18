dyn.load("loglik-call2_32.so")
source("function_library.r")
library(gtools)
library(doMC) # for parallel computing
registerDoMC()
require(methods)
MC <- detectCores()
parallel <- F
N <- 10 # N is the number of times that the optimization with different inital values is repeated. If the variable "parallel" is TRUE, then the repeated optimization is done using parallel computing R package doMC. If it is not possible, let "parallel" be FALSE, then the optimization is repeated using for loop. We recommend N to be at least 10.

Eventtable=read.table("rectal_table",header=T)
groupgene=read.table("CoreGroup2Gene.txt",header = T)

pathway.database=pathwaylist.fun(groupgene) #Generate the list of genes included for all pathways.
driver.geneset=names(pathway.database) #The vectors of names of all pathways analyzed
mutated.genes=unique(Eventtable[,2]) #All mutated genes
genesets=redefine.geneset.fun2(pathway.database,mutated.genes) #Generate the list of redefined genesets based on sets relationships
result=mut.table.fun(Eventtable,redefined.geneset=genesets,T) #Generate mutation table and functional impact score list


########estimation for one pair of pathways##########
pathway.database=pathway.database[c(2,9)]
driver.geneset=names(pathway.database)
mutated.genes=unique(Eventtable[,2])
genesets=redefine.geneset.fun2(pathway.database,mutated.genes)
result=mut.table.fun(Eventtable,redefined.geneset=genesets,T)
MAX <- 3 # assume the same distribution after MAX events. In this case, assume P_{k,i}=P_{4,i} for k>4
nozero.samp=which(rowSums(result$mut.table)>0)
mut.table <- result$mut.table[nozero.samp, ]# remove samples with no mutations
score.list=result$score.list[nozero.samp]# remove samples with no mutations
mle <- order_estimate(mut.table,N,parallel,score.list) # mle of P_{k,i}, the prob of kth mutation occurring in gene i
row.names(mle)=colnames(mut.table)

########simulations for the pair of pathways based on above estimation##########
repeats <- 100
samplesize=50
event.score.list=generate.event.score(genesets,Eventtable) #generate the list of functional impact scores for all mutations in each redefined mutually exclusive genesets, the list will be used in the simulation step
prob.nondam=probnon(mut.table,score.list) #the probabilities of a random non-functional mutations occurs in each pathways/datasets.
distance=simu.pairprob.fun(repeats=repeats,MAX=MAX,mut.table=mut.table,score.list=score.list,mle.result=mle,event.score.list=event.score.list,samplesize=samplesize[i],genesets=genesets,Eventtable=Eventtable,cutoff=0.1,max_unfix=10,prob.nondam=prob.nondam) #the mean distances of probabilities of estimated values of P(A<B), P(A>B) and P(A=B) from true values


