source("function_library.r")# a library of all functions used for estimating the order of mutations
dyn.load("loglik-call2_32.so")
library(gtools)
library(alabama)
library(doMC) # for parallel computing
registerDoMC()
require(methods)
MC <- detectCores()
parallel <- F
N <- 10 # N is the number of times that the optimization with different inital values is repeated. If the variable "parallel" is TRUE, then the repeated optimization is done using parallel computing R package doMC. If it is not possible, let "parallel" be FALSE, then the optimization is repeated using for loop. We recommend N to be at least 10.

Eventtable=read.table("rectal_table",header=T)
groupgene=read.table("CoreGroup2Gene.txt",header = T)

pathway.database=pathwaylist.fun(groupgene)
driver.geneset=names(pathway.database)
mutated.genes=unique(Eventtable[,2])
genesets=redefine.geneset.fun2(pathway.database,mutated.genes)
result=mut.table.fun(Eventtable,redefined.geneset=genesets,T)  


####MLE#######
MAX <- 3 # assume the same distribution after MAX events. In this case, assume P_{k,i}=P_{5,i} for k>5
nozero.samp=which(rowSums(result$mut.table)>0)
mut.table <- result$mut.table[nozero.samp, ]# remove samples with no mutations
score.list=result$score.list[nozero.samp]
mle <- order_estimate(mut.table,N,parallel,score.list) # mle of P_{k,i}, the prob of kth mutation occurring in gene i
row.names(mle)=colnames(mut.table)
pair.order.prob=pair.order.prob.fun(driver.geneset,mle)
pair.order.prob[,3:5]=round(pair.order.prob[,3:5],3)
write.table(pair.order.prob,"pairprob.txt",row.names=F,col.names=T)


########a simulations for one pair of pathways##########
pathway.database=pathway.database[c(2,9)]
driver.geneset=names(pathway.database)
mutated.genes=unique(Eventtable[,2])
genesets=redefine.geneset.fun2(pathway.database,mutated.genes)
result=mut.table.fun(Eventtable,redefined.geneset=genesets,T)  
MAX <- 3 # assume the same distribution after MAX events. In this case, assume P_{k,i}=P_{5,i} for k>5
nozero.samp=which(rowSums(result$mut.table)>0)
mut.table <- result$mut.table[nozero.samp, ]# remove samples with no mutations
score.list=result$score.list[nozero.samp]
mle <- order_estimate(mut.table,N,parallel,score.list) # mle of P_{k,i}, the prob of kth mutation occurring in gene i
row.names(mle)=colnames(mut.table)

repeats <- 100
samplesize=50
event.score.list=generate.event.score(genesets,Eventtable)
prob.nondam=probnon(mut.table,score.list)
distance=simu.pairprob.fun(repeats=repeats,MAX=MAX,mut.table=mut.table,score.list=score.list,mle.result=mle.result,event.score.list=event.score.list,samplesize=samplesize[i],genesets=genesets,Eventtable=Eventtable,cutoff=0.1,max_unfix=10,prob.nondam=prob.nondam)


