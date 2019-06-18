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

pathway.database=pathwaylist.fun(groupgene) #Generate the list of genes included for all pathways

#########pairwise estimation#########
for(i in 1:(length(pathway.database)-1)){
  for(j in (i+1):length(pathway.database)){
    pathway.database.pairwise=pathway.database[c(i,j)]
    driver.geneset=names(pathway.database.pairwise)
    mutated.genes=unique(Eventtable[,2])
    genesets=redefine.geneset.fun2(pathway.database.pairwise,mutated.genes)
    result=mut.table.fun(Eventtable,redefined.geneset=genesets,T)
    MAX <- 3 # assume the same distribution after MAX events. In this case, assume P_{k,i}=P_{4,i} for k>4
    nozero.samp=which(rowSums(result$mut.table)>0)
    mut.table <- result$mut.table[nozero.samp, ]# remove samples with no mutations
    score.list=result$score.list[nozero.samp]# remove samples with no mutations
    mle <- order_estimate(mut.table,N,parallel,score.list) # mle of P_{k,i}, the probability of kth mutation occurring in pathway/geneset i
    row.names(mle)=colnames(mut.table)
    pair.order.prob=pair.order.prob.fun(driver.geneset,mle) #extract the probabilities of A<B, A>B and A=B
    pair.order.prob[,3:5]=round(pair.order.prob[,3:5],3)
    filename=paste(i,j,"pairprob.txt",sep="")
    write.table(pair.order.prob,filename,row.names=F,col.names=T)
  }
}



