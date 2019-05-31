##########generate pathwaylist######
pathwaylist.fun=function(groupgene){
  pathway=unique(groupgene[,1])
  pathwaylist=list()
  for(i in 1:length(pathway)){
    pathwaylist[[i]]=groupgene[which(groupgene[,1]==pathway[i]),2]   
  }
  names(pathwaylist)=pathway
  return(pathwaylist)
}

###########redefine geneset######


redefine.geneset.fun2=function(driver.database,mutated.genes){
  driver.geneset=names(driver.database)
  no.dri=length(driver.geneset)
  #  if(all=F){
  if(no.dri==2){
    name=c(driver.geneset[1],driver.geneset[2],paste(driver.geneset[1],driver.geneset[2],sep=","))
    genes.per.genesets=list()
    genes.per.genesets[[3]]=intersect(driver.database[[1]],driver.database[[2]])
    genes.per.genesets[[1]]=setdiff(driver.database[[1]],genes.per.genesets[[3]])
    genes.per.genesets[[2]]=setdiff(driver.database[[2]],genes.per.genesets[[3]])
    names(genes.per.genesets)=name
  }else{
    indicate=sapply(1:no.dri,function(j) mutated.genes %in% driver.database[[j]] )
    indicate2=t(sapply(1:length(mutated.genes),function(i) indicate[i,]*(1:no.dri)))
    gene.score=sapply(1:length(mutated.genes),function(i) sum(2^indicate2[i,which(indicate2[i,]!=0)]))
    unique.score=sort(unique(gene.score))
    
    unique.score=unique.score[which(unique.score!=0)]
    unique.score.index=sapply(1:length(unique.score),function(i) which(gene.score==unique.score[i]))
    genes.per.genesets=sapply(1:length(unique.score),function(i) mutated.genes[unique.score.index[[i]]])
    pathway.per.genesets.index=t(sapply(1:length(unique.score),function(i) indicate[unique.score.index[[i]][1],]))
    pathway.per.genesets=sapply(1:length(unique.score),function(i) driver.geneset[which(pathway.per.genesets.index[i,]==1)])
    names(genes.per.genesets)=sapply(1:length(pathway.per.genesets), function(i) paste(pathway.per.genesets[[i]],collapse = ","))
  }
  #  }else{
  #    for(i in 1:2^no.dri){
  
  #    }
  #  }
  return(genes.per.genesets)
}






###generate mut table (and probability score) for redefined geneset###



mut.table.fun=function(Eventtable,redefined.geneset,scorelist_TF){
  if(is.numeric(Eventtable[,3])==F){
    Eventtable[,3]=as.numeric(as.character(Eventtable[,3]))
  }
  
  samples=unique(Eventtable[,1])
  
  
  #  samples.id=sample(1:length(samples),50)
  #  samples=samples[samples.id]
  
  mut.table=matrix(0,nrow=length(samples),ncol=length(redefined.geneset))
  score.list=vector(mode="list",length=length(samples))
  row.names(mut.table)=samples
  colnames(mut.table)=names(redefined.geneset)
  
  for(i in 1:length(samples)){
    sample.id=which(Eventtable[,1]==samples[i])
    sample.gene.mut=Eventtable[sample.id,2]
    if(scorelist_TF==T){
      score.gene.mut=Eventtable[sample.id,3]
    }
    for(j in 1:length(redefined.geneset)){
      for(k in 1:length(sample.gene.mut)){
        if(sample.gene.mut[k] %in% redefined.geneset[[j]]){
          mut.table[i,j]=mut.table[i,j]+1
          if(scorelist_TF==T){
            score.list[[i]]=c(score.list[[i]],score.gene.mut[k])}
        }
      }
    }
  score.list[[i]]=as.numeric(score.list[[i]])
  }
  result=list(mut.table=mut.table,score.list=score.list)
  return(result)
}

####Caculate the probility of order of two pathways###
pair.order.prob.fun=function(sort.geneset,mle.result){ 
  a=mle.result
  no.dri=length(sort.geneset)
  order.prob=matrix(NA,no.dri*(no.dri-1)/2,5)
  order.prob=as.data.frame(order.prob)
  colnames(order.prob)=c("A","B","P(A=B)","P(A<B)","P(A>B)")
  no.dri=length(sort.geneset)
  geneset.names=strsplit(row.names(mle.result),split=",")
#  geneset.names=row.names(mle.result)
  index=0
  

  
  for(i in 1:(no.dri-1)){
    for(j in (i+1):no.dri){
      index = index + 1
      A=sort.geneset[i]
      B=sort.geneset[j]
      logic.a=sapply(1:length(geneset.names),function(i) A %in% geneset.names[[i]])
      logic.b=sapply(1:length(geneset.names),function(i) B %in% geneset.names[[i]])
      pab=sapply(1:ncol(a),function(i) sum(a[,i]*logic.a*logic.b))
      pa=sapply(1:ncol(a),function(i) sum(a[,i]*logic.a*(!logic.b)))
      pb=sapply(1:ncol(a),function(i) sum(a[,i]*(!logic.a)*logic.b))
      pn=sapply(1:ncol(a),function(i) sum(a[,i]*(!logic.a)*(!logic.b)))
      pall=rbind(pab,pa,pb,pn)
      order.prob[index,3:5]=pall[1:3,1]
      
      for(l in 3:5){
        for(k in 2:ncol(a)){
          order.prob[index,l] =  order.prob[index,l]+pall[l-2,k]*prod(pall[4,1:(k-1)])
          
        }
        
      }
      order.prob[index,1]=A
      order.prob[index,2]=B
      order.prob[index,3:5]=order.prob[index,3:5]/sum(order.prob[index,3:5]) ###the prob of i=j,i>j,i<j

    }
  }
  return(order.prob)
}

###caculate the score of each pathway###
#order.score.fun=function(order.prob.pair,driver.geneset){
#  no.dri=length(driver.geneset)
#  order.score=rep(NA,no.dri)
#  names(order.score)=driver.geneset
#  for(i in 1:no.dri){
#    order.score[i]=sum(order.prob.pair[which(order.prob.pair[,1]==driver.geneset[i]),4])+sum(order.prob.pair[which(order.prob.pair[,2]==driver.geneset[i]),5])#-sum(order.prob.pair[which(order.prob.pair[,1]==driver.geneset[i]),5])-sum(order.prob.pair[which(order.prob.pair[,2]==driver.geneset[i]),4])
#  }
#  sort.geneset=names(sort(order.score,decreasing = TRUE))
#}

#prob4heatmap.fun=function(order.prob,sort.geneset) {
#  no.dri=length(sort.geneset)
#  order.prob.heatmap=matrix(NA,no.dri,no.dri)###prob for heatmap
#  colnames(order.prob.heatmap)=sort.geneset
#  rownames(order.prob.heatmap)=sort.geneset
#   for(i in 1:(no.dri-1)){
#    for(j in (i+1):no.dri){
#      index1=which((order.prob[,1]==sort.geneset[i]) & (order.prob[,2]==sort.geneset[j]))
#      if(length(index1)==1){
#        order.prob.heatmap[i,j]=order.prob[index1,4]
#        order.prob.heatmap[j,i]=order.prob[index1,5]
#      }else{
#        index2=which((order.prob[,2]==sort.geneset[i]) & (order.prob[,1]==sort.geneset[j]))
#        order.prob.heatmap[i,j]=order.prob[index2,5]
#        order.prob.heatmap[j,i]=order.prob[index2,4]
#      }
#    }
#  }
#  return(order.prob.heatmap)
#}


#######function for simulation######
# simulate.mut.table=function(mut.table,sample.size,mle.result){ #generate simulation table for a given sample size
#   no.mut.sample=apply(mut.table,1,sum)
#   no.sample=length(no.mut.sample)
#   no.genesets=nrow(prob)
#   simu.sample.index=sample(1:no.sample,sample.size,replace = T)
#   simu.no.mut=no.mut.sample[simu.sample.index]
#   mut.table=matrix(0,nrow=sample.size,ncol=no.genesets)
#   colnames(mut.table)=rownames(mle.result)
#   for(i in 1:sample.size){
#     for(j in 1:simu.no.mut[i]){
#       simu.genesets=sample(1:no.gensets,1,prob=mle.result[,j])
#       mut.table[i,simu.genesets] = mut.table[i,simu.genesets]+1
#     }
#   }
#   return(mut.table)
# }
# 
# simu.pairprob.fun=function(repeats,MAX,mle.result,no.mut.sample,sample.size,no.dri){
#   a=mle.result
#   simu.mle=array(NA,c(nrow(a),MAX+1,repeats))
#   simu.pair.prob=array(NA,c(no.dri*(no.dri-1)/2,3,repeats))
#   
#   for(i in 1:repeats){
#     simu.mut.table=simulate.mut.table(no.mut.sample,sample.size, a)
#     simu.mle[,,i]=order_estimate(simu.mut.table,N,parallel)[,1:5]
#     simu.pair.prob[,,i]=as.matrix(pair.order.prob.fun(sort.geneset,simu.mle[,,i]))
# #    print(i)
#   }
#    result=simu.mle
# #  simu.mle.mean=apply(simu.mle,c(1,2),mean)
# #  distance1=apply(abs(simu.mle.mean-a),2,mean)
# #  simu.pair.prob.mean=apply(simu.pair.prob,c(1,2),mean)
# #  pair.prob.true=as.matrix(pair.order.prob.fun(sort.geneset,a))
# #  distance2=abs(pair.prob.true-simu.pair.prob.mean)
# #  result=list(mle=simu.mle,distance1=distance1,distance2=distance2)
# #  return(result)
# }

##############id2name
#id2name.fun=function(id){
#  name=rep(NA,length(id))
#  for(i in 1:length(id)){
#    if(length(which(names(EntrezID2Name)==id[i]))==1){
#      index=which(names(EntrezID2Name)==id[i])
#      name[i]=EntrezID2Name[index]
#    }
#  }
#  return(name)
#}


####function for simulation for score####
generate.event.score=function(genesets,Eventtable) {
  event.score.list=list()
  index.list=sapply(1:length(genesets),function(i) which(Eventtable[,2]%in%genesets[[i]]))
  for(i in 1:length(genesets)){
    event.score.list[[i]]=as.numeric(Eventtable[index.list[[i]],3])
  }
  return(event.score.list)
}

simulate.mut.table=function(mut.table,score.list,mle.result,samplesize,event.score.list,prob.nondam){
  no.gene=ncol(mut.table)
  index=sample(1:nrow(mut.table),samplesize,replace = T)
  simu.mut.table=matrix(0,ncol=no.gene,nrow=samplesize)
  # print(samplesize)
  simu.score.list=vector("list",length=samplesize)
  for(i in 1:samplesize){
    no.mut=sum(mut.table[index[i],])
    score=score.list[[index[i]]]
    no.dmut=sum(sapply(1:no.mut,function(j) sample(c(1,0),1,prob=c(score[j],1-score[j]))))
    #  print("1")
    
    if(no.dmut>0){
      #print("if")
      #print(no.dmut)
      #print(no.gene)
      #print(mle.result)
      simu.dmut=sapply(1:no.dmut,function(j) sample(1:no.gene,1,prob=mle.result[,j]))
      #print("12")
      simu.dmut.count=sapply(1:no.gene,function(j) sum(simu.dmut==j))
      #print("13")
    }else{
      #print("else")
      simu.dmut.count=rep(0,no.gene)
      #print("14")
    }
    # print("2")
    if((no.mut-no.dmut)>0){
      simu.nmut=sample(1:no.gene,(no.mut-no.dmut),prob=prob.nondam,replace = T)
      simu.nmut.count=sapply(1:no.gene,function(j) sum(simu.nmut==j))
    }else{
      simu.nmut.count=rep(0,no.gene)
    }
    #print("3")
    simu.mut.table[i,]=simu.dmut.count+simu.nmut.count
    for(j in 1:no.gene){
      simu.score.list[[i]]=c(simu.score.list[[i]],sample(event.score.list[[j]],simu.dmut.count[j],replace = T,prob=event.score.list[[j]]/sum(event.score.list[[j]])))
      simu.score.list[[i]]=c(simu.score.list[[i]],sample(event.score.list[[j]],simu.nmut.count[j],replace = T,prob=(1-event.score.list[[j]])/sum(1-event.score.list[[j]])))
    }
    #print("4")
  }
  simu=list(simu.mut.table=simu.mut.table,simu.score.list=simu.score.list)
}

simu.pairprob.fun=function(repeats,MAX,mut.table,score.list,mle.result,event.score.list,samplesize,genesets,Eventtable,no.dri=2,cutoff,max_unfix,prob.nondam){
  a=mle.result
  simu.mle=array(NA,c(nrow(a),MAX+1,repeats))
  simu.pair.prob=array(NA,c(no.dri*(no.dri-1)/2,3,repeats))
  distance=matrix(NA,repeats,3)
  for(i in 1:repeats){
    simu=simulate.mut.table(mut.table,score.list,mle.result,samplesize,event.score.list,prob.nondam)
    # print(i)
    simu.mut.table=simu$simu.mut.table
    #  colnames(simu.mut.table)=colnames(mut.table)
    simu.score.list=simu$simu.score.list
    
    mle=order_estimate(simu.mut.table,N,parallel,simu.score.list,cutoff,max_unfix)
    if(ncol(mle)>MAX){
      simu.mle[,,i]=mle[,1:(MAX+1)]
    }else{
      #     mle=cbind(mle,mle[,ncol(mle)]
      simu.mle[,1:ncol(mle),i]=mle
    } 
    temp.mle=simu.mle[,,i]
    row.names(temp.mle)=colnames(mut.table)
    simu.pair.prob[,,i]=as.matrix(pair.order.prob.fun(driver.geneset,temp.mle)[,3:5])
    distance[i,]=abs(simu.mle[,1,i]-mle.result[,1])
  }
  mean.distance=apply(distance,2,mean)
  names(mean.distance)=c("P(A<B)","P(A>B)","P(A=B)")
  result=list(mle=simu.mle,distance=distance)
  return(result)#(mean.distance)
}


####function to change score in a specified pathway####
changescore=function(Eventtable,driver.database,changeset,change){
  for(i in 1:nrow(Eventtable)){
    if((Eventtable[i,2]) %in% driver.database[[changeset]]){
      Eventtable[i,3] = as.numeric(Eventtable[i,3])+change
      Eventtable[i,3] =max(Eventtable[i,3],0)
      Eventtable[i,3] =min(Eventtable[i,3],1)
    }
  }
  return(Eventtable)
}


########function for individual prediction#######

#AllPesti is the estimated probability matrix of the ith mutation occurs in A,B or A intersect B, Yj is the observated mutation in the jth individual, prob is the function impact score for the individual.#
predict <- function(Yj,AllPesti,a,prob,MAX)  
{
  storage.mode(AllPesti) <- "double"
  storage.mode(Yj) <- "integer"
  storage.mode(prob) <- "double"
  for (i in 1:length(a))
    storage.mode(a[[i]]) <- "integer"
  storage.mode(MAX) <- "integer"
  .Call("predict", Yj,AllPesti,a,prob,MAX)
}
