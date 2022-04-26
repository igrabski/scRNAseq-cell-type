library(sads)
library(BiocParallel)

# Get gene-specific distributions
load('params_EM_81020.rda')
pi.all.g <- params[[1]]
g.on.all.g <- params[[2]]
g.off.all.g <- params[[3]]
a.all.g <- params[[4]]
sigma.all.g <- params[[5]]
mu.g<-c(-12.73153,-8.301349)
gaps <- mu.g[2]+g.on.all.g-mu.g[1]-g.off.all.g
discrim.g<-rownames(pi.all.g)[which(gaps>=1&rowSums(pi.all.g)>=0.05&rowSums(pi.all.g)<=0.95)]

default <- registered()
register(MulticoreParam(workers = 6), default = TRUE)

# Helper function (from old snipEM package)
sumlog=function(x,lower=-745,upper=709) {
 n=length(x)
 # log(a + b) = log(a) + log (1 + exp(log(b) - log(a)))

 index=rep(1,n)
 i=1
 while(i<n & (x[i]<lower | x[i]>upper)) {
   i=i+1
 }
 if(i>n) {stop("Not Summable!")}
 if(i<=n) {
   s=x[i]
   index[i]=0
 }

 while(sum(index)>0) {
  idx=which(index==1)
  l=length(idx)
  i = 1 
  while(i<l & ((x[idx[i]]-s)<lower | (x[idx[i]]-s)>upper)) {
    i=i+1
  }
  if(i>l) {stop("Not Summable!")}
  if(i<=l) {
    s=s+log (1 + exp(x[idx[i]]-s))
    index[idx[i]]=0
  }
 }
 
 s
}


# Train reference data
trainReference <- function(ref,pi.all=pi.all.g,mu=mu.g,a=a.all.g,sigma.all=sigma.all.g,
                           g.on.all=g.on.all.g,g.off.all=g.off.all.g,discrim_only=TRUE) {
  # Subset to genes present in both reference data and knowledge base
  common <- intersect(names(ref),rownames(pi.all))
  if (discrim_only) {
    common <- intersect(common,discrim.g)
  }
  N <- sum(ref)
  ref <- ref[common]
  pi.all2 <- pi.all[common,]
  a.all2 <- a[common]
  sigma.all2 <- sigma.all[common,]
  g.on.all2 <- g.on.all[common]
  g.off.all2 <- g.off.all[common]
  
  # Compute mixing probability for each component
  prob.exp <- pi.all2[,1]*dnbinom(ref,1,1/((N/(a.all2))+1))
  prob.ln1 <- sapply(common,function(j) 
    pi.all2[j,2]*sads::dpoilog(ref[j],mu[1]+g.off.all2[j]+log(N),sigma.all2[j,1]))
  prob.ln2 <- sapply(common,function(j) 
    (1-pi.all2[j,1]-pi.all2[j,2])*sads::dpoilog(ref[j],mu[2]+g.on.all2[j]+log(N),sigma.all2[j,2])) 
  return(data.frame(rate=ref/N,exp=prob.exp/(prob.exp+prob.ln1+prob.ln2),ln1=prob.ln1/(prob.exp+prob.ln1+prob.ln2),
                    ln2=prob.ln2/(prob.exp+prob.ln1+prob.ln2)))
}

# Train all reference
trainAllReference <- function(data,labels,pi.all=pi.all.g,mu=mu.g,a=a.all.g,sigma.all=sigma.all.g,
                              g.on.all=g.on.all.g,g.off.all=g.off.all.g,discrim_only=TRUE) {
  cell.types <- unique(labels)
  dataSums <- lapply(cell.types,function(c) rowSums(data[,labels==c]))
  d.list <- bplapply(dataSums,trainReference,pi.all=pi.all,mu=mu,a=a,sigma.all=sigma.all,
                     g.on.all=g.on.all,g.off.all=g.off.all,discrim_only=discrim_only)
  names(d.list) <- cell.types
  return(d.list)
}

# Classify target cells (fast)
classifyTarget <- function(target,d.list,other=FALSE,pi.all=pi.all.g,mu=mu.g,a=a.all.g,sigma.all=sigma.all.g,
                           g.off.all=g.off.all.g,g.on.all=g.on.all.g,discrim=discrim.g,return.probs=FALSE) {
  genes.s <- intersect(rownames(target),discrim)
  genes.s <- intersect(genes.s,rownames(d.list[[1]]))
  n <- colSums(target)
  target <- target[genes.s,]
  pi.all2 <- pi.all[genes.s,]
  a2 <- a[genes.s]
  sigma.all2 <- sigma.all[genes.s,]
  g.on.all2 <- g.on.all[genes.s]
  g.off.all2 <- g.off.all[genes.s]
  if (other) {
    d.list[['other']] <- data.frame(rate=rep(0,nrow(d.list[[1]])),exp=pi.all.g[rownames(d.list[[1]]),1],
                                     ln1=pi.all.g[rownames(d.list[[1]]),2],
                                     ln2=1-pi.all.g[rownames(d.list[[1]]),1]-pi.all.g[rownames(d.list[[1]]),2])
  }
  d.list2 <- lapply(1:length(d.list),function(j) log(d.list[[j]][genes.s,2:4]))
  names(d.list2) <- names(d.list)
  
  # Store probabilities
  gene.components <- array(0,dim=c(length(genes.s),ncol(target),3))
  dimnames(gene.components)[[1]] <- genes.s
  probs <- array(0,dim=c(ncol(target),length(d.list2)))
  
  mu.matrix1 <- array(1/a2,dim=c(length(genes.s),1))%*%array(n,dim=c(1,length(n)))
  gene.components[,,1] <- dnbinom(target,1,1/(mu.matrix1+1),log=T)
  mu.matrix2a <- array(exp(2*(mu[1]+g.off.all2+sigma.all2[,1]^2/2))/
                         ((exp(sigma.all2[,1]^2)-1)*(exp(2*mu[1]+2*g.off.all2+sigma.all2[,1]^2))),
                       dim=c(length(genes.s),1))%*%array(1,dim=c(1,length(n)))
  mu.matrix2b <- 1/(1+(1/array(exp((mu[1]+g.off.all2+sigma.all2[,1]^2/2))/
                         ((exp(sigma.all2[,1]^2)-1)*(exp(2*mu[1]+2*g.off.all2+sigma.all2[,1]^2))),
                       dim=c(length(genes.s),1))%*%array(1/n,dim=c(1,length(n)))))
  gene.components[,,2] <- dnbinom(target,mu.matrix2a,mu.matrix2b,log=T)
  mu.matrix3a <- array(exp(2*(mu[2]+g.on.all2+sigma.all2[,2]^2/2))/
                         ((exp(sigma.all2[,2]^2)-1)*(exp(2*mu[2]+2*g.on.all2+sigma.all2[,2]^2))),
                       dim=c(length(genes.s),1))%*%array(1,dim=c(1,length(n)))
  mu.matrix3b <- 1/(1+(1/array(exp((mu[2]+g.on.all2+sigma.all2[,2]^2/2))/
                                 ((exp(sigma.all2[,2]^2)-1)*(exp(2*mu[2]+2*g.on.all2+sigma.all2[,2]^2))),
                               dim=c(length(genes.s),1))%*%array(1/n,dim=c(1,length(n)))))
  gene.components[,,3] <- dnbinom(target,mu.matrix3a,mu.matrix3b,log=T)

  # For each cell-type...
  for (t in 1:length(d.list2)) {
    prob <- exp(sweep(t(gene.components[,,1]),2,d.list2[[t]][,1],'+'))+
                               exp(sweep(t(gene.components[,,2]),2,d.list2[[t]][,2],'+'))+
                               exp(sweep(t(gene.components[,,3]),2,d.list2[[t]][,3],'+'))
    flag <- which(prob==0,arr.ind=T)
    prob <- log(prob)
    for (f in seq_len(nrow(flag))) {
      prob[flag[f,1],flag[f,2]] <- sumlog(gene.components[flag[f,2],flag[f,1],]+as.numeric(d.list2[[t]][flag[f,2],]))
    }
    probs[,t] <- rowSums(prob)
  }

  if (return.probs) {
    colnames(probs) <- names(d.list2)
    probs <- exp(probs-sapply(1:nrow(probs),function(x) max(probs[x,])))
    return(sweep(probs,1,rowSums(probs),'/'))
  }
  
  return(names(d.list2)[sapply(1:ncol(target),function(x) which.max(probs[x,]))])
}

# Return likelihoods
likelihoodTarget <- function(target,d.list,other=FALSE,pi.all=pi.all.g,mu=mu.g,a=a.all.g,sigma.all=sigma.all.g,
                           g.off.all=g.off.all.g,g.on.all=g.on.all.g,discrim=discrim.g) {
  genes.s <- intersect(rownames(target),discrim)
  genes.s <- intersect(genes.s,rownames(d.list[[1]]))
  n <- colSums(target)
  target <- target[genes.s,]
  pi.all2 <- pi.all[genes.s,]
  a2 <- a[genes.s]
  sigma.all2 <- sigma.all[genes.s,]
  g.on.all2 <- g.on.all[genes.s]
  g.off.all2 <- g.off.all[genes.s]
  if (other) {
    d.list[['other']] <- data.frame(rate=rep(0,nrow(d.list[[1]])),exp=pi.all.g[rownames(d.list[[1]]),1],
                                    ln1=pi.all.g[rownames(d.list[[1]]),2],
                                    ln2=1-pi.all.g[rownames(d.list[[1]]),1]-pi.all.g[rownames(d.list[[1]]),2])
  }
  d.list2 <- lapply(1:length(d.list),function(j) log(d.list[[j]][genes.s,2:4]))
  names(d.list2) <- names(d.list)
  
  
  # Store probabilities
  gene.components <- array(0,dim=c(length(genes.s),ncol(target),3))
  dimnames(gene.components)[[1]] <- genes.s
  probs <- array(0,dim=c(ncol(target),length(d.list2)))
  
  mu.matrix1 <- array(1/a2,dim=c(length(genes.s),1))%*%array(n,dim=c(1,length(n)))
  gene.components[,,1] <- dnbinom(target,1,1/(mu.matrix1+1),log=T)
  mu.matrix2a <- array(exp(2*(mu[1]+g.off.all2+sigma.all2[,1]^2/2))/
                         ((exp(sigma.all2[,1]^2)-1)*(exp(2*mu[1]+2*g.off.all2+sigma.all2[,1]^2))),
                       dim=c(length(genes.s),1))%*%array(1,dim=c(1,length(n)))
  mu.matrix2b <- 1/(1+(1/array(exp((mu[1]+g.off.all2+sigma.all2[,1]^2/2))/
                                 ((exp(sigma.all2[,1]^2)-1)*(exp(2*mu[1]+2*g.off.all2+sigma.all2[,1]^2))),
                               dim=c(length(genes.s),1))%*%array(1/n,dim=c(1,length(n)))))
  gene.components[,,2] <- dnbinom(target,mu.matrix2a,mu.matrix2b,log=T)
  mu.matrix3a <- array(exp(2*(mu[2]+g.on.all2+sigma.all2[,2]^2/2))/
                         ((exp(sigma.all2[,2]^2)-1)*(exp(2*mu[2]+2*g.on.all2+sigma.all2[,2]^2))),
                       dim=c(length(genes.s),1))%*%array(1,dim=c(1,length(n)))
  mu.matrix3b <- 1/(1+(1/array(exp((mu[2]+g.on.all2+sigma.all2[,2]^2/2))/
                                 ((exp(sigma.all2[,2]^2)-1)*(exp(2*mu[2]+2*g.on.all2+sigma.all2[,2]^2))),
                               dim=c(length(genes.s),1))%*%array(1/n,dim=c(1,length(n)))))
  gene.components[,,3] <- dnbinom(target,mu.matrix3a,mu.matrix3b,log=T)
  
  # For each cell-type...
  for (t in 1:length(d.list2)) {
    prob <- exp(sweep(t(gene.components[,,1]),2,d.list2[[t]][,1],'+'))+
      exp(sweep(t(gene.components[,,2]),2,d.list2[[t]][,2],'+'))+
      exp(sweep(t(gene.components[,,3]),2,d.list2[[t]][,3],'+'))
    flag <- which(prob==0,arr.ind=T)
    prob <- log(prob)
    for (f in seq_len(nrow(flag))) {
      prob[flag[f,1],flag[f,2]] <- sumlog(gene.components[flag[f,2],flag[f,1],]+as.numeric(d.list2[[t]][flag[f,2],]))
    }
    probs[,t] <- rowSums(prob)
  }
  
  return(probs)
}

# Get barcode (probability each gene is on) from trained reference
getBarcode <- function(d.list) {
  barcode <- sapply(d.list,function(x) 1-x[,2]-x[,3])
  barcode[barcode<0] <- 0
  colnames(barcode) <- names(d.list)
  rownames(barcode) <- rownames(d.list[[1]])
  return(barcode)
}

# Find markers for a given cell_type from a barcode object from getBarcode
# where the probability of being on is >thresh_up in that cell_type and 
# <thresh_below in others. Set relative=T if these markers should be relative 
# to the other cell-types in the barcode object, and relative=F if these 
# markers should be compared against the global probability parameters.
findMarkers <- function(cell_type,barcodes,thresh_up,thresh_below,relative=T,pi.all=pi.all.g) {
  if (relative) {
    return(barcodes[barcodes[,cell_type]>thresh_up&
            sapply(1:nrow(barcodes),function(x)all(barcodes[x,!colnames(barcodes)%in%c(cell_type)]<thresh_below)),])
  } else {
    pi.all2 <- pi.all[rownames(barcodes),]
    pi.on <- 1-pi.all2[,1]-pi.all2[,2]
    return(barcodes[barcodes[,cell_type]>thresh_up&pi.on<thresh_below,cell_type,drop=F])
  }
}


