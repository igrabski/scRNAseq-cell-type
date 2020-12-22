library(sads)

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

# Train reference data
trainReference <- function(ref,pi.all=pi.all.g,mu=mu.g,a=a.all.g,sigma.all=sigma.all.g,
                           g.on.all=g.on.all.g,g.off.all=g.off.all.g) {
  # Subset to genes present in both reference data and knowledge base
  common <- intersect(names(ref),rownames(pi.all))
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
    (1-rowSums(pi.all2)[j])*sads::dpoilog(ref[j],mu[2]+g.on.all2[j]+log(N),sigma.all2[j,2]))
  return(data.frame(rate=ref/N,exp=prob.exp/(prob.exp+prob.ln1+prob.ln2),ln1=prob.ln1/(prob.exp+prob.ln1+prob.ln2)))
}

# Train all reference
trainAllReference <- function(data,labels,pi.all=pi.all.g,mu=mu.g,a=a.all.g,sigma.all=sigma.all.g,
                              g.on.all=g.on.all.g,g.off.all=g.off.all.g) {
  d.list <- list()
  cell.types <- unique(labels)
  for (c in cell.types) {
    ref <- rowSums(data[,labels==c])
    d.list[[c]] <- trainReference(ref,pi.all=pi.all,mu=mu,a=a,sigma.all=sigma.all,
                                  g.on.all=g.on.all,g.off.all=g.off.all)
  }
  return(d.list)
}

# Classify target cells (fast)
classifyTarget <- function(target,d.list,pi.all=pi.all.g,mu=mu.g,a=a.all.g,sigma.all=sigma.all.g,
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
  d.list2 <- lapply(1:length(d.list),function(j) d.list[[j]][genes.s,])
  names(d.list2) <- names(d.list)
  
  # Store probabilities
  gene.components <- array(0,dim=c(length(genes.s),ncol(target),3))
  dimnames(gene.components)[[1]] <- genes.s
  probs <- array(0,dim=c(ncol(target),length(d.list2)))
  
  # For each gene...
  for (j in genes.s) {
    gene.components[j,,1] <- dnbinom(target[j,],1,1/((n/(a2[j]))+1))
    gene.components[j,,2] <- sapply(1:ncol(target),function(x)
      sads::dpoilog(target[j,x],mu[1]+g.off.all2[j]+log(n[x]),sigma.all2[j,1]))
    gene.components[j,,3] <- sapply(1:ncol(target),function(x)
      sads::dpoilog(target[j,x],mu[2]+g.on.all2[j]+log(n[x]),sigma.all2[j,2]))
  }
  
  # For each cell-type...
  for (t in 1:length(d.list2)) {
    prob.on2 <- 1-rowSums(d.list2[[t]][,2:3])
    prob.on2[prob.on2<0] <- 0
    probs[,t] <- rowSums(log(sweep(t(gene.components[,,1]),2,d.list2[[t]][,2],'*')+
                               sweep(t(gene.components[,,2]),2,d.list2[[t]][,3],'*')+
                               sweep(t(gene.components[,,3]),2,prob.on2,'*')))
  }
  
  return(names(d.list2)[sapply(1:ncol(target),function(x) which.max(probs[x,]))])
}

# Return likelihoods
likelihoodTarget <- function(target,d.list,pi.all=pi.all.g,mu=mu.g,a=a.all.g,sigma.all=sigma.all.g,
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
  d.list2 <- lapply(1:length(d.list),function(j) d.list[[j]][genes.s,])
  names(d.list2) <- names(d.list)
  
  # Store probabilities
  gene.components <- array(0,dim=c(length(genes.s),ncol(target),3))
  dimnames(gene.components)[[1]] <- genes.s
  probs <- array(0,dim=c(ncol(target),length(d.list2)))
  
  # For each gene...
  for (j in genes.s) {
    gene.components[j,,1] <- dnbinom(target[j,],1,1/((n/(a2[j]))+1))
    gene.components[j,,2] <- sapply(1:ncol(target),function(x)
      sads::dpoilog(target[j,x],mu[1]+g.off.all2[j]+log(n[x]),sigma.all2[j,1]))
    gene.components[j,,3] <- sapply(1:ncol(target),function(x)
      sads::dpoilog(target[j,x],mu[2]+g.on.all2[j]+log(n[x]),sigma.all2[j,2]))
  }
  
  # For each cell-type...
  for (t in 1:length(d.list2)) {
    prob.on2 <- 1-rowSums(d.list2[[t]][,2:3])
    prob.on2[prob.on2<0] <- 0
    probs[,t] <- rowSums(log(sweep(t(gene.components[,,1]),2,d.list2[[t]][,2],'*')+
                               sweep(t(gene.components[,,2]),2,d.list2[[t]][,3],'*')+
                               sweep(t(gene.components[,,3]),2,prob.on2,'*')))
  }
  
  return(probs)
}
