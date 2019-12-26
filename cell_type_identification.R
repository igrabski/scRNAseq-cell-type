library(dplyr)
library(Seurat)
library(Biobase)
library(reshape2)
library(rjags)

## Process reference data
process_singlecell <- function(data) {
  num_zero <- sum(rowSums(data)==0)
  data.f <- data[which(rowSums(data)!=0),]
  data.z <- data[which(rowSums(data)==0),]
  return(list(data=data.f,num_zero,data.z))
}

## Rjags model string to fit reference data
model_string <- "model{
  # Likelihood
  for (i in 1:n) {
    mu1[i] <- lambda1[i]
    mu3[i] <- N*lambda3[i]
    mu2[i] <- N*lambda2[i]
    lambda1[i] ~ dexp(a)
    lambda2[i] ~ dlnorm(c,sqrt(1/d))
    lambda3[i] ~ dlnorm(b,sqrt(1/f))
    m[i] ~ dcat(mprior[])
    mu[i] <- equals(m[i],1)*mu1[i] + equals(m[i],2)*mu2[i] + equals(m[i],3)*mu3[i]
    Y[i] ~ dpois(mu[i])
  }
  for (j in 1:nk) {
    zz[j] ~ dcat(zp[])
    mu_z[j] ~ dexp(a)
    mu_zf[j] <- equals(zz[j],1)*0 + equals(zz[j],2)*mu_z[j]
    Zeros[j] ~ dpois(mu_zf[j])
  }
  # Prior
  mprior[1] <- (1/3)
  mprior[2] <- (1/3)
  mprior[3] <- (1/3)
  a ~ dgamma(1,1)
  b ~ dnorm(0,1)
  c ~ dnorm(0,1)I(b,mx)
  d ~ dbeta(1,1)
  f ~ dbeta(1,1)
  zp[1] <- 0.5
  zp[2] <- 0.5
}"

## Fit processed reference data with a default of 5000 burn-in iterations
fit_mix <- function(process,iters=5000) {
  data <- rowSums(process[[1]])
  m <- ifelse(data>=quantile(data,0.75),2,3)
  m <- ifelse(data<=quantile(data,0.25),1,m)
  model <- jags.model(textConnection(model_string),inits=list(m=m),
                      data = list(Y=data,N=sum(data),n=length(data),Zeros=rep(0,process[[2]]),
                                  nk=process[[2]],mx=max(data)))
  update(model,iters)
  samp <- coda.samples(model,
                       variable.names=c("a","b","c","d","f","m","mu","mu_zf","zz"),
                       n.iter=2000)
  return(samp[[1]])
}

## Prepare estimated parameters for reference data
prep_df <- function(process,mix) {
  data <- rowSums(process[[1]])
  mus <- colMeans(mix)[(6+length(data)):(6+length(data)*2-1)]
  df <- as.data.frame(cbind(data,mus))
  mus <- colMeans(mix)[(6+length(data)*2):(6+length(data)*2+process[[2]]-1)]
  data <- rowSums(process[[3]])
  df <- rbind(df,as.data.frame(cbind(data,mus)))
  return(df)
}

## Compute probability of cells in target dataframe belonging to cell type, given estimated parameters
compute_prob <- function(df,target_df) {
  N <- sum(df$data)
  genes <- intersect(rownames(df),rownames(target_df))
  df <- df[genes,]
  target_df <- target_df[genes,]
  return(sapply(1:dim(target_df)[2],function(x) 
    sum(dpois(target_df[,x],lambda=(df$mus)*sum(target_df[,x])/N,log=TRUE))))
}

## Evaluate accuracy on target dataset
predict_type <- function(reference_dfs,target_df,labels,types) {
  num_cells <- length(labels)
  num_refs <- length(reference_dfs)
  probs <- as.data.frame(lapply(1:num_refs,function(y) compute_prob(reference_dfs[[y]],target_df)))
  chosen <- sapply(1:num_cells,function(x) types[which(probs[x,]==max(probs[x,]))])
  print(sum(chosen==labels)/num_cells)
  print(table(chosen,labels))
}
