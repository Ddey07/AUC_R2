require(readxl)
require(dplyr)
require(pROC)
library(glmnet)
library(modeest)
library(combinat)
library(GGally)
library(PRROC)
library(doParallel)
library(foreach)
library(tmvtnorm)
library(pcaPP)
require('foreign')
library(readxl)
library(plyr)
library(dplyr)
library(fCopulae)
library(corrplot)
library(gamlss)
library(huge)
library(MASS)
library(matrixcalc)
library(ggplot2)
#library(cowplot)
require(Matrix)
require(tibble)
require(edgebundleR)
require(pcaPP)
require("forecast")
require("ROCR")
require(DescTools)
require(glmgraph)
require(hierNet)
library(gridExtra)
library(DescTools)
library(RColorBrewer)

# October 2nd, 2017
# Author: Irina Gaynanova, parts of the code in estimateR_binary are taken from Yang Ning

# Given a binary matrix X, n x p, returns a p x p matrix of kendall tau values
Kendall_binary <- function(X){
  n <- nrow(X)
  p <- ncol(X)
  tau <- 2*t(X)%*%(diag(n) - matrix(1/n,n,n))%*%X/(n-1)
  return(tau)
}

# Given a a mixed (binry and continuous) matrix X, n x p, returns a p x p Kendall's tau matrix
ties=function(x){
  count=table(x)[which(table(x)>1)]
  return(sum(choose(count,2)))
}

Kendall_mixed=function(X){
  n=nrow(X)
  p=ncol(X)
  
  hat_tau=matrix(1,p,p)
  
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      time.ind=Sys.time()
      X_c=X[complete.cases(X[,c(i,j)]),c(i,j)]
      k=nrow(X_c)
      t=cor.fk(X_c)[1,2]
      k_1=ties(X_c[,1])
      k_2=ties(X_c[,2])
      k_0=choose(k,2)
      m=sqrt((k_0-k_1)*(k_0-k_2))/k_0
      hat_tau[i,j]=hat_tau[j,i]=t*m
      #print(c(i,j))
      #print(Sys.time()-time.ind)
    }}
  return(hat_tau)
  
}

# Given a a mixed (binary and continuous) matrix X, n x p, returns a p x p correlation matrix R which is an unbiased estimate of latent correlation matrix.
# tol - numerical tolerance for finding an inverse image of bridge function
# Part of the code is taken from Yang Ning, one of the authors of JRSS B paper
estimateR_mixed <- function(X, tol = 1e-6){
  #n sample size, total dimension p
  n = nrow(X)
  p = ncol(X)
  
  #Identify binary columns
  bin_id=which(apply(X,2,function(x) { all(na.omit(x) %in% 0:1) }))
  
  #bridge function F(t,\delta_j,\delta_k)
  bridgeF_bb = function(t,de1,de2){as.numeric(2*(pnorm2d(de1,de2,rho=t) - pnorm(de1)*pnorm(de2)))}
  bridgeF_cc = function(t){as.numeric(sin((pi/2)*t))}
  bridgeF_bc= function(t,de){as.numeric(4*(pnorm2d(de,0,rho=t/sqrt(2))) - 2*pnorm(de))}
  
  
  # Calculate hat delta for each column
  hatdelta=qnorm(1-colMeans(X,na.rm=TRUE))
  
  # Calculate Kendall matrix
  #K = Kendall_binary(X)
  K=Kendall_mixed(X)
  #print(is.positive.definite(K))
  #print(is.positive.semi.definite(K))
  #K=cor(X,method="kendall",use="pairwise")
  # calculate \hat R
  hatR=matrix(1,p,p)
  
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      time.ind=Sys.time()
      
      #Binary-Binary
      if(i %in% bin_id & j %in% bin_id){
        # Form a function to be set to zero as quadratic around K[i,j]
        f1 = function(t) (bridgeF_bb(t,hatdelta[i],hatdelta[j])-K[i,j])^2
        # Find the inverse
        op = tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = tol)[1], error = function(e) 100)
      }
      
      #Binary-Continuous
      else if(i %in% bin_id & !j %in% bin_id){
        # Form a function to be set to zero as quadratic around K[i,j]
        f1 = function(t) (bridgeF_bc(t,hatdelta[i])-K[i,j])^2
        # Find the inverse
        op = tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = tol)[1], error = function(e) 100)
      }
      
      else if(!i %in% bin_id & j %in% bin_id){
        # Form a function to be set to zero as quadratic around K[i,j]
        f1 = function(t) (bridgeF_bc(t,hatdelta[j])-K[i,j])^2
        # Find the inverse
        op = tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = tol)[1], error = function(e) 100)
      }
      
      #Continuous-Continuous
      else {
        op = bridgeF_cc(K[i,j])
      }
      
      if (op == 100){
        # This catches possible errors
        hatR[i,j]=hatR[j,i]=0
      }else {
        hatR[i,j]=hatR[j,i]=unlist(op)
      }
      #print(c(i,j))
      #print(Sys.time()-time.ind)
    }
  }  
  return(hatR)
}

CorrectCM = function(CM){
  n <- dim(var(CM))[1L]
  E <- eigen(CM)
  CM1 <- E$vectors %*% tcrossprod(diag(pmax(E$values, 0), n), E$vectors)
  Balance <- diag(1/sqrt(diag(CM1)))
  CM2 <- Balance %*% CM1 %*% Balance  
  return(CM2)
}

###########estimation of f in nonparanormal copula############
npn.f.est <- function(x){
  x <- na.omit(x)
  En <- ecdf(x)
  n <- length(x)
  Dn <- 1/(4*n^(1/4)*sqrt(pi * log(n)))
  Fn <- function(x){
    if(En(x) > Dn & En(x) < 1-Dn){
      En(x)
    }
    else if(En(x) <= Dn){
      Dn
    }
    else{
      1-Dn
    }
  }
  
  mu <- mean(x)
  sigma <- sqrt(var(x))
  
  
  f.hat <- function(x){
    qnorm(Fn(x))
  }
  return(f.hat)
}


glm.auc=function(x, T=Test){
  pred <- predict(x, newdata=T, type ="response")
  return(pROC::auc(pROC::roc(T$mortstat ~ pred)))
  
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- (cor(x, y, method="spearman"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  #text(0.5, 0.5, txt, cex = cex.cor * abs(r))
  text(0.5, 0.5, txt)
}

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}
