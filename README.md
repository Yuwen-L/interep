
<!-- README.md is generated from README.Rmd. Please edit that file -->

# interep

> **Inter**action Analysis of **Rep**eated Measure Data

<!-- badges: start -->

<!-- [![CRAN](https://www.r-pkg.org/badges/version/interep)](https://cran.r-project.org/package=interep) -->

<!-- [![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/interep)](http://www.r-pkg.org/pkg/interep) -->

[![Travis build
status](https://travis-ci.org/feizhoustat/interep.svg?branch=master)](https://travis-ci.org/feizhoustat/interep)
[![CRAN
status](https://www.r-pkg.org/badges/version/interep)](https://CRAN.R-project.org/package=interep)
[![Codecov test
coverage](https://codecov.io/gh/feizhoustat/interep/branch/master/graph/badge.svg)](https://codecov.io/gh/feizhoustat/interep?branch=master)
<!-- badges: end -->

Extensive penalized variable selection methods have been developed in the past two decades for analyzing high dimensional omics data, such as gene expressions, single nucleotide polymorphisms (SNPs), copy number variations (CNVs) and others. However, lipidomics data have been rarely investigated by using high dimensional variable selection methods. This package incorporates our recently developed penalization procedures to conduct interaction analysis for high dimensional lipidomics data with repeated measurements. The core module of this package is developed in C++. The development of this software package and the associated statistical methods have been partially supported by an Innovative Research Award from Johnson Cancer Research Center, Kansas State University.

## How to install
    - To install from Github, run these two lines of code in R

<!-- end list -->

    install.packages("devtools")
    devtools::install_github("Yuwen-L/interep")
    
  - Released versions of interep are available on CRAN
    [(link)](https://cran.r-project.org/package=interep), and can be
    installed within R via

<!-- end list -->

    install.packages("interep")

## Example

library(interep)
library(MASS)
set.seed(1000)

Data1 <- function(n,p,k,q,rho){
  # n: sample size; p: number of G factors; k: number of time points; q: number of E factors
  y = matrix(rep(0,n*k),n,k)  
  sig = matrix(0,p,p)  
  for (i in 1: p) {
    for (j in 1: p) { sig[i,j] = 0.5^abs(i-j)  }          
  }
  x = mvrnorm(n,rep(0,p),sig)
  g = x
  
  # generate binary variables
  dummy0 <- as.numeric(x[,1] <= -0.5)
  #dummy0 = dummy0+1
  dummy1 <- as.numeric(x[,1] > -0.5 & x[,1] <= 0)
  #dummy1 = dummy1*2+1
  dummy2 <- as.numeric(x[,1] > 0 & x[,1] <= 0.5)
  #dummy2 = dummy2*3+1
  
  # generate environment factors 
  e = cbind(dummy0,dummy1,dummy2)
  
  # set up the design matrix for the interaction model
  x=cbind(dummy0,dummy1,dummy2,x)
  for (i in (q+1):(p+q)) {
    for (j in 1:q) {
      x=cbind(x,x[,j]*x[,i])  }
  }
  
  x=scale(x)
  
  ll=0.4
  ul=0.8
  coef1=runif(q,ll,ul) # for interaction effects
  coef2=runif(q,ll,ul) # for interaction effects
  coef3=runif(q,ll,ul) # for interaction effects
  coef4=runif(7,ll,ul) # for E and G main effects
  coef=c(coef4,coef1,coef2,coef3)
  mat=x[,c(1,2,3,5,7,10,15,(p+q+1):(p+q+3),(p+5*q+1):(p+5*q+3),(p+10*q+1):(p+10*q+3))]
  
  for(u in 1:k){
    y[,u] =  0.6+rowSums(coef*mat) }
  sig1 = matrix(0,k,k) # AR(1) correlation
  diag(sig1)=1
  for (i in 1: k)  {
    for (j in 1: k)   { sig1[i,j] = rho^abs(i-j) }    
  }
  error = mvrnorm(n,rep(0,k),sig1)  
  y = y + error  
  index=1+c(1,2,3,5,7,10,15,(p+q+1):(p+q+3),(p+5*q+1):(p+5*q+3),(p+10*q+1):(p+10*q+3))
  dat = list(y=y,e=e,g=g,index=index,coef=coef)
  return(dat)    
}

n=250;p=75;k=5;q=3;rho=0.5
dat=Data1(n,p,k,q,rho)
e=dat$e
g=dat$g
y=dat$y

dim(e)

dim(g)

dim(y)

dat$coef

dat$index

index.true=dat $index
beta0=rep (0.1,1+q+p+p*q)

l1=c(0.45,0.5)
l2=c(0.8,1)
cv <- cv.interep(e,g,y,beta0,lambda1=l1,lambda2=l2,nfolds=5,corre="e",pmethod="mixed",maxits=30)

lambda1=cv$lam1
lambda2=cv$lam2
beta.est=interep(e, g, y, beta0, corre="e", pmethod="mixed", lam1=lambda1, lam2=lambda2, maxits=30)
beta.est[abs(beta.est)<0.05]=0
index.est=which (beta.est!=0)[-1]
tp=length(intersect(index.true, index.est))
fp=length(index.est)-tp

tp

fp



## Methods

This package provides implementation for methods proposed in

  - Zhou, F., Ren,J., Li, G., Jiang, Y., Li, X., Wang, W. and Wu, C. (2019). Penalized Variable Selection for Lipid–        Environment Interactions in a Longitudinal Lipidomics Study. Genes. 10(12), 1002
