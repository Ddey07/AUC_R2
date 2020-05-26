## clearing your environment and default packages

rm(list=ls())

detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
  
}

detachAllPackages()

### Loading necessary packages 
library(mvtnorm)
library(pcaPP)
library(Matrix)
library(pROC)
library(wCorr)
library(fMultivar)
library(dplyr)
library(spatstat)
library(Hmisc)
library(tidyr)
library(sampling)
library("nhanesdata")
library(rnhanesdata)
library("survey")
#library("gamlss")
library("SDMTools")
library("dplyr")
library(doParallel)
library(foreach)

####### Simulations ############### 

## necessary functions to load 

### changing default weighted.rank function in HMisc
wtd.rank <- function (x, weights = NULL, normwt = FALSE, na.rm = TRUE,...) 
{
  if (!length(weights)) 
    return(rank(x, na.last = if (na.rm) NA else TRUE))
  tab <- wtd.table(x, weights, normwt = normwt, na.rm = na.rm)
  freqs <- tab$sum.of.weights
  r <- cumsum(freqs) - 0.5 * (freqs - 1)
  approx(tab$x, r, xout = x,...)$y
}

assignInNamespace("wtd.rank","wtd.rank","Hmisc")

### Bridging functions
Gq=function(t,delta,m){
  if(m==0){
    return(pmvnorm(lower=c(delta,0),corr=matrix(c(1,t,t,1),ncol=2)) -  pmvnorm(lower=c(delta,-Inf),upper=c(Inf,0),corr=matrix(c(1,t,t,1),ncol=2)))
  } else {
    return(pmvnorm(upper=c(delta,0),corr=matrix(c(1,t,t,1),ncol=2)) -  pmvnorm(lower=c(-Inf,0),upper=c(delta,Inf),corr=matrix(c(1,t,t,1),ncol=2)))
  }
}

Gs=function(t,delta){
  p=pnorm(delta)
  SIG <- matrix(c(1,t/sqrt(2),t/sqrt(2),1),ncol=2)
  gamma <- pmvnorm(upper=c(0,-delta),sigma=SIG) + p*pmvnorm(lower=c(-Inf,-delta),upper=c(0,Inf),sigma=SIG)
  return(gamma)
}

Gk=function(t,delta){
  return(4*pnorm2d(delta,0,rho=t/sqrt(2)) - 2*pnorm(delta))
}

### simulation setting

seed <- seq(2,200,len=100)
rho.seq <- seq(0.005,0.995,len=50)
out.prop <- c(0,0.05,0.15)
conf.set <- c(0,1,2)
out.type <- c(1,2,3)

## the matrix setting contains all possible simulation scenarios to consider and id tells which situation to take
setting  <- expand.grid(rho=rho.seq,seed=seed,conf=conf.set,out.prop=out.prop)
#id <- as.numeric(Sys.getenv("SGE_TASK_ID"))
id <- 1 ## id can be varied from 1 to nrow(setting)

H <- 10
U <- 10
u <- 2
n <- u*H*(sum(rep(seq(0.2,1,length=H/2)))*2)*5
N <- n*100
tol <- 1e-6
r <- 0.5
d <- 0.2

seed0 = setting$seed[id]
r = setting$rho[id]

r1 <- r

### true AUC
p0 <- 1-pnorm(d)
At=abs(Gk(t=r,delta=d))/(4*p0*(1-p0)) + 1/2

### strata-specific AUCs
auc.h <- At + seq(-0.1,0.1,length.out = 10)
r.h <- numeric(length(auc.h))
for(h in 1:H){
  a=auc.h[h]
  fa=function(t){
    (abs(Gk(t,delta=d))/(4*p0*(1-p0)) + 1/2 - a)^2
  }
  r.h[h] = tryCatch(optimize(fa, lower = -0.99, upper = 0.99, tol = tol)[1], error = function(e) 100)$minimum
  if(sign(r.h[h])!=sign(r)){r.h[h]=-r.h[h]}
}

set.seed(seed0)

### simulating finite population data 
data <- NULL

for (h in 1:H) {
  
  data.h <- NULL
  data.h <- rmvnorm(N/H, sigma=matrix(c(1,r.h[h],r.h[h],1),ncol=2))
  data.h <- cbind(data.h,h)
  data <- rbind(data, data.h)
}


XY <- data[,c(1,2)]
strata <- data[,3]
psu <- lapply(table(strata), function(x){sample(c(1:U),x,replace=TRUE, prob= c(1:U))})
S <- data.frame(id = 1:N, strata = strata, psu = unlist(psu))

### informing strata with AUC (depending on stratifiedness)
auc.order <- rank(auc.h)
thg <- pp <- rep(NA,H)

quantile(1:H)[2]

confid= setting$conf[id]
if(confid==0){
  conf.prop <- rep(0.1,H)
} else if(confid==1){
  conf.prop <- rep(0.07,H)
  conf.prop[c(1:floor(quantile(1:H)[2]))] <- 0.06
  conf.prop[c((floor(quantile(1:H)[4])+1):H)] <- 0.18
} else {
  conf.prop <- rep(0.07,H)
  conf.prop[c(1:floor(quantile(1:H)[2]))] <- 0.02
  conf.prop[c((floor(quantile(1:H)[4])+1):H)] <- 0.22
}

pp <- conf.prop
sz <- round(n*pp)


### Two-stage stratified cluster sampling
psu.sample <- lapply(psu,function(x){sample(1:10,u,prob=as.numeric(table(x)))})
s <- unlist(lapply(1:H,function(x){sample((S %>% filter(strata == x & (psu %in% psu.sample[[x]])) %>% select(id))[,1], size=sz[x])}))
S$sampled <- as.numeric(S$id %in% s)


### calculating individual and joint selection probability for our simulation setting
Kh <- rep(U,H)
kh <- rep(u,H)
T <- t(lapply(psu,function(x){table(x)}) %>% bind_rows())
t <- matrix(nrow=H,ncol=U)
for(i in 1:H){
  for(j in 1:U){
    t[i,j] = sum(S$sampled[S$strata==i & S$psu==j])
  }
}

Pi2 <- matrix(nrow=n,ncol=n)

for(i in 1:(n-1)){
  for(j in (i+1):n){
    i0 <- S$id[which(S$sampled>0)[i]]
    j0 <- S$id[which(S$sampled>0)[j]]
    h <- S$strata[i0]
    hprime <- S$strata[j0]
    g <- S$psu[i0]
    gprime <- S$psu[j0]
    if(h == hprime & g==gprime){
      Pi2[i,j]=(kh[h]/Kh[h]) * (t[h,g]*(t[h,g]-1)) *(1/(T[h,g]*(T[h,g]-1)))
      Pi2[j,i]=Pi2[i,j]
    } else if(h == hprime & g!= gprime){
      Pi2[i,j]=(kh[h]/Kh[h]) * ((kh[h]-1)/(Kh[h]-1))* (t[h,g]*(t[h,gprime])) *(1/(T[h,g]*(T[h,gprime])))
      Pi2[j,i]=Pi2[i,j]
    } else {
      Pi2[i,j]=(kh[h]/Kh[h]) * ((kh[hprime])/(Kh[hprime]))* (t[h,g]*(t[hprime,gprime])) *(1/(T[h,g]*(T[hprime,gprime])))
      Pi2[j,i]=Pi2[i,j]
    }
  }
}

Pi <- numeric(n)
for(i in 1:n){
  i0 <- S$id[which(S$sampled>0)[i]]
  Pi[i] <- (kh[S$strata[i0]]/Kh[S$strata[i0]])*(t[S$strata[i0],S$psu[i0]]/T[S$strata[i0],S$psu[i0]])
}

W0 <- Pi
Wdash <- Pi2

smp <- which(S$sampled >0)

### Fixing outlyingness 

outp= setting$out.prop[id]
out.type= 1
if(outp >0){
  K <- 4
  outs <- sample(1:n,size=n*outp,replace = TRUE)
  mix <- sample(c(0,1),size=n*outp,replace = TRUE)
  if(out.type==1){
    XY[smp,][outs,] <- mix*rmvnorm(n*outp,mean=c(-K,K),sigma = 0.01^2*diag(2)) + (1-mix)*rmvnorm(n*outp,mean=c(K,-K),sigma = 0.01^2*diag(2)) 
  } else if (out.type==2){
    XY[smp,][outs,][,1] <- (mix*rmvnorm(n*outp,mean=c(-K,K),sigma = 0.01^2*diag(2)) + (1-mix)*rmvnorm(n*outp,mean=c(K,-K),sigma = 0.01^2*diag(2)))[,1]
  } else { 
    XY[smp,][outs,][,2] <- mix*rmvnorm(n*outp,mean=c(-K,K),sigma = 0.01^2*diag(2)) + (1-mix)*rmvnorm(n*outp,mean=c(K,-K),sigma = 0.01^2*diag(2))[,2] 
  }
}

### picking a finite sample from the finite population
X <- XY[smp,][,2]
Y <- as.numeric((XY[smp,])[,1] > d)

Z1=X[Y==1]
Z2=X[Y==0]

n1 <- which(Y==1)
n2 <- which(Y==0)

W1 = (1/W0)[Y==1]
W2= (1/W0)[Y==0]


### Estimation

### Calculating rank statistics

tag0=0
tag=0
tag2=0
wtsum=0
tagk=0
tagkw=0
tagkwt=0
kw=0
kwt=0

for(i in 1:(n-1)){
  for(j in (i+1):n){
    tagk=tagk + sign(X[i]-X[j])*sign(Y[i]-Y[j])
    tagkw=tagkw + (1/W0[i])*(1/W0[j])*sign(X[i]-X[j])*sign(Y[i]-Y[j])
    kw= kw + (1/W0[i])*(1/W0[j])
    tagkwt=tagkwt + (1/Wdash[i,j])*sign(X[i]-X[j])*sign(Y[i]-Y[j])
    kwt= kwt + (1/Wdash[i,j])
  }
}

kt=tagk*2/(n*(n-1))
ktw=tagkw/kw
ktwt= tagkwt/kwt

rqw=sum((1/W0)*sign((X-wtd.quantile(X,weights=1/W0,probs=0.5,normwt=TRUE))*(Y-wtd.quantile(Y,weights=1/W0,probs=0.5,normwt=TRUE))))/sum(1/W0)
rsw=weighted.mean(ewcdf(X,weights=(1/W0)/sum(1/W0))(X)*ewcdf(Y,weights = (1/W0)/sum(1/W0))(Y),1/W0)
rsw=weighted.mean((wtd.rank(X,weights = 1/W0,rule=2)/(sum(1/W0)+1))*ewcdf(Y,weights = (1/W0)/sum(1/W0))(Y),1/W0)

H=ewcdf(X,weights=(1/W0)/sum(1/W0))
What=abs((weighted.mean(H(X[Y==1]),W1) - weighted.mean(H(X[Y==0]),W2)))+0.5

hatdelta=qnorm(1-weighted.mean(Y,1/W0))
my=wtd.quantile(Y,weights=1/W0,probs=0.5,normwt=TRUE)

### Inverting bridging function to calculate latent correlation

fq = function(t) (Gq(t,delta=hatdelta,m=my)-rqw)^2
opq = tryCatch(optimize(fq, lower = -0.99, upper = 0.99, tol = tol)[1], error = function(e) 100)$minimum

fs=function(t) (Gs(t,delta=hatdelta)-rsw)^2
ops = tryCatch(optimize(fs, lower = -0.99, upper = 0.99, tol = tol)[1], error = function(e) 100)$minimum

fk=function(t) (Gk(t,delta=hatdelta)-ktw)^2
opk = tryCatch(optimize(fk, lower = -0.99, upper = 0.99, tol = tol)[1], error = function(e) 100)$minimum


p=weighted.mean(Y,1/W0)
puw=mean(Y)

### AUC estimates from bridging
Aq=abs(Gk(t=opq,delta=hatdelta))/(4*p*(1-p)) + 1/2
As=abs(Gk(t=ops,delta=hatdelta))/(4*p*(1-p)) + 1/2
Ak=abs(kt)/(4*puw*(1-puw)) + 1/2
Akw=abs(ktw)/(4*p*(1-p)) + 1/2
Akwt= abs(ktwt)/(4*p*(1-p)) + 1/2

### True AUC
p0 <- 1-pnorm(d)
At=abs(Gk(t=r,delta=d))/(4*p0*(1-p0)) + 1/2

## Latent R^2
lrq=opq^2
lrs=ops^2
lrk=opk^2

Res <- data.frame(rt=r,out.type=out.type,out.prop=outp,infsamp=confid,true=At,puw=puw,p=p,quad=Aq,spear=As,lrq=lrq,lrs=lrs,lrk =lrk ,W=What, kendall=Ak,kendallw=Akw,kendallwt=Akwt)

set.seed(seed0)
results <- Res

Res
### One has to run the id through 1 to 45000 to explore all the settings of our simulations which will also take care of different seeds. 

######### Data Example ######### 
data("lookup_full")

## Getting variable data
cohorts <- c("2003-2004","2005-2006")
var1 <- "URXUMA"
var2 <- "BMXBMI"
agerange <- c(0,85)

n <- length(cohorts)

#loading first variable
data1name <- tolower(lookup.full$`Data File Name`[which(lookup.full$`Variable Name` == var1 & lookup.full$Year %in% cohorts)])
data2name <- tolower(lookup.full$`Data File Name`[which(lookup.full$`Variable Name`=="RIDAGEYR" & lookup.full$Year %in% cohorts)])
if(identical(data1name,character(0))==TRUE | identical(data2name,character(0))==TRUE){
  print("Problem with data at cohort 1")
  return(0)
}
A <- dplyr::select(get(data1name),c("SEQN",var1))
B <- dplyr::select(get(data2name),c("SEQN","RIDAGEYR","RIAGENDR","WTMEC2YR","SDMVPSU","SDMVSTRA","WTINT2YR"))
dataset1 <- left_join(B,A,by="SEQN")
dataset1$Year <- cohorts[1]

if(n>=2){
  for(i in 2:n){
    data1name <- tolower(lookup.full$`Data File Name`[which(lookup.full$`Variable Name`==var1 & lookup.full$Year==cohorts[i])])
    data2name <- tolower(lookup.full$`Data File Name`[which(lookup.full$`Variable Name`=="RIDAGEYR" & lookup.full$Year==cohorts[i])])
    if(identical(data1name,character(0))==TRUE | identical(data2name,character(0))==TRUE){
      print(paste("Problem with data at cohort", i))
      return(0)
    }
    A <- dplyr::select(get(data1name),c("SEQN",var1))
    B <- dplyr::select(get(data2name),c("SEQN","RIDAGEYR","RIAGENDR","WTMEC2YR","SDMVPSU","SDMVSTRA","WTINT2YR"))
    dataset_temp <- left_join(B,A,by="SEQN")
    dataset_temp$Year <- cohorts[i]
    dataset1 <- rbind.data.frame(dataset1,dataset_temp)
  }
}
rm(A,B,dataset_temp)
dataset1 <- na.omit(dataset1)
dataset1 <- dplyr::filter(dataset1, RIDAGEYR %in% c(agerange[1]:agerange[2]))
dataset1 <- dataset1 %>% dplyr::select(SEQN,URXUMA)

### Getting physical activity variables
load("pavars2003_04_2005_06.RData")


### Read mortality data and other variables
##Sourcing package and functions
source('Functions.R')

#Loading uber dataset
load("Data_full_NHANES.RData")

#Systloic and Diastolic average
data_full$BPXSY=rowMeans(cbind(data_full$BPXSY1,data_full$BPXSY2,data_full$BPXSY3,data_full$BPXSY4),na.rm = TRUE)
data_full$BPXDI=rowMeans(cbind(data_full$BPXDI1,data_full$BPXDI2,data_full$BPXDI3,data_full$BPXDI4),na.rm = TRUE)

data_full = data_full %>% dplyr::filter(waveID %in% c("C","D")) %>% dplyr::select(SEQN,BMXBMI,BPXSY,SDDSRVYR,SDMVPSU,SDMVSTRA,WTINT2YR,WTMEC2YR,RIDAGEYR,RIAGENDR)
data_full$SEQN <- as.integer(data_full$SEQN)


### Mortality data
mort_ls <- rnhanesdata::process_mort()
mort <- rbind(((mort_ls$Mortality_2011_C) %>% dplyr::select(SEQN,mortstat,ucod_leading,permth_exm)),((mort_ls$Mortality_2011_D) %>% dplyr::select(SEQN,mortstat,ucod_leading,permth_exm)))

covar_ls <- process_covar()

## re-code gender for the both the 2003-2004 and 2005-2006 waves
covar_ls$Covariate_C$Gender <- factor(covar_ls$Covariate_C$RIAGENDR, levels=1:2,
                                      labels=c("Male","Female"), ordered=FALSE)
covar_ls$Covariate_D$Gender <- factor(covar_ls$Covariate_D$RIAGENDR, levels=1:2,
                                      labels=c("Male","Female"), ordered=FALSE)

andrewfulldata <- rbind(Covariate_C,Covariate_D)

### Joining data frames by sequence numbers for 2003-2004 and 2005-2006 cohort
ourdata <- data_full %>% full_join(mort) %>% full_join(dataset1) %>% full_join(data_analysis) %>% full_join(andrewfulldata)

## subsetting data to focus on our intended subgroup and redefining mortality outcome as 5-year mortality
ourdatafull= subset(ourdata,!is.na(TAC) & (is.na(ucod_leading) | ucod_leading !="004") & !is.na(RIDAGEYR) & RIDAGEYR >=50 & RIDAGEYR <85 & (!is.na(mortstat)) & !(mortstat==0 & permth_exm < 60) & !is.na(URXUMA) &!is.na(MVPA) &!is.na(ASTP) &!is.na(BPXSY))
ourdatafull$mortstat = as.numeric(ourdatafull$permth_exm <= 60 & ourdatafull$mortstat==1)

### variable list under consideration
varlist <- c("RIDAGEYR","URXUMA","TAC","MVPA","ASTP","BPXSY")

### function to calculate estiamte and bootstrap confidence intervals
se.bootstrap=function(var,b=5){
  ind=which(colnames(ourdatafull)==var)
  ourdata = ourdatafull[!is.na(ourdatafull[,ind]),]
  
  ourdata= reweight_accel(ourdata)
  
  nhanes_design <- 
    svydesign(
      id = ~SDMVPSU , 
      strata = ~SDMVSTRA ,
      nest = TRUE ,
      weights = ~wtmec4yr_adj,
      data = ourdata
    )
  
  aucse <- function(w,data){
    XY <- data.frame(data[,ind],data$mortstat,data$wtmec4yr_adj)
    XY <- XY[complete.cases(XY),]
    X <- XY[,1]
    Y <- XY[,2]
    W0 <- 1/w
    n <- length(X)
    Yord <- Y[order(X)]
    
    # fit <- glm(Y ~ X, family=binomial)
    # pred <- predict(fit)
    # 
    Z1=X[Y==1]
    Z2=X[Y==0]
    
    n1 <- which(Y==1)
    n2 <- which(Y==0)
    
    W1 = (1/W0)[Y==1]
    W2= (1/W0)[Y==0]
    
    # if(rho>=0){Z1=X[Y==1]
    # Z2=X[Y==0]} else {
    #   Z1=X[Y==0]
    #   Z2=X[Y==1]
    # }
    tag0=0
    tag=0
    tag2=0
    wtsum=0
    
    tagk=0
    tagkw=0
    tagkwt=0
    kw=0
    kwt=0
    sn=0
    
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        tagk=tagk + sign(X[i]-X[j])*sign(Y[i]-Y[j])
        tagkw=tagkw + (1/W0[i])*(1/W0[j])*sign(X[i]-X[j])*sign(Y[i]-Y[j])
        kw= kw + (1/W0[i])*(1/W0[j])
        #sn = sn + (j-i)*(Yord[i] > Yord[j])
      }
    }
    
    kt=tagk*2/(n*(n-1))
    ktw=tagkw/kw
    #rsn = 1-((12*sn)/(n*(n^2-1)))
    
    for(i in 1:length(Z1)){
      for(j in 1:length(Z2)){
        tag=tag + W1[i]*W2[j]*as.numeric(Z1[i] > Z2[j])
        #tag2=tag2 + (1/Wdash[n1[i],n2[j]])*as.numeric(Z1[i] > Z2[j])
        #wtsum=wtsum+(1/Wdash[n1[i],n2[j]])
        tag0=tag0 +as.numeric(Z1[i] > Z2[j])
      }
    }
    
    Aw = tag/ (sum(W1)*sum(W2))
    Aw=max(Aw,1-Aw)
    Awt=tag2/wtsum
    Awt=max(Awt,1-Awt)
    Auwt=tag0/(length(Z1)*length(Z2))
    Auwt=max(Auwt,1-Auwt)
    
    # rqw= sum(W0*sign((X-median(X))*(Y-median(Y))))/sum(W0)
    # rsw= weighted.mean(ecdf(X)(X)*ecdf(Y)(Y),W0)
    # 
    
    rqw=sum((1/W0)*sign((X-wtd.quantile(X,weights=1/W0,probs=0.5,normwt=TRUE))*(Y-wtd.quantile(Y,weights=1/W0,probs=0.5,normwt=TRUE))))/sum(1/W0)
    rsw=weighted.mean(ewcdf(X,weights=(1/W0)/sum(1/W0))(X)*ewcdf(Y,weights = (1/W0)/sum(1/W0))(Y),1/W0)
    rs = mean(ecdf(X)(X)*ecdf(Y)(Y))
    
    
    rs=mean((rank(X)/(n+1)) * rank(Y,ties.method = "max")/(n+1))
    rsw=weighted.mean((wtd.rank(X,weights = 1/W0,rule=2)/(sum(1/W0)+1))*ewcdf(Y,weights = (1/W0)/sum(1/W0))(Y),1/W0)
    
    
    ###Calculating W
    H=ewcdf(X,weights=(1/W0)/sum(1/W0))
    What=abs(weighted.mean(H(X[Y==1]),W1) - weighted.mean(H(X[Y==0]),W2))+0.5
    
    hatdelta=qnorm(1-weighted.mean(Y,1/W0))
    #my=median(Y)
    my=wtd.quantile(Y,weights=1/W0,probs=0.5,normwt=TRUE)
    
    ### Inverting bridging functions
    fq = function(t) (Gq(t,delta=hatdelta,m=my)-rqw)^2
    opq = tryCatch(optimize(fq, lower = -0.99, upper = 0.99, tol = tol)[1], error = function(e) 100)$minimum
    
    fs = function(t) (Gs(t,delta=hatdelta)-rsw)^2
    ops = tryCatch(optimize(fs, lower = -0.99, upper = 0.99, tol = tol)[1], error = function(e) 100)$minimum
    
    fk=function(t) (Gk(t,delta=hatdelta)-ktw)^2
    opkw = tryCatch(optimize(fk, lower = -0.99, upper = 0.99, tol = tol)[1], error = function(e) 100)$minimum
    
    fk=function(t) (Gk(t,delta=hatdelta)-kt)^2
    opk = tryCatch(optimize(fk, lower = -0.99, upper = 0.99, tol = tol)[1], error = function(e) 100)$minimum
    
    p=weighted.mean(Y,1/W0)
    
    ### AUC estimates from bridging
    Aq=abs(Gk(t=opq,delta=hatdelta))/(4*p*(1-p)) + 1/2
    #Asw=abs(Gk(t=opsw,delta=hatdelta))/(4*p*(1-p)) + 1/2
    As=abs(Gk(t=ops,delta=hatdelta))/(4*p*(1-p)) + 1/2
    Ak=abs(kt)/(4*p*(1-p)) + 1/2
    Akw=abs(ktw)/(4*p*(1-p)) + 1/2
    #Akwt= abs(ktwt)/(4*p*(1-p)) + 1/2
    
    ## Latent R^2
    lrq=opq^2
    lrs=ops^2
    lrk = opk^2
    lrkw=opkw^2
    
    Res <- c(wt=Aw,kendalluw=Ak,kendall=Akw,W=What,spear=As,quad=Aq,lrq=lrq,lrs=lrs,lrkuw=lrk,lrkpw=lrkw)
    return(Res)
  }
  
  nhanes_design_bs <- as.svrepdesign(nhanes_design,type="bootstrap",replicates=b)
  res <- withReplicates(nhanes_design_bs,aucse,return.replicates = TRUE)
  
  Est <- data.frame(var=var,coef(res))
  Bootse <- data.frame(var=var,SE(res))
  bootci <- apply(res$replicates,2,function(x){quantile(x,c(0.025,0.975))})
  
  out1 <- cbind(Est,Bootse)
  l <- rownames(out1)
  l[1] <- c("wt")
  l[3] <- c("kendall")
  l1 <- rep(l,each="4")
  l1[seq(2,38,by=4)] <- paste0(l,".se")
  l1[seq(3,39,by=4)] <- paste0(l,".upr")
  l1[seq(4,40,by=4)] <- paste0(l,".lwr")
  
  out0 <- data.frame(t(as.numeric(t(cbind(out1[,2],out1[,4],bootci[1,],bootci[2,]))))) 
  names(out0) <- l1
  out0$var <- var
  # Est <- data.frame(var=var,wt=coef(Awrep)[1],quad=coef(Aqrep)[2],kendall=coef(Akrep),spear=coef(Asrep),lrq=coef(lrqrep),lrs=coef(lrsrep),lrk=coef(lrkrep),W=coef(Wrep),unwtd=coef(Auwtrep))
  # SE <- data.frame(var=var,wt=SE(Awrep),quad=SE(Aqrep),kendall=SE(Akrep),spear=SE(Asrep),lrq=SE(lrqrep),lrs=SE(lrsrep),lrk=SE(lrkrep),W=SE(Wrep),unwtd=SE(Auwtrep))
  
  return(out0)
}

### computing 2 replicate bootstrap SE interval without parallelizing
Res_with_se <- lapply(varlist,function(x){se.bootstrap(x,b=2)}) %>% bind_rows()

### This is a time-consuming step for computing 100 replicate bootstrap SE interval, so, we parallelize
Res_with_se <- foreach(x=varlist,.combine = rbind) %dopar%{
  se.bootstrap(x,b=100)
}

save(Res_with_se,file="NHANES_results.RData")

