rm(list=ls())
set.seed(1234)
library(pROC)
library(igraph)
#library(ggplot2)
library(plyr)
library(coda)
library(Metrics)
library(mltools)
#library(ggpubr)
#theme_set(theme_minimal())
#setwd("~/Dropbox/Bayesian Graphical Models/BMIM/phd-Andrea/Out")
#setwd("C:/Users/aleav/Dropbox/Bayesian Graphical Models/BMIM/phd-Andrea/Out")

X <- 4
p <- 10
edges <- p*(p-1)/2
sims <- 10
K2<-200000
K<-50000
n<-100
burn<-n
alpha<--1
beta<-1.5

###### SCENARIO B with scale-free graphs ######
g1 <- make_empty_graph(p)
g2 <- make_empty_graph(p)
g3 <- make_empty_graph(p)
g4 <- make_empty_graph(p)
# Use the Barabási-Albert model to grow the network
g1 <- barabasi.game(p, power = 1, m = 1, out.pref = TRUE)
g2 <- barabasi.game(p, power = 1, m = 1, out.pref = TRUE)
g3 <- barabasi.game(p, power = 1, m = 1, out.pref = TRUE)
g4 <- barabasi.game(p, power = 1, m = 1, out.pref = TRUE)

# Plot the network (optional)
plot(g1, layout = layout.fruchterman.reingold)
plot(g2, layout = layout.fruchterman.reingold)
plot(g3, layout = layout.fruchterman.reingold)
plot(g4, layout = layout.fruchterman.reingold)

g1
g2
g3
g4

G0 <- as.matrix(g1[1:10])
G0[upper.tri(G0)] = t(G0)[upper.tri(G0)]

G1 <- as.matrix(g2[1:10])
G1[upper.tri(G1)] = t(G1)[upper.tri(G1)]

G2<-as.matrix(g3[1:10])
G2[upper.tri(G2)] = t(G2)[upper.tri(G2)]

G3<-as.matrix(g4[1:10])
G3[upper.tri(G3)] = t(G3)[upper.tri(G3)]


DATA<-list()


RBEL<-matrix(0,sims,2)

predBEL<-list()

AUCBEL<-c()

OutNuKBEL <- list()
OutThetaKBEL <- list()

for (x in 1:X){
  
  predBEL[[x]]<-matrix(0,sims,edges)
  
  
}



D<-matrix(0,K,edges)

ciccio<-list()

PREDBEL<-list()


for (t in 1:10){
  
  ciccio[[t]]<-D
  
}

for (x in 1:X){
  
  PREDBEL[[x]]<-ciccio
  
  
  
}

t00 <- Sys.time()
for (P in 1:sims){ # INIZIO
  
  
  #########################################################################
  #################            G0-G1-G2-G3             ####################
  #########################################################################
  dataf<-matrix(0,K2,p)
  
  dataf[1,]<-rep(1,p)
  
  for (k in 2:K2){
    
    for (r in 1:p){
      
      s<-sum(dataf[k,which(G0[r,]==1)])
      
      odds<-exp(beta*s+alpha)
      
      prob<-odds/(1+odds)
      
      dataf[k,r]<-rbinom(1,1,prob)
      
    }
    if(k%%1000==0){print(paste("data0 K",k," P",P))}
  }
  
  data<-dataf[(nrow(dataf)-(n-1)):nrow(dataf),]
  
  data0<-cbind(data,0)
  
  
  dataf<-matrix(0,K2,p)
  
  dataf[1,]<-rep(1,p)
  
  for (k in 2:K2){
    
    for (r in 1:p){
      
      s<-sum(dataf[k,which(G1[r,]==1)])
      
      odds<-exp(beta*s+alpha)
      
      prob<-odds/(1+odds)
      
      dataf[k,r]<-rbinom(1,1,prob)
      
    }
    if(k%%1000==0){print(paste("data1 K",k," P",P))}
  }
  
  data<-dataf[(nrow(dataf)-(n-1)):nrow(dataf),]
  
  data1<-cbind(data,1)
  
  
  
  
  
  dataf<-matrix(0,K2,p)
  
  dataf[1,]<-rep(1,p)
  
  for (k in 2:K2){
    
    for (r in 1:p){
      
      s<-sum(dataf[k,which(G2[r,]==1)])
      
      odds<-exp(beta*s+alpha)
      
      prob<-odds/(1+odds)
      
      dataf[k,r]<-rbinom(1,1,prob)
      
    }
    if(k%%1000==0){print(paste("data2 K",k," P",P))}
  }
  
  data<-dataf[(nrow(dataf)-(n-1)):nrow(dataf),]
  
  data2<-cbind(data,2)
  
  
  
  
  dataf<-matrix(0,K2,p)
  
  dataf[1,]<-rep(1,p)
  
  for (k in 2:K2){
    
    for (r in 1:p){
      
      s<-sum(dataf[k,which(G3[r,]==1)])
      
      odds<-exp(beta*s+alpha)
      
      prob<-odds/(1+odds)
      
      dataf[k,r]<-rbinom(1,1,prob)
      
    }
    if(k%%1000==0){print(paste("data3 K",k," P",P))}
  }
  
  data<-dataf[(nrow(dataf)-(n-1)):nrow(dataf),]
  
  data3<-cbind(data,3)
  
  
  
  ####### FINE DATA
  
  ###### INIZIO REGISTRAZIONE
  
  dataX<-list()
  
  dataX[[1]]<-data0[,1:p]
  dataX[[2]]<-data1[,1:p]
  dataX[[3]]<-data2[,1:p]
  dataX[[4]]<-data3[,1:p]
  
  
  
  DATA[[P]]<-dataX
  
  
}
t0 <- Sys.time()
t0-t00 #Time difference of 6.079958 mins

GX<-list()
GX[[1]]<-G0
GX[[2]]<-G1
GX[[3]]<-G2
GX[[4]]<-G3

MX.test<-list()


MX.test[[1]]<-as.vector(t(G0[upper.tri(G0)]))
MX.test[[2]]<-as.vector(t(G1[upper.tri(G1)]))
MX.test[[3]]<-as.vector(t(G2[upper.tri(G2)]))
MX.test[[4]]<-as.vector(t(G3[upper.tri(G3)]))


######FINE REGISTRAZIONE
t1 <-Sys.time()

for (P in 1:sims){
  
  
  dataX<-list()
  
  dataX[[1]]<-DATA[[P]][[1]]
  dataX[[2]]<-DATA[[P]][[2]]
  dataX[[3]]<-DATA[[P]][[3]]
  dataX[[4]]<-DATA[[P]][[4]]
  
  
  ####################### Bayesian exact linking ######################################
  
  ####################### MODEL SELECTION ALGORITHM ######################################
  
  
  V<-p+1
  
  # MOEBIUS MATRIX  
  
  
  mobkron<- function(V){
    
    h<-matrix(c(1,1,0,1),2,2)
    
    mat<-list(c())
    
    mat[[1]]<-h
    
    for (i in 2:V){
      
      mat[[i]]<- kronecker(mat[[i-1]],h)
      
    }
    
    return(t(mat[[i]]))
    
    
  }
  
  K<-50000
  
  nuK<-matrix(0,K,(p*(p-1)/2))
  
  
  
  aa<-1
  bb<-1
  
  aap<-aa
  bbp<-bb
  
  
  theta<-matrix(0,X,X)
  diag(theta)<-0
  w<-0.6
  
  alpha<-1
  chi<-2
  
  de0<-1
  
  epsilon<-matrix(0,X,X)
  diag(epsilon)<-0
  
  
  epsilonK<-list()
  thetaK<-list()
  
  epsilonK[[1]]<-epsilon
  thetaK[[1]]<-theta
  
  
  
  
  
  pdelta<-function(i,x,nu,the,lrx2){
    
    
    for (j in 1:length(lrx2[[x]][i,])){
      
      
      LrXrest<-c()
      
      for (z in 1:X){
        
        LrXrest[z]<-lrx2[[z]][i,j]
        
      }
      
      
      odds<-exp( lrx2[[x]][i,j]%*%(nu[j]+the[x,]%*%LrXrest ) )
      
      prob<-odds/(1+odds)
      
      q[j]<-prob
      
      
    }
    
    PROB<-prod(q)
    
    return(PROB)
    
  }
  
  
  
  LX<-list()
  L<-list(c())
  
  
  for (r in 1:p){
    L[[r]]<-matrix(0,K,p)
  }
  
  for (x in 1:X){
    
    LX[[x]]<-L
    
  }
  
  LcX<-list()
  
  aX<- list()
  
  for (x in 1:X){
    
    aX[[x]]<-1
    
  }
  
  
  ## START
  
  
  c.cX<-list()
  yX<-list()
  sX<-list()
  LrX<-list()
  lrX<-list()
  log.M.LX<-list()
  ldifX<-list()
  i.freeX<-list()
  la.prX<-list()
  start.prX<-list()
  pX<-list()
  SX<-list()
  H.prX<-list()
  log.k.prX<-list()
  H.poX<-list()
  log.k.poX<-list()
  FITX<-list()
  fit.lmX<-list()
  fit.l2X<-c()
  log.M.LX<-c()
  ldifX<-c()
  la.poX<-list()
  l.post.p<-list()
  
  PdeltaX<-list()
  
  for (x in 1:X){
    
    PdeltaX[[x]]<-0
    
  }
  
  for (x in 1:X){
    
    l.post.p[[x]]<-0
    
  }
  
  
  #B.llm.S<-function (V,dataX,K,a){ 
  
  # V   = total nodes + X
  # K   = iterations number
  # c.c = vector cell counts
  # a   = vector prior cell counts
  
  a<-0.02
  
  Z.k<-mobkron(V-1) # matrice Z  (V-1 x V-1)
  Z.kt<-mobkron(V)  # matrice Z (V x V)
  M.k<-solve(Z.k)   # matrice M                                          
  
  It<-2^V           # spazio degli stati con X
  I<-2^(V-1)        # spazio degli stati senza X
  p.p<-rep(a/I,I)   # probabilit? a priori (uniforme)
  c.c.p<-p.p*a      # conteggi a priori (uniforme)
  
  # genero conteggi da dataset
  N <- ncol(dataX[[x]])       # Number of nodes
  nSample <- nrow(dataX[[x]]) # Number of samples
  Resp <- c(0L,1L)
  n<- nrow(dataX[[x]])
  
  
  
  # genero vettore selezionatore  
  L<-rep(8,(2^(V-1))-1)  
  tZ<-t(Z.k)[-1,-1]
  tZt<-kronecker(matrix(c(1,0,0,1),2,2),tZ)
  
  
  
  for(i in 1:(2^(V-1)-1)){
    
    if ( sum(tZ[i,]) == 1)
      
      L[i]<-1
    
  } 
  
  
  p.LM<-which(L==1) # posizione effetti principali
  
  
  
  
  
  for(i in 1:(2^(V-1)-1)){
    
    if ( sum(tZ[i,]) == 3)
      
      L[i]<-2
    
  } 
  
  
  p.L2<-which(L==2) # posizione interazioni doppie
  
  
  p.LZ<-which(L==8) # posizione interazioni tiple e oltre
  
  
  
  ch<-choose(V-1,2)
  
  l<-L
  
  
  LrX2<-list()
  
  for (x in 1:X){
    
    LrX2[[x]]<-matrix(0,K,p*(p-1)/2)
    
  }
  
  #####################################
  
  for (x in 1:X){
    
    
    AllStates <- do.call(expand.grid,lapply(1:N,function(x)Resp))
    c.cX[[x]] <- apply(AllStates,1,function(z)sum(colSums(t(dataX[[x]]) == z)==N)) # conteggi tabella di contingenza
    m.c.c<-Z.k%*%c.cX[[x]]      
    m.c.c.p<-Z.k%*%c.c.p
    yX[[x]]<-m.c.c[-1]          # conteggi marginali osservati
    sX[[x]]<-m.c.c.p[-1]        # conteggi marginali a priori
    
    
    
    
    log.M.LX[[x]]<-1
    ldifX[[x]]<-1
    
    LX[[x]]<-matrix(0,K,length(l))
    LX[[x]][1,p.LM]<-rep(1,length(p.LM))
    LX[[x]][1,p.LZ]<-8
    lrX[[x]]<-LX[[x]][1,-p.LZ]
    LrX[[x]]<-matrix(0,K,length(lrX[[x]]))
    i.freeX[[x]]<-which(LX[[x]][1,]==1) # given a model for 0 we identify the parameter that are proposed free
    
    LrX2[[x]][1,]<-rep(0,p*(p-1)/2)
    
    # INITIALIZING
    
    
    
    
    # Prior parameters estimates for saturated model 
    
    
    
    la.prX[[x]]<-t(M.k)%*%log(p.p) # check log-linear parameters
    
    la.prX[[x]]<-la.prX[[x]][-1]
    
    start.prX[[x]]<-la.prX[[x]][i.freeX[[x]]]
    
    pX[[x]]<-1
    
    
    
    
    
    
    
    
    
    H.prX[[x]]<- - (a / (1+t(rep(1,I-1)) %*% exp(tZ[,i.freeX[[x]]] %*% la.prX[[x]][i.freeX[[x]]]) ))[1,1] *
      
      (t(tZ[,i.freeX[[x]]]) %*% diag(as.vector(exp(tZ[,i.freeX[[x]]]%*%la.prX[[x]][i.freeX[[x]]]))) %*% tZ[,i.freeX[[x]]]) +
      
      (t(tZ[,i.freeX[[x]]]) %*% diag(as.vector(exp(tZ[,i.freeX[[x]]]%*%la.prX[[x]][i.freeX[[x]]]))) %*% rep(1,I-1) %*%
         
         (a / (1+t(rep(1,I-1)) %*% exp(tZ[,i.freeX[[x]]] %*%la.prX[[x]][i.freeX[[x]]]) )^2)[1,1] %*%
         
         t(t(tZ[,i.freeX[[x]]]) %*% diag(as.vector(exp(tZ[,i.freeX[[x]]]%*%la.prX[[x]][i.freeX[[x]]])))  %*% rep(1,I-1)))
    
    
    
    
    
    
    log.k.prX[[x]]<- sX[[x]][i.freeX[[x]]]%*%la.prX[[x]][i.freeX[[x]]]-a*log(1+rep(1,I-1)%*%exp(tZ[,i.freeX[[x]]]%*%la.prX[[x]][i.freeX[[x]]]))
    
    
    
    # Posterior parameters estimates for saturated model 
    
    
    
    pX[[x]]<-(c.cX[[x]]+10^(-10))/n
    
    i=1
    
    S<-matrix(rep(0,(V-1)^2),V-1,V-1)
    
    
    
    
    # Prior parameters estimates for saturated model 
    
    
    
    SX[[x]]<-matrix(rep(0,(V-1)^2),V-1,V-1)
    
    
    
    SX[[x]][upper.tri(SX[[x]])] = LX[[x]][i,p.L2]
    
    SX[[x]]<-t(SX[[x]])
    
    SX[[x]][upper.tri(SX[[x]])] = LX[[x]][i,p.L2]
    
    
    
    
    FITX[[x]]<-matrix(0,V-1,V-1)
    
    fit.lmX[[x]]<-1
    fit.l2X[[x]]<-1
    
    
    for (j in 1:(V-1)) {
      
      if (sum(SX[[x]][j,])!=0) {
        
        fit<-glm(dataX[[x]][,j]~dataX[[x]][,which(SX[[x]][j,]==1)],family = binomial)
        
        FITX[[x]][j,which(SX[[x]][j,]==1)]<-fit$coefficients[-1]
        
        fit.lmX[[x]][j]<-fit$coefficients[1]
        
      } else {
        
        fit.lmX[[x]][j]<- log(pX[[x]][1+p.LM[j]]/pX[[x]][1]) 
        
      }
      
      if (j>1){
        
        fit.l2X[[x]]<-c(fit.l2X[[x]], FITX[[x]][j,(1:(j-1))])
        
      }
      
    }
    
    
    fit.l2X[[x]]<-fit.l2X[[x]][-1]
    
    
    
    
    
    
    # Posterior parameters estimates given X=0
    
    
    
    
    
    
    la.poX[[x]]<-rep(0,2^(V-1)-1)
    
    
    la.poX[[x]][p.LM]<-fit.lmX[[x]]
    
    
    la.poX[[x]][p.L2]<-fit.l2X[[x]]
    
    
    la.poX[[x]]<-la.poX[[x]][i.freeX[[x]]]
    
    
    
    
    
    
    
    
    H.poX[[x]]<- - ((a+n) / (1+t(rep(1,I-1)) %*% exp(tZ[,i.freeX[[x]]] %*% la.poX[[x]]) ))[1,1] *
      
      (t(tZ[,i.freeX[[x]]]) %*% diag(as.vector(exp(tZ[,i.freeX[[x]]]%*%la.poX[[x]]))) %*% tZ[,i.freeX[[x]]]) +
      
      (t(tZ[,i.freeX[[x]]]) %*% diag(as.vector(exp(tZ[,i.freeX[[x]]]%*%la.poX[[x]]))) %*% rep(1,I-1) %*%
         
         ((a+n) / (1+t(rep(1,I-1)) %*% exp(tZ[,i.freeX[[x]]] %*% la.poX[[x]]) )^2)[1,1] %*%
         
         t(t(tZ[,i.freeX[[x]]]) %*% diag(as.vector(exp(tZ[,i.freeX[[x]]]%*%la.poX[[x]])))  %*% rep(1,I-1)))
    
    
    
    
    
    log.k.poX[[x]]<- (sX[[x]][i.freeX[[x]]]+yX[[x]][i.freeX[[x]]])%*%la.poX[[x]]-(a+n)*log(1+rep(1,I-1)%*%exp(tZ[,i.freeX[[x]]]%*%la.poX[[x]]))
    
    
    
    
    
    # Marginal log-likelihood
    
    
    log.M.LX[[x]]<-1
    
    log.M.LX[[x]][1]<-log.k.poX[[x]]-log.k.prX[[x]]+log(sqrt(det(-H.prX[[x]])))-log(sqrt(det(-H.poX[[x]])))  
    
    l.post.p[[x]][1]<-log.M.LX[[x]][i]+log(PdeltaX[[x]][1])
    
    
  }
  
  ####################   START THE LOOP !!!!  ###############################
  
  
  for (i in 2:K){
    
    for (x in 1:X){
      
      # MODEL SETTING
      
      
      
      # We begin setting the model equal to the previous and useful objects
      
      thetaK[[i]]<-thetaK[[i-1]]
      epsilonK[[i]]<-epsilonK[[i-1]]
      nuK[i,]<-nuK[i-1,]
      
      
      LX[[x]][i,]<-LX[[x]][i-1,]
      
      LrX2[[x]][i,]<-LrX2[[x]][i-1,]
      
      
      
      # We will add or remove randomly from 0 and 1
      
      
      # IN SOSTANZA PROPONGO UN NUOVO VETTORE SELEZIONATORE IN MANIERA CASUALE SWITCHANDO UN SOLO PARAMETRO ALLA VOLTA
      
      
      q<-sample(1:length(LX[[x]][i,p.L2]),1) # adding or removing from 0
      
      if (LX[[x]][i,p.L2[q]]==1){
        
        LX[[x]][i,p.L2[q]]<-0
        
      } else {
        
        LX[[x]][i,p.L2[q]]<-1
        
      }
      
      l<-LX[[x]][i,] # vettore selezionatore completo
      
      
      
      LrX[[x]][i,]<-LX[[x]][i,-p.LZ] # vettore selezionatore di interesse
      LrX2[[x]][i,]<-LX[[x]][i,p.L2]
      
      
      
      
      
      # MARGINAL LOG-LIKELIHOOD COMPUTATION
      
      
      
      # STARTING
      
      
      
      i.freeX[[x]]<-which(LX[[x]][i,]==1) # given a model for 0 we identify the parameter that are proposed free
      
      
      
      # Prior parameters estimates for saturated model
      
      
      
      
      H.prX[[x]]<- - (a / (1+t(rep(1,I-1)) %*% exp(tZ[,i.freeX[[x]]] %*% la.prX[[x]][i.freeX[[x]]]) ))[1,1] *
        
        (t(tZ[,i.freeX[[x]]]) %*% diag(as.vector(exp(tZ[,i.freeX[[x]]]%*%la.prX[[x]][i.freeX[[x]]]))) %*% tZ[,i.freeX[[x]]]) +
        
        (t(tZ[,i.freeX[[x]]]) %*% diag(as.vector(exp(tZ[,i.freeX[[x]]]%*%la.prX[[x]][i.freeX[[x]]]))) %*% rep(1,I-1) %*%
           
           (a / (1+t(rep(1,I-1)) %*% exp(tZ[,i.freeX[[x]]] %*%la.prX[[x]][i.freeX[[x]]]) )^2)[1,1] %*%
           
           t(t(tZ[,i.freeX[[x]]]) %*% diag(as.vector(exp(tZ[,i.freeX[[x]]]%*%la.prX[[x]][i.freeX[[x]]])))  %*% rep(1,I-1)))
      
      
      log.k.prX[[x]]<- sX[[x]][i.freeX[[x]]]%*%la.prX[[x]][i.freeX[[x]]]-a*log(1+rep(1,I-1)%*%exp(tZ[,i.freeX[[x]]]%*%la.prX[[x]][i.freeX[[x]]]))
      
      
      
      # Posterior parameters estimates for saturated model 
      
      
      SX[[x]]<-matrix(rep(0,(V-1)^2),V-1,V-1)
      
      
      SX[[x]][upper.tri(SX[[x]])] = LX[[x]][i,p.L2]
      
      SX[[x]]<-t(SX[[x]])
      
      SX[[x]][upper.tri(SX[[x]])] = LX[[x]][i,p.L2]
      
      
      FITX[[x]]<-matrix(0,V-1,V-1)
      
      
      
      fit.lmX[[x]]<-1
      
      fit.l2X[[x]]<-1
      
      
      
      
      for (j in 1:(V-1)) {
        
        if (sum(SX[[x]][j,])!=0) {
          
          fit<-glm(dataX[[x]][,j]~dataX[[x]][,which(SX[[x]][j,]==1)],family = binomial)
          
          FITX[[x]][j,which(SX[[x]][j,]==1)]<-fit$coefficients[-1]
          
          fit.lmX[[x]][j]<-fit$coefficients[1]
          
        } else {
          
          fit.lmX[[x]][j]<- log(pX[[x]][1+p.LM[j]]/pX[[x]][1]) 
          
        }
        
        if (j>1){
          
          fit.l2X[[x]]<-c(fit.l2X[[x]], FITX[[x]][j,(1:(j-1))])
          
        }
        
      }
      
      
      
      fit.l2X[[x]]<-fit.l2X[[x]][-1]
      
      
      
      
      # Posterior parameters estimates given X=0
      
      
      
      
      
      
      la.poX[[x]]<-rep(0,2^(V-1)-1)
      
      
      
      la.poX[[x]][p.LM]<-fit.lmX[[x]]
      
      
      
      la.poX[[x]][p.L2]<-fit.l2X[[x]]
      
      
      
      la.poX[[x]]<-la.poX[[x]][i.freeX[[x]]]
      
      
      
      
      
      
      
      H.poX[[x]]<- - ((a+n) / (1+t(rep(1,I-1)) %*% exp(tZ[,i.freeX[[x]]] %*% la.poX[[x]]) ))[1,1] *
        
        (t(tZ[,i.freeX[[x]]]) %*% diag(as.vector(exp(tZ[,i.freeX[[x]]]%*%la.poX[[x]]))) %*% tZ[,i.freeX[[x]]]) +
        
        (t(tZ[,i.freeX[[x]]]) %*% diag(as.vector(exp(tZ[,i.freeX[[x]]]%*%la.poX[[x]]))) %*% rep(1,I-1) %*%
           
           ((a+n) / (1+t(rep(1,I-1)) %*% exp(tZ[,i.freeX[[x]]] %*% la.poX[[x]]) )^2)[1,1] %*%
           
           t(t(tZ[,i.freeX[[x]]]) %*% diag(as.vector(exp(tZ[,i.freeX[[x]]]%*%la.poX[[x]])))  %*% rep(1,I-1)))
      
      
      
      
      
      
      
      log.k.poX[[x]]<- (sX[[x]][i.freeX[[x]]]+yX[[x]][i.freeX[[x]]])%*%la.poX[[x]]-(a+n)*log(1+rep(1,I-1)%*%exp(tZ[,i.freeX[[x]]]%*%la.poX[[x]]))
      
      
      
      
      
      
      
      # Marginal log-likelihood
      
      # APPROSSIMAZIONE DI LAPLACE DELLA COSTANTE DI NORMALIZZAZIONE (MARGINAL LIKELIHOOD)
      
      log.M.LX[[x]][i]<-log.k.poX[[x]]-log.k.prX[[x]]+log(sqrt(det(-H.prX[[x]])))-log(sqrt(det(-H.poX[[x]])))   
      
      
      
      # Metropolis-Hastings step
      
      # (FIRSTLY WE SET ALL THE OBJECTS AS ACCEPTED)
      
      PdeltaX[[x]][i]<-pdelta(i,x,nuK[i,],thetaK[[i]],LrX2)
      
      
      
      ldifX[[x]][i]<-(log.M.LX[[x]][i]+log(PdeltaX[[x]][i]))-(log.M.LX[[x]][i-1]+log(PdeltaX[[x]][i-1])) 
      
      l.post.p[[x]][i]<-log.M.LX[[x]][i]+log(PdeltaX[[x]][i])
      
      
      r<-min(0,ldifX[[x]][i])
      
      if ( r < log(runif(1,0,1)) ){ # we evaluate if we reject the actual state
        
        
        LX[[x]][i,]<-LX[[x]][i-1,]
        
        log.M.LX[[x]][i]<-log.M.LX[[x]][i-1]
        
        l.post.p[[x]][i]<-l.post.p[[x]][i-1]
        
        PdeltaX[[x]][i]<-PdeltaX[[x]][i-1]
        
      } else {
        
        LX[[x]][i,]<-l
        
      }
      
      
      
      LrX[[x]][i,]<-LX[[x]][i,-p.LZ]
      LrX2[[x]][i,]<-LX[[x]][i,p.L2]
      
      
      PREDBEL[[x]][[P]][i,]<-LrX2[[x]][i,]
      
      
      
      print(c(i,1,P))
      
      
    }
    
    
    
    
    
    
    # THETA UPDATING
    
    Resp <- c(0L,1L)
    comb <- do.call(expand.grid,lapply(1:X,function(x)Resp)) # All posible states:
    comb<-as.matrix(comb)
    
    
    
    for (u in 1:X){ # coppie indici X
      
      for (h in 1:X){
        
        if (h > u){ # per ogni possibile coppia di livelli di X
          
          
          
          
          if (epsilonK[[i]][u,h]==0) {
            
            
            eps<-1
            the<-rgamma(1,alpha,chi)
            
            Eps<-epsilonK[[i]]
            The<-thetaK[[i]]
            
            Eps[u,h]<-eps
            Eps[h,u]<-eps
            
            The[u,h]<-the
            The[h,u]<-the
            
            
            
            Cost<-c()
            Costp<-c()
            
            
            lp.pr.epsilon<-log(dbinom(epsilonK[[i]][u,h],1,w))
            lp.pr.theta<-log((1-epsilonK[[i]][u,h])*de0+epsilonK[[i]][u,h]*dgamma(thetaK[[i]][u,h],shape=alpha,scale=chi))
            lq.km<-log(dgamma(thetaK[[i]][u,h],shape=alpha,scale=chi))
            
            
            lp.pr.epsilonp<-log(dbinom(eps,1,w))
            lp.pr.thetap<-log((1-eps)*de0+eps*dgamma(the,shape=alpha,scale=chi))
            lq.kmp<-log(dgamma(the,shape=alpha,scale=chi))
            
            prod.ij<- c() 
            prod.ijp<- c()  
            
            
            # produttoria sugli ij
            
            for (g in 1:length(LrX2[[u]][i,])){ # coppie indici variabili
              
              
              ker<-  exp(2*epsilonK[[i]][u,h]*LrX2[[u]][i,g]*LrX2[[h]][i,g])
              
              kerp<- exp(2*eps*LrX2[[u]][i,g]*LrX2[[h]][i,g])
              
              
              for (l in 1:(2^X)){
                
                Cost[l]<-exp(nuK[i,g]*t(rep(1,X))%*%comb[l,]+t(comb[l,])%*%thetaK[[i]]%*%comb[l,])
                
                Costp[l]<-exp(nuK[i,g]*t(rep(1,X))%*%comb[l,]+t(comb[l,])%*%The%*%comb[l,])
                
              }
              
              Cost<-sum(Cost)
              Costp<-sum(Costp)
              
              
              prod.ij[g]<-  log(ker / Cost) 
              prod.ijp[g]<- log(kerp / Costp) 
              
              
            }
            
            
            prod.ij<-sum(prod.ij)
            
            lp.post<- prod.ij + lp.pr.theta + lp.pr.epsilon
            
            prod.ijp<-sum(prod.ijp)
            
            lp.postp<- prod.ijp + lp.pr.thetap + lp.pr.epsilonp
            
            
            
            rrr<- lp.postp + lq.km - lp.post
            
            acc <- min(0,rrr)
            
            
            if ( acc > log(runif(1,0,1))){
              
              epsilonK[[i]][u,h]<-eps
              epsilonK[[i]][h,u]<-eps
              
              thetaK[[i]][u,h]<-the
              thetaK[[i]][h,u]<-the
              
            }
            
            
            
          } else { # epsilon = 1
            
            
            eps<-0
            the<-0
            
            Eps<-epsilonK[[i]]
            The<-thetaK[[i]]
            
            Eps[u,h]<-eps
            Eps[h,u]<-eps
            
            The[u,h]<-the
            The[h,u]<-the
            
            
            
            Cost<-c()
            Costp<-c()
            
            
            lp.pr.epsilon<-log(dbinom(epsilonK[[i]][u,h],1,w))
            lp.pr.theta<-log((1-epsilonK[[i]][u,h])*de0+epsilonK[[i]][u,h]*dgamma(thetaK[[i]][u,h],shape=alpha,scale=chi))
            lq.km<-log(dgamma(thetaK[[i]][u,h],shape=alpha,scale=chi))
            
            
            lp.pr.epsilonp<-log(dbinom(eps,1,w))
            lp.pr.thetap<-log((1-eps)*de0+eps*dgamma(the,shape=alpha,scale=chi))
            lq.kmp<-log(dgamma(the,shape=alpha,scale=chi))
            
            prod.ij<- c()  
            prod.ijp<- c()  
            
            
            # produttoria sugli ij
            
            for (g in 1:length(LrX2[[u]][i,])){ # coppie indici variabili
              
              
              
              ker<-  exp(2*epsilonK[[i]][u,h]*LrX2[[u]][i,g]*LrX2[[h]][i,g])
              
              kerp<- exp(2*eps*LrX2[[u]][i,g]*LrX2[[h]][i,g])
              
              
              for (l in 1:(2^X)){
                
                Cost[l]<-exp(nuK[i,g]*t(rep(1,X))%*%comb[l,]+t(comb[l,])%*%thetaK[[i]]%*%comb[l,])
                
                Costp[l]<-exp(nuK[i,g]*t(rep(1,X))%*%comb[l,]+t(comb[l,])%*%The%*%comb[l,])
                
              }
              
              Cost<-sum(Cost)
              Costp<-sum(Costp)
              
              
              prod.ij[g]<-  log(ker / Cost) 
              prod.ijp[g]<- log(kerp / Costp) 
              
              
            }
            
            
            prod.ij<-sum(prod.ij)
            
            lp.post<- prod.ij + lp.pr.theta + lp.pr.epsilon
            
            prod.ijp<-sum(prod.ijp)
            
            lp.postp<- prod.ijp + lp.pr.thetap + lp.pr.epsilonp
            
            
            
            rrr<- lp.postp + lq.km - lp.post
            
            acc <- min(0,rrr)
            
            
            if ( acc > log(runif(1,0,1))){
              
              epsilonK[[i]][u,h]<-eps
              epsilonK[[i]][h,u]<-eps
              
              thetaK[[i]][u,h]<-the
              thetaK[[i]][h,u]<-the
              
            }
            
          }
          
          
          
          if (epsilonK[[i]][u,h]==1){ #  se epsilon=1 propongo nuovo valore di theta
            
            
            the<-rgamma(1,alpha,chi)
            
            The<-thetaK[[i]]
            
            The[u,h]<-the
            The[h,u]<-the
            
            
            
            Cost<-c()
            Costp<-c()
            
            
            lp.pr.epsilon<-log(dbinom(epsilonK[[i]][u,h],1,w))
            lp.pr.theta<-log((1-epsilonK[[i]][u,h])*de0+epsilonK[[i]][u,h]*dgamma(thetaK[[i]][u,h],shape=alpha,scale=chi))
            lq.km<-log(dgamma(thetaK[[i]][u,h],shape=alpha,scale=chi))
            
            
            lp.pr.epsilonp<-log(dbinom(epsilonK[[i]][u,h],1,w))
            lp.pr.thetap<-log((1-epsilonK[[i]][u,h])*de0+epsilonK[[i]][u,h]*dgamma(the,shape=alpha,scale=chi))
            lq.kmp<-log(dgamma(the,shape=alpha,scale=chi))
            
            prod.ij<- c() 
            prod.ijp<- c()
            
            
            # produttoria sugli ij
            
            for (g in 1:length(LrX2[[u]][i,])){ # coppie indici variabili
              
              
              ker<-  exp(2*epsilonK[[i]][u,h]*LrX2[[u]][i,g]*LrX2[[h]][i,g])
              
              kerp<- exp(2*epsilonK[[i]][u,h]*LrX2[[u]][i,g]*LrX2[[h]][i,g])
              
              
              for (l in 1:(2^X)){
                
                Cost[l]<-exp(nuK[i,g]*t(rep(1,X))%*%comb[l,]+t(comb[l,])%*%thetaK[[i]]%*%comb[l,])
                
                Costp[l]<-exp(nuK[i,g]*t(rep(1,X))%*%comb[l,]+t(comb[l,])%*%The%*%comb[l,])
                
              }
              
              Cost<-sum(Cost)
              Costp<-sum(Costp)
              
              
              prod.ij[g]<-  log(ker / Cost) 
              prod.ijp[g]<- log(kerp / Costp) 
              
              
              
            }
            
            
            prod.ij<-sum(prod.ij)
            
            lp.post<- prod.ij + lp.pr.theta + lp.pr.epsilon
            
            prod.ijp<-sum(prod.ijp)
            
            lp.postp<- prod.ijp + lp.pr.thetap + lp.pr.epsilonp
            
            
            
            rrr<- lp.postp + lq.km  - lp.post - lq.kmp
            
            acc <- min(0,rrr)
            
            
            if ( acc > log(runif(1,0,1))){
              
              thetaK[[i]][u,h]<-the
              thetaK[[i]][h,u]<-the
              
            }
            
          }
          
          
          
          
        } 
      }
    } 
    
    
    
    
    ############# UPDATING NUij
    
    
    nup<- rbeta(1,aap,bbp)
    
    nup<-log(nup/(1-nup))
    
    # produttoria sugli ij
    
    for (g in 1:length(LrX2[[u]][i,])){ # coppie indici variabili
      
      deltaij<-c()
      
      for(t in 1:X){
        
        deltaij[t]<-LrX2[[t]][i,g]
        
      }
      
      
      for (l in 1:(2^X)){
        
        Cost[l]<-exp(nuK[i,g]*t(rep(1,X))%*%comb[l,]+t(comb[l,])%*%thetaK[[i]]%*%comb[l,])
        
        Costp[l]<-exp(nup*t(rep(1,X))%*%comb[l,]+t(comb[l,])%*%thetaK[[i]]%*%comb[l,])
        
      }
      
      Cost<-sum(Cost)
      Costp<-sum(Costp)
      
      
      pnu.ij<- exp(aap*nuK[i,g])/((1+exp(nuK[i,g]))^(aa+bb)) *  Cost^-1 * exp(nuK[i,g]*t(rep(1,X))%*%deltaij)
      pnu.ijp<- exp(aap*nup)/((1+exp(nup))^(aap+bbp)) * Costp^-1 * exp(nup*t(rep(1,X))%*%deltaij)
      
      q.nuij <- beta(aa,bb)^-1 * exp(aa*nuK[i,g])/((1+exp(nuK[i,g]))^(aa+bb))
      q.nuijp<- beta(aap,bbp)^-1 * exp(aap*nup)/((1+exp(nup))^(aap+bbp))
      
      
      r<-log(pnu.ijp)+log(q.nuij)-log(pnu.ij)-log(q.nuijp)
      
      acc <- min(0,r)
      
      
      if ( acc > log(runif(1,0,1))){
        
        nuK[i,g]<-nup
        nuK[i,g]<-nup
        
        
      }
      
    }
    
    
    
    
    
    
    
    
    # K cycle
    #print(thetaK[[i]])
    
    
    
  }
  
  OutThetaKBEL[[P]] <- thetaK
  OutNuKBEL[[P]] <- nuK

  
} #  P CYCLE


t2 <- Sys.time()
t2 - t1 #Time difference of 8.262904 days p=10
#Time difference of 6.190764 days p =10
