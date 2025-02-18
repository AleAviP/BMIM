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



######### METHODS COMPARISON , SIMULATION STUDY ################

# G0-G1-G2-G3

DATA<-list()

L.testXBAL<-list()

OutNuKBAL <- list()
OutThetaKBAL <- list()

for (P in 1:sims){
  
  L.testXBAL[[P]]<-list()
  
  for (k in 1:K){
    
    L.testXBAL[[P]][[k]]<-list()
    
    for (x in 1:X) {
      
      L.testXBAL[[P]][[k]][[x]]<-matrix(0,p,p)
    }
    
  }
  
}


t00 <- Sys.time()
for (P in 1:sims){ # INIZIO
  
  
  #########################################################################
  #################            G0-G1-G2-G3             ####################
  #########################################################################
  
  
  
  # G0<-matrix(0,p,p)
  # 
  # 
  # for (i in 1:(p-1)){
  #   
  #   G0[i,i+1]<-1
  #   
  # }
  # 
  # 
  # G0[lower.tri(G0)] = t(G0)[lower.tri(G0)]
  # 
  # G1<-G0
  # 
  # G2<-G0
  # 
  # G3<-G0
  # 
  
  # INIZIO DATA
  
  
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
t0-t00 #Time difference of 9.319643 mins

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

t1 <- Sys.time()

for (P in 1:sims){
  
  
  dataX<-list()
  
  dataX[[1]]<-DATA[[P]][[1]]
  dataX[[2]]<-DATA[[P]][[2]]
  dataX[[3]]<-DATA[[P]][[3]]
  dataX[[4]]<-DATA[[P]][[4]]
  
  
  
  ##### Bayesian aprox link
  
  print(paste("BAL P",P))
  
  
  ### Bayesian log-linear model selection for multiple graphs ###
  
  library(MASS)
  
  logl.r<- function(data,delta.r,L.dr,L.drc,r){ # quasi-likelihood nodo r-esimo
    L.dr<-diag(delta.r)%*%(L.dr+L.drc)
    L.drc<-diag(1-delta.r)%*%(L.dr+L.drc)
    datacn<-data
    datacn[,r]<-1
    logli.ri <- diag(data[,r])%*%(data %*% L.dr) - log(1+exp(datacn %*% L.dr))
    logli.r<-sum(logli.ri)
    return(logli.r)
  }
  
  
  grad.h.r<- function(data,delta.r,L.dr,L.drc,r){ # gradiente funzione h nodo r-esimo
    L.dr<-diag(delta.r)%*%(L.dr+L.drc)
    L.drc<-diag(1-delta.r)%*%(L.dr+L.drc) 
    grad<-rep(0,p)
    datacn<-data
    datacn[,r]<-1
    datacn1<-datacn
    datacn1[,which(delta.r==0)]<-0
    datacn1<-t(datacn1)
    data1<-data
    data1[,which(delta.r==0)]<-0
    data1<-diag(data[,r])%*%data1
    grad<- apply(data1,2,sum) - datacn1%*%(diag(as.vector(exp(datacn %*% L.dr)))%*%as.vector(( 1/ (1+exp(datacn %*% L.dr)) )))
    grad<-grad-(L.dr/rho)-(L.drc/gamma) 
    return(grad)
  }
  
  
  h.r<-function(data,rho,gamma,delta.r,L.dr,L.drc,r){ # funzione h_r
    L.dr<-diag(delta.r)%*%(L.dr+L.drc)
    L.drc<-diag(1-delta.r)%*%(L.dr+L.drc)
    h<-logl.r(data,delta.r,L.dr,L.drc,r)-(sum(L.dr^2)/(2*rho)) - (sum(L.drc^2)/(2*gamma))
    return(h)
  }
  
  
  
  K<-50000
  RES<-matrix(0,K,2)
  
  ### theta and nu start values
  
  nu<-matrix(0,p,p)
  nuK<-list()
  nuK[[1]]<-nu
  
  aa<-1
  bb<-2
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
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  pdelta.r<-function(k,x,r,nu,the,DX){
    
    
    for (j in 1:p){
      
      DXrest<-c()
      
      for (z in 1:X){
        
        DXrest[z]<-DX[[z]][[r]][k,j]
        
      }
      
      odds<-exp( DX[[x]][[r]][k,j]%*%(nu[r,j]+2*the[x,]%*%DXrest ) )
      
      prob<-odds/(1+odds)
      
      q[j]<-prob
      
    }
    
    PROB<-prod(q)
    
    return(PROB)
    
  }
  
  
  
  
  
  
  lpost.r<-function(data,omega.r,gamma,rho,delta.r,L.dr,L.drc,r){ # log quasi-posterior nodo r
    L.dr<-diag(delta.r)%*%(L.dr+L.drc)
    L.drc<-diag(1-delta.r)%*%(L.dr+L.drc)
    lpost<- log(omega.r) + log( (sqrt(gamma/rho))^sum(delta.r) ) + logl.r(data,delta.r,L.dr,L.drc,r)- ( sum(L.dr^2)/(2*rho) ) -
      ( sum(L.drc^2)/(2*gamma) )
    return(lpost)
  }
  
  
  G.r <- function(data,delta.r,L.dr,L.drc,r){ # funzione G del gradiente troncato h_r
    L.dr<-diag(delta.r)%*%(L.dr+L.drc)
    L.drc<-diag(1-delta.r)%*%(L.dr+L.drc)
    c<- 1
    G <- ( c/( max(c,sqrt(sum(grad.h.r(data,delta.r,L.dr,L.drc,r)^2))) ) ) * grad.h.r(data,delta.r,L.dr,L.drc,r) # gradiente troncato
    return(G)
  }
  
  
  
  
  
  #u<-1.5
  #q<- p^-(1+u)
  #c0<- 1
  #c1<- 1
  #gamma<-c0/max(n,p)
  #rho<-c1*sqrt(n/log(p))
  
  
  q<-0.5
  gamma<-0.5
  rho<-2
  sigma<-2 # costante step size?
  D.r<-matrix(0,K,p)
  La.dr<-matrix(0,K,p)
  La.drc<-matrix(0,K,p)
  delta.r<-rep(0,p)
  delta.r[1:(p/2)]<-1
  L.dr<-delta.r
  L.drc<-1-L.dr
  
  
  
  
  D<-list(c())
  L<-list(c())
  Lc<-list(c())
  
  
  for (r in 1:p){
    D[[r]]<-matrix(0,K,p)
    L[[r]]<-matrix(0,K,p)
    Lc[[r]]<-matrix(0,K,p)
    
  }
  
  
  
  
  
  DX<-list()
  
  for (x in 1:X){
    
    DX[[x]]<-D
    
  }
  
  for (x in 1:X){
    
    for (r in 1:p){
      
      DX[[x]][[r]][1,r]<-1
      
    }
    
  }
  
  LX<-list()
  
  for (x in 1:X){
    
    LX[[x]]<-L
    
  }
  
  LcX<-list()
  
  for (x in 1:X){
    
    LcX[[x]]<-Lc
    
  }
  
  L.trueX<-list()
  
  #for (x in 1:X){
  
  #  L.trueX[[x]]<-L.true
  
  #}
  
  L.trueX[[1]]<-G0
  L.trueX[[2]]<-G1
  L.trueX[[3]]<-G2
  L.trueX[[4]]<-G3
  #L.trueX[[5]]<-G4
  #L.trueX[[6]]<-G5
  #L.trueX[[7]]<-G6
  #L.trueX[[8]]<-G7
  #L.trueX[[9]]<-G8
  #L.trueX[[10]]<-G9
  
  
  
  AX<-list()
  RESX<-list()
  
  for (x in 1:X) {
    
    AX[[x]]<-matrix(0,p,p)
    RESX[[x]]<-matrix(0,K,2)
  } 
  
  
  
  
  OMEGA.R<-c()
  
  
  
  
  for (k in 2:K){
    
    thetaK[[k]]<-thetaK[[k-1]]
    epsilonK[[k]]<-epsilonK[[k-1]]
    nuK[[k]]<-nuK[[k-1]]
    
    
    for (r in 1:p){
      
      for (x in 1:X) {
        
        data<-dataX[[x]]   
        L.dr<-LX[[x]][[r]][k-1,]
        L.drc<-LcX[[x]][[r]][k-1,]
        delta.r<-DX[[x]][[r]][k-1,]
        omega.r<-pdelta.r(k,x,r,nuK[[k]],thetaK[[k]],DX)
        
        OMEGA.R[k]<-omega.r
        
        for (j in 1:p){ # per ogni parametro
          
          if (delta.r[j]==1){ # se il parametro ? attivo
            
            L.drj.p<-rnorm(1,L.dr[j],sigma) # estraggo un valore dalla distribuzione g
            L.drp<-L.dr
            L.drp[j]<-L.drj.p
            acc <- min(0, dnorm(L.dr[j],L.drp[j]+(sigma/2)*grad.h.r(data,delta.r,L.drp,L.drc,r)[j],sigma,log=TRUE) - dnorm(L.drp[j],L.dr[j]+(sigma/2)*grad.h.r(data,delta.r,L.dr,L.drc,r)[j],sigma, log=TRUE)  +
                         lpost.r(data,omega.r,gamma,rho,delta.r,L.drp,L.drc,r) - lpost.r(data,omega.r,gamma,rho,delta.r,L.dr,L.drc,r) ) 
            
            # + sigma/2)*grad.h.r(data,delta.r,L.drp,L.drc,r)[j]
            # +(sigma/2)*grad.h.r(data,delta.r,L.dr,L.drc,r)[j]
            
            
            if (acc > log(runif(1,0,1))){ # accettazione/rifiuto
              L.dr[j]<- L.drj.p
            }
          }
        }
        
        if (sum(delta.r)!=p){
          
          ############################# AGGIORNAMENTO PARAMETRI INATTIVI ####################################
          
          L.drcp<-mvrnorm(1,rep(0,sum(1-delta.r)),gamma*diag(sum(1-delta.r)))
          L.drc[which((1-delta.r)==1)]<-L.drcp
          
          ############################# VETTORE SELEZIONATORE ##############################################
        }
        
        for (j in 1:p){ # per ogni parametro
          if (j!=r){
            delta.rp<-delta.r
            delta.rp[j]<-1-delta.r[j] # switcho il parametro selezionatore del parametro j
            
            DX[[x]][[r]][k,]<-delta.rp
            
            omega.rp<-pdelta.r(k,x,r,nuK[[k]],thetaK[[k]],DX)
            
            DX[[x]][[r]][k,]<-delta.r
            
            
            b <- lpost.r(data,omega.rp,gamma,rho,delta.rp,L.dr,L.drc,r) - lpost.r(data,omega.r,gamma,rho,delta.r,L.dr,L.drc,r) 
            tau <- min(0,b)
            # calcolo la probabilit? di accettare lo switch
            if (tau > log(runif(1,0,1))){ # se estraggo 1 allora switcho il parametro altrimenti no
              delta.r[j]<- 1-delta.r[j]
            }
          }
        }
        # SALVATAGGIO
        DX[[x]][[r]][k,]<-delta.r
        LX[[x]][[r]][k,]<-diag(delta.r)%*%L.dr
        LcX[[x]][[r]][k,]<-diag(1-delta.r)%*%L.drc
        
      } # x cycle
      
    } # r cycle
    
    ############   
    
    for (x in 1:X) {
      
      for (i in 1:p){
        
        L.testXBAL[[P]][[k]][[x]][i,]<-DX[[x]][[i]][k,]
        
      }
      
    }  
    
    ############  
    
    
    print(paste("BAL K",k," P",P))
    
    
    
    # THETA UPDATING
    
    Resp <- c(0L,1L)
    comb <- do.call(expand.grid,lapply(1:X,function(x)Resp)) # All posible states:
    comb<-as.matrix(comb)
    
    
    
    for (u in 1:X){ # coppie indici X
      
      for (h in 1:X){
        
        if (h != u){ # per ogni possibile coppia di livelli di X
          
          
          
          
          if (epsilonK[[k]][u,h]==0) {
            
            
            eps<-1
            the<-rgamma(1,alpha,chi)
            
            Eps<-epsilonK[[k]]
            The<-thetaK[[k]]
            
            Eps[u,h]<-eps
            
            The[u,h]<-the
            
            
            
            Cost<-c()
            Costp<-c()
            
            
            lp.pr.epsilon<-log(dbinom(epsilonK[[k]][u,h],1,w))
            lp.pr.theta<-log((1-epsilonK[[k]][u,h])*de0+epsilonK[[k]][u,h]*dgamma(thetaK[[k]][u,h],shape=alpha,scale=chi))
            lq.km<-log(dgamma(thetaK[[k]][u,h],shape=alpha,scale=chi))
            
            
            lp.pr.epsilonp<-log(dbinom(eps,1,w))
            lp.pr.thetap<-log((1-eps)*de0+eps*dgamma(the,shape=alpha,scale=chi))
            lq.kmp<-log(dgamma(the,shape=alpha,scale=chi))
            
            prod.ij<- matrix(0,p,p)  
            prod.ijp<- matrix(0,p,p)  
            
            
            # produttoria sugli ij
            
            for (i in 1:p){ # coppie indici variabili
              
              for (j in 1:p){
                
                if (j != i){                      
                  
                  ker<-  exp(2*epsilonK[[k]][u,h]*L.testXBAL[[P]][[k]][[u]][i,j]*L.testXBAL[[P]][[k]][[h]][i,j])
                  
                  kerp<- exp(2*eps*L.testXBAL[[P]][[k]][[u]][i,j]*L.testXBAL[[P]][[k]][[h]][i,j])
                  
                  
                  for (l in 1:(2^X)){
                    
                    Cost[l]<-exp(nuK[[k]][i,j]*t(rep(1,X))%*%comb[l,]+t(comb[l,])%*%thetaK[[k]]%*%comb[l,])
                    
                    Costp[l]<-exp(nuK[[k]][i,j]*t(rep(1,X))%*%comb[l,]+t(comb[l,])%*%The%*%comb[l,])
                    
                  }
                  
                  Cost<-sum(Cost)
                  Costp<-sum(Costp)
                  
                  
                  prod.ij[i,j]<-  log(ker / Cost) 
                  prod.ijp[i,j]<- log(kerp / Costp) 
                  
                  
                }
                
              }
              
            }
            
            
            prod.ij<-sum(prod.ij)
            
            lp.post<- prod.ij + lp.pr.theta + lp.pr.epsilon
            
            prod.ijp<-sum(prod.ijp)
            
            lp.postp<- prod.ijp + lp.pr.thetap + lp.pr.epsilonp
            
            
            
            rrr<- lp.postp + lq.km - lp.post
            
            acc <- min(0,rrr)
            
            
            if ( acc > log(runif(1,0,1))){
              
              epsilonK[[k]][u,h]<-eps
              
              thetaK[[k]][u,h]<-the
              
            }
            
            
            
          } else { # epsilon = 1
            
            
            eps<-0
            the<-0
            
            Eps<-epsilonK[[k]]
            The<-thetaK[[k]]
            
            Eps[u,h]<-eps
            
            The[u,h]<-the
            
            
            
            Cost<-c()
            Costp<-c()
            
            
            lp.pr.epsilon<-log(dbinom(epsilonK[[k]][u,h],1,w))
            lp.pr.theta<-log((1-epsilonK[[k]][u,h])*de0+epsilonK[[k]][u,h]*dgamma(thetaK[[k]][u,h],shape=alpha,scale=chi))
            lq.km<-log(dgamma(thetaK[[k]][u,h],shape=alpha,scale=chi))
            
            
            lp.pr.epsilonp<-log(dbinom(eps,1,w))
            lp.pr.thetap<-log((1-eps)*de0+eps*dgamma(the,shape=alpha,scale=chi))
            lq.kmp<-log(dgamma(the,shape=alpha,scale=chi))
            
            prod.ij<- matrix(0,p,p)  
            prod.ijp<- matrix(0,p,p)  
            
            
            # produttoria sugli ij
            
            for (i in 1:p){ # coppie indici variabili
              
              for (j in 1:p){
                
                if (j != i){
                  
                  ker<-  exp(2*epsilonK[[k]][u,h]*L.testXBAL[[P]][[k]][[u]][i,j]*L.testXBAL[[P]][[k]][[h]][i,j])
                  
                  kerp<- exp(2*eps*L.testXBAL[[P]][[k]][[u]][i,j]*L.testXBAL[[P]][[k]][[h]][i,j])
                  
                  
                  for (l in 1:(2^X)){
                    
                    Cost[l]<-exp(nuK[[k]][i,j]*t(rep(1,X))%*%comb[l,]+t(comb[l,])%*%thetaK[[k]]%*%comb[l,])
                    
                    Costp[l]<-exp(nuK[[k]][i,j]*t(rep(1,X))%*%comb[l,]+t(comb[l,])%*%The%*%comb[l,])
                    
                  }
                  
                  Cost<-sum(Cost)
                  Costp<-sum(Costp)
                  
                  
                  prod.ij[i,j]<-  log(ker / Cost) 
                  prod.ijp[i,j]<- log(kerp / Costp) 
                  
                  
                }
                
              }
              
            }
            
            
            prod.ij<-sum(prod.ij)
            
            lp.post<- prod.ij + lp.pr.theta + lp.pr.epsilon
            
            prod.ijp<-sum(prod.ijp)
            
            lp.postp<- prod.ijp + lp.pr.thetap + lp.pr.epsilonp
            
            
            
            rrr<- lp.postp + lq.km - lp.post
            
            acc <- min(0,rrr)
            
            
            if ( acc > log(runif(1,0,1))){
              
              epsilonK[[k]][u,h]<-eps
              
              thetaK[[k]][u,h]<-the
              
            }
            
          }
          
          
          
          if (epsilonK[[k]][u,h]==1){ #  se epsilon=1 propongo nuovo valore di theta
            
            
            the<-rgamma(1,alpha,chi)
            
            The<-thetaK[[k]]
            
            The[u,h]<-the
            
            
            
            Cost<-c()
            Costp<-c()
            
            
            lp.pr.epsilon<-log(dbinom(epsilonK[[k]][u,h],1,w))
            lp.pr.theta<-log((1-epsilonK[[k]][u,h])*de0+epsilonK[[k]][u,h]*dgamma(thetaK[[k]][u,h],shape=alpha,scale=chi))
            lq.km<-log(dgamma(thetaK[[k]][u,h],shape=alpha,scale=chi))
            
            
            lp.pr.epsilonp<-log(dbinom(epsilonK[[k]][u,h],1,w))
            lp.pr.thetap<-log((1-epsilonK[[k]][u,h])*de0+epsilonK[[k]][u,h]*dgamma(the,shape=alpha,scale=chi))
            lq.kmp<-log(dgamma(the,shape=alpha,scale=chi))
            
            prod.ij<- matrix(0,p,p)  
            prod.ijp<- matrix(0,p,p)  
            
            
            # produttoria sugli ij
            
            for (i in 1:p){ # coppie indici variabili
              
              for (j in 1:p){
                
                if (j != i){
                  
                  ker<-  exp(2*epsilonK[[k]][u,h]*L.testXBAL[[P]][[k]][[u]][i,j]*L.testXBAL[[P]][[k]][[h]][i,j])
                  
                  kerp<- exp(2*epsilonK[[k]][u,h]*L.testXBAL[[P]][[k]][[u]][i,j]*L.testXBAL[[P]][[k]][[h]][i,j])
                  
                  
                  for (l in 1:(2^X)){
                    
                    Cost[l]<-exp(nuK[[k]][i,j]*t(rep(1,X))%*%comb[l,]+t(comb[l,])%*%thetaK[[k]]%*%comb[l,])
                    
                    Costp[l]<-exp(nuK[[k]][i,j]*t(rep(1,X))%*%comb[l,]+t(comb[l,])%*%The%*%comb[l,])
                    
                  }
                  
                  Cost<-sum(Cost)
                  Costp<-sum(Costp)
                  
                  
                  prod.ij[i,j]<-  log(ker / Cost) 
                  prod.ijp[i,j]<- log(kerp / Costp) 
                  
                  
                }
                
              }
              
            }
            
            
            prod.ij<-sum(prod.ij)
            
            lp.post<- prod.ij + lp.pr.theta + lp.pr.epsilon
            
            prod.ijp<-sum(prod.ijp)
            
            lp.postp<- prod.ijp + lp.pr.thetap + lp.pr.epsilonp
            
            
            
            rrr<- lp.postp + lq.km  - lp.post - lq.kmp
            
            acc <- min(0,rrr)
            
            
            if ( acc > log(runif(1,0,1))){
              
              thetaK[[k]][u,h]<-the
            }
            
          }
          
          
          
          
        } 
      }
    } 
    
    
    
    
    
    ############# UPDATING NUij
    
    
    nup<- rbeta(1,aap,bbp)
    
    nup<-log(nup/(1-nup))
    
    # produttoria sugli ij
    
    for (i in 1:p){ # coppie indici variabili
      
      for (j in 1:p){
        
        if (j != i){
          
          deltaij<-c()
          
          for(t in 1:X){
            
            deltaij[t]<-L.testXBAL[[P]][[k]][[t]][i,j]
            
          }
          
          
          for (l in 1:(2^X)){
            
            Cost[l]<-exp(nuK[[k]][i,j]*t(rep(1,X))%*%comb[l,]+t(comb[l,])%*%thetaK[[k]]%*%comb[l,])
            
            Costp[l]<-exp(nup*t(rep(1,X))%*%comb[l,]+t(comb[l,])%*%thetaK[[k]]%*%comb[l,])
            
          }
          
          Cost<-sum(Cost)
          Costp<-sum(Costp)
          
          
          pnu.ij<- exp(aap*nuK[[k]][i,j])/((1+exp(nuK[[k]][i,j]))^(aa+bb)) *  Cost^-1 * exp(nuK[[k]][i,j]*t(rep(1,X))%*%deltaij)
          pnu.ijp<- exp(aap*nup)/((1+exp(nup))^(aap+bbp)) * Costp^-1 * exp(nup*t(rep(1,X))%*%deltaij)
          
          q.nuij <- beta(aa,bb)^-1 * exp(aa*nuK[[k]][i,j])/((1+exp(nuK[[k]][i,j]))^(aa+bb))
          q.nuijp<- beta(aap,bbp)^-1 * exp(aap*nup)/((1+exp(nup))^(aap+bbp))
          
          
          r<-log(pnu.ijp)+log(q.nuij)-log(pnu.ij)-log(q.nuijp)
          
          acc <- min(0,r)
          
          
          if ( acc > log(runif(1,0,1))){
            
            nuK[[k]][i,j]<-nup
            
          }
        }
        
      }
      
    }
    
    
    
  }# K cycle
  
  OutThetaKBAL[[P]] <- thetaK
  OutNuKBAL[[P]] <- nuK
  
  
} #  P CYCLE

t2 <- Sys.time()
t2 - t1 #Time difference of 4.187562 days p =10