rm(list=ls())

library(pROC)
library(igraph)
library(ggm)
library(glmnet)

load("C:/Users/aleav/Dropbox/Bayesian Graphical Models/BMIM/phd-Andrea/Data/appl18paper_november.RData")
#load("~/Dropbox/Bayesian Graphical Models/BMIM/phd-Andrea/Data/appl10paper_november.RData")

# APPLICATION


#DATA<-list()
set.seed(123)

RSL<-matrix(0,18,1)
RDSSL<-matrix(0,18,1)
RBA<-matrix(0,18,1)
RBAL<-matrix(0,18,1)
RBE<-matrix(0,18,1)

RBEL<-matrix(0,18,1)

predSL<-list()
predDSSL<-list()
predBA<-list()
predBAL<-list()
predBE<-list()

predBEL<-list()



AUCSL<-c()
AUCDSSL<-c()
AUCBA<-c()
AUCBAL<-c()
AUCBE<-c()

AUCBEL<-c()

X<-3
p<-18
n<-c(263,250,255)
K<-5000

for (x in 1:X){
  
  predSL[[x]]<-matrix(0,18,153)
  predDSSL[[x]]<-matrix(0,18,153)
  predBA[[x]]<-matrix(0,18,153)
  predBAL[[x]]<-matrix(0,18,153)
  predBE[[x]]<-matrix(0,18,153)
  predBEL[[x]]<-matrix(0,18,153)
  
  
}



D<-matrix(0,5000,153)

ciccio<-list()

PREDBAL<-list()
PREDBA<-list()
PREDBE<-list()
PREDBEL<-list()


for (t in 1:1){
  
  ciccio[[t]]<-D
  
}

for (x in 1:X){
  
  PREDBAL[[x]]<-ciccio
  
  PREDBA[[x]]<-ciccio
  
  PREDBE[[x]]<-ciccio
  
  PREDBEL[[x]]<-ciccio
  
  
  
}


K<-5000

L.testXBAL<-list()

for (P in 1:1){
  
  L.testXBAL[[P]]<-list()
  
  for (k in 1:K){
    
    L.testXBAL[[P]][[k]]<-list()
    
    for (x in 1:X) {
      
      L.testXBAL[[P]][[k]][[x]]<-matrix(0,p,p)
    }
    
  }
  
}




L.testXBA<-list()

for (P in 1:1){
  
  L.testXBA[[P]]<-list()
  
  for (k in 1:K){
    
    L.testXBA[[P]][[k]]<-list()
    
    for (x in 1:X) {
      
      L.testXBA[[P]][[k]][[x]]<-matrix(0,p,p)
    }
    
  }
  
}

P<-1

dataX<-list()

dataX[[1]] <-as.matrix(datax[which(datax$profile=="0"),-19])
dataX[[2]] <-as.matrix(datax[which(datax$profile=="1"),-19])
dataX[[3]] <-as.matrix(datax[which(datax$profile=="2"),-19])




# dataX[[1]]<-DATA[[P]][[1]]
# dataX[[2]]<-DATA[[P]][[2]]
# dataX[[3]]<-DATA[[P]][[3]]
# dataX[[4]]<-DATA[[P]][[4]]




######FINE REGISTRAZIONE

#for (P in 1:1){


# dataX<-list()
# 
# dataX[[1]]<-DATA[[P]][[1]]
# dataX[[2]]<-DATA[[P]][[2]]
# dataX[[3]]<-DATA[[P]][[3]]
# dataX[[4]]<-DATA[[P]][[4]]






library(igraph)

library(ggm)

library(glmnet)






########################################## start sl *************************************


#### SEPlogit ####

print("SL"); print(P)


err<-matrix(0,1,2)

ppar<-list()

for(x in 1:X){
  
  ppar[[x]]<-matrix(0,p,p)
  
  a<-1
  
  for (i in 1:p) {
    
    
    
    # Perform 10-fold cross-validation to select lambda ---------------------------
    lambdas_to_try <- seq(0.01,0.5, length.out = 1000)
    
    
    # Use information criteria to select lambda -----------------------------------
    aic <- c()
    bic <- c()
    aic2<- c()
    bic2<- c()
    for (lambda in seq(lambdas_to_try)) {
      # Run model
      mod<-glmnet(dataX[[x]][,-i],dataX[[x]][,i],lambda = lambdas_to_try[lambda],family = "binomial")
      
      # Extract coefficients and residuals (remove first row for the intercept)
      #betas <- as.vector((as.matrix(coef(mod))[-1, ]))
      #resid <- dataX[[x]][,i] - (dataX[[x]][,-i] %*% betas)
      
      # Compute information criteria
      #aic[lambda] <- nrow(dataX[[x]][,-i]) * log((t(resid) %*% resid)/nrow(dataX[[x]][,-i])) + 2 * (sum(betas!=0))
      #bic[lambda] <- nrow(dataX[[x]][,-i]) * log((t(resid) %*% resid)/nrow(dataX[[x]][,-i])) + log(nrow(dataX[[x]][,-i])) * (sum(betas!=0))
      
      
      tLL <- - deviance(mod)
      kk <- mod$df
      nn <- mod$nobs
      aic2[lambda] <- -tLL+2*kk
      bic2[lambda]<- -tLL+log(nn)*kk
      
      
    }
    
    
    # Optimal lambdas according to both criteria
    lambda_aic <- lambdas_to_try[which.min(aic)]
    lambda_bic <- lambdas_to_try[which.min(bic)]
    lambda_aic2 <- lambdas_to_try[which.min(aic2)]
    lambda_bic2 <- lambdas_to_try[which.min(bic2)]
    
    mod<-glmnet(dataX[[x]][,-i],dataX[[x]][,i],lambda = lambda_bic2,family = "binomial")
    
    par<-coef(mod,lambda_bic2)
    
    ppar[[x]][i,i]<-par[1]
    
    par<-as.vector(par[-1])
    
    indic<-0
    
    for (z in 1:(p-1)){
      
      if (ppar[[x]][i,z]==0 & indic==0){
        
        ppar[[x]][i,z]<-par[z]
        
      } else {
        
        ppar[[x]][i,z+1]<-par[z]
        
        indic<-1
      }
      
    }
    
    for (j in 1:p){
      
      if (ppar[[x]][i,j]!=0){
        
        ppar[[x]][i,j]<-1
        
      }
      
      if (i==j){
        
        ppar[[x]][i,j]<-0
        
      }
      
    }
    
  }
  
  for (i in 1:p){
    for (j in 1:p){
      if (j>i){
        
        if (ppar[[x]][i,j]+ppar[[x]][j,i]<2){
          
          ppar[[x]][i,j]<-0
          ppar[[x]][j,i]<-0
          
          
          
        }
        
      }
      
    }
    
  }
  
}

#predSL[[x]][P,]<-as.vector(t(ppar[upper.tri(ppar)]))
# ER<-ppar-GX[[x]]
# 
# FP<-(length(which(ER==1))/2)
# FN<-(length(which(ER==-1))/2)
# 
# err<-rbind(err,c(FN,FP))
# 

ppar2<-ppar[[1]]+ppar[[2]]+ppar[[3]]

ppar3<-ppar

for (x in 1:X){
  
  for (i in 1:p){
    
    for (j in 1:p){
      
      if (ppar2[i,j]==3 & i!=j){
        
        ppar3[[x]][i,j]<-10
        
      } else if (ppar[[x]][i,j]!=0 & ppar2[i,j]!=3 & i!=j) {
        
        ppar3[[x]][i,j]<-100
        
        
      }
      
    }
    
  }
  
  rownames(ppar3[[x]])<-c("welf","high","sec","trans","park","chi","sci","ern","for",
                          "mil","bla","spa","env","hea","cit","crime","drug","edu")
  colnames(ppar3[[x]])<-c("welf","high","sec","trans","park","chi","sci","ern","for",
                          "mil","bla","spa","env","hea","cit","crime","drug","edu")
  
  graph<-graph_from_adjacency_matrix(ppar3[[x]])
  
  plotGraph(graph,layout = layout_in_circle,vc= "white",nodesize=22,dashed=TRUE,coloth="red",colbid ="black",cex=1)
  
}

############################################### end sl


# err<-err[-1,]
# 
# 
# 
# 
# 
# RSL[P,]<-apply(err,2,mean)
# 
# aucSL<-c()
# 
# for (x in 1:X){
#   
#   aucc<-roc(predSL[[x]][P,],as.vector(t(GX[[x]][upper.tri(GX[[x]])])))
#   
#   aucSL[x]<-aucc$auc
#   
# }
# 
# AUCSL[P]<-mean(aucSL)

#####################

#### DataShared-SEPlogit ####

############################################# start dssl ****************************************



print("DSSL"); print(P)

X<-3

p<-18


RES<-matrix(0,1,2)


par<-matrix(0,(p-1)*(X+1),p)

for(i in 1:p){
  
  Y<-p-1
  
  for (x in 1:X){
    
    Y<-c(Y,dataX[[x]][,i])
    
  }
  
  Y<-Y[-1] # vettore risposta ordinato per profilo
  
  
  
  a<-1
  
  
  XX1<-rep(0,p-1)
  
  
  datastX<-list()
  
  datastX[[1]]<-dataX[[1]]
  datastX[[2]]<-dataX[[2]]
  datastX[[3]]<-dataX[[3]]
  
  
  
  for(x in 1:X){
    
    for(r in 1:p){
      
      datastX[[x]][,r]<-scale(dataX[[x]][,r]) # dataset scalato
      
    }
    
  }
  
  
  
  for (v in 1:X){
    
    data2<-datastX[[v]][,-i]
    
    XX1<-rbind(XX1,data2)
    
  }
  
  XX1<-XX1[-1,]
  
  
  # XX2<-list()
  # 
  # for (v in 1:X){
  #   
  #   data2<-datastX[[v]][,-i]
  #   
  #   d<-rep(0,X)
  #   
  #   d[v]<-1
  #   
  #   XX2[[v]]<-kronecker(diag(d),as.matrix(data2))
  #   
  # }
  # 
  # XXb<-matrix(0,sum(n),(p-1)*X)
  # 
  # for (v in 1:X){
  #   
  #   XXb<-XXb+XX2[[v]]
  #   
  # }
  
  XXX1<-matrix(0,sum(n),p-1)
  
  XXX2<-matrix(0,sum(n),p-1)
  
  XXX3<-matrix(0,sum(n),p-1)
  
  
  XXX1[1:n[1],]<-datastX[[1]][,-i]
  
  XXX2[(n[1]+1):(n[1]+n[2]),]<-datastX[[2]][,-i]
  
  XXX3[(n[1]+n[2]+1):(sum(n)),]<-datastX[[3]][,-i]
  
  
  XXX<-cbind(XXX1,XXX2,XXX3)
  
  dim(XXX)
  
  XXb<-XXX
  
  XXb<-XXb*(1/sqrt(X))
  
  XX<-cbind(XX1,XXb)
  
  
  
  
  # Perform 10-fold cross-validation to select lambda ---------------------------
  lambdas_to_try <- seq(0.01,0.5, length.out = 1000)
  
  
  # Use information criteria to select lambda -----------------------------------
  aic <- c()
  bic <- c()
  aic2<- c()
  bic2<- c()
  for (lambda in seq(lambdas_to_try)) {
    # Run model
    mod<-glmnet(XX,Y,lambda = lambdas_to_try[lambda],family = "binomial")
    
    
    # Extract coefficients and residuals (remove first row for the intercept)
    #betas <- as.vector((as.matrix(coef(mod))[-1, ]))
    #resid <- Y - (XX %*% betas)
    # Compute information criteria
    #aic[lambda] <- nrow(XX) * log((t(resid) %*% resid)/nrow(XX)) + 2 * (sum(betas!=0))
    #bic[lambda] <- nrow(XX) * log((t(resid) %*% resid)/nrow(XX)) + log(nrow(XX)) * (sum(betas!=0))
    
    tLL <- - deviance(mod)
    kk <- mod$df
    nn <- mod$nobs
    aic2[lambda] <- -tLL+2*kk
    bic2[lambda]<- -tLL+log(nn)*kk
    
  }
  
  
  # Optimal lambdas according to both criteria
  #lambda_aic <- lambdas_to_try[which.min(aic)]
  #lambda_bic <- lambdas_to_try[which.min(bic)]
  
  lambda_aic2 <- lambdas_to_try[which.min(aic2)]
  lambda_bic2 <- lambdas_to_try[which.min(bic2)]
  
  
  mod<-glmnet(XX,Y,lambda = lambda_bic2,family = "binomial")
  
  par[,i]<-coef(mod,lambda_bic2)[2:((p-1)*(X+1)+1)]
  
  
  
}



TEST<-list()

for (d in 1:X){
  
  par2<-par
  
  PAR<-list()
  
  for (x in 0:X){
    
    PAR[[(x+1)]]<-par2[(1+(p-1)*x):((p-1)+(p-1)*x),]
    
  }
  
  PAR2<-list()
  
  
  for (x in 1:X){
    
    PAR2[[(x+1)]]<-PAR[[1]]+ PAR[[(x+1)]]
    
  }
  
  
  mat<-matrix(0,p,p)
  
  
  
  for (i in 1:p) {
    
    for (j in 1:p) {
      
      
      if (j<i){
        
        mat[j,i]<-PAR2[[d+1]][j,i]
        
      } else if (j==i){
        
        mat[j,i]<-0
        
        
      } else if (j>i){
        
        mat[j,i]<-PAR2[[d+1]][j-1,i]
        
        
      }
      
      
    }
  }
  
  
  
  mat2<-matrix(0,p,p)
  
  
  for (i in 1:p) {
    
    for (j in 1:(p-1)) {
      
      
      if (mat[j,i]!=0 & mat[i,j]!=0){
        
        mat2[j,i]<-1
        mat2[i,j]<-1
        
        
      } else {
        
        mat2[j,i]<-0
        mat2[i,j]<-0
        
      }
      
      
    }
  }
  
  TEST[[d]]<-mat2
  
}  




ppar2<-TEST[[1]]+TEST[[2]]+TEST[[3]]

ppar3<-TEST

for (x in 1:X){
  
  for (i in 1:p){
    
    for (j in 1:p){
      
      if (ppar2[i,j]==3 & i!=j){
        
        ppar3[[x]][i,j]<-10
        
      } else if (TEST[[x]][i,j]!=0 & ppar2[i,j]!=3 & i!=j) {
        
        ppar3[[x]][i,j]<-100
        
        
      }
      
    }
    
  }
  
  rownames(ppar3[[x]])<-c("welf","high","sec","trans","park","chi","sci","ern","for",
                          "mil","bla","spa","env","hea","cit","crime","drug","edu")
  colnames(ppar3[[x]])<-c("welf","high","sec","trans","park","chi","sci","ern","for",
                          "mil","bla","spa","env","hea","cit","crime","drug","edu")
  
  
  graph<-graph_from_adjacency_matrix(ppar3[[x]])
  
  plotGraph(graph,layout = layout_in_circle,vc= "white",nodesize=22,dashed=TRUE,coloth="red",colbid ="black",cex=1)
  
}




# library(igraph)
# 
# rownames(TEST[[1]])<-c("tv","press","exec","edu","rel","bank","court","sci","milit","congr")
# colnames(TEST[[1]])<-c("tv","press","exec","edu","rel","bank","court","sci","milit","congr")
# rownames(TEST[[2]])<-c("tv","press","exec","edu","rel","bank","court","sci","milit","congr")
# colnames(TEST[[2]])<-c("tv","press","exec","edu","rel","bank","court","sci","milit","congr")
# rownames(TEST[[3]])<-c("tv","press","exec","edu","rel","bank","court","sci","milit","congr")
# colnames(TEST[[3]])<-c("tv","press","exec","edu","rel","bank","court","sci","milit","congr")
# 
# graph1<-graph_from_adjacency_matrix(TEST[[1]],mode="undirected")
# graph2<-graph_from_adjacency_matrix(TEST[[2]],mode="undirected")
# graph3<-graph_from_adjacency_matrix(TEST[[3]],mode="undirected")
# 
# 
# plotGraph(graph1,layout = layout_in_circle,vc= "white",nodesize=30)
# plotGraph(graph2,layout = layout_in_circle,vc= "white",nodesize=30)
# plotGraph(graph3,layout = layout_in_circle,vc= "white",nodesize=30)
# 



################################### end dssl



# AX<-list()
# ER<-matrix(0,X,2)
# 
# for (x in 1:X) {
#   
#   AX[[x]]<-TEST[[x]]-GX[[x]]
#   
#   predDSSL[[x]][P,]<-as.vector(t(TEST[[x]][upper.tri(TEST[[x]])]))
#   
#   
#   ER[x,1]<-sum(AX[[x]]==-1)/2
#   ER[[x,2]]<-(sum(AX[[x]]==1))/2
# } 
# 
# 
# 
# RES<-apply(ER,2,mean)
# 
# 
# 
# 
# aucDSSL<-c()
# 
# for (x in 1:X){
#   
#   aucc<-roc(predDSSL[[x]][P,],as.vector(t(GX[[x]][upper.tri(GX[[x]])])))
#   
#   aucDSSL[x]<-aucc$auc
#   
# }
# 
# AUCDSSL[P]<-mean(aucDSSL)
# 
# 
# 
# RDSSL[P,]<- RES



####################################################





###################################################### start bal *********************************



##### BAL Bayesian Approximate-Likelihood with linking prior

print("BAL"); print(P)


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



K<-5000
RES<-matrix(0,K,2)

### theta and nu start values

nu<-matrix(0,p,p)
nuK<-list()
nuK[[1]]<-nu

aa<-1
bb<-3
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
      
      data<-as.matrix(dataX[[x]])   
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
  
  
  print(c(k,1,P))
  
  
  
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

sumbal<-list()

for(x in 1:X){
  
  sumbal[[x]]<-matrix(0,p,p)
  
  
}



for(x in 1:X){
  
  for(k in 4001:5000){
    
    sumbal[[x]]<-sumbal[[x]]+L.testXBAL[[P]][[k]][[x]]
    
    
    
  }
  
}


for(x in 1:X){
  
  
  sumbal[[x]]<-sumbal[[x]]/1000
  
  
  
}




mod<-list()


for(x in 1:X){
  
  mod[[x]]<-(sumbal[[x]]+t(sumbal[[x]]))/2
  
  mod[[x]]<-round(mod[[x]])
  
  diag(mod[[x]])<-0
  
  # rr<-roc(mod[[x]],MX.test[[x]])
  
  # aucBAL[x]<-rr$auc
  
}

library(igraph)

ppar2<-mod[[1]]+mod[[2]]+mod[[3]]

ppar3<-mod

for (x in 1:X){
  
  for (i in 1:p){
    
    for (j in 1:p){
      
      if (ppar2[i,j]==3 & i!=j){
        
        ppar3[[x]][i,j]<-10
        
      } else if (mod[[x]][i,j]!=0 & ppar2[i,j]!=3 & i!=j) {
        
        ppar3[[x]][i,j]<-100
        
        
      }
      
    }
    
  }
  
  rownames(ppar3[[x]])<-c("welf","high","sec","trans","park","chi","sci","ern","for",
                          "mil","bla","spa","env","hea","cit","crime","drug","edu")
  colnames(ppar3[[x]])<-c("welf","high","sec","trans","park","chi","sci","ern","for",
                          "mil","bla","spa","env","hea","cit","crime","drug","edu")
  
  
  graph<-graph_from_adjacency_matrix(ppar3[[x]])
  
  plotGraph(graph,layout = layout_in_circle,vc= "white",nodesize=22,dashed=TRUE,coloth="red",colbid ="black",cex=1)
  
}


# shared edges mat


sh_mat<-matrix(0,3,3)

for (i in 1:X){
  
  for (j in 1:X){
    
    if (i==j){
      
      sh_mat[i,j]<-(length(which(ppar3[[i]]!=0)))/2
      
      
    } else if (i<j) {
      
      
      sh_mat[i,j]<-(length(which(ppar3[[i]]+ppar3[[j]]==20 | ppar3[[i]]+ppar3[[j]]==200)))/2
      
      
      
    }
    
    
    
  } 
  
}

sh_mat.18.BAL<-sh_mat

# PPI

ssssum<-matrix(0,3,3)

meanssssum<-matrix(0,3,3)

for (k in 4001:5000){
  
  
  ssssum1<-matrix(0,3,3)
  
  
  ssssum1[which(thetaK[[k]]!=0)]<-1
  
  meanssssum<-meanssssum+thetaK[[k]]
  
  ssssum<-ssssum+ssssum1
  
}

ssssum<-ssssum/1000

meanssssum<-meanssssum/1000

theta.BAL.18<-thetaK

meanssssum.18.BAL<-meanssssum

round(ssssum,3)


# rownames(mod[[1]])<-c("tv","press","exec","edu","rel","bank","court","sci","milit","congr")
# colnames(mod[[1]])<-c("tv","press","exec","edu","rel","bank","court","sci","milit","congr")
# rownames(mod[[2]])<-c("tv","press","exec","edu","rel","bank","court","sci","milit","congr")
# colnames(mod[[2]])<-c("tv","press","exec","edu","rel","bank","court","sci","milit","congr")
# rownames(mod[[3]])<-c("tv","press","exec","edu","rel","bank","court","sci","milit","congr")
# colnames(mod[[3]])<-c("tv","press","exec","edu","rel","bank","court","sci","milit","congr")
# 
# graph1<-graph_from_adjacency_matrix(mod[[1]],mode="undirected")
# graph2<-graph_from_adjacency_matrix(mod[[2]],mode="undirected")
# graph3<-graph_from_adjacency_matrix(mod[[3]],mode="undirected")
# 
# 
# plotGraph(mod[[1]],layout = layout_in_circle,vc= "white",nodesize=30)
# plotGraph(mod[[2]],layout = layout_in_circle,vc= "white",nodesize=30)
# plotGraph(mod[[3]],layout = layout_in_circle,vc= "white",nodesize=30)
# 




######################################### end bal 



######### BA BAyesian Exact-likelihood without linking prior


print("BA"); print(P)


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


lpost.r<-function(data,gamma,rho,delta.r,L.dr,L.drc,r){ # log quasi-posterior nodo r
  L.dr<-diag(delta.r)%*%(L.dr+L.drc)
  L.drc<-diag(1-delta.r)%*%(L.dr+L.drc)
  omega.r<-(q^(sum(delta.r)))*((1-q)^(p-sum(delta.r))) # calcolo la probabilit? del vettore selezionatore generato
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

p<-18
q<-0.35
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


for (k in 2:K){
  
  for (r in 1:p){
    
    for (x in 1:X) {
      
      data<-as.matrix(dataX[[x]])   
      L.dr<-LX[[x]][[r]][k-1,]
      L.drc<-LcX[[x]][[r]][k-1,]
      delta.r<-DX[[x]][[r]][k-1,]
      
      for (j in 1:p){ # per ogni parametro
        
        if (delta.r[j]==1){ # se il parametro ? attivo
          
          L.drj.p<-rnorm(1,L.dr[j],sigma) # estraggo un valore dalla distribuzione g
          L.drp<-L.dr
          L.drp[j]<-L.drj.p
          acc <- min(0, dnorm(L.dr[j],L.drp[j]+(sigma/2)*grad.h.r(data,delta.r,L.drp,L.drc,r)[j],sigma,log=TRUE) - dnorm(L.drp[j],L.dr[j]+(sigma/2)*grad.h.r(data,delta.r,L.dr,L.drc,r)[j],sigma, log=TRUE)  +
                       lpost.r(data,gamma,rho,delta.r,L.drp,L.drc,r) - lpost.r(data,gamma,rho,delta.r,L.dr,L.drc,r) ) 
          
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
      
      
      
      ############################# VETTORE SELEZIONATORE ##############################################
      for (j in 1:p){ # per ogni parametro
        if (j!=r){
          delta.rp<-delta.r
          delta.rp[j]<-1-delta.r[j] # switcho il parametro selezionatore del parametro j
          b<- log( (q/(1-q))^(sum(delta.rp)) ) + log( (gamma/rho)^( (sum(delta.rp))/2 ) ) + h.r(data,rho,gamma,delta.rp,L.dr,L.drc,r) -
            log( (q/(1-q))^(sum(delta.r))  ) - log( (gamma/rho)^( (sum(delta.r))/2 ) )  - h.r(data,rho,gamma,delta.r,L.dr,L.drc,r)
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
      
      L.testXBA[[P]][[k]][[x]][i,]<-DX[[x]][[i]][k,]
      
    }
    
  } 
  
  ############  
  
  
  print(c(k,1,P))
  
} # k cycle




sumba<-list()

for(x in 1:X){
  
  sumba[[x]]<-matrix(0,p,p)
  
  
}



for(x in 1:X){
  
  for(k in 4001:5000){
    
    sumba[[x]]<-sumba[[x]]+L.testXBA[[P]][[k]][[x]]
    
    
    
  }
  
}


for(x in 1:X){
  
  
  sumba[[x]]<-sumba[[x]]/1000
  
  
  
}




mod<-list()


for(x in 1:X){
  
  mod[[x]]<-(sumba[[x]]+t(sumba[[x]]))/2
  
  mod[[x]]<-round(mod[[x]])
  
  diag(mod[[x]])<-0
  
  # rr<-roc(mod[[x]],MX.test[[x]])
  
  # aucBAL[x]<-rr$auc
  
}

library(igraph)

# rownames(mod[[1]])<-c("tv","press","exec","edu","rel","bank","court","sci","milit","congr")
# colnames(mod[[1]])<-c("tv","press","exec","edu","rel","bank","court","sci","milit","congr")
# rownames(mod[[2]])<-c("tv","press","exec","edu","rel","bank","court","sci","milit","congr")
# colnames(mod[[2]])<-c("tv","press","exec","edu","rel","bank","court","sci","milit","congr")
# rownames(mod[[3]])<-c("tv","press","exec","edu","rel","bank","court","sci","milit","congr")
# colnames(mod[[3]])<-c("tv","press","exec","edu","rel","bank","court","sci","milit","congr")
# 
# graph1<-graph_from_adjacency_matrix(mod[[1]],mode="undirected")
# graph2<-graph_from_adjacency_matrix(mod[[2]],mode="undirected")
# graph3<-graph_from_adjacency_matrix(mod[[3]],mode="undirected")
# 
# 
# plotGraph(mod[[1]],layout = layout_in_circle,vc= "white",nodesize=30)
# plotGraph(mod[[2]],layout = layout_in_circle,vc= "white",nodesize=30)
# plotGraph(mod[[3]],layout = layout_in_circle,vc= "white",nodesize=30)
# 
# 

sumba<-list()

for(x in 1:X){
  
  sumba[[x]]<-matrix(0,p,p)
  
  
}



for(x in 1:X){
  
  for(k in 4001:5000){
    
    sumba[[x]]<-sumba[[x]]+L.testXBA[[P]][[k]][[x]]
    
    
    
  }
  
}


for(x in 1:X){
  
  
  sumba[[x]]<-sumba[[x]]/1000
  
  
  
}




modba<-list()


for(x in 1:X){
  
  modba[[x]]<-(sumba[[x]]+t(sumba[[x]]))/2
  
  modba[[x]]<-round(modba[[x]])
  
  diag(modba[[x]])<-0
  
  # rr<-roc(mod[[x]],MX.test[[x]])
  
  # aucBAL[x]<-rr$auc
  
}

library(igraph)

ppar2<-modba[[1]]+modba[[2]]+modba[[3]]

ppar3<-modba

for (x in 1:X){
  
  for (i in 1:p){
    
    for (j in 1:p){
      
      if (ppar2[i,j]==3 & i!=j){
        
        ppar3[[x]][i,j]<-10
        
      } else if (modba[[x]][i,j]!=0 & ppar2[i,j]!=3 & i!=j) {
        
        ppar3[[x]][i,j]<-100
        
        
      }
      
    }
    
  }
  
  rownames(ppar3[[x]])<-c("welf","high","sec","trans","park","chi","sci","ern","for",
                          "mil","bla","spa","env","hea","cit","crime","drug","edu")
  colnames(ppar3[[x]])<-c("welf","high","sec","trans","park","chi","sci","ern","for",
                          "mil","bla","spa","env","hea","cit","crime","drug","edu")
  
  
  graph<-graph_from_adjacency_matrix(ppar3[[x]])
  
  plotGraph(graph,layout = layout_in_circle,vc= "white",nodesize=22,dashed=TRUE,coloth="red",colbid ="black",cex=1)
  
}


# shared edges mat


sh_mat<-matrix(0,3,3)

for (i in 1:X){
  
  for (j in 1:X){
    
    if (i==j){
      
      sh_mat[i,j]<-(length(which(ppar3[[i]]!=0)))/2
      
      
    } else if (i<j) {
      
      
      sh_mat[i,j]<-(length(which(ppar3[[i]]+ppar3[[j]]==20 | ppar3[[i]]+ppar3[[j]]==200)))/2
      
      
      
    }
    
    
    
  } 
  
}















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

K<-2000

nuK<-matrix(0,K,(p*(p-1)/2))



aa<-1
bb<-3

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


TEST<-list()

library(igraph)

for(x in 1:X){
  
  
  sssum<-apply(PREDBEL[[x]][[P]][1501:2000,],2,sum)/500
  
  predd<-rep(0,45)
  
  predd[which(sssum>=0.5)]<-1
  
  preddmat<-matrix(0,p,p)
  
  preddmat[upper.tri(preddmat)]<-predd
  
  TEST[[x]]<-preddmat+t(preddmat)
  
  
}






ppar2<-TEST[[1]]+TEST[[2]]+TEST[[3]]

ppar3<-TEST

for (x in 1:X){
  
  for (i in 1:p){
    
    for (j in 1:p){
      
      if (ppar2[i,j]==3 & i!=j){
        
        ppar3[[x]][i,j]<-10
        
      } else if (TEST[[x]][i,j]!=0 & ppar2[i,j]!=3 & i!=j) {
        
        ppar3[[x]][i,j]<-100
        
        
      }
      
    }
    
  }
  
  rownames(ppar3[[x]])<-c("tv","press","lab","exec","edu","rel","comp","bank","court","congr")
  colnames(ppar3[[x]])<-c("tv","press","lab","exec","edu","rel","comp","bank","court","congr")
  
  
  graph<-graph_from_adjacency_matrix(ppar3[[x]])
  
  plotGraph(graph,layout = layout_in_circle,vc= "white",nodesize=30,dashed=TRUE,coloth="red",colbid ="black")
  
}








# FP[x]<-length(which(predd-MX.test[[x]]==1))
# FN[x]<-length(which(predd-MX.test[[x]]==-1))
# 
# 
# rr<-roc(predd,MX.test[[x]])
# 
# aucBEL[x]<-rr$auc

}



############################# start be **********************************************

################################################################



print("BE"); print(P)



#############

RES<-matrix(0,4,2)





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






LX<-list()

for (x in 1:X){
  
  LX[[x]]<-L
  
}

LcX<-list()

aX<- list()

for (x in 1:X){
  
  aX[[x]]<-1
  
}


LrX2<-list()

for (x in 1:X){
  
  LrX2[[x]]<-matrix(0,K,p*(p-1)/2)
  
}



PdeltaX<-list()

for (x in 1:X){
  
  PdeltaX[[x]]<-0
  
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

a<-0.02

qq<-0.2


# V   = total nodes + X
# K   = iterations number
# c.c = vector cell counts
# a   = vector prior cell counts


Z.k<-mobkron(V-1) # matrice Z  (V-1 x V-1)
Z.kt<-mobkron(V)  # matrice Z (V x V)
M.k<-solve(Z.k)   # matrice M                                          

It<-2^V           # spazio degli stati con X
I<-2^(V-1)        # spazio degli stati senza X
p.p<-rep(a/I,I)   # probabilit? a priori (uniforme)
c.c.p<-p.p*a      # conteggi a priori (uniforme)

# genero conteggi da dataset
N <- ncol(dataX[[1]])       # Number of nodes
nSample <- nrow(dataX[[1]]) # Number of samples
Resp <- c(0L,1L)
n<- nrow(dataX[[1]])



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
  
  ldifX[[x]]<--10^6
  
}

####################   START THE LOOP !!!!  ###############################


for (i in 2:K){
  
  for (x in 1:X){
    
    # MODEL SETTING
    
    
    
    # We begin setting the model equal to the previous and useful objects
    
    
    
    LX[[x]][i,]<-LX[[x]][i-1,]
    
    
    
    
    
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
    
    PdeltaX[[x]][i]<-(qq/(1-qq))^(sum(LrX2[[x]][i,])) 
    
    
    ldifX[[x]][i]<-(log(PdeltaX[[x]][i])+log.M.LX[[x]][i])-(log(PdeltaX[[x]][i-1])+log.M.LX[[x]][i-1]) 
    
    
    
    
    
    if ( ldifX[[x]][i] < log(runif(1,0,1)) ){ # we evaluate if we reject the actual state
      
      
      LX[[x]][i,]<-LX[[x]][i-1,]
      
      log.M.LX[[x]][i]<-log.M.LX[[x]][i-1]
      
      PdeltaX[[x]][i]<- PdeltaX[[x]][i-1]
      
    } else {
      
      LX[[x]][i,]<-l
      
    }
    
    
    
    LrX[[x]][i,]<-LX[[x]][i,-p.LZ]
    LrX2[[x]][i,]<-LX[[x]][i,p.L2]
    
    
    PREDBE[[x]][[P]][i,]<- LrX2[[x]][i,]
    
    
    print(c(i,2,P))
    
    
    
  }
  
}



library(igraph)

for(x in 1:X){
  
  
  sssum<-apply(PREDBE[[x]][[P]][1501:2000,],2,sum)/500
  
  predd<-rep(0,45)
  
  predd[which(sssum>=0.5)]<-1
  
  preddmat<-matrix(0,p,p)
  
  preddmat[upper.tri(preddmat)]<-predd
  
  preddmat<-preddmat+t(preddmat)
  
  
  rownames(preddmat)<-c("tv","press","exec","edu","rel","bank","court","sci","milit","congr")
  colnames(preddmat)<-c("tv","press","exec","edu","rel","bank","court","sci","milit","congr")
  
  
  graph1<-graph_from_adjacency_matrix(preddmat,mode="undirected")
  
  plotGraph(preddmat,layout = layout_in_circle,vc= "white",nodesize=30)
  
  # FP[x]<-length(which(predd-MX.test[[x]]==1))
  # FN[x]<-length(which(predd-MX.test[[x]]==-1))
  # 
  # 
  # rr<-roc(predd,MX.test[[x]])
  # 
  # aucBEL[x]<-rr$auc
  
}









RBEL<-matrix(0,10,2)

aucBEL<-c()

FP<-c()
FN<<-c()

sumbal<-list()

for(x in 1:X){
  
  sumbal[[x]]<-matrix(0,p,p)
  
  
}



for(x in 1:X){
  
  for(k in 20001:25000){
    
    sumbal[[x]]<-sumbal[[x]]+L.testXBAL[[P]][[k]][[x]]
    
    
    
  }
  
}


MCCBEL<-c()


for (P in 1:10){
  
  MCCBEL[P]<- ((9-RBEL[P,1])*(36-RBEL[P,2])-RBEL[P,2]*RBEL[P,1])/sqrt( ((9-RBEL[P,1])+RBEL[P,2]) * ((9-RBEL[P,1])+RBEL[P,1]) * ((36-RBEL[P,2])+RBEL[P,2]) * ((36-RBEL[P,2])+RBEL[P,1]) )
  
}





RBE<-matrix(0,4,2)

aucBE<-c()

FP<-c()
FN<-c()

for(P in 1:4){
  
  for(x in 1:X){
    
    
    sssum<-apply(PREDBE[[x]][[P]][2000:2500,],2,sum)/501
    
    predd<-rep(0,45)
    
    predd[which(sssum>=0.5)]<-1
    
    
    FP[x]<-length(which(predd-MX.test[[x]]==1))
    FN[x]<-length(which(predd-MX.test[[x]]==-1))
    
    
    rr<-roc(predd,MX.test[[x]])
    
    aucBE[x]<-rr$auc
    
  }
  
  AUCBE[P]<-mean(aucBE)
  RBE[P,]<-c(mean(FN),mean(FP))
}



MCCBE<-c()


for (P in 1:4){
  
  MCCBE[P]<- ((9-RBE[P,1])*(36-RBE[P,2])-RBE[P,2]*RBE[P,1])/sqrt( ((9-RBE[P,1])+RBE[P,2]) * ((9-RBE[P,1])+RBE[P,1]) * ((36-RBE[P,2])+RBE[P,2]) * ((36-RBE[P,2])+RBE[P,1]) )
  
}





RBAL<-matrix(0,1,2)

aucBAL<-c()

FN<-c()
FP<-c()

for(P in 1:1){
  
  for(x in 1:X){
    
    
    sssum<-apply(PREDBAL[[x]][[P]][20000:25000,],2,sum)/5001
    
    predd<-rep(0,45)
    
    predd[which(sssum>0.5)]<-1
    
    
    FP[x]<-length(which(predd-MX.test[[x]]==1))
    FN[x]<-length(which(predd-MX.test[[x]]==-1))
    
    
    rr<-roc(predd,MX.test[[x]])
    
    aucBAL[x]<-rr$auc
    
  }
  
  AUCBAL[P]<-mean(aucBAL)
  
  RBAL[P,]<-c(mean(FN),mean(FP))
}







RBA<-matrix(0,4,2)

aucBA<-c()

for(P in 1:1){
  
  for(x in 1:X){
    
    
    sssum<-apply(PREDBA[[x]][[P]][20000:25000,],2,sum)/5001
    
    predd<-rep(0,45)
    
    predd[which(sssum>0.5)]<-1
    
    FP[x]<-length(which(predd-MX.test[[x]]==1))
    FN[x]<-length(which(predd-MX.test[[x]]==-1))
    
    
    rr<-roc(predd,MX.test[[x]])
    
    aucBA[x]<-rr$auc
    
  }
  
  AUCBA[P]<-mean(aucBA)
  
  RBA[P,]<-c(mean(FN),mean(FP))
}



MCCBAL<-c()


for (P in 1:2){
  
  MCCBAL[P]<- ((9-RBAL[P,1])*(36-RBAL[P,2])-RBAL[P,2]*RBAL[P,1])/sqrt( ((9-RBAL[P,1])+RBAL[P,2]) * ((9-RBAL[P,1])+RBAL[P,1]) * ((36-RBAL[P,2])+RBAL[P,2]) * ((36-RBAL[P,2])+RBAL[P,1]) )
  
}


MCCBA<-c()


for (P in 1:2){
  
  MCCBA[P]<- ((9-RBA[P,1])*(36-RBA[P,2])-RBA[P,2]*RBA[P,1])/sqrt( ((9-RBA[P,1])+RBA[P,2]) * ((9-RBA[P,1])+RBA[P,1]) * ((36-RBA[P,2])+RBA[P,2]) * ((36-RBA[P,2])+RBA[P,1]) )
  
}






boxplot(MCCSL,MCCDSSL,MCCBA,MCCBAL,ylim=c(0.5,1),ylab="MCC",las=2,names=c("SL","DSSL","BA","BAL"),main="(A)")

+boxplot(AUCSL,AUCDSSL,AUCBA,AUCBAL,ylim=c(0.5,1),ylab="AUC",las=2,names=c("SL","DSSL","BA","BAL"),main="(A)")



boxplot(RBE[1:10,1],RBEL[1:10,1],ylim=c(0,6),ylab="FN",las=2,names=c("BE","BEL"),main="(A)")

boxplot(RBE[1:10,2],RBEL[1:4,2],ylim=c(0,6),ylab="FP",las=2,names=c("BE","BEL"),main="(A)")



boxplot(RBA[1:4,1],RBAL[1:4,1],ylim=c(0,6),ylab="FN",las=2,names=c("BA","BAL"),main="(A)")

boxplot(RBA[1:4,2],RBAL[1:4,2],ylim=c(0,5),ylab="FP",las=2,names=c("BA","BAL"),main="(A)")


apply(RBA,2,mean)
apply(RBAL,2,mean)


apply(RBE[1:6,],2,mean)
apply(RBEL[1:6,],2,mean)

apply(RSL,2,mean)



ssssum<-matrix(0,4,4)

for (k in 20000:25000){
  
  ssssum1<-matrix(0,4,4)
  
  ssssum1[which(thetaK[[k]]!=0)]<-1
  
  ssssum<-ssssum+ssssum1
  
}

ssssum<-ssssum/5001

round(ssssum)


boxplot((9-RSL[1:10,1])/9,(9-RDSSL[1:10,1])/9,(9-RBA[1:10,1])/9,(9-RBAL[1:10,1])/9,(9-RBE[1:10,1])/9,(9-RBEL[1:10,1])/9,ylim=c(0.5,1),ylab="TPR",las=2,names=c("SL","DSSL","BA","BAL","BE","BEL"),main="(A)")




boxplot((RSL[1:10,2])/36,(RDSSL[1:10,2])/36,(RBA[1:10,2])/36,(RBAL[1:10,2])/36,(RBE[1:10,2])/36,(RBEL[1:10,2])/36,ylim=c(0,0.2),ylab="FPR",las=2,names=c("SL","DSSL","BA","BAL","BE","BEL"),main="(A)")



RBAL<-matrix(0,2,2)


for(P in 1:1){
  
  sumbal<-list()
  
  for(x in 1:X){
    
    sumbal[[x]]<-matrix(0,p,p)
    
    
  }
  
  
  
  for(x in 1:X){
    
    for(k in 13001:15000){
      
      sumbal[[x]]<-sumbal[[x]]+L.testXBAL[[P]][[k]][[x]]
      
      
      
    }
    
  }
  
  
  for(x in 1:X){
    
    
    sumbal[[x]]<-sumbal[[x]]/2000
    
    
    
  }
  
  
  
  
  mod<-list()
  
  for(x in 1:X){
    
    mod[[x]]<-rep(0,45)
    
  }
  
  for(x in 1:X){
    
    mod[[x]]<-(as.vector(t(sumbal[[x]][upper.tri(sumbal[[x]])]))+as.vector(t(sumbal[[x]])[upper.tri(sumbal[[x]])]))/2
    
    pos<-which(mod[[x]]>0.5)
    
    mod[[x]]<-rep(0,45)
    
    mod[[x]][pos]<-1
    
    for(x in 1:X){
      
      
      sssum<-apply(PREDBA[[x]][[P]][20000:25000,],2,sum)/5001
      
      predd<-rep(0,45)
      
      predd[which(sssum>0.5)]<-1
      
      FP[x]<-length(which(predd-MX.test[[x]]==1))
      FN[x]<-length(which(predd-MX.test[[x]]==-1))
      
      
      rr<-roc(predd,MX.test[[x]])
      
      aucBA[x]<-rr$auc
      
    }
    
    AUCBA[P]<-mean(aucBA)
    
  }
  
  FP<-c()
  FN<-c()
  
  for(x in 1:X){
    
    
    FP[x]<-length(which(mod[[x]]-MX.test[[x]]==1))
    FN[x]<-length(which(mod[[x]]-MX.test[[x]]==-1))
    
  }
  
  
  RBAL[P,]<-c(mean(FN),mean(FP))
  
}


RBA<-matrix(0,2,2)


for(P in 1:2){
  
  sumba<-list()
  
  for(x in 1:X){
    
    sumba[[x]]<-matrix(0,p,p)
    
    
  }
  
  
  
  for(x in 1:X){
    
    for(k in 20001:25000){
      
      sumba[[x]]<-sumba[[x]]+L.testXBA[[P]][[k]][[x]]
      
      
      
    }
    
  }
  
  
  for(x in 1:X){
    
    
    sumba[[x]]<-sumba[[x]]/5000
    
    
    
  }
  
  
  
  
  mod<-list()
  
  for(x in 1:X){
    
    mod[[x]]<-rep(0,45)
    
  }
  
  for(x in 1:X){
    
    mod[[x]]<-(as.vector(t(sumba[[x]][upper.tri(sumba[[x]])]))+as.vector(t(sumba[[x]])[upper.tri(sumba[[x]])]))/2
    
    pos<-which(mod[[x]]>0.5)
    
    mod[[x]]<-rep(0,45)
    
    mod[[x]][pos]<-1
    
    rr<-roc(mod[[x]],MX.test[[x]])
    
    aucBA[x]<-rr$auc
    
  }
  
  AUCBA[P]<-mean(aucBA)
  
  FP<-c()
  FN<-c()
  
  for(x in 1:X){
    
    
    FP[x]<-length(which(mod[[x]]-MX.test[[x]]==1))
    FN[x]<-length(which(mod[[x]]-MX.test[[x]]==-1))
    
  }
  
  
  RBA[P,]<-c(mean(FN),mean(FP))
  
}


















ssssum<-matrix(0,4,4)

for (k in 40001:50000){
  
  ssssum1<-thetaK[[k]]
  
  #ssssum1[which(thetaK[[k]]!=0)]<-1
  
  ssssum<-ssssum+ssssum1
  
}

ssssum<-ssssum/10000

round(ssssum)






ssssum<-matrix(0,4,4)

for (k in 5001:6000){
  
  ssssum1<-matrix(0,4,4)
  
  ssssum1[which(thetaK[[k]]!=0)]<-1
  
  ssssum<-ssssum+ssssum1
  
}

ssssum<-ssssum/1000

round(ssssum)




sumba<-list()

for(x in 1:X){
  
  sumba[[x]]<-matrix(0,p,p)
  
  
}



for(x in 1:X){
  
  for(k in 20001:25000){
    
    sumba[[x]]<-sumba[[x]]+L.testXBA[[P]][[k]][[x]]
    
    
    
  }
  
}


for(x in 1:X){
  
  
  sumba[[x]]<-sumba[[x]]/5000
  
  
  
}

}


mod<-list()

for(x in 1:X){
  
  mod[[x]]<-rep(0,45)
  
}

for(x in 1:X){
  
  mod[[x]]<-(as.vector(t(sumba[[x]][upper.tri(sumba[[x]])]))+as.vector(t(sumba[[x]])[upper.tri(sumba[[x]])]))/2
  
  pos<-which(mod[[x]]>0.5)
  
  mod[[x]]<-rep(0,45)
  
  mod[[x]][pos]<-1
  
}

FP<-c()
FN<-c()

for(x in 1:X){
  
  
  FP[x]<-length(which(mod[[x]]-MX.test[[x]]==1))
  FN[x]<-length(which(mod[[x]]-MX.test[[x]]==-1))
  
}




boxplot(0.79,0.98,0.82,0.87,0.88,0.91,ylim=c(0.5,1),ylab="MCC",las=2,names=c("SL","DSSL","BA","BAL","BE","BEL"),main="(A)")





ssssum<-matrix(0,4,4)

for (k in 40001:50000){
  
  ssssum1<-matrix(0,4,4)
  
  ssssum1[which(thetaK[[k]]!=0)]<-1
  
  ssssum<-ssssum+ssssum1
  
}

ssssum<-ssssum/10000

round(ssssum)




ssssum<-matrix(0,4,4)

for (k in 500:1000){
  
  ssssum1<-matrix(0,4,4)
  
  ssssum1[which(thetaK[[k]]!=0)]<-1
  
  ssssum<-ssssum+ssssum1
  
}

ssssum<-ssssum/500

round(ssssum)









################### OK


MCCSL<-c()


for (P in 1:6){
  
  MCCSL[P]<- ((9-RSL[P,1])*(36-RSL[P,2])-RSL[P,2]*RSL[P,1])/sqrt( ((9-RSL[P,1])+RSL[P,2]) * ((9-RSL[P,1])+RSL[P,1]) * ((36-RSL[P,2])+RSL[P,2]) * ((36-RSL[P,2])+RSL[P,1]) )
  
}

MCCDSSL<-c()


for (P in 1:6){
  
  MCCDSSL[P]<- ((9-RDSSL[P,1])*(36-RDSSL[P,2])-RDSSL[P,2]*RDSSL[P,1])/sqrt( ((9-RDSSL[P,1])+RDSSL[P,2]) * ((9-RDSSL[P,1])+RDSSL[P,1]) * ((36-RDSSL[P,2])+RDSSL[P,2]) * ((36-RDSSL[P,2])+RDSSL[P,1]) )
  
}


RBAL<-matrix(0,7,2)

ssssum<-matrix(0,3,3)

for (k in 1501:2000){
  
  ssssum1<-matrix(0,3,3)
  
  ssssum1[which(thetaK[[k]]!=0)]<-1
  
  ssssum<-ssssum+ssssum1
  
}

ssssum<-ssssum/500


for(P in 1:1){
  
  sumbal<-list()
  
  for(x in 1:X){
    
    sumbal[[x]]<-matrix(0,p,p)
    
    
  }
  
  
  
  for(x in 1:X){
    
    for(k in 8001:10000){
      
      sumbal[[x]]<-sumbal[[x]]+L.testXBAL[[P]][[k]][[x]]
      
      
      
    }
    
  }
  
  
  for(x in 1:X){
    
    
    sumbal[[x]]<-sumbal[[x]]/2000
    
    
    
  }
  
  
  
  
  mod<-list()
  
  for(x in 1:X){
    
    mod[[x]]<-rep(0,21)
    
  }
  
  for(x in 1:X){
    
    mod[[x]]<-(as.vector(t(sumbal[[x]][upper.tri(sumbal[[x]])]))+as.vector(t(sumbal[[x]])[upper.tri(sumbal[[x]])]))/2
    
    pos<-which(mod[[x]]>0.5)
    
    mod[[x]]<-rep(0,21)
    
    mod[[x]][pos]<-1
    
    # rr<-roc(mod[[x]],MX.test[[x]])
    
    # aucBAL[x]<-rr$auc
    
  }
  
  AUCBAL[P]<-mean(aucBAL)
  
  
  FP<-c()
  FN<-c()
  
  for(x in 1:X){
    
    
    FP[x]<-length(which(mod[[x]]-MX.test[[x]]==1))
    FN[x]<-length(which(mod[[x]]-MX.test[[x]]==-1))
    
  }
  
  
  RBAL[P,]<-c(mean(FN),mean(FP))
  
}







sumba<-list()

for(x in 1:X){
  
  sumba[[x]]<-matrix(0,p,p)
  
  
}



for(x in 1:X){
  
  for(k in 801:1000){
    
    sumba[[x]]<-sumba[[x]]+L.testXBA[[P]][[k]][[x]]
    
    
    
  }
  
}


for(x in 1:X){
  
  
  sumba[[x]]<-sumba[[x]]/200
  
  
  
}




mod<-list()

for(x in 1:X){
  
  mod[[x]]<-rep(0,21)
  
}

for(x in 1:X){
  
  mod[[x]]<-(as.vector(t(sumba[[x]][upper.tri(sumba[[x]])]))+as.vector(t(sumba[[x]])[upper.tri(sumba[[x]])]))/2
  
  pos<-which(mod[[x]]>0.5)
  
  mod[[x]]<-rep(0,21)
  
  mod[[x]][pos]<-1
  
  # rr<-roc(mod[[x]],MX.test[[x]])
  
  # aucBAL[x]<-rr$auc
  
}


















RBA<-matrix(0,10,2)


for(P in 1:6){
  
  sumba<-list()
  
  for(x in 1:X){
    
    sumba[[x]]<-matrix(0,p,p)
    
    
  }
  
  
  
  for(x in 1:X){
    
    for(k in 40001:50000){
      
      sumba[[x]]<-sumba[[x]]+L.testXBA[[P]][[k]][[x]]
      
      
      
    }
    
  }
  
  
  for(x in 1:X){
    
    
    sumba[[x]]<-sumba[[x]]/10000
    
    
    
  }
  
  
  
  
  mod<-list()
  
  for(x in 1:X){
    
    mod[[x]]<-rep(0,45)
    
  }
  
  for(x in 1:X){
    
    mod[[x]]<-(as.vector(t(sumba[[x]][upper.tri(sumba[[x]])]))+as.vector(t(sumba[[x]])[upper.tri(sumba[[x]])]))/2
    
    pos<-which(mod[[x]]>0.5)
    
    mod[[x]]<-rep(0,45)
    
    mod[[x]][pos]<-1
    
    rr<-roc(mod[[x]],MX.test[[x]])
    
    aucBA[x]<-rr$auc
    
  }
  
  AUCBA[P]<-mean(aucBA)
  
  FP<-c()
  FN<-c()
  
  for(x in 1:X){
    
    
    FP[x]<-length(which(mod[[x]]-MX.test[[x]]==1))
    FN[x]<-length(which(mod[[x]]-MX.test[[x]]==-1))
    
  }
  
  
  RBA[P,]<-c(mean(FN),mean(FP))
  
}



MCCBAL<-c()


for (P in 1:5){
  
  MCCBAL[P]<- ((9-RBAL[P,1])*(36-RBAL[P,2])-RBAL[P,2]*RBAL[P,1])/sqrt( ((9-RBAL[P,1])+RBAL[P,2]) * ((9-RBAL[P,1])+RBAL[P,1]) * ((36-RBAL[P,2])+RBAL[P,2]) * ((36-RBAL[P,2])+RBAL[P,1]) )
  
}


MCCBA<-c()


for (P in 1:5){
  
  MCCBA[P]<- ((9-RBA[P,1])*(36-RBA[P,2])-RBA[P,2]*RBA[P,1])/sqrt( ((9-RBA[P,1])+RBA[P,2]) * ((9-RBA[P,1])+RBA[P,1]) * ((36-RBA[P,2])+RBA[P,2]) * ((36-RBA[P,2])+RBA[P,1]) )
  
}



boxplot(MCCSL,MCCDSSL,MCCBA,MCCBAL,ylim=c(0.5,1),ylab="MCC",las=2,names=c("SL","DSSL","MCCBA","MCCBAL"),main="(A)")

boxplot(AUCSL,AUCDSSL,AUCBA,AUCBAL,ylim=c(0.5,1),ylab="MCC",las=2,names=c("SL","DSSL","MCCBA","MCCBAL"),main="(A)")


mean(MCCSL)
mean(MCCDSSL)
mean(MCCBA)
mean(MCCBAL)

apply(RBA[1:6,],2,mean)
apply(RBAL[1:6,],2,mean)


mean(AUCSL[1:6])
mean(AUCDSSL[1:6])
mean(AUCBA[1:6])
mean(AUCBAL[1:6])





RBEL<-matrix(0,10,2)

aucBEL<-c()

FP<-c()
FN<-c()

for(P in 1:4){
  
  for(x in 1:X){
    
    
    sssum<-apply(PREDBEL[[x]][[P]][8001:10000,],2,sum)/2000
    
    predd<-rep(0,21)
    
    predd[which(sssum>=0.5)]<-1
    
    
    # FP[x]<-length(which(predd-MX.test[[x]]==1))
    # FN[x]<-length(which(predd-MX.test[[x]]==-1))
    # 
    # 
    # rr<-roc(predd,MX.test[[x]])
    # 
    # aucBEL[x]<-rr$auc
    
  }
  
  AUCBEL[P]<-mean(aucBEL)
  RBEL[P,]<-c(mean(FN),mean(FP))
}









########################## COUPLE-COMPARISON PLOTS

ppar3<-matrix(0,p,p)

couple_sharing<-c()

B<-matrix(c(1,1,2,2,3,3),3,2)

for (t in 1:3){

  ppar2<-mod[[B[t,1]]]+mod[[B[t,2]]]
  
  for (i in 1:p){
    
    for (j in 1:p){
      
      if (ppar2[i,j]==2 & i!=j){
        
        ppar3[i,j]<-10
        
      } else if (ppar2[i,j]==1 & i!=j) {
        
        ppar3[i,j]<-100
        
        
      }
      
    }
    
  }
  
  rownames(ppar3)<-c("welf","high","sec","trans","park","chi","sci","ern","for",
                          "mil","bla","spa","env","hea","cit","crime","drug","edu")
  colnames(ppar3)<-c("welf","high","sec","trans","park","chi","sci","ern","for",
                          "mil","bla","spa","env","hea","cit","crime","drug","edu")
  
  graph<-graph_from_adjacency_matrix(ppar3)
  
  plotGraph(graph,layout = layout_in_circle,vc= "white",nodesize=22,dashed=TRUE,coloth="blue",colbid ="black",cex=1)
  
  couple_sharing[t]<- round(sum(ppar3==10)/(p*(p-1)),2)

}


couple_sharingBA<-c(0.02, 0.05, 0.04)











