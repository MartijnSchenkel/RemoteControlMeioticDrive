library(tidyverse)
library(viridis)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


SA_function <- function(G, kD, s, t, hA, z, r) {
  # G is the number of generations
  # kD is the skew in favor of the Y chromosome
  # s is the selection in males
  # t is the selection in females
  # hA is the dominance of SA locus
  # z is the fitness cost of the driver
  # r is the recombination rate
  
  DfreqVector<-numeric(G)
  AfreqVector<-numeric(G)
  SexRatioVector<-numeric(G)
  DfreqVector<-numeric(G)
  AfreqVector<-numeric(G)
  SexRatioVector<-numeric(G)
  
  DA2D2Vector<-numeric(G)
  DA2YVector<-numeric(G)
  DD2YVector<-numeric(G)
  
  #BURN IN 
  
  u11f <- 1
  u12f <- 1
  u21f <- 1
  u22f <- 1
  
  u11m <- 1
  u12m <- 1-z
  u21m <- 1-z
  u22m <- 1-z
  
  v11f <- 1
  v12f <- 1-hA*t
  v21f <- 1-hA*t
  v22f <- 1-t
  
  v11m <- 1-s
  v12m <- 1-(1-hA)*s
  v21m <- 1-(1-hA)*s
  v22m <- 1
  
  w11m <- u11m*v11m
  w12m <- u11m*v12m
  w13m <- u12m*v11m
  w14m <- u12m*v12m
  w21m <- u11m*v21m
  w22m <- u11m*v22m
  w23m <- u12m*v21m
  w24m <- u12m*v22m
  w31m <- u21m*v11m
  w32m <- u21m*v12m
  w33m <- u22m*v11m
  w34m <- u22m*v12m
  w41m <- u21m*v21m
  w42m <- u21m*v22m
  w43m <- u22m*v21m
  w44m <- u22m*v22m
  
  w11f <- u11f*v11f
  w12f <- u11f*v12f
  w13f <- u12f*v11f
  w14f <- u12f*v12f
  w21f <- u11f*v21f
  w22f <- u11f*v22f
  w23f <- u12f*v21f
  w24f <- u12f*v22f
  w31f <- u21f*v11f
  w32f <- u21f*v12f
  w33f <- u22f*v11f
  w34f <- u22f*v12f
  w41f <- u21f*v21f
  w42f <- u21f*v22f
  w43f <- u22f*v21f
  w44f <- u22f*v22f
  
  g1e <- 0.5
  g2e <- 0.5
  g3e <- 0.0
  g4e <- 0.0
  
  g1Xs <- 0.25
  g1Ys <- 0.25
  g2Xs <- 0.25
  g2Ys <- 0.25
  
  g3Xs <- 0.0
  g3Ys <- 0.0
  g4Xs <- 0.0
  g4Ys <- 0.0
  
  for (d in 1:500){
    wfbar<-((g1e*g1Xs*w11f)+(g1e*g2Xs*w12f)+(g1e*g3Xs*w13f)+(g1e*g4Xs*w14f)+(g2e*g1Xs*w21f)+(g2e*g2Xs*w22f)+(g2e*g3Xs*w23f)+(g2e*g4Xs*w24f)+(g3e*g1Xs*w31f)+(g3e*g2Xs*w32f)+(g3e*g3Xs*w33f)+(g3e*g4Xs*w34f)+(g4e*g1Xs*w41f)+(g4e*g2Xs*w42f)+(g4e*g3Xs*w43f)+(g4e*g4Xs*w44f))
    g1e_pr<- (1/wfbar)*(1*g1e*g1Xs*w11f + 0.5* g1e*g2Xs*w12f + 0.5* g1e*g3Xs*w13f + 0.5*(1-r)* g1e*g4Xs*w14f + 0.5* g2e*g1Xs*w21f + 0.5*(r)* g2e*g3Xs*w23f + 0.5* g3e*g1Xs*w31f + 0.5*(r)* g3e*g2Xs*w32f + 0.5*(1-r)* g4e*g1Xs*w41f) 
    g2e_pr<- (1/wfbar)*(0.5* g1e*g2Xs*w12f + 0.5*(r)* g1e*g4Xs*w14f + 0.5* g2e*g1Xs*w21f + 1* g2e*g2Xs*w22f + 0.5*(1-r)* g2e*g3Xs*w23f + 0.5* g2e*g4Xs*w24f + 0.5*(1-r)* g3e*g2Xs*w32f + 0.5*(r)* g4e*g1Xs*w41f + 0.5* g4e*g2Xs*w42f)
    g3e_pr<- (1/wfbar)*(0.5* g1e*g3Xs*w13f + 0.5*(r)* g1e*g4Xs*w14f + 0.5*(1-r)* g2e*g3Xs*w23f + 0.5* g3e*g1Xs*w31f + 0.5*(1-r)* g3e*g2Xs*w32f + 1* g3e*g3Xs*w33f + 0.5* g3e*g4Xs*w34f + 0.5*(r)* g4e*g1Xs*w41f + 0.5* g4e*g3Xs*w43f)
    g4e_pr<- (1/wfbar)*(0.5*(1-r)* g1e*g4Xs*w14f + 0.5*(r)* g2e*g3Xs*w23f + 0.5* g2e*g4Xs*w24f + 0.5*(r)* g3e*g2Xs*w32f + 0.5* g3e*g4Xs*w34f + 0.5*(1-r)* g4e*g1Xs*w41f + 0.5* g4e*g2Xs*w42f + 0.5* g4e*g3Xs*w43f + 1* g4e*g4Xs*w44f)
    
    #Sperm haplotypes (8 total)
    
    wmbar<-((g1e*g1Ys*w11m)+(g1e*g2Ys*w12m)+(g1e*g3Ys*w13m)+(g1e*g4Ys*w14m)+(g2e*g1Ys*w21m)+(g2e*g2Ys*w22m)+(g2e*g3Ys*w23m)+(g2e*g4Ys*w24m)+(g3e*g1Ys*w31m)+(g3e*g2Ys*w32m)+(g3e*g3Ys*w33m)+(g3e*g4Ys*w34m)+(g4e*g1Ys*w41m)+(g4e*g2Ys*w42m)+(g4e*g3Ys*w43m)+(g4e*g4Ys*w44m))
    g1Xs_pr<- (1/wmbar)*(0.5*1* g1e*g1Ys*w11m + 0.5*0.5* g1e*g2Ys*w12m + (0.5-kD)*0.5* g1e*g3Ys*w13m + (0.5-kD)*0.5*(1-r)* g1e*g4Ys*w14m + 0.5*0.5* g2e*g1Ys*w21m + (0.5-kD)*0.5*(r)* g2e*g3Ys*w23m + (0.5-kD)*0.5* g3e*g1Ys*w31m + (0.5-kD)*0.5*(r)* g3e*g2Ys*w32m + (0.5-kD)*0.5*(1-r)* g4e*g1Ys*w41m)
    g2Xs_pr<- (1/wmbar)*(0.5*0.5* g1e*g2Ys*w12m + (0.5-kD)*0.5*(r)* g1e*g4Ys*w14m + 0.5*0.5* g2e*g1Ys*w21m + 0.5*1* g2e*g2Ys*w22m + (0.5-kD)*0.5*(1-r)* g2e*g3Ys*w23m + (0.5-kD)*0.5* g2e*g4Ys*w24m + (0.5-kD)*0.5*(1-r)* g3e*g2Ys*w32m + (0.5-kD)*0.5*(r)* g4e*g1Ys*w41m + (0.5-kD)*0.5* g4e*g2Ys*w42m)
    g3Xs_pr<- (1/wmbar)*((0.5-kD)*0.5* g1e*g3Ys*w13m + (0.5-kD)*0.5*(r)* g1e*g4Ys*w14m + (0.5-kD)*0.5*(1-r)* g2e*g3Ys*w23m + (0.5-kD)*0.5* g3e*g1Ys*w31m + (0.5-kD)*0.5*(1-r)* g3e*g2Ys*w32m + (0.5-kD)*1* g3e*g3Ys*w33m + (0.5-kD)*0.5* g3e*g4Ys*w34m + (0.5-kD)*0.5*(r)* g4e*g1Ys*w41m + (0.5-kD)*0.5* g4e*g3Ys*w43m)
    g4Xs_pr<- (1/wmbar)*((0.5-kD)*0.5*(1-r)* g1e*g4Ys*w14m + (0.5-kD)*0.5*(r)* g2e*g3Ys*w23m + (0.5-kD)*0.5* g2e*g4Ys*w24m + (0.5-kD)*0.5*(r)* g3e*g2Ys*w32m + (0.5-kD)*0.5* g3e*g4Ys*w34m + (0.5-kD)*0.5*(1-r)* g4e*g1Ys*w41m + (0.5-kD)*0.5* g4e*g2Ys*w42m + (0.5-kD)*0.5* g4e*g3Ys*w43m + (0.5-kD)*1* g4e*g4Ys*w44m) 
    
    g1Ys_pr<- (1/wmbar)*(0.5*1* g1e*g1Ys*w11m + 0.5*0.5* g1e*g2Ys*w12m + (0.5+kD)*0.5* g1e*g3Ys*w13m + (0.5+kD)*0.5*(1-r)* g1e*g4Ys*w14m + 0.5*0.5* g2e*g1Ys*w21m + (0.5+kD)*0.5*(r)* g2e*g3Ys*w23m + (0.5+kD)*0.5* g3e*g1Ys*w31m + (0.5+kD)*0.5*(r)* g3e*g2Ys*w32m + (0.5+kD)*0.5*(1-r)* g4e*g1Ys*w41m)
    g2Ys_pr<- (1/wmbar)*(0.5*0.5* g1e*g2Ys*w12m + (0.5+kD)*0.5*(r)* g1e*g4Ys*w14m + 0.5*0.5* g2e*g1Ys*w21m + 0.5*1* g2e*g2Ys*w22m + (0.5+kD)*0.5*(1-r)* g2e*g3Ys*w23m + (0.5+kD)*0.5* g2e*g4Ys*w24m + (0.5+kD)*0.5*(1-r)* g3e*g2Ys*w32m + (0.5+kD)*0.5*(r)* g4e*g1Ys*w41m + (0.5+kD)*0.5* g4e*g2Ys*w42m)
    g3Ys_pr<- (1/wmbar)*((0.5+kD)*0.5* g1e*g3Ys*w13m + (0.5+kD)*0.5*(r)* g1e*g4Ys*w14m + (0.5+kD)*0.5*(1-r)* g2e*g3Ys*w23m + (0.5+kD)*0.5* g3e*g1Ys*w31m + (0.5+kD)*0.5*(1-r)* g3e*g2Ys*w32m + (0.5+kD)*1* g3e*g3Ys*w33m + (0.5+kD)*0.5* g3e*g4Ys*w34m + (0.5+kD)*0.5*(r)* g4e*g1Ys*w41m + (0.5+kD)*0.5* g4e*g3Ys*w43m)
    g4Ys_pr<- (1/wmbar)*((0.5+kD)*0.5*(1-r)* g1e*g4Ys*w14m + (0.5+kD)*0.5*(r)* g2e*g3Ys*w23m + (0.5+kD)*0.5* g2e*g4Ys*w24m + (0.5+kD)*0.5*(r)* g3e*g2Ys*w32m + (0.5+kD)*0.5* g3e*g4Ys*w34m + (0.5+kD)*0.5*(1-r)* g4e*g1Ys*w41m + (0.5+kD)*0.5* g4e*g2Ys*w42m + (0.5+kD)*0.5* g4e*g3Ys*w43m + (0.5+kD)*1* g4e*g4Ys*w44m)
    
    g1e <- g1e_pr 
    g2e <- g2e_pr 
    g3e <- g3e_pr 
    g4e <- g4e_pr 
    
    g1Xs <- g1Xs_pr 
    g2Xs <- g2Xs_pr
    g3Xs <- g3Xs_pr
    g4Xs <- g4Xs_pr
    g1Ys <- g1Ys_pr 
    g2Ys <- g2Ys_pr
    g3Ys <- g3Ys_pr
    g4Ys <- g4Ys_pr
    
  }
  
  g1e <- g1e_pr 
  g2e <- g2e_pr 
  g3e <- g3e_pr 
  g4e <- g4e_pr 
  
  g1Xs <- g1Xs_pr 
  g2Xs <- g2Xs_pr
  g3Xs <- g3Xs_pr
  g4Xs <- g4Xs_pr
  
  g1Ys <- g1Ys_pr-0.0001
  g2Ys <- g2Ys_pr
  g3Ys <- g3Ys_pr
  g4Ys <- g4Ys_pr+0.0001
  
  #SIMULATION 
  
  #u is the trans drive locus (1 is neutral, 2 drives)
  u11f <- 1
  u12f <- 1
  u21f <- 1
  u22f <- 1
  
  u11m <- 1
  u12m <- 1-z
  u21m <- 1-z
  u22m <- 1-z
  
  # v is the female beneficial locus (1 is female beneficial, 2 is male beneficial) 
  v11f <- 1
  v12f <- 1-hA*t
  v21f <- 1-hA*t
  v22f <- 1-t
  
  v11m <- 1-s
  v12m <- 1-(1-hA)*s
  v21m <- 1-(1-hA)*s
  v22m <- 1
  
  w11m <- u11m*v11m
  w12m <- u11m*v12m
  w13m <- u12m*v11m
  w14m <- u12m*v12m
  w21m <- u11m*v21m
  w22m <- u11m*v22m
  w23m <- u12m*v21m
  w24m <- u12m*v22m
  w31m <- u21m*v11m
  w32m <- u21m*v12m
  w33m <- u22m*v11m
  w34m <- u22m*v12m
  w41m <- u21m*v21m
  w42m <- u21m*v22m
  w43m <- u22m*v21m
  w44m <- u22m*v22m
  
  w11f <- u11f*v11f
  w12f <- u11f*v12f
  w13f <- u12f*v11f
  w14f <- u12f*v12f
  w21f <- u11f*v21f
  w22f <- u11f*v22f
  w23f <- u12f*v21f
  w24f <- u12f*v22f
  w31f <- u21f*v11f
  w32f <- u21f*v12f
  w33f <- u22f*v11f
  w34f <- u22f*v12f
  w41f <- u21f*v21f
  w42f <- u21f*v22f
  w43f <- u22f*v21f
  w44f <- u22f*v22f
  
  SexRatioVector[d]<-(g1Ys + g2Ys + g3Ys + g4Ys)
  
  #Haplotype Frequency Recursions
  
  for (d in 1:G){
    wfbar<-((g1e*g1Xs*w11f)+(g1e*g2Xs*w12f)+(g1e*g3Xs*w13f)+(g1e*g4Xs*w14f)+(g2e*g1Xs*w21f)+(g2e*g2Xs*w22f)+(g2e*g3Xs*w23f)+(g2e*g4Xs*w24f)+(g3e*g1Xs*w31f)+(g3e*g2Xs*w32f)+(g3e*g3Xs*w33f)+(g3e*g4Xs*w34f)+(g4e*g1Xs*w41f)+(g4e*g2Xs*w42f)+(g4e*g3Xs*w43f)+(g4e*g4Xs*w44f))
    g1e_pr<- (1/wfbar)*(1*g1e*g1Xs*w11f + 0.5* g1e*g2Xs*w12f + 0.5* g1e*g3Xs*w13f + 0.5*(1-r)* g1e*g4Xs*w14f + 0.5* g2e*g1Xs*w21f + 0.5*(r)* g2e*g3Xs*w23f + 0.5* g3e*g1Xs*w31f + 0.5*(r)* g3e*g2Xs*w32f + 0.5*(1-r)* g4e*g1Xs*w41f) 
    g2e_pr<- (1/wfbar)*(0.5* g1e*g2Xs*w12f + 0.5*(r)* g1e*g4Xs*w14f + 0.5* g2e*g1Xs*w21f + 1* g2e*g2Xs*w22f + 0.5*(1-r)* g2e*g3Xs*w23f + 0.5* g2e*g4Xs*w24f + 0.5*(1-r)* g3e*g2Xs*w32f + 0.5*(r)* g4e*g1Xs*w41f + 0.5* g4e*g2Xs*w42f)
    g3e_pr<- (1/wfbar)*(0.5* g1e*g3Xs*w13f + 0.5*(r)* g1e*g4Xs*w14f + 0.5*(1-r)* g2e*g3Xs*w23f + 0.5* g3e*g1Xs*w31f + 0.5*(1-r)* g3e*g2Xs*w32f + 1* g3e*g3Xs*w33f + 0.5* g3e*g4Xs*w34f + 0.5*(r)* g4e*g1Xs*w41f + 0.5* g4e*g3Xs*w43f)
    g4e_pr<- (1/wfbar)*(0.5*(1-r)* g1e*g4Xs*w14f + 0.5*(r)* g2e*g3Xs*w23f + 0.5* g2e*g4Xs*w24f + 0.5*(r)* g3e*g2Xs*w32f + 0.5* g3e*g4Xs*w34f + 0.5*(1-r)* g4e*g1Xs*w41f + 0.5* g4e*g2Xs*w42f + 0.5* g4e*g3Xs*w43f + 1* g4e*g4Xs*w44f)
    
    #Sperm haplotypes (8 total)
    
    wmbar<-((g1e*g1Ys*w11m)+(g1e*g2Ys*w12m)+(g1e*g3Ys*w13m)+(g1e*g4Ys*w14m)+(g2e*g1Ys*w21m)+(g2e*g2Ys*w22m)+(g2e*g3Ys*w23m)+(g2e*g4Ys*w24m)+(g3e*g1Ys*w31m)+(g3e*g2Ys*w32m)+(g3e*g3Ys*w33m)+(g3e*g4Ys*w34m)+(g4e*g1Ys*w41m)+(g4e*g2Ys*w42m)+(g4e*g3Ys*w43m)+(g4e*g4Ys*w44m))
    g1Xs_pr<- (1/wmbar)*(0.5*1* g1e*g1Ys*w11m + 0.5*0.5* g1e*g2Ys*w12m + (0.5-kD)*0.5* g1e*g3Ys*w13m + (0.5-kD)*0.5*(1-r)* g1e*g4Ys*w14m + 0.5*0.5* g2e*g1Ys*w21m + (0.5-kD)*0.5*(r)* g2e*g3Ys*w23m + (0.5-kD)*0.5* g3e*g1Ys*w31m + (0.5-kD)*0.5*(r)* g3e*g2Ys*w32m + (0.5-kD)*0.5*(1-r)* g4e*g1Ys*w41m)
    g2Xs_pr<- (1/wmbar)*(0.5*0.5* g1e*g2Ys*w12m + (0.5-kD)*0.5*(r)* g1e*g4Ys*w14m + 0.5*0.5* g2e*g1Ys*w21m + 0.5*1* g2e*g2Ys*w22m + (0.5-kD)*0.5*(1-r)* g2e*g3Ys*w23m + (0.5-kD)*0.5* g2e*g4Ys*w24m + (0.5-kD)*0.5*(1-r)* g3e*g2Ys*w32m + (0.5-kD)*0.5*(r)* g4e*g1Ys*w41m + (0.5-kD)*0.5* g4e*g2Ys*w42m)
    g3Xs_pr<- (1/wmbar)*((0.5-kD)*0.5* g1e*g3Ys*w13m + (0.5-kD)*0.5*(r)* g1e*g4Ys*w14m + (0.5-kD)*0.5*(1-r)* g2e*g3Ys*w23m + (0.5-kD)*0.5* g3e*g1Ys*w31m + (0.5-kD)*0.5*(1-r)* g3e*g2Ys*w32m + (0.5-kD)*1* g3e*g3Ys*w33m + (0.5-kD)*0.5* g3e*g4Ys*w34m + (0.5-kD)*0.5*(r)* g4e*g1Ys*w41m + (0.5-kD)*0.5* g4e*g3Ys*w43m)
    g4Xs_pr<- (1/wmbar)*((0.5-kD)*0.5*(1-r)* g1e*g4Ys*w14m + (0.5-kD)*0.5*(r)* g2e*g3Ys*w23m + (0.5-kD)*0.5* g2e*g4Ys*w24m + (0.5-kD)*0.5*(r)* g3e*g2Ys*w32m + (0.5-kD)*0.5* g3e*g4Ys*w34m + (0.5-kD)*0.5*(1-r)* g4e*g1Ys*w41m + (0.5-kD)*0.5* g4e*g2Ys*w42m + (0.5-kD)*0.5* g4e*g3Ys*w43m + (0.5-kD)*1* g4e*g4Ys*w44m) 
    
    g1Ys_pr<- (1/wmbar)*(0.5*1* g1e*g1Ys*w11m + 0.5*0.5* g1e*g2Ys*w12m + (0.5+kD)*0.5* g1e*g3Ys*w13m + (0.5+kD)*0.5*(1-r)* g1e*g4Ys*w14m + 0.5*0.5* g2e*g1Ys*w21m + (0.5+kD)*0.5*(r)* g2e*g3Ys*w23m + (0.5+kD)*0.5* g3e*g1Ys*w31m + (0.5+kD)*0.5*(r)* g3e*g2Ys*w32m + (0.5+kD)*0.5*(1-r)* g4e*g1Ys*w41m)
    g2Ys_pr<- (1/wmbar)*(0.5*0.5* g1e*g2Ys*w12m + (0.5+kD)*0.5*(r)* g1e*g4Ys*w14m + 0.5*0.5* g2e*g1Ys*w21m + 0.5*1* g2e*g2Ys*w22m + (0.5+kD)*0.5*(1-r)* g2e*g3Ys*w23m + (0.5+kD)*0.5* g2e*g4Ys*w24m + (0.5+kD)*0.5*(1-r)* g3e*g2Ys*w32m + (0.5+kD)*0.5*(r)* g4e*g1Ys*w41m + (0.5+kD)*0.5* g4e*g2Ys*w42m)
    g3Ys_pr<- (1/wmbar)*((0.5+kD)*0.5* g1e*g3Ys*w13m + (0.5+kD)*0.5*(r)* g1e*g4Ys*w14m + (0.5+kD)*0.5*(1-r)* g2e*g3Ys*w23m + (0.5+kD)*0.5* g3e*g1Ys*w31m + (0.5+kD)*0.5*(1-r)* g3e*g2Ys*w32m + (0.5+kD)*1* g3e*g3Ys*w33m + (0.5+kD)*0.5* g3e*g4Ys*w34m + (0.5+kD)*0.5*(r)* g4e*g1Ys*w41m + (0.5+kD)*0.5* g4e*g3Ys*w43m)
    g4Ys_pr<- (1/wmbar)*((0.5+kD)*0.5*(1-r)* g1e*g4Ys*w14m + (0.5+kD)*0.5*(r)* g2e*g3Ys*w23m + (0.5+kD)*0.5* g2e*g4Ys*w24m + (0.5+kD)*0.5*(r)* g3e*g2Ys*w32m + (0.5+kD)*0.5* g3e*g4Ys*w34m + (0.5+kD)*0.5*(1-r)* g4e*g1Ys*w41m + (0.5+kD)*0.5* g4e*g2Ys*w42m + (0.5+kD)*0.5* g4e*g3Ys*w43m + (0.5+kD)*1* g4e*g4Ys*w44m)
    
    g1e <- g1e_pr 
    g2e <- g2e_pr 
    g3e <- g3e_pr 
    g4e <- g4e_pr 
    
    g1Xs <- g1Xs_pr 
    g2Xs <- g2Xs_pr
    g3Xs <- g3Xs_pr
    g4Xs <- g4Xs_pr
    g1Ys <- g1Ys_pr 
    g2Ys <- g2Ys_pr
    g3Ys <- g3Ys_pr
    g4Ys <- g4Ys_pr
    
    SexRatioVector[d]<-(g1Ys + g2Ys + g3Ys + g4Ys)
  }
  return(SexRatioVector[d])
}

  svector<-c(50)
  kDvector<-c(50)
  z1vector<-c(2500)
  z2vector<-c(2500)
  z3vector<-c(2500)
  z4vector<-c(2500)
  z5vector<-c(2500)
  z6vector<-c(2500)
  z7vector<-c(2500)
  z8vector<-c(2500)
  z9vector<-c(2500)
  i<-1
  
  for (a in 1:50){
    s <- 0.02*a
    t<-s
    for (b in 1:50){
      kD <- 0.01*b
      svector[[i]]<-s
      kDvector[[i]]<-kD
      z1vector[i]<- SA_function(1000,kD,s,t,0.5,0.1,0.1)
      z2vector[i]<- SA_function(1000,kD,s,t,0.5,0.1,0.01)
      z3vector[i]<- SA_function(1000,kD,s,t,0.5,0.1,0.001)
      z4vector[i]<- SA_function(1000,kD,s,t,0.5,0.01,0.1)
      z5vector[i]<- SA_function(1000,kD,s,t,0.5,0.01,0.01)
      z6vector[i]<- SA_function(1000,kD,s,t,0.5,0.01,0.001)
      z7vector[i]<- SA_function(1000,kD,s,t,0.5,0.001,0.1)
      z8vector[i]<- SA_function(1000,kD,s,t,0.5,0.001,0.01)
      z9vector[i]<- SA_function(1000,kD,s,t,0.5,0.001,0.001)
      i<-i+1

    }
    print(a)
  }



d <- tibble(kD = kDvector, s = svector, z1vector, z2vector, z3vector, z4vector, z5vector, 
            z6vector, z7vector, z8vector, z9vector)
write.table(d, file = "data_S1.txt", row.names = F, col.names = T)


d2 <- d %>% gather(z1vector, z2vector, z3vector, z4vector, z5vector, z6vector, z7vector, z8vector, z9vector, 
                   key = "Vector", value = "SexRatio")



d3 <- d2 %>% mutate(z = as.numeric(str_sub(Vector, 2,2)),
              Xcat = ifelse(z %in% c(1,2,3), "italic(r)==0.1", 
                            ifelse(z %in% c(4,5,6), "italic(r)==0.01", 
                                  "italic(r)==0.001")),
              Ycat = ifelse(z %in% c(1,4,7),"italic(z)==0.1", 
                            ifelse(z %in% c(2,5,8), "italic(z)==0.01", 
                                   "italic(z)==0.001"))) %>% 
  mutate(Ycat = factor(Ycat, c("italic(z)==0.1", "italic(z)==0.01", "italic(z)==0.001")))

  
p <- d3 %>% ggplot(aes(kD, s, z = SexRatio)) + 
  geom_contour_filled(color = "black", bins=11) + 
  scale_fill_viridis(discrete = T, option = "mako", labels = seq(0.5, 0.68, 0.02)) +
  facet_grid(Ycat ~ Xcat, labeller = label_parsed) + 
  scale_x_continuous(limits = c(min(d3$kD), max(d3$kD)), expand = c(0,0)) + 
  scale_y_continuous(limits = c(min(d3$s), max(d3$s)), expand = c(0,0)) + 
  theme(strip.background = element_rect(fill = rgb(0.1, 0.1, 0.1),
                                        color = NULL),
        strip.text = element_text(color = "white"),
        legend.key.size = unit(0.5, "cm"),
        legend.key.height = unit(0.35, "cm")) +
  guides(fill = guide_legend(reverse = TRUE)) + 
  labs(x = bquote(italic(trans)*" drive strength ("*italic(k[D])*")"), 
       y = bquote("Selection strength ("*italic('s,t')*")"),
       fill = "Sex ratio")

pdf(file = "2024_12_30_Figure_S1.pdf", width = 7, height = 6)
p
dev.off()

png(file = "2024_12_30_Figure_S1.png", width = 7, height = 6, res = 1000, units = "in")
p
dev.off()