#1= D1A1 (neutral, neutral)
#2= D1A2 (neutral, cosegregates with Y)
#3= D2A1 (driver, neutral)
#4= D2A2 (driver, cosegregates with Y)
library(tidyverse)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#parameters
G <- 1000 #this is the number of generations 
kD <- 0.0 #proportion males
c <- 0.75  #cosegregation rate 
hA <- 0.0 #dominance of SA locus #if relevant 
s <- 0.01 #fitness cost of cosegregator allele 
z <- 0.00 #fitness cost of the driver 
r <- 0.000 #recombination 
l<- 1 #fitness of ones with A2D2 

DfreqVector<-numeric(G)
AfreqVector<-numeric(G)
SexRatioVector<-numeric(G)
DzVector<-numeric(G)
DA2D2Vector<-numeric(G)
DA2YVector<-numeric(G)
DD2YVector<-numeric(G)


XFreq <-numeric(G)
YFreq <-numeric(G)
g1eFreq<-numeric(G)
g2eFreq<-numeric(G)
g3eFreq<-numeric(G)
g4eFreq<-numeric(G)
g1XsFreq<-numeric(G)
g2XsFreq<-numeric(G)
g3XsFreq<-numeric(G)
g4XsFreq<-numeric(G)
g1YsFreq<-numeric(G)
g2YsFreq<-numeric(G)
g3YsFreq<-numeric(G)
g4YsFreq<-numeric(G)

#u is the trans drive locus (1 is neutral, 2 drives)
u11f <- 1
u12f <- 1
u21f <- 1
u22f <- 1

u11m <- 1
u12m <- 1-z
u21m <- 1-z
u22m <- 1-z

# v is the cosegregation locus 
v11f <- 1
v12f <- 1
v21f <- 1
v22f <- 1

v11m <- 1
v12m <- 1-s
v21m <- 1-s
v22m <- 1-s 

w11m <- u11m*v11m
w12m <- u11m*v12m
w13m <- u12m*v11m
w14m <- u12m*v12m #l
w21m <- u11m*v21m
w22m <- u11m*v22m
w23m <- u12m*v21m
w24m <- u12m*v22m #l
w31m <- u21m*v11m
w32m <- u21m*v12m
w33m <- u22m*v11m
w34m <- u22m*v12m #l
w41m <- u21m*v21m #l
w42m <- u21m*v22m #l
w43m <- u22m*v21m #l
w44m <- u22m*v22m #l

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

# for loop
for (d in 1:500){

wfbar<-((g1e*g1Xs*w11f)+(g1e*g2Xs*w12f)+(g1e*g3Xs*w13f)+(g1e*g4Xs*w14f)+(g2e*g1Xs*w21f)+(g2e*g2Xs*w22f)+(g2e*g3Xs*w23f)+(g2e*g4Xs*w24f)+(g3e*g1Xs*w31f)+(g3e*g2Xs*w32f)+(g3e*g3Xs*w33f)+(g3e*g4Xs*w34f)+(g4e*g1Xs*w41f)+(g4e*g2Xs*w42f)+(g4e*g3Xs*w43f)+(g4e*g4Xs*w44f))
g1e_pr<- (1/wfbar)*(1*g1e*g1Xs*w11f + 0.5* g1e*g2Xs*w12f + 0.5* g1e*g3Xs*w13f + 0.5*(1-r)* g1e*g4Xs*w14f + 0.5* g2e*g1Xs*w21f + 0.5*(r)* g2e*g3Xs*w23f + 0.5* g3e*g1Xs*w31f + 0.5*(r)* g3e*g2Xs*w32f + 0.5*(1-r)* g4e*g1Xs*w41f) 
g2e_pr<- (1/wfbar)*(0.5* g1e*g2Xs*w12f + 0.5*(r)* g1e*g4Xs*w14f + 0.5* g2e*g1Xs*w21f + 1* g2e*g2Xs*w22f + 0.5*(1-r)* g2e*g3Xs*w23f + 0.5* g2e*g4Xs*w24f + 0.5*(1-r)* g3e*g2Xs*w32f + 0.5*(r)* g4e*g1Xs*w41f + 0.5* g4e*g2Xs*w42f)
g3e_pr<- (1/wfbar)*(0.5* g1e*g3Xs*w13f + 0.5*(r)* g1e*g4Xs*w14f + 0.5*(1-r)* g2e*g3Xs*w23f + 0.5* g3e*g1Xs*w31f + 0.5*(1-r)* g3e*g2Xs*w32f + 1* g3e*g3Xs*w33f + 0.5* g3e*g4Xs*w34f + 0.5*(r)* g4e*g1Xs*w41f + 0.5* g4e*g3Xs*w43f)
g4e_pr<- (1/wfbar)*(0.5*(1-r)* g1e*g4Xs*w14f + 0.5*(r)* g2e*g3Xs*w23f + 0.5* g2e*g4Xs*w24f + 0.5*(r)* g3e*g2Xs*w32f + 0.5* g3e*g4Xs*w34f + 0.5*(1-r)* g4e*g1Xs*w41f + 0.5* g4e*g2Xs*w42f + 0.5* g4e*g3Xs*w43f + 1* g4e*g4Xs*w44f)

wmbar<-((g1e*g1Ys*w11m)+(g1e*g2Ys*w12m)+(g1e*g3Ys*w13m)+(g1e*g4Ys*w14m)+(g2e*g1Ys*w21m)+(g2e*g2Ys*w22m)+(g2e*g3Ys*w23m)+(g2e*g4Ys*w24m)+(g3e*g1Ys*w31m)+(g3e*g2Ys*w32m)+(g3e*g3Ys*w33m)+(g3e*g4Ys*w34m)+(g4e*g1Ys*w41m)+(g4e*g2Ys*w42m)+(g4e*g3Ys*w43m)+(g4e*g4Ys*w44m))
g1Xs_pr<- (1/wmbar)*(0.5*1* g1e*g1Ys*w11m + c*0.5* g1e*g2Ys*w12m + (0.5-kD)*0.5* g1e*g3Ys*w13m + c*(0.5-kD)*(1-r)*g1e*g4Ys*w14m + c*0.5*g2e*g1Ys*w21m + c*(0.5-kD)*(r)*g2e*g3Ys*w23m + (0.5-kD)*0.5* g3e*g1Ys*w31m + c*(0.5-kD)*(r)*g3e*g2Ys*w32m + c*(0.5-kD)*(1-r)* g4e*g1Ys*w41m)
g2Xs_pr<- (1/wmbar)*((1-c)*0.5* g1e*g2Ys*w12m + (1-c)*(0.5-kD)*(r)* g1e*g4Ys*w14m + (1-c)*0.5* g2e*g1Ys*w21m + 0.5*1* g2e*g2Ys*w22m + (1-c)*(0.5-kD)*(1-r)* g2e*g3Ys*w23m + (0.5-kD)*0.5* g2e*g4Ys*w24m + (1-c)*(0.5-kD)*(1-r)* g3e*g2Ys*w32m + (1-c)*(0.5-kD)*(r)* g4e*g1Ys*w41m + (0.5-kD)*0.5* g4e*g2Ys*w42m)
g3Xs_pr<- (1/wmbar)*((0.5-kD)*0.5* g1e*g3Ys*w13m + c*(0.5-kD)*(r)*g1e*g4Ys*w14m + c*(0.5-kD)*(1-r)* g2e*g3Ys*w23m + (0.5-kD)*0.5* g3e*g1Ys*w31m + c*(0.5-kD)*(1-r)* g3e*g2Ys*w32m + (0.5-kD)*1* g3e*g3Ys*w33m + c*(0.5-kD)* g3e*g4Ys*w34m + c*(0.5-kD)*(r)* g4e*g1Ys*w41m + c*(0.5-kD)*g4e*g3Ys*w43m)
g4Xs_pr<- (1/wmbar)*((1-c)*(0.5-kD)*(1-r)* g1e*g4Ys*w14m + (1-c)*(0.5-kD)*(r)* g2e*g3Ys*w23m + (0.5-kD)*0.5* g2e*g4Ys*w24m + (1-c)*(0.5-kD)*(r)* g3e*g2Ys*w32m + (1-c)*(0.5-kD)* g3e*g4Ys*w34m + (1-c)*(0.5-kD)*(1-r)* g4e*g1Ys*w41m + (0.5-kD)*0.5* g4e*g2Ys*w42m +(1-c)* (0.5-kD)*g4e*g3Ys*w43m + (0.5-kD)*1* g4e*g4Ys*w44m)

g1Ys_pr<- (1/wmbar)*(0.5*1* g1e*g1Ys*w11m + (1-c)*0.5* g1e*g2Ys*w12m + (0.5+kD)*0.5* g1e*g3Ys*w13m + (1-c)*(0.5+kD)*(1-r)*g1e*g4Ys*w14m + (1-c)*0.5*g2e*g1Ys*w21m + (1-c)*(0.5+kD)*(r)*g2e*g3Ys*w23m + (0.5+kD)*0.5* g3e*g1Ys*w31m + (1-c)*(0.5+kD)*(r)* g3e*g2Ys*w32m + (1-c)*(0.5+kD)*(1-r)* g4e*g1Ys*w41m)
g2Ys_pr<- (1/wmbar)*(c*0.5*g1e*g2Ys*w12m + c*(0.5+kD)*(r)* g1e*g4Ys*w14m + c*0.5*g2e*g1Ys*w21m + 0.5*1* g2e*g2Ys*w22m + c*(0.5+kD)*(1-r)* g2e*g3Ys*w23m + (0.5+kD)*0.5* g2e*g4Ys*w24m + c*(0.5+kD)*(1-r)* g3e*g2Ys*w32m + c*(0.5+kD)*(r)* g4e*g1Ys*w41m + (0.5+kD)*0.5* g4e*g2Ys*w42m)
g3Ys_pr<- (1/wmbar)*((0.5+kD)*0.5* g1e*g3Ys*w13m + (1-c)*(0.5+kD)*(r)* g1e*g4Ys*w14m + (1-c)*(0.5+kD)*(1-r)* g2e*g3Ys*w23m + (0.5+kD)*0.5* g3e*g1Ys*w31m + (1-c)*(0.5+kD)*(1-r)* g3e*g2Ys*w32m + (0.5+kD)*1* g3e*g3Ys*w33m + (1-c)*(0.5+kD)*g3e*g4Ys*w34m + (1-c)*(0.5+kD)*(r)* g4e*g1Ys*w41m + (1-c)*(0.5+kD)*g4e*g3Ys*w43m)
g4Ys_pr<- (1/wmbar)*(c*(0.5+kD)*(1-r)*g1e*g4Ys*w14m + c*(0.5+kD)*(r)* g2e*g3Ys*w23m + (0.5+kD)*0.5* g2e*g4Ys*w24m + c*(0.5+kD)*(r)* g3e*g2Ys*w32m + c*(0.5+kD)*g3e*g4Ys*w34m + c*(0.5+kD)*(1-r)* g4e*g1Ys*w41m + (0.5+kD)*0.5* g4e*g2Ys*w42m + c*(0.5+kD)*g4e*g3Ys*w43m + (0.5+kD)*1* g4e*g4Ys*w44m)

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

#RECURSION 
g1e <- g1e_pr 
g2e <- g2e_pr 
g3e <- g3e_pr 
g4e <- g4e_pr 

g1Xs <- g1Xs_pr 
g2Xs <- g2Xs_pr
g3Xs <- g3Xs_pr
g4Xs <- g4Xs_pr


g3Ys <- g1Ys_pr*0.001
g4Ys <- g2Ys_pr*0.001
g1Ys <- g1Ys_pr-g1Ys_pr*0.001
g2Ys <- g2Ys_pr-g2Ys_pr*0.001


G <- 1000 #this is the number of generations 
kD <- 0.25 #proportion males
c <- 0.75  #cosegregation rate 
#h <- 0.5 #dominance of SA locus #if relevant 
s <- 0.01 #fitness cost of cosegregator allele 
z <- 0.01 #fitness cost of the driver 
r <- 0.001 #recombination 
l<- 1 ###fitness of ones with A2D2 

DfreqVector<-numeric(G)
AfreqVector<-numeric(G)
SexRatioVector<-numeric(G)
DA2D2Vector<-numeric(G)
DA2YVector<-numeric(G)
DD2YVector<-numeric(G)

#u is the trans drive locus (1 is neutral, 2 drives)
u11f <- 1
u12f <- 1
u21f <- 1
u22f <- 1

u11m <- 1
u12m <- 1-z
u21m <- 1-z
u22m <- 1-z

##WHAT IMPACT WOULD THIS HAVE ON FITNESS OF THE COSEGREGATION LOCUS? 
# v is the cosegregation locus (1 is female beneficial, 2 is male beneficial) 
v11f <- 1
v12f <- 1
v21f <- 1
v22f <- 1

v11m <- 1
v12m <- 1-s
v21m <- 1-s
v22m <- 1-s 

w11m <- u11m*v11m
w12m <- u11m*v12m
w13m <- u12m*v11m
w14m <- u12m*v12m #l
w21m <- u11m*v21m
w22m <- u11m*v22m
w23m <- u12m*v21m
w24m <- u12m*v22m #l
w31m <- u21m*v11m
w32m <- u21m*v12m
w33m <- u22m*v11m
w34m <- u22m*v12m #l
w41m <- u21m*v21m #l
w42m <- u21m*v22m #l
w43m <- u22m*v21m #l
w44m <- u22m*v22m #l

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

SexRatioVector[1]<-g1Ys+g2Ys+g3Ys+g4Ys
DfreqVector[1]<-g3Xs+g4Xs+g3Ys+g4Ys
AfreqVector[1]<-g2Xs+g4Xs+g2Ys+g4Ys

# for loop

for (d in 1:G){
  
  wfbar<-((g1e*g1Xs*w11f)+(g1e*g2Xs*w12f)+(g1e*g3Xs*w13f)+(g1e*g4Xs*w14f)+(g2e*g1Xs*w21f)+(g2e*g2Xs*w22f)+(g2e*g3Xs*w23f)+(g2e*g4Xs*w24f)+(g3e*g1Xs*w31f)+(g3e*g2Xs*w32f)+(g3e*g3Xs*w33f)+(g3e*g4Xs*w34f)+(g4e*g1Xs*w41f)+(g4e*g2Xs*w42f)+(g4e*g3Xs*w43f)+(g4e*g4Xs*w44f))
  g1e_pr<- (1/wfbar)*(1*g1e*g1Xs*w11f + 0.5* g1e*g2Xs*w12f + 0.5* g1e*g3Xs*w13f + 0.5*(1-r)* g1e*g4Xs*w14f + 0.5* g2e*g1Xs*w21f + 0.5*(r)* g2e*g3Xs*w23f + 0.5* g3e*g1Xs*w31f + 0.5*(r)* g3e*g2Xs*w32f + 0.5*(1-r)* g4e*g1Xs*w41f) 
  g2e_pr<- (1/wfbar)*(0.5* g1e*g2Xs*w12f + 0.5*(r)* g1e*g4Xs*w14f + 0.5* g2e*g1Xs*w21f + 1* g2e*g2Xs*w22f + 0.5*(1-r)* g2e*g3Xs*w23f + 0.5* g2e*g4Xs*w24f + 0.5*(1-r)* g3e*g2Xs*w32f + 0.5*(r)* g4e*g1Xs*w41f + 0.5* g4e*g2Xs*w42f)
  g3e_pr<- (1/wfbar)*(0.5* g1e*g3Xs*w13f + 0.5*(r)* g1e*g4Xs*w14f + 0.5*(1-r)* g2e*g3Xs*w23f + 0.5* g3e*g1Xs*w31f + 0.5*(1-r)* g3e*g2Xs*w32f + 1* g3e*g3Xs*w33f + 0.5* g3e*g4Xs*w34f + 0.5*(r)* g4e*g1Xs*w41f + 0.5* g4e*g3Xs*w43f)
  g4e_pr<- (1/wfbar)*(0.5*(1-r)* g1e*g4Xs*w14f + 0.5*(r)* g2e*g3Xs*w23f + 0.5* g2e*g4Xs*w24f + 0.5*(r)* g3e*g2Xs*w32f + 0.5* g3e*g4Xs*w34f + 0.5*(1-r)* g4e*g1Xs*w41f + 0.5* g4e*g2Xs*w42f + 0.5* g4e*g3Xs*w43f + 1* g4e*g4Xs*w44f)
  
  wmbar<-((g1e*g1Ys*w11m)+(g1e*g2Ys*w12m)+(g1e*g3Ys*w13m)+(g1e*g4Ys*w14m)+(g2e*g1Ys*w21m)+(g2e*g2Ys*w22m)+(g2e*g3Ys*w23m)+(g2e*g4Ys*w24m)+(g3e*g1Ys*w31m)+(g3e*g2Ys*w32m)+(g3e*g3Ys*w33m)+(g3e*g4Ys*w34m)+(g4e*g1Ys*w41m)+(g4e*g2Ys*w42m)+(g4e*g3Ys*w43m)+(g4e*g4Ys*w44m))
  g1Xs_pr<- (1/wmbar)*(0.5*1* g1e*g1Ys*w11m + c*0.5* g1e*g2Ys*w12m + (0.5-kD)*0.5* g1e*g3Ys*w13m + c*(0.5-kD)*(1-r)*g1e*g4Ys*w14m + c*0.5*g2e*g1Ys*w21m + c*(0.5-kD)*(r)*g2e*g3Ys*w23m + (0.5-kD)*0.5* g3e*g1Ys*w31m + c*(0.5-kD)*(r)*g3e*g2Ys*w32m + c*(0.5-kD)*(1-r)* g4e*g1Ys*w41m)
  g2Xs_pr<- (1/wmbar)*((1-c)*0.5* g1e*g2Ys*w12m + (1-c)*(0.5-kD)*(r)* g1e*g4Ys*w14m + (1-c)*0.5* g2e*g1Ys*w21m + 0.5*1* g2e*g2Ys*w22m + (1-c)*(0.5-kD)*(1-r)* g2e*g3Ys*w23m + (0.5-kD)*0.5* g2e*g4Ys*w24m + (1-c)*(0.5-kD)*(1-r)* g3e*g2Ys*w32m + (1-c)*(0.5-kD)*(r)* g4e*g1Ys*w41m + (0.5-kD)*0.5* g4e*g2Ys*w42m)
  g3Xs_pr<- (1/wmbar)*((0.5-kD)*0.5* g1e*g3Ys*w13m + c*(0.5-kD)*(r)*g1e*g4Ys*w14m + c*(0.5-kD)*(1-r)* g2e*g3Ys*w23m + (0.5-kD)*0.5* g3e*g1Ys*w31m + c*(0.5-kD)*(1-r)* g3e*g2Ys*w32m + (0.5-kD)*1* g3e*g3Ys*w33m + c*(0.5-kD)* g3e*g4Ys*w34m + c*(0.5-kD)*(r)* g4e*g1Ys*w41m + c*(0.5-kD)*g4e*g3Ys*w43m)
  g4Xs_pr<- (1/wmbar)*((1-c)*(0.5-kD)*(1-r)* g1e*g4Ys*w14m + (1-c)*(0.5-kD)*(r)* g2e*g3Ys*w23m + (0.5-kD)*0.5* g2e*g4Ys*w24m + (1-c)*(0.5-kD)*(r)* g3e*g2Ys*w32m + (1-c)*(0.5-kD)* g3e*g4Ys*w34m + (1-c)*(0.5-kD)*(1-r)* g4e*g1Ys*w41m + (0.5-kD)*0.5* g4e*g2Ys*w42m +(1-c)* (0.5-kD)*g4e*g3Ys*w43m + (0.5-kD)*1* g4e*g4Ys*w44m)
  
  g1Ys_pr<- (1/wmbar)*(0.5*1* g1e*g1Ys*w11m + (1-c)*0.5* g1e*g2Ys*w12m + (0.5+kD)*0.5* g1e*g3Ys*w13m + (1-c)*(0.5+kD)*(1-r)*g1e*g4Ys*w14m + (1-c)*0.5*g2e*g1Ys*w21m + (1-c)*(0.5+kD)*(r)*g2e*g3Ys*w23m + (0.5+kD)*0.5* g3e*g1Ys*w31m + (1-c)*(0.5+kD)*(r)* g3e*g2Ys*w32m + (1-c)*(0.5+kD)*(1-r)* g4e*g1Ys*w41m)
  g2Ys_pr<- (1/wmbar)*(c*0.5*g1e*g2Ys*w12m + c*(0.5+kD)*(r)* g1e*g4Ys*w14m + c*0.5* g2e*g1Ys*w21m + 0.5*1* g2e*g2Ys*w22m + c*(0.5+kD)*(1-r)* g2e*g3Ys*w23m + (0.5+kD)*0.5* g2e*g4Ys*w24m + c*(0.5+kD)*(1-r)* g3e*g2Ys*w32m + c*(0.5+kD)*(r)* g4e*g1Ys*w41m + (0.5+kD)*0.5* g4e*g2Ys*w42m)
  g3Ys_pr<- (1/wmbar)*((0.5+kD)*0.5* g1e*g3Ys*w13m + (1-c)*(0.5+kD)*(r)* g1e*g4Ys*w14m + (1-c)*(0.5+kD)*(1-r)* g2e*g3Ys*w23m + (0.5+kD)*0.5* g3e*g1Ys*w31m + (1-c)*(0.5+kD)*(1-r)* g3e*g2Ys*w32m + (0.5+kD)*1* g3e*g3Ys*w33m + (1-c)*(0.5+kD)*g3e*g4Ys*w34m + (1-c)*(0.5+kD)*(r)* g4e*g1Ys*w41m + (1-c)*(0.5+kD)*g4e*g3Ys*w43m)
  g4Ys_pr<- (1/wmbar)*(c*(0.5+kD)*(1-r)*g1e*g4Ys*w14m + c*(0.5+kD)*(r)* g2e*g3Ys*w23m + (0.5+kD)*0.5* g2e*g4Ys*w24m + c*(0.5+kD)*(r)* g3e*g2Ys*w32m + c*(0.5+kD)*g3e*g4Ys*w34m + c*(0.5+kD)*(1-r)* g4e*g1Ys*w41m + (0.5+kD)*0.5* g4e*g2Ys*w42m + c*(0.5+kD)*g4e*g3Ys*w43m + (0.5+kD)*1* g4e*g4Ys*w44m)
  
  
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


  SexRatioVector[d]<-g1Ys+g2Ys+g3Ys+g4Ys
  DfreqVector[d]<-g3Xs+g4Xs+g3Ys+g4Ys
  AfreqVector[d]<-g2Xs+g4Xs+g2Ys+g4Ys
  
 
  DA2D2Vector[d]<-((g1Xs+g1Ys)*(g4Ys+g4Xs))-((g2Xs+g2Ys)*(g3Xs+g3Ys))
  DA2YVector[d]<- (g2Ys+g4Ys)*(g1Xs+g3Xs)-(g1Ys+g3Ys)*(g2Xs+g4Xs)
  DD2YVector[d]<-((g3Ys+g4Ys)*(g1Xs+g2Xs))-((g3Xs+g4Xs)*(g1Ys+g2Ys))
  YFreq[d]<-g1Ys+g2Ys+g3Ys+g4Ys
  XFreq[d]<-g1Xs+g2Xs+g3Xs+g4Xs
  
  g1eFreq[d]<-g1e
  g2eFreq[d]<-g2e
  g3eFreq[d]<-g3e
  g4eFreq[d]<-g4e
  g1XsFreq[d]<-g1Xs
  g2XsFreq[d]<-g2Xs
  g3XsFreq[d]<-g3Xs
  g4XsFreq[d]<-g4Xs
  g1YsFreq[d]<-g1Ys
  g2YsFreq[d]<-g2Ys
  g3YsFreq[d]<-g3Ys
  g4YsFreq[d]<-g4Ys
  
  DA2D2Vector[d]<-((g1Xs+g1Ys)*(g4Ys+g4Xs))-((g2Xs+g2Ys)*(g3Xs+g3Ys))
  DA2YVector[d]<- (g2Ys+g4Ys)*(g1Xs+g3Xs)-(g1Ys+g3Ys)*(g2Xs+g4Xs)
  DD2YVector[d]<-((g3Ys+g4Ys)*(g1Xs+g2Xs))-((g3Xs+g4Xs)*(g1Ys+g2Ys))
}

layout(matrix(c(1, 3, 2, 3, 4, 5), 3, 2, byrow = TRUE))
x <- 1:G
y1 <- XFreq
y2 <- YFreq
y3 <- AfreqVector
y4 <- DfreqVector
y5 <- 0.5



d <- tibble(G = x, X = y1, Y = y2, A = y3, D = y4)
write.table(x = d, file = "data_S4A.txt", row.names = F, col.names = T)

# #### LINKAGE DISEQUILIBRIUM
# x <- 1:G
y1 <- DA2D2Vector
y2 <- DA2YVector
y3 <- DD2YVector
y4 <- 0.0


d <- tibble(G = x, DA2D2 = y1, DA2Y = y2, DD2Y = y3)
write.table(x = d, file = "data_S4B.txt", row.names = F, col.names = T)

d <- tibble(G = x, 
            A1D1X = g1XsFreq,
            A2D1X = g2XsFreq,
            A1D1Y = g1YsFreq,
            A2D1Y = g2YsFreq, 
            A1D2X = g3XsFreq,
            A2D2X = g4XsFreq,
            A1D2Y = g3YsFreq,
            A2D2Y = g4YsFreq)
write.table(x = d, file = "data_S4C.txt", row.names = F, col.names = T)


CoFunction<-function(G,kD,c,z,r,l){
  
  original_kD <- kD
  original_z <- z
  z<- 0 
  kD <-0 
  
  DfreqVector<-numeric(G)
  AfreqVector<-numeric(G)
  SexRatioVector<-numeric(G)
  DA2D2Vector<-numeric(G)
  DA2YVector<-numeric(G)
  DD2YVector<-numeric(G)
  
  u11f <- 1
  u12f <- 1
  u21f <- 1
  u22f <- 1
  
  u11m <- 1
  u12m <- 1-z
  u21m <- 1-z
  u22m <- 1-z
  
  # v is the cosegregation locus (1 is female beneficial, 2 is male beneficial) 
  v11f <- 1
  v12f <- 1
  v21f <- 1
  v22f <- 1
  
  v11m <- 1
  v12m <- 1-s
  v21m <- 1-s
  v22m <- 1-s 
  
  w11m <- u11m*v11m
  w12m <- u11m*v12m
  w13m <- u12m*v11m
  w14m <- l#u12m*v12m
  w21m <- u11m*v21m
  w22m <- u11m*v22m
  w23m <- u12m*v21m
  w24m <- l #u12m*v22m
  w31m <- u21m*v11m
  w32m <- u21m*v12m
  w33m <- u22m*v11m
  w34m <- l#u22m*v12m
  w41m <- l #u21m*v21m
  w42m <- l#u21m*v22m
  w43m <- l#u22m*v21m
  w44m <- l#u22m*v22m
  
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
  
  # for loop
  for (d in 1:500){
    
    wfbar<-((g1e*g1Xs*w11f)+(g1e*g2Xs*w12f)+(g1e*g3Xs*w13f)+(g1e*g4Xs*w14f)+(g2e*g1Xs*w21f)+(g2e*g2Xs*w22f)+(g2e*g3Xs*w23f)+(g2e*g4Xs*w24f)+(g3e*g1Xs*w31f)+(g3e*g2Xs*w32f)+(g3e*g3Xs*w33f)+(g3e*g4Xs*w34f)+(g4e*g1Xs*w41f)+(g4e*g2Xs*w42f)+(g4e*g3Xs*w43f)+(g4e*g4Xs*w44f))
    g1e_pr<- (1/wfbar)*(1*g1e*g1Xs*w11f + 0.5* g1e*g2Xs*w12f + 0.5* g1e*g3Xs*w13f + 0.5*(1-r)* g1e*g4Xs*w14f + 0.5* g2e*g1Xs*w21f + 0.5*(r)* g2e*g3Xs*w23f + 0.5* g3e*g1Xs*w31f + 0.5*(r)* g3e*g2Xs*w32f + 0.5*(1-r)* g4e*g1Xs*w41f) 
    g2e_pr<- (1/wfbar)*(0.5* g1e*g2Xs*w12f + 0.5*(r)* g1e*g4Xs*w14f + 0.5* g2e*g1Xs*w21f + 1* g2e*g2Xs*w22f + 0.5*(1-r)* g2e*g3Xs*w23f + 0.5* g2e*g4Xs*w24f + 0.5*(1-r)* g3e*g2Xs*w32f + 0.5*(r)* g4e*g1Xs*w41f + 0.5* g4e*g2Xs*w42f)
    g3e_pr<- (1/wfbar)*(0.5* g1e*g3Xs*w13f + 0.5*(r)* g1e*g4Xs*w14f + 0.5*(1-r)* g2e*g3Xs*w23f + 0.5* g3e*g1Xs*w31f + 0.5*(1-r)* g3e*g2Xs*w32f + 1* g3e*g3Xs*w33f + 0.5* g3e*g4Xs*w34f + 0.5*(r)* g4e*g1Xs*w41f + 0.5* g4e*g3Xs*w43f)
    g4e_pr<- (1/wfbar)*(0.5*(1-r)* g1e*g4Xs*w14f + 0.5*(r)* g2e*g3Xs*w23f + 0.5* g2e*g4Xs*w24f + 0.5*(r)* g3e*g2Xs*w32f + 0.5* g3e*g4Xs*w34f + 0.5*(1-r)* g4e*g1Xs*w41f + 0.5* g4e*g2Xs*w42f + 0.5* g4e*g3Xs*w43f + 1* g4e*g4Xs*w44f)
    
    wmbar<-((g1e*g1Ys*w11m)+(g1e*g2Ys*w12m)+(g1e*g3Ys*w13m)+(g1e*g4Ys*w14m)+(g2e*g1Ys*w21m)+(g2e*g2Ys*w22m)+(g2e*g3Ys*w23m)+(g2e*g4Ys*w24m)+(g3e*g1Ys*w31m)+(g3e*g2Ys*w32m)+(g3e*g3Ys*w33m)+(g3e*g4Ys*w34m)+(g4e*g1Ys*w41m)+(g4e*g2Ys*w42m)+(g4e*g3Ys*w43m)+(g4e*g4Ys*w44m))
    g1Xs_pr<- (1/wmbar)*(0.5*1* g1e*g1Ys*w11m + c*0.5* g1e*g2Ys*w12m + (0.5-kD)*0.5* g1e*g3Ys*w13m + c*(0.5-kD)*(1-r)*g1e*g4Ys*w14m + c*0.5*g2e*g1Ys*w21m + c*(0.5-kD)*(r)*g2e*g3Ys*w23m + (0.5-kD)*0.5* g3e*g1Ys*w31m + c*(0.5-kD)*(r)*g3e*g2Ys*w32m + c*(0.5-kD)*(1-r)* g4e*g1Ys*w41m)
    g2Xs_pr<- (1/wmbar)*((1-c)*0.5* g1e*g2Ys*w12m + (1-c)*(0.5-kD)*(r)* g1e*g4Ys*w14m + (1-c)*0.5* g2e*g1Ys*w21m + 0.5*1* g2e*g2Ys*w22m + (1-c)*(0.5-kD)*(1-r)* g2e*g3Ys*w23m + (0.5-kD)*0.5* g2e*g4Ys*w24m + (1-c)*(0.5-kD)*(1-r)* g3e*g2Ys*w32m + (1-c)*(0.5-kD)*(r)* g4e*g1Ys*w41m + (0.5-kD)*0.5* g4e*g2Ys*w42m)
    g3Xs_pr<- (1/wmbar)*((0.5-kD)*0.5* g1e*g3Ys*w13m + c*(0.5-kD)*(r)*g1e*g4Ys*w14m + c*(0.5-kD)*(1-r)* g2e*g3Ys*w23m + (0.5-kD)*0.5* g3e*g1Ys*w31m + c*(0.5-kD)*(1-r)* g3e*g2Ys*w32m + (0.5-kD)*1* g3e*g3Ys*w33m + c*(0.5-kD)* g3e*g4Ys*w34m + c*(0.5-kD)*(r)* g4e*g1Ys*w41m + c*(0.5-kD)*g4e*g3Ys*w43m)
    g4Xs_pr<- (1/wmbar)*((1-c)*(0.5-kD)*(1-r)* g1e*g4Ys*w14m + (1-c)*(0.5-kD)*(r)* g2e*g3Ys*w23m + (0.5-kD)*0.5* g2e*g4Ys*w24m + (1-c)*(0.5-kD)*(r)* g3e*g2Ys*w32m + (1-c)*(0.5-kD)* g3e*g4Ys*w34m + (1-c)*(0.5-kD)*(1-r)* g4e*g1Ys*w41m + (0.5-kD)*0.5* g4e*g2Ys*w42m +(1-c)* (0.5-kD)*g4e*g3Ys*w43m + (0.5-kD)*1* g4e*g4Ys*w44m)
    
    g1Ys_pr<- (1/wmbar)*(0.5*1* g1e*g1Ys*w11m + (1-c)*0.5* g1e*g2Ys*w12m + (0.5+kD)*0.5* g1e*g3Ys*w13m + (1-c)*(0.5+kD)*(1-r)*g1e*g4Ys*w14m + (1-c)*0.5*g2e*g1Ys*w21m + (1-c)*(0.5+kD)*(r)*g2e*g3Ys*w23m + (0.5+kD)*0.5* g3e*g1Ys*w31m + (1-c)*(0.5+kD)*(r)* g3e*g2Ys*w32m + (1-c)*(0.5+kD)*(1-r)* g4e*g1Ys*w41m)
    g2Ys_pr<- (1/wmbar)*(c*0.5*g1e*g2Ys*w12m + c*(0.5+kD)*(r)* g1e*g4Ys*w14m + c*0.5*g2e*g1Ys*w21m + 0.5*1* g2e*g2Ys*w22m + c*(0.5+kD)*(1-r)* g2e*g3Ys*w23m + (0.5+kD)*0.5* g2e*g4Ys*w24m + c*(0.5+kD)*(1-r)* g3e*g2Ys*w32m + c*(0.5+kD)*(r)* g4e*g1Ys*w41m + (0.5+kD)*0.5* g4e*g2Ys*w42m)
    g3Ys_pr<- (1/wmbar)*((0.5+kD)*0.5* g1e*g3Ys*w13m + (1-c)*(0.5+kD)*(r)* g1e*g4Ys*w14m + (1-c)*(0.5+kD)*(1-r)* g2e*g3Ys*w23m + (0.5+kD)*0.5* g3e*g1Ys*w31m + (1-c)*(0.5+kD)*(1-r)* g3e*g2Ys*w32m + (0.5+kD)*1* g3e*g3Ys*w33m + (1-c)*(0.5+kD)*g3e*g4Ys*w34m + (1-c)*(0.5+kD)*(r)* g4e*g1Ys*w41m + (1-c)*(0.5+kD)*g4e*g3Ys*w43m)
    g4Ys_pr<- (1/wmbar)*(c*(0.5+kD)*(1-r)*g1e*g4Ys*w14m + c*(0.5+kD)*(r)* g2e*g3Ys*w23m + (0.5+kD)*0.5* g2e*g4Ys*w24m + c*(0.5+kD)*(r)* g3e*g2Ys*w32m + c*(0.5+kD)*g3e*g4Ys*w34m + c*(0.5+kD)*(1-r)* g4e*g1Ys*w41m + (0.5+kD)*0.5* g4e*g2Ys*w42m + c*(0.5+kD)*g4e*g3Ys*w43m + (0.5+kD)*1* g4e*g4Ys*w44m)
    
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
  
  #RECURSION 
  
  
  kD <- original_kD
  z <- original_z
  g1e <- g1e_pr 
  g2e <- g2e_pr 
  g3e <- g3e_pr 
  g4e <- g4e_pr 
  
  g1Xs <- g1Xs_pr 
  g2Xs <- g2Xs_pr
  g3Xs <- g3Xs_pr
  g4Xs <- g4Xs_pr
  
  
  g3Ys <- g1Ys_pr*0.001
  g4Ys <- g2Ys_pr*0.001
  g1Ys <- g1Ys_pr-g1Ys_pr*0.001
  g2Ys <- g2Ys_pr-g2Ys_pr*0.001
  
  u11f <- 1
  u12f <- 1
  u21f <- 1
  u22f <- 1
  
  u11m <- 1
  u12m <- 1-z
  u21m <- 1-z
  u22m <- 1-z

  # v is the cosegregation locus (1 is female beneficial, 2 is male beneficial) 
  v11f <- 1
  v12f <- 1
  v21f <- 1
  v22f <- 1
  
  v11m <- 1
  v12m <- 1-s
  v21m <- 1-s
  v22m <- 1-s 
  
  w11m <- u11m*v11m
  w12m <- u11m*v12m
  w13m <- u12m*v11m
  w14m <- l#u12m*v12m
  w21m <- u11m*v21m
  w22m <- u11m*v22m
  w23m <- u12m*v21m
  w24m <- l #u12m*v22m
  w31m <- u21m*v11m
  w32m <- u21m*v12m
  w33m <- u22m*v11m
  w34m <- l#u22m*v12m
  w41m <- l #u21m*v21m
  w42m <- l#u21m*v22m
  w43m <- l#u22m*v21m
  w44m <- l#u22m*v22m
  
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
  
  SexRatioVector[1]<-g1Ys+g2Ys+g3Ys+g4Ys
  DfreqVector[1]<-g3Xs+g4Xs+g3Ys+g4Ys
  AfreqVector[1]<-g2Xs+g4Xs+g2Ys+g4Ys
  
  # for loop
  for (d in 1:G){
    #Egg haplotypes (4 total)
    
    wfbar<-((g1e*g1Xs*w11f)+(g1e*g2Xs*w12f)+(g1e*g3Xs*w13f)+(g1e*g4Xs*w14f)+(g2e*g1Xs*w21f)+(g2e*g2Xs*w22f)+(g2e*g3Xs*w23f)+(g2e*g4Xs*w24f)+(g3e*g1Xs*w31f)+(g3e*g2Xs*w32f)+(g3e*g3Xs*w33f)+(g3e*g4Xs*w34f)+(g4e*g1Xs*w41f)+(g4e*g2Xs*w42f)+(g4e*g3Xs*w43f)+(g4e*g4Xs*w44f))
    g1e_pr<- (1/wfbar)*(1*g1e*g1Xs*w11f + 0.5* g1e*g2Xs*w12f + 0.5* g1e*g3Xs*w13f + 0.5*(1-r)* g1e*g4Xs*w14f + 0.5* g2e*g1Xs*w21f + 0.5*(r)* g2e*g3Xs*w23f + 0.5* g3e*g1Xs*w31f + 0.5*(r)* g3e*g2Xs*w32f + 0.5*(1-r)* g4e*g1Xs*w41f) 
    g2e_pr<- (1/wfbar)*(0.5* g1e*g2Xs*w12f + 0.5*(r)* g1e*g4Xs*w14f + 0.5* g2e*g1Xs*w21f + 1* g2e*g2Xs*w22f + 0.5*(1-r)* g2e*g3Xs*w23f + 0.5* g2e*g4Xs*w24f + 0.5*(1-r)* g3e*g2Xs*w32f + 0.5*(r)* g4e*g1Xs*w41f + 0.5* g4e*g2Xs*w42f)
    g3e_pr<- (1/wfbar)*(0.5* g1e*g3Xs*w13f + 0.5*(r)* g1e*g4Xs*w14f + 0.5*(1-r)* g2e*g3Xs*w23f + 0.5* g3e*g1Xs*w31f + 0.5*(1-r)* g3e*g2Xs*w32f + 1* g3e*g3Xs*w33f + 0.5* g3e*g4Xs*w34f + 0.5*(r)* g4e*g1Xs*w41f + 0.5* g4e*g3Xs*w43f)
    g4e_pr<- (1/wfbar)*(0.5*(1-r)* g1e*g4Xs*w14f + 0.5*(r)* g2e*g3Xs*w23f + 0.5* g2e*g4Xs*w24f + 0.5*(r)* g3e*g2Xs*w32f + 0.5* g3e*g4Xs*w34f + 0.5*(1-r)* g4e*g1Xs*w41f + 0.5* g4e*g2Xs*w42f + 0.5* g4e*g3Xs*w43f + 1* g4e*g4Xs*w44f)
    
    wmbar<-((g1e*g1Ys*w11m)+(g1e*g2Ys*w12m)+(g1e*g3Ys*w13m)+(g1e*g4Ys*w14m)+(g2e*g1Ys*w21m)+(g2e*g2Ys*w22m)+(g2e*g3Ys*w23m)+(g2e*g4Ys*w24m)+(g3e*g1Ys*w31m)+(g3e*g2Ys*w32m)+(g3e*g3Ys*w33m)+(g3e*g4Ys*w34m)+(g4e*g1Ys*w41m)+(g4e*g2Ys*w42m)+(g4e*g3Ys*w43m)+(g4e*g4Ys*w44m))
    g1Xs_pr<- (1/wmbar)*(0.5*1* g1e*g1Ys*w11m + c*0.5* g1e*g2Ys*w12m + (0.5-kD)*0.5* g1e*g3Ys*w13m + c*(0.5-kD)*(1-r)*g1e*g4Ys*w14m + c*0.5*g2e*g1Ys*w21m + c*(0.5-kD)*(r)*g2e*g3Ys*w23m + (0.5-kD)*0.5* g3e*g1Ys*w31m + c*(0.5-kD)*(r)*g3e*g2Ys*w32m + c*(0.5-kD)*(1-r)* g4e*g1Ys*w41m)
    g2Xs_pr<- (1/wmbar)*((1-c)*0.5* g1e*g2Ys*w12m + (1-c)*(0.5-kD)*(r)* g1e*g4Ys*w14m + (1-c)*0.5* g2e*g1Ys*w21m + 0.5*1* g2e*g2Ys*w22m + (1-c)*(0.5-kD)*(1-r)* g2e*g3Ys*w23m + (0.5-kD)*0.5* g2e*g4Ys*w24m + (1-c)*(0.5-kD)*(1-r)* g3e*g2Ys*w32m + (1-c)*(0.5-kD)*(r)* g4e*g1Ys*w41m + (0.5-kD)*0.5* g4e*g2Ys*w42m)
    g3Xs_pr<- (1/wmbar)*((0.5-kD)*0.5* g1e*g3Ys*w13m + c*(0.5-kD)*(r)*g1e*g4Ys*w14m + c*(0.5-kD)*(1-r)* g2e*g3Ys*w23m + (0.5-kD)*0.5* g3e*g1Ys*w31m + c*(0.5-kD)*(1-r)* g3e*g2Ys*w32m + (0.5-kD)*1* g3e*g3Ys*w33m + c*(0.5-kD)* g3e*g4Ys*w34m + c*(0.5-kD)*(r)* g4e*g1Ys*w41m + c*(0.5-kD)*g4e*g3Ys*w43m)
    g4Xs_pr<- (1/wmbar)*((1-c)*(0.5-kD)*(1-r)* g1e*g4Ys*w14m + (1-c)*(0.5-kD)*(r)* g2e*g3Ys*w23m + (0.5-kD)*0.5* g2e*g4Ys*w24m + (1-c)*(0.5-kD)*(r)* g3e*g2Ys*w32m + (1-c)*(0.5-kD)* g3e*g4Ys*w34m + (1-c)*(0.5-kD)*(1-r)* g4e*g1Ys*w41m + (0.5-kD)*0.5* g4e*g2Ys*w42m +(1-c)* (0.5-kD)*g4e*g3Ys*w43m + (0.5-kD)*1* g4e*g4Ys*w44m)
    
    g1Ys_pr<- (1/wmbar)*(0.5*1* g1e*g1Ys*w11m + (1-c)*0.5* g1e*g2Ys*w12m + (0.5+kD)*0.5* g1e*g3Ys*w13m + (1-c)*(0.5+kD)*(1-r)*g1e*g4Ys*w14m + (1-c)*0.5*g2e*g1Ys*w21m + (1-c)*(0.5+kD)*(r)*g2e*g3Ys*w23m + (0.5+kD)*0.5* g3e*g1Ys*w31m + (1-c)*(0.5+kD)*(r)* g3e*g2Ys*w32m + (1-c)*(0.5+kD)*(1-r)* g4e*g1Ys*w41m)
    g2Ys_pr<- (1/wmbar)*(c*0.5*g1e*g2Ys*w12m + c*(0.5+kD)*(r)* g1e*g4Ys*w14m + c*0.5*g2e*g1Ys*w21m + 0.5*1* g2e*g2Ys*w22m + c*(0.5+kD)*(1-r)* g2e*g3Ys*w23m + (0.5+kD)*0.5* g2e*g4Ys*w24m + c*(0.5+kD)*(1-r)* g3e*g2Ys*w32m + c*(0.5+kD)*(r)* g4e*g1Ys*w41m + (0.5+kD)*0.5* g4e*g2Ys*w42m)
    g3Ys_pr<- (1/wmbar)*((0.5+kD)*0.5* g1e*g3Ys*w13m + (1-c)*(0.5+kD)*(r)* g1e*g4Ys*w14m + (1-c)*(0.5+kD)*(1-r)* g2e*g3Ys*w23m + (0.5+kD)*0.5* g3e*g1Ys*w31m + (1-c)*(0.5+kD)*(1-r)* g3e*g2Ys*w32m + (0.5+kD)*1* g3e*g3Ys*w33m + (1-c)*(0.5+kD)*g3e*g4Ys*w34m + (1-c)*(0.5+kD)*(r)* g4e*g1Ys*w41m + (1-c)*(0.5+kD)*g4e*g3Ys*w43m)
    g4Ys_pr<- (1/wmbar)*(c*(0.5+kD)*(1-r)*g1e*g4Ys*w14m + c*(0.5+kD)*(r)* g2e*g3Ys*w23m + (0.5+kD)*0.5* g2e*g4Ys*w24m + c*(0.5+kD)*(r)* g3e*g2Ys*w32m + c*(0.5+kD)*g3e*g4Ys*w34m + c*(0.5+kD)*(1-r)* g4e*g1Ys*w41m + (0.5+kD)*0.5* g4e*g2Ys*w42m + c*(0.5+kD)*g4e*g3Ys*w43m + (0.5+kD)*1* g4e*g4Ys*w44m)
    
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
    
    DA2D2Vector[d]<-((g1Xs+g1Ys)*(g4Ys+g4Xs))-((g2Xs+g2Ys)*(g3Xs+g3Ys))
    DA2YVector[d]<- (g2Ys+g4Ys)*(g1Xs+g3Xs)-(g1Ys+g3Ys)*(g2Xs+g4Xs)
    DD2YVector[d]<-((g3Ys+g4Ys)*(g1Xs+g2Xs))-((g3Xs+g4Xs)*(g1Ys+g2Ys))
    
    SexRatioVector[d]<-g1Ys+g2Ys+g3Ys+g4Ys
    DfreqVector[d]<-g3Xs+g4Xs+g3Ys+g4Ys
    AfreqVector[d]<-g2Xs+g4Xs+g2Ys+g4Ys
  }
  return(c(DfreqVector[d],SexRatioVector[d]))
}

#A_function <- function(G, k, s, t, h, f, z, r) {
cvector<-c(50)
kDvector<-c(50)
SRvector<-c(50)
DFvector<-c(50)

i<-1


i<-1 
n <- 100
for (a in 1:n){
  c <- 0.5+a/(2*n)
  for (b in 1:n){
    kD <- b/(2*n)
    t <- CoFunction(1000,kD, c, 0.01, 0.001, 1)

    DF<- t[1]
    SR<- t[2]
    # DF<-CoFunction(1000,kD,c,0.01,0.001,1)
    # print(DF)
    cvector[i]<-c
    kDvector[i]<-kD
    SRvector[i]<-SR
    DFvector[i]<-DF
    i<-i+1
  }
  print(a)
}


DFvector[1]
DF <- rep(NA, length(DFvector))
SR <- rep(NA, length(SRvector))
for(i in 1:length(DF))
{
  DF[i] <- DFvector[[i]]
  SR[i] <- SRvector[[i]]
}

library(tidyverse)
d <- tibble(k = kDvector, c = cvector, DF = DFvector)
write.table(x = d, file = "data_S4D.txt", row.names = F, col.names = T)

library(tidyverse)
d <- tibble(k = kDvector, c = cvector, SR = SRvector)
write.table(x = d, file = "data_S4E.txt", row.names = F, col.names = T)

