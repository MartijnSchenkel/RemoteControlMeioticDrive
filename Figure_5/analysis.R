# EVOLUTION of XO from XY; using a male-limited  drive system 
#1= D1A1 (X-shredder, male beneficial)
#2= D1A2 (X-shredder, female beneficial)
#3= D2A1 (neutral, male beneficial)
#4= D2A2 (neutral, female beneficial)

#parameters
G <- 500
z <- 0.00  # Fitness cost of the trans driver
s <- 0.50  # Fitness cost of cis driver
hA <- 0.0   # Dominance of cis driver
kD <- 0.0   # Skew in favor of Y chromosome (trans drive) in males
kA <- 0.25   # Skew in favor of A2 (cis drive) in males
r <- 0.001 # Recombination rate
w <-0

###BURN IN 
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
v12f <- 1-hA*s
v21f <- 1-hA*s
v22f <- 1-s

v11m <- 1
v12m <- 1-hA*s
v21m <- 1-hA*s
v22m <- 1-s

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

w11Ym <- u11m*v11m
w12Ym <- u11m*v12m
w13Ym <- u12m*v11m
w14Ym <- u12m*v12m
w21Ym <- u11m*v21m
w22Ym <- u11m*v22m
w23Ym <- u12m*v21m
w24Ym <- u12m*v22m
w31Ym <- u21m*v11m
w32Ym <- u21m*v12m
w33Ym <- u22m*v11m
w34Ym <- u22m*v12m
w41Ym <- u21m*v21m
w42Ym <- u21m*v22m
w43Ym <- u22m*v21m
w44Ym <- u22m*v22m

w11Om <- u11m*v11m*(1-w)
w12Om <- u11m*v12m*(1-w)
w13Om <- u12m*v11m*(1-w)
w14Om <- u12m*v12m*(1-w)
w21Om <- u11m*v21m*(1-w)
w22Om <- u11m*v22m*(1-w)
w23Om <- u12m*v21m*(1-w)
w24Om <- u12m*v22m*(1-w)
w31Om <- u21m*v11m*(1-w)
w32Om <- u21m*v12m*(1-w)
w33Om <- u22m*v11m*(1-w)
w34Om <- u22m*v12m*(1-w)
w41Om <- u21m*v21m*(1-w)
w42Om <- u21m*v22m*(1-w)
w43Om <- u22m*v21m*(1-w)
w44Om <- u22m*v22m*(1-w)

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

g3Xs <- 0.0
g3Ys <- 0.0
g4Xs <- 0.0
g4Ys <- 0.0

g1Os <- 0.0
g2Os <- 0.0
g3Os <- 0.0
g4Os <- 0.0

# in the process of changin "M" to either "Om" or "Xm" for male fitness (so that we can differentiate between these and add fitness costs!)

for (d in 1:500){
  
  wfbar<-(g1e*g1Xs*w11f)+(g1e*g2Xs*w12f)+(g1e*g3Xs*w13f)+(g1e*g4Xs*w14f)+(g2e*g1Xs*w21f)+(g2e*g2Xs*w22f)+(g2e*g3Xs*w23f)+(g2e*g4Xs*w24f)+(g3e*g1Xs*w31f)+(g3e*g2Xs*w32f)+(g3e*g3Xs*w33f)+(g3e*g4Xs*w34f)+(g4e*g1Xs*w41f)+(g4e*g2Xs*w42f)+(g4e*g3Xs*w43f)+(g4e*g4Xs*w44f)
  
  g1e_pr<- (1/wfbar)*(1*g1e*g1Xs*w11f + 0.5* g1e*g2Xs*w12f + 0.5* g1e*g3Xs*w13f + 0.5*(1-r)* g1e*g4Xs*w14f + 0.5* g2e*g1Xs*w21f + 0.5*(r)* g2e*g3Xs*w23f + 0.5* g3e*g1Xs*w31f + 0.5*(r)* g3e*g2Xs*w32f + 0.5*(1-r)* g4e*g1Xs*w41f)
  g2e_pr<- (1/wfbar)*(0.5* g1e*g2Xs*w12f + 0.5*(r)* g1e*g4Xs*w14f + 0.5* g2e*g1Xs*w21f + 1* g2e*g2Xs*w22f + 0.5*(1-r)* g2e*g3Xs*w23f + 0.5* g2e*g4Xs*w24f + 0.5*(1-r)* g3e*g2Xs*w32f + 0.5*(r)* g4e*g1Xs*w41f + 0.5* g4e*g2Xs*w42f) 
  g3e_pr<- (1/wfbar)*(0.5* g1e*g3Xs*w13f + 0.5*(r)* g1e*g4Xs*w14f + 0.5*(1-r)* g2e*g3Xs*w23f + 0.5* g3e*g1Xs*w31f + 0.5*(1-r)* g3e*g2Xs*w32f + 1* g3e*g3Xs*w33f + 0.5* g3e*g4Xs*w34f + 0.5*(r)* g4e*g1Xs*w41f + 0.5* g4e*g3Xs*w43f) 
  g4e_pr<- (1/wfbar)*(0.5*(1-r)* g1e*g4Xs*w14f + 0.5*(r)* g2e*g3Xs*w23f + 0.5* g2e*g4Xs*w24f + 0.5*(r)* g3e*g2Xs*w32f + 0.5* g3e*g4Xs*w34f + 0.5*(1-r)* g4e*g1Xs*w41f + 0.5* g4e*g2Xs*w42f + 0.5* g4e*g3Xs*w43f + 1* g4e*g4Xs*w44f) 
  
  wmbar<-(g1e*g1Ys*w11Ym)+(g1e*g2Ys*w12Ym)+(g1e*g3Ys*w13Ym)+(g1e*g4Ys*w14Ym)+(g2e*g1Ys*w21Ym)+(g2e*g2Ys*w22Ym)+(g2e*g3Ys*w23Ym)+(g2e*g4Ys*w24Ym)+(g3e*g1Ys*w31Ym)+(g3e*g2Ys*w32Ym)+(g3e*g3Ys*w33Ym)+(g3e*g4Ys*w34Ym)+(g4e*g1Ys*w41Ym)+(g4e*g2Ys*w42Ym)+(g4e*g3Ys*w43Ym)+(g4e*g4Ys*w44Ym)+(g1e*g1Os*w11Om)+(g1e*g2Os*w12Om)+(g1e*g3Os*w13Om)+(g1e*g4Os*w14Om)+(g2e*g1Os*w21Om)+(g2e*g2Os*w22Om)+(g2e*g3Os*w23Om)+(g2e*g4Os*w24Om)+(g3e*g1Os*w31Om)+(g3e*g2Os*w32Om)+(g3e*g3Os*w33Om)+(g3e*g4Os*w34Om)+(g4e*g1Os*w41Om)+(g4e*g2Os*w42Om)+(g4e*g3Os*w43Om)+(g4e*g4Os*w44Om)
  g1Xs_pr<- (1/wmbar)*(0.5*1*g1e*g1Ys*w11Ym + 0.5*(0.5-kA)*g1e*g2Ys*w12Ym + (0.5-kD)*0.5*g1e*g3Ys*w13Ym + (0.5-kD)*(0.5-kA)*(1-r)*g1e*g4Ys*w14Ym + 0.5*(0.5-kA)*g2e*g1Ys*w21Ym + (0.5-kD)*(0.5-kA)*(r)*g2e*g3Ys*w23Ym + (0.5-kD)*0.5*g3e*g1Ys*w31Ym + (0.5-kD)*(0.5-kA)*(r)*g3e*g2Ys*w32Ym + (0.5-kD)*(0.5-kA)*(1-r)*g4e*g1Ys*w41Ym +0.5*1*g1e*g1Os*w11Om + 0.5*(0.5-kA)*g1e*g2Os*w12Om + (0.5-kD)*0.5*g1e*g3Os*w13Om + (0.5-kD)*(0.5-kA)*(1-r)*g1e*g4Os*w14Om + 0.5*(0.5-kA)*g2e*g1Os*w21Om + (0.5-kD)*(0.5-kA)*(r)*g2e*g3Os*w23Om + (0.5-kD)*0.5*g3e*g1Os*w31Om + (0.5-kD)*(0.5-kA)*(r)*g3e*g2Os*w32Om + (0.5-kD)*(0.5-kA)*(1-r)*g4e*g1Os*w41Om)
  g2Xs_pr<- (1/wmbar)*(0.5*(0.5+kA)*g1e*g2Ys*w12Ym + (0.5-kD)*(0.5+kA)*(r)*g1e*g4Ys*w14Ym + 0.5*(0.5+kA)*g2e*g1Ys*w21Ym + 0.5*1*g2e*g2Ys*w22Ym + (0.5-kD)*(0.5+kA)*(1-r)*g2e*g3Ys*w23Ym + (0.5-kD)*0.5*g2e*g4Ys*w24Ym + (0.5-kD)*(0.5+kA)*(1-r)*g3e*g2Ys*w32Ym + (0.5-kD)*(0.5+kA)*(r)*g4e*g1Ys*w41Ym + (0.5-kD)*0.5*g4e*g2Ys*w42Ym + 0.5*(0.5+kA)*g1e*g2Os*w12Om + (0.5-kD)*(0.5+kA)*(r)*g1e*g4Os*w14Om + (0.5)*(0.5+kA)*g2e*g1Os*w21Om + 0.5*1*g2e*g2Os*w22Om + (0.5-kD)*(0.5+kA)*(1-r)*g2e*g3Os*w23Om + (0.5-kD)*0.5*g2e*g4Os*w24Om + (0.5-kD)*(0.5+kA)*(1-r)*g3e*g2Os*w32Om + (0.5-kD)*(0.5+kA)*(r)*g4e*g1Os*w41Om + (0.5-kD)*0.5*g4e*g2Os*w42Om)
  g3Xs_pr<- (1/wmbar)*((0.5-kD)*0.5*g1e*g3Ys*w13Ym + (0.5-kD)*(0.5-kA)*(r)*g1e*g4Ys*w14Ym + (0.5-kD)*(0.5-kA)*(1-r)*g2e*g3Ys*w23Ym + (0.5-kD)*0.5*g3e*g1Ys*w31Ym + (0.5-kD)*(0.5-kA)*(1-r)*g3e*g2Ys*w32Ym + (0.5-kD)*1*g3e*g3Ys*w33Ym + (0.5-kD)*(0.5-kA)*g3e*g4Ys*w34Ym + (0.5-kD)*(0.5-kA)*(r)*g4e*g1Ys*w41Ym + (0.5-kD)*(0.5-kA)*g4e*g3Ys*w43Ym+(0.5-kD)*0.5*g1e*g3Os*w13Om + (0.5-kD)*(0.5-kA)*(r)*g1e*g4Os*w14Om + (0.5-kD)*(0.5-kA)*(1-r)*g2e*g3Os*w23Om + (0.5-kD)*0.5*g3e*g1Os*w31Om + (0.5-kD)*(0.5-kA)*(1-r)*g3e*g2Os*w32Om + (0.5-kD)*1*g3e*g3Os*w33Om + (0.5-kD)*(0.5-kA)*g3e*g4Os*w34Om + (0.5-kD)*(0.5-kA)*(r)*g4e*g1Os*w41Om + (0.5-kD)*(0.5-kA)*g4e*g3Os*w43Om)
  g4Xs_pr<- (1/wmbar)*((0.5-kD)*(0.5+kA)*(1-r)*g1e*g4Ys*w14Ym + (0.5-kD)*(0.5+kA)*(r)*g2e*g3Ys*w23Ym + (0.5-kD)*0.5*g2e*g4Ys*w24Ym + (0.5-kD)*(0.5+kA)*(r)*g3e*g2Ys*w32Ym + (0.5-kD)*(0.5+kA)*g3e*g4Ys*w34Ym + (0.5-kD)*(0.5+kA)*(1-r)*g4e*g1Ys*w41Ym + (0.5-kD)*0.5*g4e*g2Ys*w42Ym + (0.5-kD)*(0.5+kA)*g4e*g3Ys*w43Ym + (0.5-kD)*1* g4e*g4Ys*w44Ym + (0.5-kD)*(0.5+kA)*(1-r)*g1e*g4Os*w14Om + (0.5-kD)*(0.5+kA)*(r)*g2e*g3Os*w23Om + (0.5-kD)*0.5*g2e*g4Os*w24Om + (0.5-kD)*(0.5+kA)*(r)*g3e*g2Os*w32Om + (0.5-kD)*(0.5+kA)*g3e*g4Os*w34Om + (0.5-kD)*(0.5+kA)*(1-r)*g4e*g1Os*w41Om + (0.5-kD)*0.5*g4e*g2Os*w42Om + (0.5-kD)*(0.5+kA)*g4e*g3Os*w43Om + (0.5-kD)*1*g4e*g4Os*w44Om)
  
  g1Ys_pr<- (1/wmbar)*(0.5*1* g1e*g1Ys*w11Ym + 0.5*(0.5-kA)* g1e*g2Ys*w12Ym + 0.5*0.5* g1e*g3Ys*w13Ym + 0.5*(0.5-kA)*(1-r)*g1e*g4Ys*w14Ym + 0.5*(0.5-kA)* g2e*g1Ys*w21Ym + 0.5*(0.5-kA)*(r)* g2e*g3Ys*w23Ym + 0.5*0.5* g3e*g1Ys*w31Ym + 0.5*(0.5-kA)*(r)* g3e*g2Ys*w32Ym + 0.5*(0.5-kA)*(1-r)* g4e*g1Ys*w41Ym)
  g2Ys_pr<- (1/wmbar)*(0.5*(0.5+kA)* g1e*g2Ys*w12Ym + 0.5*(0.5+kA)*(r)* g1e*g4Ys*w14Ym + 0.5*(0.5+kA)* g2e*g1Ys*w21Ym + 0.5*1* g2e*g2Ys*w22Ym + 0.5*(0.5+kA)*(1-r)* g2e*g3Ys*w23Ym + 0.5*0.5* g2e*g4Ys*w24Ym + 0.5*(0.5+kA)*(1-r)* g3e*g2Ys*w32Ym + 0.5*(0.5+kA)*(r)* g4e*g1Ys*w41Ym + 0.5*0.5* g4e*g2Ys*w42Ym)
  g3Ys_pr<- (1/wmbar)*(0.5*0.5* g1e*g3Ys*w13Ym + 0.5*(0.5-kA)*(r)* g1e*g4Ys*w14Ym + 0.5*(0.5-kA)*(1-r)* g2e*g3Ys*w23Ym + 0.5*0.5* g3e*g1Ys*w31Ym + 0.5*(0.5-kA)*(1-r)* g3e*g2Ys*w32Ym + 0.5*1* g3e*g3Ys*w33Ym + 0.5*(0.5-kA)* g3e*g4Ys*w34Ym + 0.5*(0.5-kA)*(r)* g4e*g1Ys*w41Ym + 0.5*(0.5-kA)* g4e*g3Ys*w43Ym)
  g4Ys_pr<- (1/wmbar)*(0.5*(0.5+kA)*(1-r)* g1e*g4Ys*w14Ym + 0.5*(0.5+kA)*(r)* g2e*g3Ys*w23Ym + 0.5*0.5* g2e*g4Ys*w24Ym + 0.5*(0.5+kA)*(r)* g3e*g2Ys*w32Ym + 0.5*(0.5+kA)* g3e*g4Ys*w34Ym + 0.5*(0.5+kA)*(1-r)* g4e*g1Ys*w41Ym + 0.5*0.5* g4e*g2Ys*w42Ym + 0.5*(0.5+kA)* g4e*g3Ys*w43Ym + 0.5*1* g4e*g4Ys*w44Ym)
  g1Os_pr<- (1/wmbar)*(0*1*g1e*g1Ys*w11Ym + 0*(0.5-kA)* g1e*g2Ys*w12Ym + kD*0.5*g1e*g3Ys*w13Ym + kD*(0.5-kA)*(1-r)*g1e*g4Ys*w14Ym + 0*(0.5-kA)* g2e*g1Ys*w21Ym + kD*(0.5-kA)*(r)* g2e*g3Ys*w23Ym + kD*0.5* g3e*g1Ys*w31Ym + kD*(0.5-kA)*(r)*g3e*g2Ys*w32Ym + kD*(0.5-kA)*(1-r)*g4e*g1Ys*w41Ym +0.5*1*g1e*g1Os*w11Om + 0.5*(0.5-kA)* g1e*g2Os*w12Om + (0.5+kD)*0.5* g1e*g3Os*w13Om + (0.5+kD)*(0.5-kA)*(1-r)* g1e*g4Os*w14Om + 0.5*(0.5-kA)*g2e*g1Os*w21Om + (0.5+kD)*(0.5-kA)*(r)* g2e*g3Os*w23Om + (0.5+kD)*0.5* g3e*g1Os*w31Om + (0.5+kD)*(0.5-kA)*(r)* g3e*g2Os*w32Om + (0.5+kD)*(0.5-kA)*(1-r)* g4e*g1Os*w41Om)
  g2Os_pr<- (1/wmbar)*(0*(0.5+kA)*g1e*g2Ys*w12Ym + kD*(0.5+kA)*(r)* g1e*g4Ys*w14Ym + 0*(0.5+kA)* g2e*g1Ys*w21Ym + 0*1* g2e*g2Ys*w22Ym + kD*(0.5+kA)*(1-r)* g2e*g3Ys*w23Ym + kD*0.5* g2e*g4Ys*w24Ym + kD*(0.5+kA)*(1-r)* g3e*g2Ys*w32Ym + kD*(0.5+kA)*(r)* g4e*g1Ys*w41Ym + kD*0.5* g4e*g2Ys*w42Ym + 0.5*(0.5+kA)* g1e*g2Os*w12Om + (0.5+kD)*(0.5+kA)*(r)* g1e*g4Os*w14Om + 0.5*(0.5+kA)* g2e*g1Os*w21Om + 0.5*1*g2e*g2Os*w22Om + (0.5+kD)*(0.5+kA)*(1-r)* g2e*g3Os*w23Om + (0.5+kD)*0.5* g2e*g4Os*w24Om + (0.5+kD)*(0.5+kA)*(1-r)* g3e*g2Os*w32Om + (0.5+kD)*(0.5+kA)*(r)* g4e*g1Os*w41Om + (0.5+kD)*0.5* g4e*g2Os*w42Om)
  g3Os_pr<- (1/wmbar)*(kD*0.5*g1e*g3Ys*w13Ym + kD*(0.5-kA)*(r)*g1e*g4Ys*w14Ym + kD*(0.5-kA)*(1-r)* g2e*g3Ys*w23Ym + kD*0.5* g3e*g1Ys*w31Ym + kD*(0.5-kA)*(1-r)*g3e*g2Ys*w32Ym + kD*1* g3e*g3Ys*w33Ym + kD*(0.5-kA)* g3e*g4Ys*w34Ym + kD*(0.5-kA)*(r)* g4e*g1Ys*w41Ym + kD*(0.5-kA)* g4e*g3Ys*w43Ym + (0.5+kD)*0.5* g1e*g3Os*w13Om + (0.5+kD)*(0.5-kA)*(r)* g1e*g4Os*w14Om + (0.5+kD)*(0.5-kA)*(1-r)* g2e*g3Os*w23Om + (0.5+kD)*0.5* g3e*g1Os*w31Om + (0.5+kD)*(0.5-kA)*(1-r)* g3e*g2Os*w32Om + (0.5+kD)*1* g3e*g3Os*w33Om + (0.5+kD)*(0.5-kA)* g3e*g4Os*w34Om + (0.5+kD)*(0.5-kA)*(r)* g4e*g1Os*w41Om + (0.5+kD)*(0.5-kA)* g4e*g3Os*w43Om)
  g4Os_pr<- (1/wmbar)*(kD*(0.5+kA)*(1-r)*g1e*g4Ys*w14Ym + kD*(0.5+kA)*(r)*g2e*g3Ys*w23Ym + kD*0.5*g2e*g4Ys*w24Ym + kD*(0.5+kA)*(r)*g3e*g2Ys*w32Ym + kD*(0.5+kA)*g3e*g4Ys*w34Ym + kD*(0.5+kA)*(1-r)*g4e*g1Ys*w41Ym + kD*0.5*g4e*g2Ys*w42Ym + kD*(0.5+kA)*g4e*g3Ys*w43Ym + kD*(0.5+kA)*1*g4e*g4Ys*w44Ym + (0.5+kD)*(0.5+kA)*(1-r)*g1e*g4Os*w14Om + (0.5+kD)*(0.5+kA)*(r)*g2e*g3Os*w23Om + (0.5+kD)*0.5*g2e*g4Os*w24Om + (0.5+kD)*(0.5+kA)*(r)*g3e*g2Os*w32Om + (0.5+kD)*(0.5+kA)*g3e*g4Os*w34Om + (0.5+kD)*(0.5+kA)*(1-r)*g4e*g1Os*w41Om + (0.5+kD)*0.5*g4e*g2Os*w42Om + (0.5+kD)*(0.5+kA)*g4e*g3Os*w43Om + (0.5+kD)*1*g4e*g4Os*w44Om)
  
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

#SIMULATION 
g1e <- g1e
g2e <- g2e
g3e <- g3e
g4e <- g4e

g1Xs <- g1Xs
g2Xs <- g2Xs
g3Xs <- g3Xs
g4Xs <- g4Xs

g1Os <- g1Os
g2Os <- g2Os
g3Os <- g3Os
g4Os <- g4Os

g3Ys <- g1Ys_pr*0.001
g4Ys <- g2Ys_pr*0.001
g1Ys <- g1Ys_pr-g1Ys_pr*0.001
g2Ys <- g2Ys_pr-g2Ys_pr*0.001

G <- 1000 #this is the number of generations
z <- 0.001 #fitness cost of the trans driver
s <- 0.5 #fitness cost of cis driver
hA <- 0 # dominance of cis driver
kD <- 0.25 #skew in favor of Y chromosome (trans driver)
kA <- 0.25 #skew in favor of A2, in males (cis driver)
r <- 0.001 #recombination
w <- 0.001

DfreqVector<-numeric(G)
AfreqVector<-numeric(G)
SexRatioVector<-numeric(G)
OFreq <-numeric(G)
YFreq<-numeric(G)
XFreq<-numeric(G)
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
g1OsFreq<-numeric(G)
g2OsFreq<-numeric(G)
g3OsFreq<-numeric(G)
g4OsFreq<-numeric(G)

DA2D2Vector<-numeric(G)
DA2YVector<-numeric(G)
DD2YVector<-numeric(G)

g1eFreq[0]<-g1e
g2eFreq[0]<-g2e
g3eFreq[0]<-g3e
g4eFreq[0]<-g4e
g1XsFreq[0]<-g1Xs
g2XsFreq[0]<-g2Xs
g3XsFreq[0]<-g3Xs
g4XsFreq[0]<-g4Xs
g1YsFreq[0]<-g1Ys
g2YsFreq[0]<-g2Ys
g3YsFreq[0]<-g3Ys
g4YsFreq[0]<-g4Ys
g1OsFreq[0]<-g1Os
g2OsFreq[0]<-g2Os
g3OsFreq[0]<-g3Os
g4OsFreq[0]<-g4Os

u11f <- 1
u12f <- 1
u21f <- 1
u22f <- 1

u11m <- 1
u12m <- 1-z
u21m <- 1-z
u22m <- 1-z

#v is the cis drive locus (1 is neutral, 2 drives)
v11f <- 1 
v12f <- 1-hA*s
v21f <- 1-hA*s
v22f <- 1-s

v11m <- 1
v12m <- 1-hA*s
v21m <- 1-hA*s
v22m <- 1-s

w11Ym <- u11m*v11m
w12Ym <- u11m*v12m
w13Ym <- u12m*v11m
w14Ym <- u12m*v12m
w21Ym <- u11m*v21m
w22Ym <- u11m*v22m
w23Ym <- u12m*v21m
w24Ym <- u12m*v22m
w31Ym <- u21m*v11m
w32Ym <- u21m*v12m
w33Ym <- u22m*v11m
w34Ym <- u22m*v12m
w41Ym <- u21m*v21m
w42Ym <- u21m*v22m
w43Ym <- u22m*v21m
w44Ym <- u22m*v22m

w11Om <- u11m*v11m*(1-w)
w12Om <- u11m*v12m*(1-w)
w13Om <- u12m*v11m*(1-w)
w14Om <- u12m*v12m*(1-w)
w21Om <- u11m*v21m*(1-w)
w22Om <- u11m*v22m*(1-w)
w23Om <- u12m*v21m*(1-w)
w24Om <- u12m*v22m*(1-w)
w31Om <- u21m*v11m*(1-w)
w32Om <- u21m*v12m*(1-w)
w33Om <- u22m*v11m*(1-w)
w34Om <- u22m*v12m*(1-w)
w41Om <- u21m*v21m*(1-w)
w42Om <- u21m*v22m*(1-w)
w43Om <- u22m*v21m*(1-w)
w44Om <- u22m*v22m*(1-w)

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

# for loop
for (d in 1:G){
  
  wfbar<-(g1e*g1Xs*w11f)+(g1e*g2Xs*w12f)+(g1e*g3Xs*w13f)+(g1e*g4Xs*w14f)+(g2e*g1Xs*w21f)+(g2e*g2Xs*w22f)+(g2e*g3Xs*w23f)+(g2e*g4Xs*w24f)+(g3e*g1Xs*w31f)+(g3e*g2Xs*w32f)+(g3e*g3Xs*w33f)+(g3e*g4Xs*w34f)+(g4e*g1Xs*w41f)+(g4e*g2Xs*w42f)+(g4e*g3Xs*w43f)+(g4e*g4Xs*w44f)
  
  g1e_pr<- (1/wfbar)*(1*g1e*g1Xs*w11f + 0.5* g1e*g2Xs*w12f + 0.5* g1e*g3Xs*w13f + 0.5*(1-r)* g1e*g4Xs*w14f + 0.5* g2e*g1Xs*w21f + 0.5*(r)* g2e*g3Xs*w23f + 0.5* g3e*g1Xs*w31f + 0.5*(r)* g3e*g2Xs*w32f + 0.5*(1-r)* g4e*g1Xs*w41f)
  g2e_pr<- (1/wfbar)*(0.5* g1e*g2Xs*w12f + 0.5*(r)* g1e*g4Xs*w14f + 0.5* g2e*g1Xs*w21f + 1* g2e*g2Xs*w22f + 0.5*(1-r)* g2e*g3Xs*w23f + 0.5* g2e*g4Xs*w24f + 0.5*(1-r)* g3e*g2Xs*w32f + 0.5*(r)* g4e*g1Xs*w41f + 0.5* g4e*g2Xs*w42f) 
  g3e_pr<- (1/wfbar)*(0.5* g1e*g3Xs*w13f + 0.5*(r)* g1e*g4Xs*w14f + 0.5*(1-r)* g2e*g3Xs*w23f + 0.5* g3e*g1Xs*w31f + 0.5*(1-r)* g3e*g2Xs*w32f + 1* g3e*g3Xs*w33f + 0.5* g3e*g4Xs*w34f + 0.5*(r)* g4e*g1Xs*w41f + 0.5* g4e*g3Xs*w43f) 
  g4e_pr<- (1/wfbar)*(0.5*(1-r)* g1e*g4Xs*w14f + 0.5*(r)* g2e*g3Xs*w23f + 0.5* g2e*g4Xs*w24f + 0.5*(r)* g3e*g2Xs*w32f + 0.5* g3e*g4Xs*w34f + 0.5*(1-r)* g4e*g1Xs*w41f + 0.5* g4e*g2Xs*w42f + 0.5* g4e*g3Xs*w43f + 1* g4e*g4Xs*w44f) 
  
  wmbar<-(g1e*g1Ys*w11Ym)+(g1e*g2Ys*w12Ym)+(g1e*g3Ys*w13Ym)+(g1e*g4Ys*w14Ym)+(g2e*g1Ys*w21Ym)+(g2e*g2Ys*w22Ym)+(g2e*g3Ys*w23Ym)+(g2e*g4Ys*w24Ym)+(g3e*g1Ys*w31Ym)+(g3e*g2Ys*w32Ym)+(g3e*g3Ys*w33Ym)+(g3e*g4Ys*w34Ym)+(g4e*g1Ys*w41Ym)+(g4e*g2Ys*w42Ym)+(g4e*g3Ys*w43Ym)+(g4e*g4Ys*w44Ym)+(g1e*g1Os*w11Om)+(g1e*g2Os*w12Om)+(g1e*g3Os*w13Om)+(g1e*g4Os*w14Om)+(g2e*g1Os*w21Om)+(g2e*g2Os*w22Om)+(g2e*g3Os*w23Om)+(g2e*g4Os*w24Om)+(g3e*g1Os*w31Om)+(g3e*g2Os*w32Om)+(g3e*g3Os*w33Om)+(g3e*g4Os*w34Om)+(g4e*g1Os*w41Om)+(g4e*g2Os*w42Om)+(g4e*g3Os*w43Om)+(g4e*g4Os*w44Om)
  g1Xs_pr<- (1/wmbar)*(0.5*1*g1e*g1Ys*w11Ym + 0.5*(0.5-kA)*g1e*g2Ys*w12Ym + (0.5-kD)*0.5*g1e*g3Ys*w13Ym + (0.5-kD)*(0.5-kA)*(1-r)*g1e*g4Ys*w14Ym + 0.5*(0.5-kA)*g2e*g1Ys*w21Ym + (0.5-kD)*(0.5-kA)*(r)*g2e*g3Ys*w23Ym + (0.5-kD)*0.5*g3e*g1Ys*w31Ym + (0.5-kD)*(0.5-kA)*(r)*g3e*g2Ys*w32Ym + (0.5-kD)*(0.5-kA)*(1-r)*g4e*g1Ys*w41Ym +0.5*1*g1e*g1Os*w11Om + 0.5*(0.5-kA)*g1e*g2Os*w12Om + (0.5-kD)*0.5*g1e*g3Os*w13Om + (0.5-kD)*(0.5-kA)*(1-r)*g1e*g4Os*w14Om + 0.5*(0.5-kA)*g2e*g1Os*w21Om + (0.5-kD)*(0.5-kA)*(r)*g2e*g3Os*w23Om + (0.5-kD)*0.5*g3e*g1Os*w31Om + (0.5-kD)*(0.5-kA)*(r)*g3e*g2Os*w32Om + (0.5-kD)*(0.5-kA)*(1-r)*g4e*g1Os*w41Om)
  g2Xs_pr<- (1/wmbar)*(0.5*(0.5+kA)*g1e*g2Ys*w12Ym + (0.5-kD)*(0.5+kA)*(r)*g1e*g4Ys*w14Ym + 0.5*(0.5+kA)*g2e*g1Ys*w21Ym + 0.5*1*g2e*g2Ys*w22Ym + (0.5-kD)*(0.5+kA)*(1-r)*g2e*g3Ys*w23Ym + (0.5-kD)*0.5*g2e*g4Ys*w24Ym + (0.5-kD)*(0.5+kA)*(1-r)*g3e*g2Ys*w32Ym + (0.5-kD)*(0.5+kA)*(r)*g4e*g1Ys*w41Ym + (0.5-kD)*0.5*g4e*g2Ys*w42Ym + 0.5*(0.5+kA)*g1e*g2Os*w12Om + (0.5-kD)*(0.5+kA)*(r)*g1e*g4Os*w14Om + (0.5)*(0.5+kA)*g2e*g1Os*w21Om + 0.5*1*g2e*g2Os*w22Om + (0.5-kD)*(0.5+kA)*(1-r)*g2e*g3Os*w23Om + (0.5-kD)*0.5*g2e*g4Os*w24Om + (0.5-kD)*(0.5+kA)*(1-r)*g3e*g2Os*w32Om + (0.5-kD)*(0.5+kA)*(r)*g4e*g1Os*w41Om + (0.5-kD)*0.5*g4e*g2Os*w42Om)
  g3Xs_pr<- (1/wmbar)*((0.5-kD)*0.5*g1e*g3Ys*w13Ym + (0.5-kD)*(0.5-kA)*(r)*g1e*g4Ys*w14Ym + (0.5-kD)*(0.5-kA)*(1-r)*g2e*g3Ys*w23Ym + (0.5-kD)*0.5*g3e*g1Ys*w31Ym + (0.5-kD)*(0.5-kA)*(1-r)*g3e*g2Ys*w32Ym + (0.5-kD)*1*g3e*g3Ys*w33Ym + (0.5-kD)*(0.5-kA)*g3e*g4Ys*w34Ym + (0.5-kD)*(0.5-kA)*(r)*g4e*g1Ys*w41Ym + (0.5-kD)*(0.5-kA)*g4e*g3Ys*w43Ym+(0.5-kD)*0.5*g1e*g3Os*w13Om + (0.5-kD)*(0.5-kA)*(r)*g1e*g4Os*w14Om + (0.5-kD)*(0.5-kA)*(1-r)*g2e*g3Os*w23Om + (0.5-kD)*0.5*g3e*g1Os*w31Om + (0.5-kD)*(0.5-kA)*(1-r)*g3e*g2Os*w32Om + (0.5-kD)*1*g3e*g3Os*w33Om + (0.5-kD)*(0.5-kA)*g3e*g4Os*w34Om + (0.5-kD)*(0.5-kA)*(r)*g4e*g1Os*w41Om + (0.5-kD)*(0.5-kA)*g4e*g3Os*w43Om)
  g4Xs_pr<- (1/wmbar)*((0.5-kD)*(0.5+kA)*(1-r)*g1e*g4Ys*w14Ym + (0.5-kD)*(0.5+kA)*(r)*g2e*g3Ys*w23Ym + (0.5-kD)*0.5*g2e*g4Ys*w24Ym + (0.5-kD)*(0.5+kA)*(r)*g3e*g2Ys*w32Ym + (0.5-kD)*(0.5+kA)*g3e*g4Ys*w34Ym + (0.5-kD)*(0.5+kA)*(1-r)*g4e*g1Ys*w41Ym + (0.5-kD)*0.5*g4e*g2Ys*w42Ym + (0.5-kD)*(0.5+kA)*g4e*g3Ys*w43Ym + (0.5-kD)*1* g4e*g4Ys*w44Ym + (0.5-kD)*(0.5+kA)*(1-r)*g1e*g4Os*w14Om + (0.5-kD)*(0.5+kA)*(r)*g2e*g3Os*w23Om + (0.5-kD)*0.5*g2e*g4Os*w24Om + (0.5-kD)*(0.5+kA)*(r)*g3e*g2Os*w32Om + (0.5-kD)*(0.5+kA)*g3e*g4Os*w34Om + (0.5-kD)*(0.5+kA)*(1-r)*g4e*g1Os*w41Om + (0.5-kD)*0.5*g4e*g2Os*w42Om + (0.5-kD)*(0.5+kA)*g4e*g3Os*w43Om + (0.5-kD)*1*g4e*g4Os*w44Om)

  g1Ys_pr<- (1/wmbar)*(0.5*1* g1e*g1Ys*w11Ym + 0.5*(0.5-kA)* g1e*g2Ys*w12Ym + 0.5*0.5* g1e*g3Ys*w13Ym + 0.5*(0.5-kA)*(1-r)*g1e*g4Ys*w14Ym + 0.5*(0.5-kA)* g2e*g1Ys*w21Ym + 0.5*(0.5-kA)*(r)* g2e*g3Ys*w23Ym + 0.5*0.5* g3e*g1Ys*w31Ym + 0.5*(0.5-kA)*(r)* g3e*g2Ys*w32Ym + 0.5*(0.5-kA)*(1-r)* g4e*g1Ys*w41Ym)
  g2Ys_pr<- (1/wmbar)*(0.5*(0.5+kA)* g1e*g2Ys*w12Ym + 0.5*(0.5+kA)*(r)* g1e*g4Ys*w14Ym + 0.5*(0.5+kA)* g2e*g1Ys*w21Ym + 0.5*1* g2e*g2Ys*w22Ym + 0.5*(0.5+kA)*(1-r)* g2e*g3Ys*w23Ym + 0.5*0.5* g2e*g4Ys*w24Ym + 0.5*(0.5+kA)*(1-r)* g3e*g2Ys*w32Ym + 0.5*(0.5+kA)*(r)* g4e*g1Ys*w41Ym + 0.5*0.5* g4e*g2Ys*w42Ym)
  g3Ys_pr<- (1/wmbar)*(0.5*0.5* g1e*g3Ys*w13Ym + 0.5*(0.5-kA)*(r)* g1e*g4Ys*w14Ym + 0.5*(0.5-kA)*(1-r)* g2e*g3Ys*w23Ym + 0.5*0.5* g3e*g1Ys*w31Ym + 0.5*(0.5-kA)*(1-r)* g3e*g2Ys*w32Ym + 0.5*1* g3e*g3Ys*w33Ym + 0.5*(0.5-kA)* g3e*g4Ys*w34Ym + 0.5*(0.5-kA)*(r)* g4e*g1Ys*w41Ym + 0.5*(0.5-kA)* g4e*g3Ys*w43Ym)
  g4Ys_pr<- (1/wmbar)*(0.5*(0.5+kA)*(1-r)* g1e*g4Ys*w14Ym + 0.5*(0.5+kA)*(r)* g2e*g3Ys*w23Ym + 0.5*0.5* g2e*g4Ys*w24Ym + 0.5*(0.5+kA)*(r)* g3e*g2Ys*w32Ym + 0.5*(0.5+kA)* g3e*g4Ys*w34Ym + 0.5*(0.5+kA)*(1-r)* g4e*g1Ys*w41Ym + 0.5*0.5* g4e*g2Ys*w42Ym + 0.5*(0.5+kA)* g4e*g3Ys*w43Ym + 0.5*1* g4e*g4Ys*w44Ym)
  g1Os_pr<- (1/wmbar)*(0*1*g1e*g1Ys*w11Ym + 0*(0.5-kA)* g1e*g2Ys*w12Ym + kD*0.5*g1e*g3Ys*w13Ym + kD*(0.5-kA)*(1-r)*g1e*g4Ys*w14Ym + 0*(0.5-kA)* g2e*g1Ys*w21Ym + kD*(0.5-kA)*(r)* g2e*g3Ys*w23Ym + kD*0.5* g3e*g1Ys*w31Ym + kD*(0.5-kA)*(r)*g3e*g2Ys*w32Ym + kD*(0.5-kA)*(1-r)*g4e*g1Ys*w41Ym +0.5*1*g1e*g1Os*w11Om + 0.5*(0.5-kA)* g1e*g2Os*w12Om + (0.5+kD)*0.5* g1e*g3Os*w13Om + (0.5+kD)*(0.5-kA)*(1-r)* g1e*g4Os*w14Om + 0.5*(0.5-kA)*g2e*g1Os*w21Om + (0.5+kD)*(0.5-kA)*(r)* g2e*g3Os*w23Om + (0.5+kD)*0.5* g3e*g1Os*w31Om + (0.5+kD)*(0.5-kA)*(r)* g3e*g2Os*w32Om + (0.5+kD)*(0.5-kA)*(1-r)* g4e*g1Os*w41Om)
  g2Os_pr<- (1/wmbar)*(0*(0.5+kA)*g1e*g2Ys*w12Ym + kD*(0.5+kA)*(r)* g1e*g4Ys*w14Ym + 0*(0.5+kA)* g2e*g1Ys*w21Ym + 0*1* g2e*g2Ys*w22Ym + kD*(0.5+kA)*(1-r)* g2e*g3Ys*w23Ym + kD*0.5* g2e*g4Ys*w24Ym + kD*(0.5+kA)*(1-r)* g3e*g2Ys*w32Ym + kD*(0.5+kA)*(r)* g4e*g1Ys*w41Ym + kD*0.5* g4e*g2Ys*w42Ym + 0.5*(0.5+kA)* g1e*g2Os*w12Om + (0.5+kD)*(0.5+kA)*(r)* g1e*g4Os*w14Om + 0.5*(0.5+kA)* g2e*g1Os*w21Om + 0.5*1*g2e*g2Os*w22Om + (0.5+kD)*(0.5+kA)*(1-r)* g2e*g3Os*w23Om + (0.5+kD)*0.5* g2e*g4Os*w24Om + (0.5+kD)*(0.5+kA)*(1-r)* g3e*g2Os*w32Om + (0.5+kD)*(0.5+kA)*(r)* g4e*g1Os*w41Om + (0.5+kD)*0.5* g4e*g2Os*w42Om)
  g3Os_pr<- (1/wmbar)*(kD*0.5*g1e*g3Ys*w13Ym + kD*(0.5-kA)*(r)*g1e*g4Ys*w14Ym + kD*(0.5-kA)*(1-r)* g2e*g3Ys*w23Ym + kD*0.5* g3e*g1Ys*w31Ym + kD*(0.5-kA)*(1-r)*g3e*g2Ys*w32Ym + kD*1* g3e*g3Ys*w33Ym + kD*(0.5-kA)* g3e*g4Ys*w34Ym + kD*(0.5-kA)*(r)* g4e*g1Ys*w41Ym + kD*(0.5-kA)* g4e*g3Ys*w43Ym + (0.5+kD)*0.5* g1e*g3Os*w13Om + (0.5+kD)*(0.5-kA)*(r)* g1e*g4Os*w14Om + (0.5+kD)*(0.5-kA)*(1-r)* g2e*g3Os*w23Om + (0.5+kD)*0.5* g3e*g1Os*w31Om + (0.5+kD)*(0.5-kA)*(1-r)* g3e*g2Os*w32Om + (0.5+kD)*1* g3e*g3Os*w33Om + (0.5+kD)*(0.5-kA)* g3e*g4Os*w34Om + (0.5+kD)*(0.5-kA)*(r)* g4e*g1Os*w41Om + (0.5+kD)*(0.5-kA)* g4e*g3Os*w43Om)
  g4Os_pr<- (1/wmbar)*(kD*(0.5+kA)*(1-r)*g1e*g4Ys*w14Ym + kD*(0.5+kA)*(r)*g2e*g3Ys*w23Ym + kD*0.5*g2e*g4Ys*w24Ym + kD*(0.5+kA)*(r)*g3e*g2Ys*w32Ym + kD*(0.5+kA)*g3e*g4Ys*w34Ym + kD*(0.5+kA)*(1-r)*g4e*g1Ys*w41Ym + kD*0.5*g4e*g2Ys*w42Ym + kD*(0.5+kA)*g4e*g3Ys*w43Ym + kD*(0.5+kA)*1*g4e*g4Ys*w44Ym + (0.5+kD)*(0.5+kA)*(1-r)*g1e*g4Os*w14Om + (0.5+kD)*(0.5+kA)*(r)*g2e*g3Os*w23Om + (0.5+kD)*0.5*g2e*g4Os*w24Om + (0.5+kD)*(0.5+kA)*(r)*g3e*g2Os*w32Om + (0.5+kD)*(0.5+kA)*g3e*g4Os*w34Om + (0.5+kD)*(0.5+kA)*(1-r)*g4e*g1Os*w41Om + (0.5+kD)*0.5*g4e*g2Os*w42Om + (0.5+kD)*(0.5+kA)*g4e*g3Os*w43Om + (0.5+kD)*1*g4e*g4Os*w44Om)
 
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
  
  g1Os <- g1Os_pr 
  g2Os <- g2Os_pr
  g3Os <- g3Os_pr
  g4Os <- g4Os_pr
  
  SexRatioVector[d]<-g1Ys+g2Ys+g3Ys+g4Ys+g1Os+g2Os+g3Os+g4Os
  DfreqVector[d]<-g3Xs+g4Xs+g3Ys+g4Ys+g3Os+g4Os
  AfreqVector[d]<-g2Xs+g4Xs+g2Ys+g4Ys+g2Os+g4Os
  OFreq[d]<-g1Os+g2Os+g3Os+g4Os
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
  g1OsFreq[d]<-g1Os
  g2OsFreq[d]<-g2Os
  g3OsFreq[d]<-g3Os
  g4OsFreq[d]<-g4Os
  
  DA2D2Vector[d]<-((g1Xs+g1Ys+g1Os)*(g4Ys+g4Os+g4Xs))-((g2Xs+g2Ys+g2Os)*(g3Xs+g3Ys+g3Os))
  DA2YVector[d]<- (g2Ys+g4Ys+g2Os+g4Os)*(g1Xs+g3Xs)-(g1Ys+g3Ys+g1Os+g3Os)*(g2Xs+g4Xs)
  DD2YVector[d]<-((g3Ys+g4Ys+g3Os+g4Os)*(g1Xs+g2Xs))-((g3Xs+g4Xs)*(g1Ys+g2Ys+g1Os+g2Os))
  
}

x <- 1:G
y1 <- XFreq
y2 <- YFreq
y3 <- OFreq
y4 <- AfreqVector
y5 <- DfreqVector
y6 <- 0.5

library(tidyverse)
d <- tibble(G = x, X = y1, Y = y2, A = y3, D = y4, O = y5)
write.table(x = d, file = "data_5A.txt", row.names = F, col.names = T)

#### LINKAGE DISEQUILIBRIUM 

x <- 1:G
y1 <- DA2D2Vector
y2 <- DA2YVector
y3 <- DD2YVector
y4 <- 0.0

d <- tibble(G = x, DA2D2 = y1, DA2Y = y2, DD2Y = y3)
write.table(x = d, file = "data_5B.txt", row.names = F, col.names = T)

df=data.frame(g1XsFreq,g2XsFreq,g1YsFreq,g2YsFreq,g1OsFreq,g2OsFreq,g3XsFreq,g4XsFreq,g3YsFreq,g4YsFreq,g3OsFreq,g4OsFreq)

d <- tibble(G = x, 
            A1D1X = g1XsFreq,
            A2D1X = g2XsFreq,
            A1D1Y = g1YsFreq,
            A2D1Y = g2YsFreq,
            A1D1O = g1OsFreq,
            A2D1O = g2OsFreq,
            A1D2X = g3XsFreq,
            A2D2X = g4XsFreq,
            A1D2Y = g3YsFreq,
            A2D2Y = g4YsFreq,
            A1D2O = g3OsFreq,
            A2D2O = g4OsFreq)
write.table(x = d, file = "data_5C.txt", row.names = F, col.names = T)


ML_function <- function(G, kD, kA, s, hA, z, r, w) {
  
  original_kD <- kD
  original_z <- z
  
  ###BURN IN 
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
  v12f <- 1-hA*s
  v21f <- 1-hA*s
  v22f <- 1-s
  
  v11m <- 1
  v12m <- 1-hA*s
  v21m <- 1-hA*s
  v22m <- 1-s
  
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
  
  w11Ym <- u11m*v11m
  w12Ym <- u11m*v12m
  w13Ym <- u12m*v11m
  w14Ym <- u12m*v12m
  w21Ym <- u11m*v21m
  w22Ym <- u11m*v22m
  w23Ym <- u12m*v21m
  w24Ym <- u12m*v22m
  w31Ym <- u21m*v11m
  w32Ym <- u21m*v12m
  w33Ym <- u22m*v11m
  w34Ym <- u22m*v12m
  w41Ym <- u21m*v21m
  w42Ym <- u21m*v22m
  w43Ym <- u22m*v21m
  w44Ym <- u22m*v22m
  
  w11Om <- u11m*v11m*(1-w)
  w12Om <- u11m*v12m*(1-w)
  w13Om <- u12m*v11m*(1-w)
  w14Om <- u12m*v12m*(1-w)
  w21Om <- u11m*v21m*(1-w)
  w22Om <- u11m*v22m*(1-w)
  w23Om <- u12m*v21m*(1-w)
  w24Om <- u12m*v22m*(1-w)
  w31Om <- u21m*v11m*(1-w)
  w32Om <- u21m*v12m*(1-w)
  w33Om <- u22m*v11m*(1-w)
  w34Om <- u22m*v12m*(1-w)
  w41Om <- u21m*v21m*(1-w)
  w42Om <- u21m*v22m*(1-w)
  w43Om <- u22m*v21m*(1-w)
  w44Om <- u22m*v22m*(1-w)
  
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
  
  g3Xs <- 0.0
  g3Ys <- 0.0
  g4Xs <- 0.0
  g4Ys <- 0.0
  
  g1Os <- 0.0
  g2Os <- 0.0
  g3Os <- 0.0
  g4Os <- 0.0
  
  # in the process of changin "M" to either "Om" or "Xm" for male fitness (so that we can differentiate between these and add fitness costs!)
  
  for (d in 1:500){
    
    wfbar<-(g1e*g1Xs*w11f)+(g1e*g2Xs*w12f)+(g1e*g3Xs*w13f)+(g1e*g4Xs*w14f)+(g2e*g1Xs*w21f)+(g2e*g2Xs*w22f)+(g2e*g3Xs*w23f)+(g2e*g4Xs*w24f)+(g3e*g1Xs*w31f)+(g3e*g2Xs*w32f)+(g3e*g3Xs*w33f)+(g3e*g4Xs*w34f)+(g4e*g1Xs*w41f)+(g4e*g2Xs*w42f)+(g4e*g3Xs*w43f)+(g4e*g4Xs*w44f)
    
    g1e_pr<- (1/wfbar)*(1*g1e*g1Xs*w11f + 0.5* g1e*g2Xs*w12f + 0.5* g1e*g3Xs*w13f + 0.5*(1-r)* g1e*g4Xs*w14f + 0.5* g2e*g1Xs*w21f + 0.5*(r)* g2e*g3Xs*w23f + 0.5* g3e*g1Xs*w31f + 0.5*(r)* g3e*g2Xs*w32f + 0.5*(1-r)* g4e*g1Xs*w41f)
    g2e_pr<- (1/wfbar)*(0.5* g1e*g2Xs*w12f + 0.5*(r)* g1e*g4Xs*w14f + 0.5* g2e*g1Xs*w21f + 1* g2e*g2Xs*w22f + 0.5*(1-r)* g2e*g3Xs*w23f + 0.5* g2e*g4Xs*w24f + 0.5*(1-r)* g3e*g2Xs*w32f + 0.5*(r)* g4e*g1Xs*w41f + 0.5* g4e*g2Xs*w42f) 
    g3e_pr<- (1/wfbar)*(0.5* g1e*g3Xs*w13f + 0.5*(r)* g1e*g4Xs*w14f + 0.5*(1-r)* g2e*g3Xs*w23f + 0.5* g3e*g1Xs*w31f + 0.5*(1-r)* g3e*g2Xs*w32f + 1* g3e*g3Xs*w33f + 0.5* g3e*g4Xs*w34f + 0.5*(r)* g4e*g1Xs*w41f + 0.5* g4e*g3Xs*w43f) 
    g4e_pr<- (1/wfbar)*(0.5*(1-r)* g1e*g4Xs*w14f + 0.5*(r)* g2e*g3Xs*w23f + 0.5* g2e*g4Xs*w24f + 0.5*(r)* g3e*g2Xs*w32f + 0.5* g3e*g4Xs*w34f + 0.5*(1-r)* g4e*g1Xs*w41f + 0.5* g4e*g2Xs*w42f + 0.5* g4e*g3Xs*w43f + 1* g4e*g4Xs*w44f) 
    
    wmbar<-(g1e*g1Ys*w11Ym)+(g1e*g2Ys*w12Ym)+(g1e*g3Ys*w13Ym)+(g1e*g4Ys*w14Ym)+(g2e*g1Ys*w21Ym)+(g2e*g2Ys*w22Ym)+(g2e*g3Ys*w23Ym)+(g2e*g4Ys*w24Ym)+(g3e*g1Ys*w31Ym)+(g3e*g2Ys*w32Ym)+(g3e*g3Ys*w33Ym)+(g3e*g4Ys*w34Ym)+(g4e*g1Ys*w41Ym)+(g4e*g2Ys*w42Ym)+(g4e*g3Ys*w43Ym)+(g4e*g4Ys*w44Ym)+(g1e*g1Os*w11Om)+(g1e*g2Os*w12Om)+(g1e*g3Os*w13Om)+(g1e*g4Os*w14Om)+(g2e*g1Os*w21Om)+(g2e*g2Os*w22Om)+(g2e*g3Os*w23Om)+(g2e*g4Os*w24Om)+(g3e*g1Os*w31Om)+(g3e*g2Os*w32Om)+(g3e*g3Os*w33Om)+(g3e*g4Os*w34Om)+(g4e*g1Os*w41Om)+(g4e*g2Os*w42Om)+(g4e*g3Os*w43Om)+(g4e*g4Os*w44Om)
    g1Xs_pr<- (1/wmbar)*(0.5*1*g1e*g1Ys*w11Ym + 0.5*(0.5-kA)*g1e*g2Ys*w12Ym + (0.5-kD)*0.5*g1e*g3Ys*w13Ym + (0.5-kD)*(0.5-kA)*(1-r)*g1e*g4Ys*w14Ym + 0.5*(0.5-kA)*g2e*g1Ys*w21Ym + (0.5-kD)*(0.5-kA)*(r)*g2e*g3Ys*w23Ym + (0.5-kD)*0.5*g3e*g1Ys*w31Ym + (0.5-kD)*(0.5-kA)*(r)*g3e*g2Ys*w32Ym + (0.5-kD)*(0.5-kA)*(1-r)*g4e*g1Ys*w41Ym +0.5*1*g1e*g1Os*w11Om + 0.5*(0.5-kA)*g1e*g2Os*w12Om + (0.5-kD)*0.5*g1e*g3Os*w13Om + (0.5-kD)*(0.5-kA)*(1-r)*g1e*g4Os*w14Om + 0.5*(0.5-kA)*g2e*g1Os*w21Om + (0.5-kD)*(0.5-kA)*(r)*g2e*g3Os*w23Om + (0.5-kD)*0.5*g3e*g1Os*w31Om + (0.5-kD)*(0.5-kA)*(r)*g3e*g2Os*w32Om + (0.5-kD)*(0.5-kA)*(1-r)*g4e*g1Os*w41Om)
    g2Xs_pr<- (1/wmbar)*(0.5*(0.5+kA)*g1e*g2Ys*w12Ym + (0.5-kD)*(0.5+kA)*(r)*g1e*g4Ys*w14Ym + 0.5*(0.5+kA)*g2e*g1Ys*w21Ym + 0.5*1*g2e*g2Ys*w22Ym + (0.5-kD)*(0.5+kA)*(1-r)*g2e*g3Ys*w23Ym + (0.5-kD)*0.5*g2e*g4Ys*w24Ym + (0.5-kD)*(0.5+kA)*(1-r)*g3e*g2Ys*w32Ym + (0.5-kD)*(0.5+kA)*(r)*g4e*g1Ys*w41Ym + (0.5-kD)*0.5*g4e*g2Ys*w42Ym + 0.5*(0.5+kA)*g1e*g2Os*w12Om + (0.5-kD)*(0.5+kA)*(r)*g1e*g4Os*w14Om + (0.5)*(0.5+kA)*g2e*g1Os*w21Om + 0.5*1*g2e*g2Os*w22Om + (0.5-kD)*(0.5+kA)*(1-r)*g2e*g3Os*w23Om + (0.5-kD)*0.5*g2e*g4Os*w24Om + (0.5-kD)*(0.5+kA)*(1-r)*g3e*g2Os*w32Om + (0.5-kD)*(0.5+kA)*(r)*g4e*g1Os*w41Om + (0.5-kD)*0.5*g4e*g2Os*w42Om)
    g3Xs_pr<- (1/wmbar)*((0.5-kD)*0.5*g1e*g3Ys*w13Ym + (0.5-kD)*(0.5-kA)*(r)*g1e*g4Ys*w14Ym + (0.5-kD)*(0.5-kA)*(1-r)*g2e*g3Ys*w23Ym + (0.5-kD)*0.5*g3e*g1Ys*w31Ym + (0.5-kD)*(0.5-kA)*(1-r)*g3e*g2Ys*w32Ym + (0.5-kD)*1*g3e*g3Ys*w33Ym + (0.5-kD)*(0.5-kA)*g3e*g4Ys*w34Ym + (0.5-kD)*(0.5-kA)*(r)*g4e*g1Ys*w41Ym + (0.5-kD)*(0.5-kA)*g4e*g3Ys*w43Ym+(0.5-kD)*0.5*g1e*g3Os*w13Om + (0.5-kD)*(0.5-kA)*(r)*g1e*g4Os*w14Om + (0.5-kD)*(0.5-kA)*(1-r)*g2e*g3Os*w23Om + (0.5-kD)*0.5*g3e*g1Os*w31Om + (0.5-kD)*(0.5-kA)*(1-r)*g3e*g2Os*w32Om + (0.5-kD)*1*g3e*g3Os*w33Om + (0.5-kD)*(0.5-kA)*g3e*g4Os*w34Om + (0.5-kD)*(0.5-kA)*(r)*g4e*g1Os*w41Om + (0.5-kD)*(0.5-kA)*g4e*g3Os*w43Om)
    g4Xs_pr<- (1/wmbar)*((0.5-kD)*(0.5+kA)*(1-r)*g1e*g4Ys*w14Ym + (0.5-kD)*(0.5+kA)*(r)*g2e*g3Ys*w23Ym + (0.5-kD)*0.5*g2e*g4Ys*w24Ym + (0.5-kD)*(0.5+kA)*(r)*g3e*g2Ys*w32Ym + (0.5-kD)*(0.5+kA)*g3e*g4Ys*w34Ym + (0.5-kD)*(0.5+kA)*(1-r)*g4e*g1Ys*w41Ym + (0.5-kD)*0.5*g4e*g2Ys*w42Ym + (0.5-kD)*(0.5+kA)*g4e*g3Ys*w43Ym + (0.5-kD)*1* g4e*g4Ys*w44Ym + (0.5-kD)*(0.5+kA)*(1-r)*g1e*g4Os*w14Om + (0.5-kD)*(0.5+kA)*(r)*g2e*g3Os*w23Om + (0.5-kD)*0.5*g2e*g4Os*w24Om + (0.5-kD)*(0.5+kA)*(r)*g3e*g2Os*w32Om + (0.5-kD)*(0.5+kA)*g3e*g4Os*w34Om + (0.5-kD)*(0.5+kA)*(1-r)*g4e*g1Os*w41Om + (0.5-kD)*0.5*g4e*g2Os*w42Om + (0.5-kD)*(0.5+kA)*g4e*g3Os*w43Om + (0.5-kD)*1*g4e*g4Os*w44Om)
    
    g1Ys_pr<- (1/wmbar)*(0.5*1* g1e*g1Ys*w11Ym + 0.5*(0.5-kA)* g1e*g2Ys*w12Ym + 0.5*0.5* g1e*g3Ys*w13Ym + 0.5*(0.5-kA)*(1-r)*g1e*g4Ys*w14Ym + 0.5*(0.5-kA)* g2e*g1Ys*w21Ym + 0.5*(0.5-kA)*(r)* g2e*g3Ys*w23Ym + 0.5*0.5* g3e*g1Ys*w31Ym + 0.5*(0.5-kA)*(r)* g3e*g2Ys*w32Ym + 0.5*(0.5-kA)*(1-r)* g4e*g1Ys*w41Ym)
    g2Ys_pr<- (1/wmbar)*(0.5*(0.5+kA)* g1e*g2Ys*w12Ym + 0.5*(0.5+kA)*(r)* g1e*g4Ys*w14Ym + 0.5*(0.5+kA)* g2e*g1Ys*w21Ym + 0.5*1* g2e*g2Ys*w22Ym + 0.5*(0.5+kA)*(1-r)* g2e*g3Ys*w23Ym + 0.5*0.5* g2e*g4Ys*w24Ym + 0.5*(0.5+kA)*(1-r)* g3e*g2Ys*w32Ym + 0.5*(0.5+kA)*(r)* g4e*g1Ys*w41Ym + 0.5*0.5* g4e*g2Ys*w42Ym)
    g3Ys_pr<- (1/wmbar)*(0.5*0.5* g1e*g3Ys*w13Ym + 0.5*(0.5-kA)*(r)* g1e*g4Ys*w14Ym + 0.5*(0.5-kA)*(1-r)* g2e*g3Ys*w23Ym + 0.5*0.5* g3e*g1Ys*w31Ym + 0.5*(0.5-kA)*(1-r)* g3e*g2Ys*w32Ym + 0.5*1* g3e*g3Ys*w33Ym + 0.5*(0.5-kA)* g3e*g4Ys*w34Ym + 0.5*(0.5-kA)*(r)* g4e*g1Ys*w41Ym + 0.5*(0.5-kA)* g4e*g3Ys*w43Ym)
    g4Ys_pr<- (1/wmbar)*(0.5*(0.5+kA)*(1-r)* g1e*g4Ys*w14Ym + 0.5*(0.5+kA)*(r)* g2e*g3Ys*w23Ym + 0.5*0.5* g2e*g4Ys*w24Ym + 0.5*(0.5+kA)*(r)* g3e*g2Ys*w32Ym + 0.5*(0.5+kA)* g3e*g4Ys*w34Ym + 0.5*(0.5+kA)*(1-r)* g4e*g1Ys*w41Ym + 0.5*0.5* g4e*g2Ys*w42Ym + 0.5*(0.5+kA)* g4e*g3Ys*w43Ym + 0.5*1* g4e*g4Ys*w44Ym)
    g1Os_pr<- (1/wmbar)*(0*1*g1e*g1Ys*w11Ym + 0*(0.5-kA)* g1e*g2Ys*w12Ym + kD*0.5*g1e*g3Ys*w13Ym + kD*(0.5-kA)*(1-r)*g1e*g4Ys*w14Ym + 0*(0.5-kA)* g2e*g1Ys*w21Ym + kD*(0.5-kA)*(r)* g2e*g3Ys*w23Ym + kD*0.5* g3e*g1Ys*w31Ym + kD*(0.5-kA)*(r)*g3e*g2Ys*w32Ym + kD*(0.5-kA)*(1-r)*g4e*g1Ys*w41Ym +0.5*1*g1e*g1Os*w11Om + 0.5*(0.5-kA)* g1e*g2Os*w12Om + (0.5+kD)*0.5* g1e*g3Os*w13Om + (0.5+kD)*(0.5-kA)*(1-r)* g1e*g4Os*w14Om + 0.5*(0.5-kA)*g2e*g1Os*w21Om + (0.5+kD)*(0.5-kA)*(r)* g2e*g3Os*w23Om + (0.5+kD)*0.5* g3e*g1Os*w31Om + (0.5+kD)*(0.5-kA)*(r)* g3e*g2Os*w32Om + (0.5+kD)*(0.5-kA)*(1-r)* g4e*g1Os*w41Om)
    g2Os_pr<- (1/wmbar)*(0*(0.5+kA)*g1e*g2Ys*w12Ym + kD*(0.5+kA)*(r)* g1e*g4Ys*w14Ym + 0*(0.5+kA)* g2e*g1Ys*w21Ym + 0*1* g2e*g2Ys*w22Ym + kD*(0.5+kA)*(1-r)* g2e*g3Ys*w23Ym + kD*0.5* g2e*g4Ys*w24Ym + kD*(0.5+kA)*(1-r)* g3e*g2Ys*w32Ym + kD*(0.5+kA)*(r)* g4e*g1Ys*w41Ym + kD*0.5* g4e*g2Ys*w42Ym + 0.5*(0.5+kA)* g1e*g2Os*w12Om + (0.5+kD)*(0.5+kA)*(r)* g1e*g4Os*w14Om + 0.5*(0.5+kA)* g2e*g1Os*w21Om + 0.5*1*g2e*g2Os*w22Om + (0.5+kD)*(0.5+kA)*(1-r)* g2e*g3Os*w23Om + (0.5+kD)*0.5* g2e*g4Os*w24Om + (0.5+kD)*(0.5+kA)*(1-r)* g3e*g2Os*w32Om + (0.5+kD)*(0.5+kA)*(r)* g4e*g1Os*w41Om + (0.5+kD)*0.5* g4e*g2Os*w42Om)
    g3Os_pr<- (1/wmbar)*(kD*0.5*g1e*g3Ys*w13Ym + kD*(0.5-kA)*(r)*g1e*g4Ys*w14Ym + kD*(0.5-kA)*(1-r)* g2e*g3Ys*w23Ym + kD*0.5* g3e*g1Ys*w31Ym + kD*(0.5-kA)*(1-r)*g3e*g2Ys*w32Ym + kD*1* g3e*g3Ys*w33Ym + kD*(0.5-kA)* g3e*g4Ys*w34Ym + kD*(0.5-kA)*(r)* g4e*g1Ys*w41Ym + kD*(0.5-kA)* g4e*g3Ys*w43Ym + (0.5+kD)*0.5* g1e*g3Os*w13Om + (0.5+kD)*(0.5-kA)*(r)* g1e*g4Os*w14Om + (0.5+kD)*(0.5-kA)*(1-r)* g2e*g3Os*w23Om + (0.5+kD)*0.5* g3e*g1Os*w31Om + (0.5+kD)*(0.5-kA)*(1-r)* g3e*g2Os*w32Om + (0.5+kD)*1* g3e*g3Os*w33Om + (0.5+kD)*(0.5-kA)* g3e*g4Os*w34Om + (0.5+kD)*(0.5-kA)*(r)* g4e*g1Os*w41Om + (0.5+kD)*(0.5-kA)* g4e*g3Os*w43Om)
    g4Os_pr<- (1/wmbar)*(kD*(0.5+kA)*(1-r)*g1e*g4Ys*w14Ym + kD*(0.5+kA)*(r)*g2e*g3Ys*w23Ym + kD*0.5*g2e*g4Ys*w24Ym + kD*(0.5+kA)*(r)*g3e*g2Ys*w32Ym + kD*(0.5+kA)*g3e*g4Ys*w34Ym + kD*(0.5+kA)*(1-r)*g4e*g1Ys*w41Ym + kD*0.5*g4e*g2Ys*w42Ym + kD*(0.5+kA)*g4e*g3Ys*w43Ym + kD*(0.5+kA)*1*g4e*g4Ys*w44Ym + (0.5+kD)*(0.5+kA)*(1-r)*g1e*g4Os*w14Om + (0.5+kD)*(0.5+kA)*(r)*g2e*g3Os*w23Om + (0.5+kD)*0.5*g2e*g4Os*w24Om + (0.5+kD)*(0.5+kA)*(r)*g3e*g2Os*w32Om + (0.5+kD)*(0.5+kA)*g3e*g4Os*w34Om + (0.5+kD)*(0.5+kA)*(1-r)*g4e*g1Os*w41Om + (0.5+kD)*0.5*g4e*g2Os*w42Om + (0.5+kD)*(0.5+kA)*g4e*g3Os*w43Om + (0.5+kD)*1*g4e*g4Os*w44Om)
    
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
  
  #SIMULATION 
  g1e <- g1e
  g2e <- g2e
  g3e <- g3e
  g4e <- g4e
  
  g1Xs <- g1Xs
  g2Xs <- g2Xs
  g3Xs <- g3Xs
  g4Xs <- g4Xs
  
  g1Os <- g1Os
  g2Os <- g2Os
  g3Os <- g3Os
  g4Os <- g4Os
  
  g3Ys <- g1Ys_pr*0.001
  g4Ys <- g2Ys_pr*0.001
  g1Ys <- g1Ys_pr-g1Ys_pr*0.001
  g2Ys <- g2Ys_pr-g2Ys_pr*0.001

  
  kD <- original_kD 
  z <- original_z
  
  DfreqVector<-numeric(G)
  AfreqVector<-numeric(G)
  SexRatioVector<-numeric(G)
  OFreq <-numeric(G)
  YFreq<-numeric(G)
  XFreq<-numeric(G)
  
  # fitness of the remote-control driving locus (1 is neutrla, 2 drives)
  u11f <- 1
  u12f <- 1
  u21f <- 1
  u22f <- 1
  
  u11m <- 1
  u12m <- 1-z
  u21m <- 1-z
  u22m <- 1-z
  
  #v is the cis drive locus (1 is neutral, 2 drives)
  v11f <- 1 
  v12f <- 1-hA*s
  v21f <- 1-hA*s
  v22f <- 1-s
  
  v11m <- 1
  v12m <- 1-hA*s
  v21m <- 1-hA*s
  v22m <- 1-s
  
  w11Ym <- u11m*v11m
  w12Ym <- u11m*v12m
  w13Ym <- u12m*v11m
  w14Ym <- u12m*v12m
  w21Ym <- u11m*v21m
  w22Ym <- u11m*v22m
  w23Ym <- u12m*v21m
  w24Ym <- u12m*v22m
  w31Ym <- u21m*v11m
  w32Ym <- u21m*v12m
  w33Ym <- u22m*v11m
  w34Ym <- u22m*v12m
  w41Ym <- u21m*v21m
  w42Ym <- u21m*v22m
  w43Ym <- u22m*v21m
  w44Ym <- u22m*v22m
  
  w11Om <- u11m*v11m*(1-w)
  w12Om <- u11m*v12m*(1-w)
  w13Om <- u12m*v11m*(1-w)
  w14Om <- u12m*v12m*(1-w)
  w21Om <- u11m*v21m*(1-w)
  w22Om <- u11m*v22m*(1-w)
  w23Om <- u12m*v21m*(1-w)
  w24Om <- u12m*v22m*(1-w)
  w31Om <- u21m*v11m*(1-w)
  w32Om <- u21m*v12m*(1-w)
  w33Om <- u22m*v11m*(1-w)
  w34Om <- u22m*v12m*(1-w)
  w41Om <- u21m*v21m*(1-w)
  w42Om <- u21m*v22m*(1-w)
  w43Om <- u22m*v21m*(1-w)
  w44Om <- u22m*v22m*(1-w)
  
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
  
  #Allele Frequency Recursions 
  for (d in 1:G){
    
    wfbar<-(g1e*g1Xs*w11f)+(g1e*g2Xs*w12f)+(g1e*g3Xs*w13f)+(g1e*g4Xs*w14f)+(g2e*g1Xs*w21f)+(g2e*g2Xs*w22f)+(g2e*g3Xs*w23f)+(g2e*g4Xs*w24f)+(g3e*g1Xs*w31f)+(g3e*g2Xs*w32f)+(g3e*g3Xs*w33f)+(g3e*g4Xs*w34f)+(g4e*g1Xs*w41f)+(g4e*g2Xs*w42f)+(g4e*g3Xs*w43f)+(g4e*g4Xs*w44f)
    
    g1e_pr<- (1/wfbar)*(1*g1e*g1Xs*w11f + 0.5* g1e*g2Xs*w12f + 0.5* g1e*g3Xs*w13f + 0.5*(1-r)* g1e*g4Xs*w14f + 0.5* g2e*g1Xs*w21f + 0.5*(r)* g2e*g3Xs*w23f + 0.5* g3e*g1Xs*w31f + 0.5*(r)* g3e*g2Xs*w32f + 0.5*(1-r)* g4e*g1Xs*w41f)
    g2e_pr<- (1/wfbar)*(0.5* g1e*g2Xs*w12f + 0.5*(r)* g1e*g4Xs*w14f + 0.5* g2e*g1Xs*w21f + 1* g2e*g2Xs*w22f + 0.5*(1-r)* g2e*g3Xs*w23f + 0.5* g2e*g4Xs*w24f + 0.5*(1-r)* g3e*g2Xs*w32f + 0.5*(r)* g4e*g1Xs*w41f + 0.5* g4e*g2Xs*w42f) 
    g3e_pr<- (1/wfbar)*(0.5* g1e*g3Xs*w13f + 0.5*(r)* g1e*g4Xs*w14f + 0.5*(1-r)* g2e*g3Xs*w23f + 0.5* g3e*g1Xs*w31f + 0.5*(1-r)* g3e*g2Xs*w32f + 1* g3e*g3Xs*w33f + 0.5* g3e*g4Xs*w34f + 0.5*(r)* g4e*g1Xs*w41f + 0.5* g4e*g3Xs*w43f) 
    g4e_pr<- (1/wfbar)*(0.5*(1-r)* g1e*g4Xs*w14f + 0.5*(r)* g2e*g3Xs*w23f + 0.5* g2e*g4Xs*w24f + 0.5*(r)* g3e*g2Xs*w32f + 0.5* g3e*g4Xs*w34f + 0.5*(1-r)* g4e*g1Xs*w41f + 0.5* g4e*g2Xs*w42f + 0.5* g4e*g3Xs*w43f + 1* g4e*g4Xs*w44f) 
    
    wmbar<-(g1e*g1Ys*w11Ym)+(g1e*g2Ys*w12Ym)+(g1e*g3Ys*w13Ym)+(g1e*g4Ys*w14Ym)+(g2e*g1Ys*w21Ym)+(g2e*g2Ys*w22Ym)+(g2e*g3Ys*w23Ym)+(g2e*g4Ys*w24Ym)+(g3e*g1Ys*w31Ym)+(g3e*g2Ys*w32Ym)+(g3e*g3Ys*w33Ym)+(g3e*g4Ys*w34Ym)+(g4e*g1Ys*w41Ym)+(g4e*g2Ys*w42Ym)+(g4e*g3Ys*w43Ym)+(g4e*g4Ys*w44Ym)+(g1e*g1Os*w11Om)+(g1e*g2Os*w12Om)+(g1e*g3Os*w13Om)+(g1e*g4Os*w14Om)+(g2e*g1Os*w21Om)+(g2e*g2Os*w22Om)+(g2e*g3Os*w23Om)+(g2e*g4Os*w24Om)+(g3e*g1Os*w31Om)+(g3e*g2Os*w32Om)+(g3e*g3Os*w33Om)+(g3e*g4Os*w34Om)+(g4e*g1Os*w41Om)+(g4e*g2Os*w42Om)+(g4e*g3Os*w43Om)+(g4e*g4Os*w44Om)
    g1Xs_pr<- (1/wmbar)*(0.5*1*g1e*g1Ys*w11Ym + 0.5*(0.5-kA)*g1e*g2Ys*w12Ym + (0.5-kD)*0.5*g1e*g3Ys*w13Ym + (0.5-kD)*(0.5-kA)*(1-r)*g1e*g4Ys*w14Ym + 0.5*(0.5-kA)*g2e*g1Ys*w21Ym + (0.5-kD)*(0.5-kA)*(r)*g2e*g3Ys*w23Ym + (0.5-kD)*0.5*g3e*g1Ys*w31Ym + (0.5-kD)*(0.5-kA)*(r)*g3e*g2Ys*w32Ym + (0.5-kD)*(0.5-kA)*(1-r)*g4e*g1Ys*w41Ym +0.5*1*g1e*g1Os*w11Om + 0.5*(0.5-kA)*g1e*g2Os*w12Om + (0.5-kD)*0.5*g1e*g3Os*w13Om + (0.5-kD)*(0.5-kA)*(1-r)*g1e*g4Os*w14Om + 0.5*(0.5-kA)*g2e*g1Os*w21Om + (0.5-kD)*(0.5-kA)*(r)*g2e*g3Os*w23Om + (0.5-kD)*0.5*g3e*g1Os*w31Om + (0.5-kD)*(0.5-kA)*(r)*g3e*g2Os*w32Om + (0.5-kD)*(0.5-kA)*(1-r)*g4e*g1Os*w41Om)
    g2Xs_pr<- (1/wmbar)*(0.5*(0.5+kA)*g1e*g2Ys*w12Ym + (0.5-kD)*(0.5+kA)*(r)*g1e*g4Ys*w14Ym + 0.5*(0.5+kA)*g2e*g1Ys*w21Ym + 0.5*1*g2e*g2Ys*w22Ym + (0.5-kD)*(0.5+kA)*(1-r)*g2e*g3Ys*w23Ym + (0.5-kD)*0.5*g2e*g4Ys*w24Ym + (0.5-kD)*(0.5+kA)*(1-r)*g3e*g2Ys*w32Ym + (0.5-kD)*(0.5+kA)*(r)*g4e*g1Ys*w41Ym + (0.5-kD)*0.5*g4e*g2Ys*w42Ym + 0.5*(0.5+kA)*g1e*g2Os*w12Om + (0.5-kD)*(0.5+kA)*(r)*g1e*g4Os*w14Om + (0.5)*(0.5+kA)*g2e*g1Os*w21Om + 0.5*1*g2e*g2Os*w22Om + (0.5-kD)*(0.5+kA)*(1-r)*g2e*g3Os*w23Om + (0.5-kD)*0.5*g2e*g4Os*w24Om + (0.5-kD)*(0.5+kA)*(1-r)*g3e*g2Os*w32Om + (0.5-kD)*(0.5+kA)*(r)*g4e*g1Os*w41Om + (0.5-kD)*0.5*g4e*g2Os*w42Om)
    g3Xs_pr<- (1/wmbar)*((0.5-kD)*0.5*g1e*g3Ys*w13Ym + (0.5-kD)*(0.5-kA)*(r)*g1e*g4Ys*w14Ym + (0.5-kD)*(0.5-kA)*(1-r)*g2e*g3Ys*w23Ym + (0.5-kD)*0.5*g3e*g1Ys*w31Ym + (0.5-kD)*(0.5-kA)*(1-r)*g3e*g2Ys*w32Ym + (0.5-kD)*1*g3e*g3Ys*w33Ym + (0.5-kD)*(0.5-kA)*g3e*g4Ys*w34Ym + (0.5-kD)*(0.5-kA)*(r)*g4e*g1Ys*w41Ym + (0.5-kD)*(0.5-kA)*g4e*g3Ys*w43Ym+(0.5-kD)*0.5*g1e*g3Os*w13Om + (0.5-kD)*(0.5-kA)*(r)*g1e*g4Os*w14Om + (0.5-kD)*(0.5-kA)*(1-r)*g2e*g3Os*w23Om + (0.5-kD)*0.5*g3e*g1Os*w31Om + (0.5-kD)*(0.5-kA)*(1-r)*g3e*g2Os*w32Om + (0.5-kD)*1*g3e*g3Os*w33Om + (0.5-kD)*(0.5-kA)*g3e*g4Os*w34Om + (0.5-kD)*(0.5-kA)*(r)*g4e*g1Os*w41Om + (0.5-kD)*(0.5-kA)*g4e*g3Os*w43Om)
    g4Xs_pr<- (1/wmbar)*((0.5-kD)*(0.5+kA)*(1-r)*g1e*g4Ys*w14Ym + (0.5-kD)*(0.5+kA)*(r)*g2e*g3Ys*w23Ym + (0.5-kD)*0.5*g2e*g4Ys*w24Ym + (0.5-kD)*(0.5+kA)*(r)*g3e*g2Ys*w32Ym + (0.5-kD)*(0.5+kA)*g3e*g4Ys*w34Ym + (0.5-kD)*(0.5+kA)*(1-r)*g4e*g1Ys*w41Ym + (0.5-kD)*0.5*g4e*g2Ys*w42Ym + (0.5-kD)*(0.5+kA)*g4e*g3Ys*w43Ym + (0.5-kD)*1* g4e*g4Ys*w44Ym + (0.5-kD)*(0.5+kA)*(1-r)*g1e*g4Os*w14Om + (0.5-kD)*(0.5+kA)*(r)*g2e*g3Os*w23Om + (0.5-kD)*0.5*g2e*g4Os*w24Om + (0.5-kD)*(0.5+kA)*(r)*g3e*g2Os*w32Om + (0.5-kD)*(0.5+kA)*g3e*g4Os*w34Om + (0.5-kD)*(0.5+kA)*(1-r)*g4e*g1Os*w41Om + (0.5-kD)*0.5*g4e*g2Os*w42Om + (0.5-kD)*(0.5+kA)*g4e*g3Os*w43Om + (0.5-kD)*1*g4e*g4Os*w44Om)
    
    g1Ys_pr<- (1/wmbar)*(0.5*1* g1e*g1Ys*w11Ym + 0.5*(0.5-kA)* g1e*g2Ys*w12Ym + 0.5*0.5* g1e*g3Ys*w13Ym + 0.5*(0.5-kA)*(1-r)*g1e*g4Ys*w14Ym + 0.5*(0.5-kA)* g2e*g1Ys*w21Ym + 0.5*(0.5-kA)*(r)* g2e*g3Ys*w23Ym + 0.5*0.5* g3e*g1Ys*w31Ym + 0.5*(0.5-kA)*(r)* g3e*g2Ys*w32Ym + 0.5*(0.5-kA)*(1-r)* g4e*g1Ys*w41Ym)
    g2Ys_pr<- (1/wmbar)*(0.5*(0.5+kA)* g1e*g2Ys*w12Ym + 0.5*(0.5+kA)*(r)* g1e*g4Ys*w14Ym + 0.5*(0.5+kA)* g2e*g1Ys*w21Ym + 0.5*1* g2e*g2Ys*w22Ym + 0.5*(0.5+kA)*(1-r)* g2e*g3Ys*w23Ym + 0.5*0.5* g2e*g4Ys*w24Ym + 0.5*(0.5+kA)*(1-r)* g3e*g2Ys*w32Ym + 0.5*(0.5+kA)*(r)* g4e*g1Ys*w41Ym + 0.5*0.5* g4e*g2Ys*w42Ym)
    g3Ys_pr<- (1/wmbar)*(0.5*0.5* g1e*g3Ys*w13Ym + 0.5*(0.5-kA)*(r)* g1e*g4Ys*w14Ym + 0.5*(0.5-kA)*(1-r)* g2e*g3Ys*w23Ym + 0.5*0.5* g3e*g1Ys*w31Ym + 0.5*(0.5-kA)*(1-r)* g3e*g2Ys*w32Ym + 0.5*1* g3e*g3Ys*w33Ym + 0.5*(0.5-kA)* g3e*g4Ys*w34Ym + 0.5*(0.5-kA)*(r)* g4e*g1Ys*w41Ym + 0.5*(0.5-kA)* g4e*g3Ys*w43Ym)
    g4Ys_pr<- (1/wmbar)*(0.5*(0.5+kA)*(1-r)* g1e*g4Ys*w14Ym + 0.5*(0.5+kA)*(r)* g2e*g3Ys*w23Ym + 0.5*0.5* g2e*g4Ys*w24Ym + 0.5*(0.5+kA)*(r)* g3e*g2Ys*w32Ym + 0.5*(0.5+kA)* g3e*g4Ys*w34Ym + 0.5*(0.5+kA)*(1-r)* g4e*g1Ys*w41Ym + 0.5*0.5* g4e*g2Ys*w42Ym + 0.5*(0.5+kA)* g4e*g3Ys*w43Ym + 0.5*1* g4e*g4Ys*w44Ym)
    g1Os_pr<- (1/wmbar)*(0*1*g1e*g1Ys*w11Ym + 0*(0.5-kA)* g1e*g2Ys*w12Ym + kD*0.5*g1e*g3Ys*w13Ym + kD*(0.5-kA)*(1-r)*g1e*g4Ys*w14Ym + 0*(0.5-kA)* g2e*g1Ys*w21Ym + kD*(0.5-kA)*(r)* g2e*g3Ys*w23Ym + kD*0.5* g3e*g1Ys*w31Ym + kD*(0.5-kA)*(r)*g3e*g2Ys*w32Ym + kD*(0.5-kA)*(1-r)*g4e*g1Ys*w41Ym +0.5*1*g1e*g1Os*w11Om + 0.5*(0.5-kA)* g1e*g2Os*w12Om + (0.5+kD)*0.5* g1e*g3Os*w13Om + (0.5+kD)*(0.5-kA)*(1-r)* g1e*g4Os*w14Om + 0.5*(0.5-kA)*g2e*g1Os*w21Om + (0.5+kD)*(0.5-kA)*(r)* g2e*g3Os*w23Om + (0.5+kD)*0.5* g3e*g1Os*w31Om + (0.5+kD)*(0.5-kA)*(r)* g3e*g2Os*w32Om + (0.5+kD)*(0.5-kA)*(1-r)* g4e*g1Os*w41Om)
    g2Os_pr<- (1/wmbar)*(0*(0.5+kA)*g1e*g2Ys*w12Ym + kD*(0.5+kA)*(r)* g1e*g4Ys*w14Ym + 0*(0.5+kA)* g2e*g1Ys*w21Ym + 0*1* g2e*g2Ys*w22Ym + kD*(0.5+kA)*(1-r)* g2e*g3Ys*w23Ym + kD*0.5* g2e*g4Ys*w24Ym + kD*(0.5+kA)*(1-r)* g3e*g2Ys*w32Ym + kD*(0.5+kA)*(r)* g4e*g1Ys*w41Ym + kD*0.5* g4e*g2Ys*w42Ym + 0.5*(0.5+kA)* g1e*g2Os*w12Om + (0.5+kD)*(0.5+kA)*(r)* g1e*g4Os*w14Om + 0.5*(0.5+kA)* g2e*g1Os*w21Om + 0.5*1*g2e*g2Os*w22Om + (0.5+kD)*(0.5+kA)*(1-r)* g2e*g3Os*w23Om + (0.5+kD)*0.5* g2e*g4Os*w24Om + (0.5+kD)*(0.5+kA)*(1-r)* g3e*g2Os*w32Om + (0.5+kD)*(0.5+kA)*(r)* g4e*g1Os*w41Om + (0.5+kD)*0.5* g4e*g2Os*w42Om)
    g3Os_pr<- (1/wmbar)*(kD*0.5*g1e*g3Ys*w13Ym + kD*(0.5-kA)*(r)*g1e*g4Ys*w14Ym + kD*(0.5-kA)*(1-r)* g2e*g3Ys*w23Ym + kD*0.5* g3e*g1Ys*w31Ym + kD*(0.5-kA)*(1-r)*g3e*g2Ys*w32Ym + kD*1* g3e*g3Ys*w33Ym + kD*(0.5-kA)* g3e*g4Ys*w34Ym + kD*(0.5-kA)*(r)* g4e*g1Ys*w41Ym + kD*(0.5-kA)* g4e*g3Ys*w43Ym + (0.5+kD)*0.5* g1e*g3Os*w13Om + (0.5+kD)*(0.5-kA)*(r)* g1e*g4Os*w14Om + (0.5+kD)*(0.5-kA)*(1-r)* g2e*g3Os*w23Om + (0.5+kD)*0.5* g3e*g1Os*w31Om + (0.5+kD)*(0.5-kA)*(1-r)* g3e*g2Os*w32Om + (0.5+kD)*1* g3e*g3Os*w33Om + (0.5+kD)*(0.5-kA)* g3e*g4Os*w34Om + (0.5+kD)*(0.5-kA)*(r)* g4e*g1Os*w41Om + (0.5+kD)*(0.5-kA)* g4e*g3Os*w43Om)
    g4Os_pr<- (1/wmbar)*(kD*(0.5+kA)*(1-r)*g1e*g4Ys*w14Ym + kD*(0.5+kA)*(r)*g2e*g3Ys*w23Ym + kD*0.5*g2e*g4Ys*w24Ym + kD*(0.5+kA)*(r)*g3e*g2Ys*w32Ym + kD*(0.5+kA)*g3e*g4Ys*w34Ym + kD*(0.5+kA)*(1-r)*g4e*g1Ys*w41Ym + kD*0.5*g4e*g2Ys*w42Ym + kD*(0.5+kA)*g4e*g3Ys*w43Ym + kD*(0.5+kA)*1*g4e*g4Ys*w44Ym + (0.5+kD)*(0.5+kA)*(1-r)*g1e*g4Os*w14Om + (0.5+kD)*(0.5+kA)*(r)*g2e*g3Os*w23Om + (0.5+kD)*0.5*g2e*g4Os*w24Om + (0.5+kD)*(0.5+kA)*(r)*g3e*g2Os*w32Om + (0.5+kD)*(0.5+kA)*g3e*g4Os*w34Om + (0.5+kD)*(0.5+kA)*(1-r)*g4e*g1Os*w41Om + (0.5+kD)*0.5*g4e*g2Os*w42Om + (0.5+kD)*(0.5+kA)*g4e*g3Os*w43Om + (0.5+kD)*1*g4e*g4Os*w44Om)
    
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
    
    g1Os <- g1Os_pr 
    g2Os <- g2Os_pr
    g3Os <- g3Os_pr
    g4Os <- g4Os_pr
    
    #add the Os to this as well! 
    SexRatioVector[d]<-g1Ys+g2Ys+g3Ys+g4Ys+g1Os+g2Os+g3Os+g4Os
    DfreqVector[d]<-g3Xs+g4Xs+g3Ys+g4Ys+g3Os+g4Os
    AfreqVector[d]<-g2Xs+g4Xs+g2Ys+g4Ys+g2Os+g4Os
    OFreq[d]<-g1Os+g2Os+g3Os+g4Os
    YFreq[d]<-g1Ys+g2Ys+g3Ys+g4Ys
    XFreq[d]<-g1Xs+g2Xs+g3Xs+g4Xs
  }
  #return(OFreq[d]/(OFreq[d]+YFreq[d]))
  return(list(DfreqVector[d],SexRatioVector[d]))
  #return(SexRatioVector[d])
}

#Visualization 
kAvector<-c(50)
kDvector<-c(50)
DFvector<-c(2500)
SRvector<-c(2500)
#equalkvector<-c(2500)
#equalwvector<-c(2500)

i<-1 
q<-1
n <- 200
for (a in 1:n){
  kA <- a/(2*n)
  for (b in 1:n){
    kD <- b/(2*n)
    #z<- MLMD_function(1000,0.01,0.5,0.6,0.0,k,j,0.001,0)
    DF<-ML_function(1000, kD, kA, 0.5, 0, 0.01, 0.001, 0.001)[[1]]
    SR<-ML_function(1000, kD, kA, 0.5, 0, 0.01, 0.001, 0.001)[[2]]
    kAvector[i]<-kA
    kDvector[i]<-kD
    DFvector[i]<-DF
    SRvector[i]<-SR 
    i<-i+1
    
  }
  
  print(a)
}
# zvector[i]<-1
library(plotly)
library(RColorBrewer)
# First plot
plot1 <- plot_ly(
  x = kDvector,
  y = kAvector,
  z = DFvector,
  #z=SRvector
  type = "contour",
  colorscale = rev("Viridis"),  # Define the colorscale here
  autocontour = TRUE,
  contours = list(
    start = 0.5,
    end = 1,
    size = 0.01,
    showlines = FALSE,
    showlabs = FALSE,
    coloring = "fill"
  ),
  colorbar = list(
    titlefont = list(size = 20), 
    title = "% male",
    thickness = 30,
    len = 0.8,  # Adjust the length of the colorbar
    xanchor = "left",  # Anchor colorbar to the left side
    y = 1,  # Place it at the center vertically
    yanchor = "top",
    tickfont = list(size = 20)  # Adjust font size of colorbar ticks
  )
) %>%
  layout(
    xaxis = list(
      #titlefont = list(size = 20), 
      #title = 'trans drive (k)', 
      tickfont = list(size = 20), 
      tickangle = 0,  # Optional: Adjust tick angle if necessary
      axisfont = list(size = 50)
    ),
    yaxis = list(
      #titlefont = list(size = 20), 
      #title = 'fitness cost (s,t)', 
      tickfont = list(size = 20),
      axisfont = list(size = 50)
    )
  )

# Show plot
plot1


DFvector[1]
DF <- rep(NA, length(DFvector))
for(i in 1:length(DF))
{
  DF[i] <- DFvector[[i]]
}

library(tidyverse)
d <- tibble(kD = kDvector, kA = kAvector, DF = DFvector)

write.table(x = d, file = "data_5D.txt", row.names = F, col.names = T)


SRvector[1]
SR <- rep(NA, length(SRvector))
for(i in 1:length(SR))
{
  SR[i] <- SRvector[[i]]
}

library(tidyverse)
d <- tibble(kD = kDvector, kA = kAvector, SR = SRvector)

write.table(x = d, file = "data_5E.txt", row.names = F, col.names = T)