
#Transition from XY to XO / Sexual Antagonism and Male-Limited Drive 

#1= D1A1 (X-shredder, male beneficial)
#2= D1A2 (X-shredder, female beneficial)
#3= D2A1 (neutral, male beneficial)
#4= D2A2 (neutral, female beneficial)

SAfunction<-function(G,kD,s,t,hA,z,r,w){
  
  original_kD <- kD
  original_z <- z
  
  z<- 0 
  kD <-0 
  
  DfreqVector<-numeric(G)
  AfreqVector<-numeric(G)
  SexRatioVector<-numeric(G)
  OFreq <-numeric(G)
  XFreq <-numeric(G)
  YFreq <-numeric(G)
  
  u11f <- 1
  u12f <- 1
  u21f <- 1
  u22f <- 1
  
  u11Ym <- 1
  u12Ym <- 1-z
  u21Ym <- 1-z
  u22Ym <- 1-z
  
  u11Om <- 1
  u12Om <- 1-z
  u21Om <- 1-z
  u22Om <- 1-z
  
  # v is the female beneficial locus (1 is female beneficial, 2 is male beneficial)
  v11f <- 1
  v12f <- 1-hA*t
  v21f <- 1-hA*t
  v22f <- 1-t
  
  v11m <- 1-s
  v12m <- 1-hA*s
  v21m <- 1-hA*s
  v22m <- 1
  
  w11Ym <- u11Ym*v11m
  w12Ym <- u11Ym*v12m
  w13Ym <- u12Ym*v11m
  w14Ym <- u12Ym*v12m
  w21Ym <- u11Ym*v21m
  w22Ym <- u11Ym*v22m
  w23Ym <- u12Ym*v21m
  w24Ym <- u12Ym*v22m
  w31Ym <- u21Ym*v11m
  w32Ym <- u21Ym*v12m
  w33Ym <- u22Ym*v11m
  w34Ym <- u22Ym*v12m
  w41Ym <- u21Ym*v21m
  w42Ym <- u21Ym*v22m
  w43Ym <- u22Ym*v21m
  w44Ym <- u22Ym*v22m
  
  w11Om <- u11Om*v11m*(1-w)
  w12Om <- u11Om*v12m*(1-w)
  w13Om <- u12Om*v11m*(1-w)
  w14Om <- u12Om*v12m*(1-w)
  w21Om <- u11Om*v21m*(1-w)
  w22Om <- u11Om*v22m*(1-w)
  w23Om <- u12Om*v21m*(1-w)
  w24Om <- u12Om*v22m*(1-w)
  w31Om <- u21Om*v11m*(1-w)
  w32Om <- u21Om*v12m*(1-w)
  w33Om <- u22Om*v11m*(1-w)
  w34Om <- u22Om*v12m*(1-w)
  w41Om <- u21Om*v21m*(1-w)
  w42Om <- u21Om*v22m*(1-w)
  w43Om <- u22Om*v21m*(1-w)
  w44Om <- u22Om*v22m*(1-w)
  
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
  
  # initialize haplotype frequencies
  g1e <- 0.5
  g2e <- 0.5
  g3e <- 0.00
  g4e <- 0.00
  
  g1Xs <- 0.25
  g1Ys <- 0.25
  g2Xs <- 0.25
  g2Ys <- 0.25
  
  g3Xs <- 0.00
  g3Ys <- 0.00
  g4Xs <- 0.00
  g4Ys <- 0.00
  
  g1Os <- 0.0
  g2Os <- 0.0
  g3Os <- 0.0
  g4Os <- 0.0
  
  for (d in 1:500){
    
    wfbar<-(g1e*g1Xs*w11f)+(g1e*g2Xs*w12f)+(g1e*g3Xs*w13f)+(g1e*g4Xs*w14f)+(g2e*g1Xs*w21f)+(g2e*g2Xs*w22f)+(g2e*g3Xs*w23f)+(g2e*g4Xs*w24f)+(g3e*g1Xs*w31f)+(g3e*g2Xs*w32f)+(g3e*g3Xs*w33f)+(g3e*g4Xs*w34f)+(g4e*g1Xs*w41f)+(g4e*g2Xs*w42f)+(g4e*g3Xs*w43f)+(g4e*g4Xs*w44f)
    g1e_pr<-(1/wfbar)*(1* g1e*g1Xs*w11f + 0.5* g1e*g2Xs*w12f + 0.5* g1e*g3Xs*w13f + 0.5*(1-r)* g1e*g4Xs*w14f + 0.5* g2e*g1Xs*w21f + 0.5*(r)* g2e*g3Xs*w23f + 0.5* g3e*g1Xs*w31f + 0.5*(r)* g3e*g2Xs*w32f + 0.5*(1-r)* g4e*g1Xs*w41f)
    g2e_pr<-(1/wfbar)*(0.5* g1e*g2Xs*w12f + 0.5*(r)* g1e*g4Xs*w14f + 0.5* g2e*g1Xs*w21f + 1* g2e*g2Xs*w22f + 0.5*(1-r)* g2e*g3Xs*w23f + 0.5* g2e*g4Xs*w24f + 0.5*(1-r)* g3e*g2Xs*w32f + 0.5*(r)* g4e*g1Xs*w41f + 0.5* g4e*g2Xs*w42f)
    g3e_pr<-(1/wfbar)*(0.5* g1e*g3Xs*w13f + 0.5*(r)* g1e*g4Xs*w14f + 0.5*(1-r)* g2e*g3Xs*w23f + 0.5* g3e*g1Xs*w31f + 0.5*(1-r)* g3e*g2Xs*w32f + 1* g3e*g3Xs*w33f + 0.5* g3e*g4Xs*w34f + 0.5*(r)* g4e*g1Xs*w41f + 0.5* g4e*g3Xs*w43f)
    g4e_pr<-(1/wfbar)*(0.5*(1-r)* g1e*g4Xs*w14f + 0.5*(r)* g2e*g3Xs*w23f + 0.5* g2e*g4Xs*w24f + 0.5*(r)* g3e*g2Xs*w32f + 0.5* g3e*g4Xs*w34f + 0.5*(1-r)* g4e*g1Xs*w41f + 0.5* g4e*g2Xs*w42f + 0.5* g4e*g3Xs*w43f + 1* g4e*g4Xs*w44f)
    
    wmbar<-(g1e*g1Ys*w11Ym)+(g1e*g2Ys*w12Ym)+(g1e*g3Ys*w13Ym)+(g1e*g4Ys*w14Ym)+(g2e*g1Ys*w21Ym)+(g2e*g2Ys*w22Ym)+(g2e*g3Ys*w23Ym)+(g2e*g4Ys*w24Ym)+(g3e*g1Ys*w31Ym)+(g3e*g2Ys*w32Ym)+(g3e*g3Ys*w33Ym)+(g3e*g4Ys*w34Ym)+(g4e*g1Ys*w41Ym)+(g4e*g2Ys*w42Ym)+(g4e*g3Ys*w43Ym)+(g4e*g4Ys*w44Ym)+(g1e*g1Os*w11Om)+(g1e*g2Os*w12Om)+(g1e*g3Os*w13Om)+(g1e*g4Os*w14Om)+(g2e*g1Os*w21Om)+(g2e*g2Os*w22Om)+(g2e*g3Os*w23Om)+(g2e*g4Os*w24Om)+(g3e*g1Os*w31Om)+(g3e*g2Os*w32Om)+(g3e*g3Os*w33Om)+(g3e*g4Os*w34Om)+(g4e*g1Os*w41Om)+(g4e*g2Os*w42Om)+(g4e*g3Os*w43Om)+(g4e*g4Os*w44Om)
    g1Xs_pr<- (1/wmbar)*(0.5*1* g1e*g1Ys*w11Ym + 0.5*0.5* g1e*g2Ys*w12Ym + (0.5-kD)*0.5* g1e*g3Ys*w13Ym + (0.5-kD)*0.5*(1-r)* g1e*g4Ys*w14Ym + 0.5*0.5* g2e*g1Ys*w21Ym + (0.5-kD)*0.5*(r)* g2e*g3Ys*w23Ym + (0.5-kD)*0.5* g3e*g1Ys*w31Ym + (0.5-kD)*0.5*(r)* g3e*g2Ys*w32Ym + (0.5-kD)*0.5*(1-r)* g4e*g1Ys*w41Ym +0.5*1* g1e*g1Os*w11Om + 0.5*0.5* g1e*g2Os*w12Om + (0.5-kD)*0.5* g1e*g3Os*w13Om + (0.5-kD)*0.5*(1-r)* g1e*g4Os*w14Om + 0.5*0.5* g2e*g1Os*w21Om + (0.5-kD)*0.5*(r)* g2e*g3Os*w23Om + (0.5-kD)*0.5* g3e*g1Os*w31Om + (0.5-kD)*0.5*(r)* g3e*g2Os*w32Om + (0.5-kD)*0.5*(1-r)* g4e*g1Os*w41Om)
    g2Xs_pr<- (1/wmbar)*(0.5*0.5* g1e*g2Ys*w12Ym + (0.5-kD)*0.5*(r)* g1e*g4Ys*w14Ym + 0.5*0.5* g2e*g1Ys*w21Ym + 0.5*1* g2e*g2Ys*w22Ym + (0.5-kD)*0.5*(1-r)* g2e*g3Ys*w23Ym + (0.5-kD)*0.5* g2e*g4Ys*w24Ym + (0.5-kD)*0.5*(1-r)* g3e*g2Ys*w32Ym + (0.5-kD)*0.5*(r)* g4e*g1Ys*w41Ym + (0.5-kD)*0.5* g4e*g2Ys*w42Ym + 0.5*0.5* g1e*g2Os*w12Om + (0.5-kD)*0.5*(r)* g1e*g4Os*w14Om + (0.5-kD)*0.5* g2e*g1Os*w21Om + 0.5*1* g2e*g2Os*w22Om + (0.5-kD)*0.5*(1-r)* g2e*g3Os*w23Om + (0.5-kD)*0.5* g2e*g4Os*w24Om + (0.5-kD)*0.5*(1-r)* g3e*g2Os*w32Om + (0.5-kD)*0.5*(r)* g4e*g1Os*w41Om + (0.5-kD)*0.5* g4e*g2Os*w42Om)
    g3Xs_pr<- (1/wmbar)*((0.5-kD)*0.5* g1e*g3Ys*w13Ym + (0.5-kD)*0.5*(r)* g1e*g4Ys*w14Ym + (0.5-kD)*0.5*(1-r)* g2e*g3Ys*w23Ym + (0.5-kD)*0.5* g3e*g1Ys*w31Ym + (0.5-kD)*0.5*(1-r)* g3e*g2Ys*w32Ym + (0.5-kD)*1* g3e*g3Ys*w33Ym + (0.5-kD)*0.5* g3e*g4Ys*w34Ym + (0.5-kD)*0.5*(r)* g4e*g1Ys*w41Ym + (0.5-kD)*0.5* g4e*g3Ys*w43Ym + (0.5-kD)*0.5* g1e*g3Os*w13Om + (0.5-kD)*0.5*(r)* g1e*g4Os*w14Om + (0.5-kD)*0.5*(1-r)* g2e*g3Os*w23Om + (0.5-kD)*0.5* g3e*g1Os*w31Om + (0.5-kD)*0.5*(1-r)* g3e*g2Os*w32Om + (0.5-kD)*1* g3e*g3Os*w33Om + (0.5-kD)*0.5* g3e*g4Os*w34Om + (0.5-kD)*0.5*(r)* g4e*g1Os*w41Om + (0.5-kD)*0.5* g4e*g3Os*w43Om)
    g4Xs_pr<- (1/wmbar)*((0.5-kD)*0.5*(1-r)* g1e*g4Ys*w14Ym + (0.5-kD)*0.5*(r)* g2e*g3Ys*w23Ym + (0.5-kD)*0.5* g2e*g4Ys*w24Ym + (0.5-kD)*0.5*(r)* g3e*g2Ys*w32Ym + (0.5-kD)*0.5* g3e*g4Ys*w34Ym + (0.5-kD)*0.5*(1-r)* g4e*g1Ys*w41Ym + (0.5-kD)*0.5* g4e*g2Ys*w42Ym + (0.5-kD)*0.5* g4e*g3Ys*w43Ym + (0.5-kD)*1* g4e*g4Ys*w44Ym + (0.5-kD)*0.5*(1-r)* g1e*g4Os*w14Om + (0.5-kD)*0.5*(r)* g2e*g3Os*w23Om + (0.5-kD)*0.5* g2e*g4Os*w24Om + (0.5-kD)*0.5*(r)* g3e*g2Os*w32Om + (0.5-kD)*0.5* g3e*g4Os*w34Om + (0.5-kD)*0.5*(1-r)* g4e*g1Os*w41Om + (0.5-kD)*0.5* g4e*g2Os*w42Om + (0.5-kD)*0.5* g4e*g3Os*w43Om + (0.5-kD)*1* g4e*g4Os*w44Om)
    g1Ys_pr<- (1/wmbar)*(0.5*1* g1e*g1Ys*w11Ym + 0.5*0.5* g1e*g2Ys*w12Ym + 0.5*0.5* g1e*g3Ys*w13Ym + 0.5*0.5*(1-r)* g1e*g4Ys*w14Ym + 0.5*0.5* g2e*g1Ys*w21Ym + 0.5*0.5*(r)* g2e*g3Ys*w23Ym + 0.5*0.5* g3e*g1Ys*w31Ym + 0.5*0.5*(r)* g3e*g2Ys*w32Ym + 0.5*0.5*(1-r)* g4e*g1Ys*w41Ym)
    g2Ys_pr<- (1/wmbar)*(0.5*0.5* g1e*g2Ys*w12Ym + 0.5*0.5*(r)* g1e*g4Ys*w14Ym + 0.5*0.5* g2e*g1Ys*w21Ym + 0.5*1* g2e*g2Ys*w22Ym + 0.5*0.5*(1-r)* g2e*g3Ys*w23Ym + 0.5*0.5* g2e*g4Ys*w24Ym + 0.5*0.5*(1-r)* g3e*g2Ys*w32Ym + 0.5*0.5*(r)* g4e*g1Ys*w41Ym + 0.5*0.5* g4e*g2Ys*w42Ym)
    g3Ys_pr<- (1/wmbar)*(0.5*0.5* g1e*g3Ys*w13Ym + 0.5*0.5*(r)* g1e*g4Ys*w14Ym + 0.5*0.5*(1-r)* g2e*g3Ys*w23Ym + 0.5*0.5* g3e*g1Ys*w31Ym + 0.5*0.5*(1-r)* g3e*g2Ys*w32Ym + 0.5*1* g3e*g3Ys*w33Ym + 0.5*0.5* g3e*g4Ys*w34Ym + 0.5*0.5*(r)* g4e*g1Ys*w41Ym + 0.5*0.5* g4e*g3Ys*w43Ym)
    g4Ys_pr<- (1/wmbar)*(0.5*0.5*(1-r)* g1e*g4Ys*w14Ym + 0.5*0.5*(r)* g2e*g3Ys*w23Ym + 0.5*0.5* g2e*g4Ys*w24Ym + 0.5*0.5*(r)* g3e*g2Ys*w32Ym + 0.5*0.5* g3e*g4Ys*w34Ym + 0.5*0.5*(1-r)* g4e*g1Ys*w41Ym + 0.5*0.5* g4e*g2Ys*w42Ym + 0.5*0.5* g4e*g3Ys*w43Ym + 0.5*1* g4e*g4Ys*w44Ym)
    
    g1Os_pr<- (1/wmbar)*(0*1* g1e*g1Ys*w11Ym + 0*0.5* g1e*g2Ys*w12Ym + kD*0.5* g1e*g3Ys*w13Ym + kD*0.5*(1-r)* g1e*g4Ys*w14Ym + 0*0.5* g2e*g1Ys*w21Ym + kD*0.5*(r)* g2e*g3Ys*w23Ym + kD*0.5* g3e*g1Ys*w31Ym + kD*0.5*(r)* g3e*g2Ys*w32Ym + kD*0.5*(1-r)* g4e*g1Ys*w41Ym +0.5*1* g1e*g1Os*w11Ym + 0.5*0.5* g1e*g2Os*w12Om + (0.5+kD)*0.5* g1e*g3Os*w13Om + (0.5+kD)*0.5*(1-r)* g1e*g4Os*w14Om + 0.5*0.5* g2e*g1Os*w21Om + (0.5+kD)*0.5*(r)* g2e*g3Os*w23Om + (0.5+kD)*0.5* g3e*g1Os*w31Om + (0.5+kD)*0.5*(r)* g3e*g2Os*w32Om + (0.5+kD)*0.5*(1-r)* g4e*g1Os*w41Om)
    g2Os_pr<- (1/wmbar)*(0*0.5* g1e*g2Ys*w12Ym + kD*0.5*(r)* g1e*g4Ys*w14Ym + 0*0.5* g2e*g1Ys*w21Ym + 0*1* g2e*g2Ys*w22Ym + kD*0.5*(1-r)* g2e*g3Ys*w23Ym + kD*0.5* g2e*g4Ys*w24Ym + kD*0.5*(1-r)* g3e*g2Ys*w32Ym + kD*0.5*(r)* g4e*g1Ys*w41Ym + kD*0.5* g4e*g2Ys*w42Ym + 0.5*0.5* g1e*g2Os*w12Ym + (0.5+kD)*0.5*(r)* g1e*g4Os*w14Om + (0.5+kD)*0.5* g2e*g1Os*w21Om + 0.5*1* g2e*g2Os*w22Om + (0.5+kD)*0.5*(1-r)* g2e*g3Os*w23Om + (0.5+kD)*0.5* g2e*g4Os*w24Om + (0.5+kD)*0.5*(1-r)* g3e*g2Os*w32Om + (0.5+kD)*0.5*(r)* g4e*g1Os*w41Om + (0.5+kD)*0.5* g4e*g2Os*w42Om)
    g3Os_pr<- (1/wmbar)*(kD*0.5* g1e*g3Ys*w13Ym + kD*0.5*(r)* g1e*g4Ys*w14Ym + kD*0.5*(1-r)* g2e*g3Ys*w23Ym + kD*0.5* g3e*g1Ys*w31Ym + kD*0.5*(1-r)* g3e*g2Ys*w32Ym + kD*1* g3e*g3Ys*w33Ym + kD*0.5* g3e*g4Ys*w34Ym + kD*0.5*(r)* g4e*g1Ys*w41Ym + kD*0.5* g4e*g3Ys*w43Ym + (0.5+kD)*0.5* g1e*g3Os*w13Om + (0.5+kD)*0.5*(r)* g1e*g4Os*w14Om + (0.5+kD)*0.5*(1-r)* g2e*g3Os*w23Om + (0.5+kD)*0.5* g3e*g1Os*w31Om + (0.5+kD)*0.5*(1-r)* g3e*g2Os*w32Om + (0.5+kD)*1* g3e*g3Os*w33Om + (0.5+kD)*0.5* g3e*g4Os*w34Om + (0.5+kD)*0.5*(r)* g4e*g1Os*w41Om + (0.5+kD)*0.5* g4e*g3Os*w43Om)
    g4Os_pr<- (1/wmbar)*(kD*0.5*(1-r)* g1e*g4Ys*w14Ym + kD*0.5*(r)* g2e*g3Ys*w23Ym + kD*0.5* g2e*g4Ys*w24Ym + kD*0.5*(r)* g3e*g2Ys*w32Ym + kD*0.5* g3e*g4Ys*w34Ym + kD*0.5*(1-r)* g4e*g1Ys*w41Ym + kD*0.5* g4e*g2Ys*w42Ym + kD*0.5* g4e*g3Ys*w43Ym + kD*1* g4e*g4Ys*w44Ym + (0.5+kD)*0.5*(1-r)* g1e*g4Os*w14Om + (0.5+kD)*0.5*(r)* g2e*g3Os*w23Om + (0.5+kD)*0.5* g2e*g4Os*w24Om + (0.5+kD)*0.5*(r)* g3e*g2Os*w32Om + (0.5+kD)*0.5* g3e*g4Os*w34Om + (0.5+kD)*0.5*(1-r)* g4e*g1Os*w41Om + (0.5+kD)*0.5* g4e*g2Os*w42Om + (0.5+kD)*0.5* g4e*g3Os*w43Om + (0.5+kD)*1* g4e*g4Os*w44Om)
    
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
  }
  
  z <- original_z
  kD <- original_kD
  
  #RECURSION 
  
  g1e <- g1e
  g2e <- g2e
  g3e <- g3e
  g4e <- g4e
  
  g1Xs <- g1Xs
  g2Xs <- g2Xs
  g3Xs <- g3Xs
  g4Xs <- g4Xs
  
  
  g3Ys <- g1Ys_pr*0.001
  g4Ys <- g2Ys_pr*0.001
  g1Ys <- g1Ys_pr-g1Ys_pr*0.001
  g2Ys <- g2Ys_pr-g2Ys_pr*0.001
  
  g1Os <- 0.0
  g2Os <- 0.0
  g3Os <- 0.0
  g4Os <- 0.0
  
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
  v12m <- 1-hA*s
  v21m <- 1-hA*s
  v22m <- 1
  
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
    # Egg Haplotypes 
    wfbar<-(g1e*g1Xs*w11f)+(g1e*g2Xs*w12f)+(g1e*g3Xs*w13f)+(g1e*g4Xs*w14f)+(g2e*g1Xs*w21f)+(g2e*g2Xs*w22f)+(g2e*g3Xs*w23f)+(g2e*g4Xs*w24f)+(g3e*g1Xs*w31f)+(g3e*g2Xs*w32f)+(g3e*g3Xs*w33f)+(g3e*g4Xs*w34f)+(g4e*g1Xs*w41f)+(g4e*g2Xs*w42f)+(g4e*g3Xs*w43f)+(g4e*g4Xs*w44f)
    g1e_pr<-(1/wfbar)*(1* g1e*g1Xs*w11f + 0.5* g1e*g2Xs*w12f + 0.5* g1e*g3Xs*w13f + 0.5*(1-r)* g1e*g4Xs*w14f + 0.5* g2e*g1Xs*w21f + 0.5*(r)* g2e*g3Xs*w23f + 0.5* g3e*g1Xs*w31f + 0.5*(r)* g3e*g2Xs*w32f + 0.5*(1-r)* g4e*g1Xs*w41f)
    g2e_pr<-(1/wfbar)*(0.5* g1e*g2Xs*w12f + 0.5*(r)* g1e*g4Xs*w14f + 0.5* g2e*g1Xs*w21f + 1* g2e*g2Xs*w22f + 0.5*(1-r)* g2e*g3Xs*w23f + 0.5* g2e*g4Xs*w24f + 0.5*(1-r)* g3e*g2Xs*w32f + 0.5*(r)* g4e*g1Xs*w41f + 0.5* g4e*g2Xs*w42f)
    g3e_pr<-(1/wfbar)*(0.5* g1e*g3Xs*w13f + 0.5*(r)* g1e*g4Xs*w14f + 0.5*(1-r)* g2e*g3Xs*w23f + 0.5* g3e*g1Xs*w31f + 0.5*(1-r)* g3e*g2Xs*w32f + 1* g3e*g3Xs*w33f + 0.5* g3e*g4Xs*w34f + 0.5*(r)* g4e*g1Xs*w41f + 0.5* g4e*g3Xs*w43f)
    g4e_pr<-(1/wfbar)*(0.5*(1-r)* g1e*g4Xs*w14f + 0.5*(r)* g2e*g3Xs*w23f + 0.5* g2e*g4Xs*w24f + 0.5*(r)* g3e*g2Xs*w32f + 0.5* g3e*g4Xs*w34f + 0.5*(1-r)* g4e*g1Xs*w41f + 0.5* g4e*g2Xs*w42f + 0.5* g4e*g3Xs*w43f + 1* g4e*g4Xs*w44f)
    
    wmbar<-(g1e*g1Ys*w11Ym)+(g1e*g2Ys*w12Ym)+(g1e*g3Ys*w13Ym)+(g1e*g4Ys*w14Ym)+(g2e*g1Ys*w21Ym)+(g2e*g2Ys*w22Ym)+(g2e*g3Ys*w23Ym)+(g2e*g4Ys*w24Ym)+(g3e*g1Ys*w31Ym)+(g3e*g2Ys*w32Ym)+(g3e*g3Ys*w33Ym)+(g3e*g4Ys*w34Ym)+(g4e*g1Ys*w41Ym)+(g4e*g2Ys*w42Ym)+(g4e*g3Ys*w43Ym)+(g4e*g4Ys*w44Ym)+(g1e*g1Os*w11Om)+(g1e*g2Os*w12Om)+(g1e*g3Os*w13Om)+(g1e*g4Os*w14Om)+(g2e*g1Os*w21Om)+(g2e*g2Os*w22Om)+(g2e*g3Os*w23Om)+(g2e*g4Os*w24Om)+(g3e*g1Os*w31Om)+(g3e*g2Os*w32Om)+(g3e*g3Os*w33Om)+(g3e*g4Os*w34Om)+(g4e*g1Os*w41Om)+(g4e*g2Os*w42Om)+(g4e*g3Os*w43Om)+(g4e*g4Os*w44Om)
    g1Xs_pr<- (1/wmbar)*(0.5*1* g1e*g1Ys*w11Ym + 0.5*0.5* g1e*g2Ys*w12Ym + (0.5-kD)*0.5* g1e*g3Ys*w13Ym + (0.5-kD)*0.5*(1-r)* g1e*g4Ys*w14Ym + 0.5*0.5* g2e*g1Ys*w21Ym + (0.5-kD)*0.5*(r)* g2e*g3Ys*w23Ym + (0.5-kD)*0.5* g3e*g1Ys*w31Ym + (0.5-kD)*0.5*(r)* g3e*g2Ys*w32Ym + (0.5-kD)*0.5*(1-r)* g4e*g1Ys*w41Ym +0.5*1* g1e*g1Os*w11Om + 0.5*0.5* g1e*g2Os*w12Om + (0.5-kD)*0.5* g1e*g3Os*w13Om + (0.5-kD)*0.5*(1-r)* g1e*g4Os*w14Om + 0.5*0.5* g2e*g1Os*w21Om + (0.5-kD)*0.5*(r)* g2e*g3Os*w23Om + (0.5-kD)*0.5* g3e*g1Os*w31Om + (0.5-kD)*0.5*(r)* g3e*g2Os*w32Om + (0.5-kD)*0.5*(1-r)* g4e*g1Os*w41Om)
    g2Xs_pr<- (1/wmbar)*(0.5*0.5* g1e*g2Ys*w12Ym + (0.5-kD)*0.5*(r)* g1e*g4Ys*w14Ym + 0.5*0.5* g2e*g1Ys*w21Ym + 0.5*1* g2e*g2Ys*w22Ym + (0.5-kD)*0.5*(1-r)* g2e*g3Ys*w23Ym + (0.5-kD)*0.5* g2e*g4Ys*w24Ym + (0.5-kD)*0.5*(1-r)* g3e*g2Ys*w32Ym + (0.5-kD)*0.5*(r)* g4e*g1Ys*w41Ym + (0.5-kD)*0.5* g4e*g2Ys*w42Ym + 0.5*0.5* g1e*g2Os*w12Om + (0.5-kD)*0.5*(r)* g1e*g4Os*w14Om + (0.5-kD)*0.5* g2e*g1Os*w21Om + 0.5*1* g2e*g2Os*w22Om + (0.5-kD)*0.5*(1-r)* g2e*g3Os*w23Om + (0.5-kD)*0.5* g2e*g4Os*w24Om + (0.5-kD)*0.5*(1-r)* g3e*g2Os*w32Om + (0.5-kD)*0.5*(r)* g4e*g1Os*w41Om + (0.5-kD)*0.5* g4e*g2Os*w42Om)
    g3Xs_pr<- (1/wmbar)*((0.5-kD)*0.5* g1e*g3Ys*w13Ym + (0.5-kD)*0.5*(r)* g1e*g4Ys*w14Ym + (0.5-kD)*0.5*(1-r)* g2e*g3Ys*w23Ym + (0.5-kD)*0.5* g3e*g1Ys*w31Ym + (0.5-kD)*0.5*(1-r)* g3e*g2Ys*w32Ym + (0.5-kD)*1* g3e*g3Ys*w33Ym + (0.5-kD)*0.5* g3e*g4Ys*w34Ym + (0.5-kD)*0.5*(r)* g4e*g1Ys*w41Ym + (0.5-kD)*0.5* g4e*g3Ys*w43Ym + (0.5-kD)*0.5* g1e*g3Os*w13Om + (0.5-kD)*0.5*(r)* g1e*g4Os*w14Om + (0.5-kD)*0.5*(1-r)* g2e*g3Os*w23Om + (0.5-kD)*0.5* g3e*g1Os*w31Om + (0.5-kD)*0.5*(1-r)* g3e*g2Os*w32Om + (0.5-kD)*1* g3e*g3Os*w33Om + (0.5-kD)*0.5* g3e*g4Os*w34Om + (0.5-kD)*0.5*(r)* g4e*g1Os*w41Om + (0.5-kD)*0.5* g4e*g3Os*w43Om)
    g4Xs_pr<- (1/wmbar)*((0.5-kD)*0.5*(1-r)* g1e*g4Ys*w14Ym + (0.5-kD)*0.5*(r)* g2e*g3Ys*w23Ym + (0.5-kD)*0.5* g2e*g4Ys*w24Ym + (0.5-kD)*0.5*(r)* g3e*g2Ys*w32Ym + (0.5-kD)*0.5* g3e*g4Ys*w34Ym + (0.5-kD)*0.5*(1-r)* g4e*g1Ys*w41Ym + (0.5-kD)*0.5* g4e*g2Ys*w42Ym + (0.5-kD)*0.5* g4e*g3Ys*w43Ym + (0.5-kD)*1* g4e*g4Ys*w44Ym + (0.5-kD)*0.5*(1-r)* g1e*g4Os*w14Om + (0.5-kD)*0.5*(r)* g2e*g3Os*w23Om + (0.5-kD)*0.5* g2e*g4Os*w24Om + (0.5-kD)*0.5*(r)* g3e*g2Os*w32Om + (0.5-kD)*0.5* g3e*g4Os*w34Om + (0.5-kD)*0.5*(1-r)* g4e*g1Os*w41Om + (0.5-kD)*0.5* g4e*g2Os*w42Om + (0.5-kD)*0.5* g4e*g3Os*w43Om + (0.5-kD)*1* g4e*g4Os*w44Om)
    g1Ys_pr<- (1/wmbar)*(0.5*1* g1e*g1Ys*w11Ym + 0.5*0.5* g1e*g2Ys*w12Ym + 0.5*0.5* g1e*g3Ys*w13Ym + 0.5*0.5*(1-r)* g1e*g4Ys*w14Ym + 0.5*0.5* g2e*g1Ys*w21Ym + 0.5*0.5*(r)* g2e*g3Ys*w23Ym + 0.5*0.5* g3e*g1Ys*w31Ym + 0.5*0.5*(r)* g3e*g2Ys*w32Ym + 0.5*0.5*(1-r)* g4e*g1Ys*w41Ym)
    g2Ys_pr<- (1/wmbar)*(0.5*0.5* g1e*g2Ys*w12Ym + 0.5*0.5*(r)* g1e*g4Ys*w14Ym + 0.5*0.5* g2e*g1Ys*w21Ym + 0.5*1* g2e*g2Ys*w22Ym + 0.5*0.5*(1-r)* g2e*g3Ys*w23Ym + 0.5*0.5* g2e*g4Ys*w24Ym + 0.5*0.5*(1-r)* g3e*g2Ys*w32Ym + 0.5*0.5*(r)* g4e*g1Ys*w41Ym + 0.5*0.5* g4e*g2Ys*w42Ym)
    g3Ys_pr<- (1/wmbar)*(0.5*0.5* g1e*g3Ys*w13Ym + 0.5*0.5*(r)* g1e*g4Ys*w14Ym + 0.5*0.5*(1-r)* g2e*g3Ys*w23Ym + 0.5*0.5* g3e*g1Ys*w31Ym + 0.5*0.5*(1-r)* g3e*g2Ys*w32Ym + 0.5*1* g3e*g3Ys*w33Ym + 0.5*0.5* g3e*g4Ys*w34Ym + 0.5*0.5*(r)* g4e*g1Ys*w41Ym + 0.5*0.5* g4e*g3Ys*w43Ym)
    g4Ys_pr<- (1/wmbar)*(0.5*0.5*(1-r)* g1e*g4Ys*w14Ym + 0.5*0.5*(r)* g2e*g3Ys*w23Ym + 0.5*0.5* g2e*g4Ys*w24Ym + 0.5*0.5*(r)* g3e*g2Ys*w32Ym + 0.5*0.5* g3e*g4Ys*w34Ym + 0.5*0.5*(1-r)* g4e*g1Ys*w41Ym + 0.5*0.5* g4e*g2Ys*w42Ym + 0.5*0.5* g4e*g3Ys*w43Ym + 0.5*1* g4e*g4Ys*w44Ym)
    
    g1Os_pr<- (1/wmbar)*(0*1* g1e*g1Ys*w11Ym + 0*0.5* g1e*g2Ys*w12Ym + kD*0.5* g1e*g3Ys*w13Ym + kD*0.5*(1-r)* g1e*g4Ys*w14Ym + 0*0.5* g2e*g1Ys*w21Ym + kD*0.5*(r)* g2e*g3Ys*w23Ym + kD*0.5* g3e*g1Ys*w31Ym + kD*0.5*(r)* g3e*g2Ys*w32Ym + kD*0.5*(1-r)* g4e*g1Ys*w41Ym +0.5*1* g1e*g1Os*w11Ym + 0.5*0.5* g1e*g2Os*w12Om + (0.5+kD)*0.5* g1e*g3Os*w13Om + (0.5+kD)*0.5*(1-r)* g1e*g4Os*w14Om + 0.5*0.5* g2e*g1Os*w21Om + (0.5+kD)*0.5*(r)* g2e*g3Os*w23Om + (0.5+kD)*0.5* g3e*g1Os*w31Om + (0.5+kD)*0.5*(r)* g3e*g2Os*w32Om + (0.5+kD)*0.5*(1-r)* g4e*g1Os*w41Om)
    g2Os_pr<- (1/wmbar)*(0*0.5* g1e*g2Ys*w12Ym + kD*0.5*(r)* g1e*g4Ys*w14Ym + 0*0.5* g2e*g1Ys*w21Ym + 0*1* g2e*g2Ys*w22Ym + kD*0.5*(1-r)* g2e*g3Ys*w23Ym + kD*0.5* g2e*g4Ys*w24Ym + kD*0.5*(1-r)* g3e*g2Ys*w32Ym + kD*0.5*(r)* g4e*g1Ys*w41Ym + kD*0.5* g4e*g2Ys*w42Ym + 0.5*0.5* g1e*g2Os*w12Ym + (0.5+kD)*0.5*(r)* g1e*g4Os*w14Om + (0.5+kD)*0.5* g2e*g1Os*w21Om + 0.5*1* g2e*g2Os*w22Om + (0.5+kD)*0.5*(1-r)* g2e*g3Os*w23Om + (0.5+kD)*0.5* g2e*g4Os*w24Om + (0.5+kD)*0.5*(1-r)* g3e*g2Os*w32Om + (0.5+kD)*0.5*(r)* g4e*g1Os*w41Om + (0.5+kD)*0.5* g4e*g2Os*w42Om)
    g3Os_pr<- (1/wmbar)*(kD*0.5* g1e*g3Ys*w13Ym + kD*0.5*(r)* g1e*g4Ys*w14Ym + kD*0.5*(1-r)* g2e*g3Ys*w23Ym + kD*0.5* g3e*g1Ys*w31Ym + kD*0.5*(1-r)* g3e*g2Ys*w32Ym + kD*1* g3e*g3Ys*w33Ym + kD*0.5* g3e*g4Ys*w34Ym + kD*0.5*(r)* g4e*g1Ys*w41Ym + kD*0.5* g4e*g3Ys*w43Ym + (0.5+kD)*0.5* g1e*g3Os*w13Om + (0.5+kD)*0.5*(r)* g1e*g4Os*w14Om + (0.5+kD)*0.5*(1-r)* g2e*g3Os*w23Om + (0.5+kD)*0.5* g3e*g1Os*w31Om + (0.5+kD)*0.5*(1-r)* g3e*g2Os*w32Om + (0.5+kD)*1* g3e*g3Os*w33Om + (0.5+kD)*0.5* g3e*g4Os*w34Om + (0.5+kD)*0.5*(r)* g4e*g1Os*w41Om + (0.5+kD)*0.5* g4e*g3Os*w43Om)
    g4Os_pr<- (1/wmbar)*(kD*0.5*(1-r)* g1e*g4Ys*w14Ym + kD*0.5*(r)* g2e*g3Ys*w23Ym + kD*0.5* g2e*g4Ys*w24Ym + kD*0.5*(r)* g3e*g2Ys*w32Ym + kD*0.5* g3e*g4Ys*w34Ym + kD*0.5*(1-r)* g4e*g1Ys*w41Ym + kD*0.5* g4e*g2Ys*w42Ym + kD*0.5* g4e*g3Ys*w43Ym + kD*1* g4e*g4Ys*w44Ym + (0.5+kD)*0.5*(1-r)* g1e*g4Os*w14Om + (0.5+kD)*0.5*(r)* g2e*g3Os*w23Om + (0.5+kD)*0.5* g2e*g4Os*w24Om + (0.5+kD)*0.5*(r)* g3e*g2Os*w32Om + (0.5+kD)*0.5* g3e*g4Os*w34Om + (0.5+kD)*0.5*(1-r)* g4e*g1Os*w41Om + (0.5+kD)*0.5* g4e*g2Os*w42Om + (0.5+kD)*0.5* g4e*g3Os*w43Om + (0.5+kD)*1* g4e*g4Os*w44Om)
    
    
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
  return(OFreq[d]/(OFreq[d]+YFreq[d]))
}

#Visualization 
n <- 200
kvector<-c(n)
wvector<-c(n)
zvector<-c(n)
equalvector<-matrix (0,n^2,2)

i<-1 

for (a in 1:n){
  kD <- a/(2*n)
  for (b in 1:n){
    w <- 0+b/(4*n)
    z<-SAfunction(1000,kD,0.5,0.5,0.5,0.01,0.001,w)
    if(0.45<= z & z <0.55){
      equalvector[i,1]<-kD
      equalvector[i,2]<-w
    } 
    kvector[i]<-kD
    wvector[i]<-w
    zvector[i]<-z
    i<-i+1
  }
  print(a)
}

library(tibble)
d <- tibble(k = kvector, w = wvector, z = zvector)
write.table(d, file = "data_S5A.txt", row.names = F, col.names = T)



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
  return(OFreq[d]/(OFreq[d]+YFreq[d]))
}

#Visualization 

n <- 200
kvector<-c(n)
wvector<-c(n)
zvector<-c(n)
equalvector<-matrix (0,n^2,2)

i<-1 

for (a in 1:n){
  kD <- a/(2*n)
  for (b in 1:n){
    w <- 0+b/(4*n)
    z<- ML_function(1000, kD, 0.25, 0.5, 0, 0.01, 0.01, w) 
    if(0.45<= z & z <0.55){
      equalvector[i,1]<-kD
      equalvector[i,2]<-w
    } 
    kvector[i]<-kD
    wvector[i]<-w
    zvector[i]<-z
    i<-i+1

  }
  print(a)
}
library(tibble)
d <- tibble(k = kvector, w = wvector, z = zvector)
write.table(d, file = "data_S5B.txt", row.names = F, col.names = T)
dA <- read.table("data_S5A.txt", T) %>% as_tibble()
dB <- read.table("data_S5B.txt", T) %>% as_tibble()

library(ggplot2)
library("viridis")
library(tidyverse)
library("cowplot")
pA <- dA %>% ggplot(aes(k, w, z = z)) + 
  geom_contour_filled(breaks = c(-0.01, 0.01, 0.99, 1.01), color = "black") +
  scale_fill_viridis(discrete = T, alpha = 0.8,
                     labels = c("Y-chromosome", "Y/0 polymorphism", "0-chromosome")) + 
  scale_x_continuous(limits = c(min(dA$k), max(dA$k)), expand = c(0,0)) + 
  scale_y_continuous(limits = c(min(dA$w), max(dA$w)), expand = c(0,0)) + 
  labs(fill = "Sex determination", 
       x = bquote(italic(trans)~"drive ("*italic(k)[D]*")"), 
       y = bquote("0-chromosome fitness cost ("~italic('w')*")")) +
  theme(panel.background = element_rect(fill = "white"))


pB <- dB %>% ggplot(aes(k, w, z = z)) + 
  geom_contour_filled(breaks = c(-0.01, 0.01, 0.99, 1.01), color = "black") +
  scale_fill_viridis(discrete = T, alpha = 0.8,
                     labels = c("Y-chromosome", "Y/0 polymorphism", "0-chromosome")) + 
  scale_x_continuous(limits = c(min(dA$k), max(dA$k)), expand = c(0,0)) + 
  scale_y_continuous(limits = c(min(dA$w), max(dA$w)), expand = c(0,0)) + 
  labs(fill = "Sex determination", 
       x = bquote(italic(trans)~"drive ("*italic(k)[D]*")"), 
       y = "") +
  theme(panel.background = element_rect(fill = "white"))
pB

legend <- get_legend(pA +
                       theme(legend.key.size = unit(0.5, "cm"),
                             legend.key.height = unit(0.35, "cm"))
)


top <- plot_grid(pA + theme(legend.position = "none"),
               pB + theme(legend.position = "none"),
               labels = LETTERS[1:2])
top
p <- plot_grid(top, legend, ncol = 1, rel_heights = c(1, 0.2))
p

pdf(file = "2025_04_25_S5.pdf", height = 4, width = 7)
p
dev.off()

png(file = "2025_04_25_S5.png", height = 4, width = 7, units = "in", res = 1000)
p
dev.off()
