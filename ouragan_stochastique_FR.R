#Ceci est un test pour un modèle avec les ouragans sans utiliser la fonction ode()

corailstochastique <- function (R0,f,lambda=1/30,impact=2,dt=1){
  #on instaure une étendue de temps
  tmax <- 100
  t <- 0
  #le temps au dernier ouragan
  to <- 0
  
  #on pose les valeurs de départ et les paramètres
  C <- 0.56
  M <- 0.1
  Y <- 0.43
  R <- R0
  P <- 1
  a<-0.1
  g<-0.8
  d<-0.44
  r<-1
  f<-f
  s<-0.49
  ca<--(3.21)
  cb<-3.65
  kmax<-17.745
  delta<-4.557
  v<-0.9877
  hG<-0.03
  hE<-0.01
  o <- 0.05
  mH <- 0.8
  
  
  #on initie les vecteurs de résultats
  vecC <- C
  vecM <- M
  vecY <- Y
  vecR <- R
  vecP <- P
  
  
  while (t<=tmax){
    
    #on effectue le tirage pour déterminer si il y a ou pas ouragan
    ouragan <- pmin(rpois(1,lambda*dt),1)
    #si oui ou non, on ajuste dq
    if (ouragan==1){
      
      #on calcule le temps depuis le dernier ouragan
      te <- t - to
      
      #on sauve le nouveau temps de dernier ouragan
      to <- t
      
      #on calcule F
      F <- pmin(0.0512 + 0.0176*te, 0.5)
      
      #différents dommage selon valeur d'impact
      if (impact==3){
        #on calcule les dommages d'un fort ouragan
        dommC <- C*(F + o*runif(1))
        dommM <- M*mH
        dommR <- (R-1)*0.5
      }
      if (impact==2){
        #on calcule les dommages d'un moyen ouragan
        dommC <- C*(F + o*runif(1))
        dommM <- M*mH
        dommR <- (R-1)*(F + o*runif(1))
      }
      if (impact==1){
        #on calcule les dommages d'un faible ouragan
        #dommC <- C*(F + o*runif(1))
        dommC <- 0
        dommM <- M*mH
        dommR <- (R-1)*0.05
      }
    }else{
      dommC <- 0
      dommM <- 0
      dommR <- 0
    }
    
    #on calcule les changements
    dM <- (a*dt)*M*C - ((1*dt)*P*M/(M+Y)) + (g*dt)*M*Y - dommM
    dC <- (r*dt)*Y*C - (d*dt)*C - (a*dt)*M*C - dommC
    dP <- s*dt*P*(1 - (P/((((ca+cb*R)/kmax)*((delta*(M+Y))/(1+v*(M+Y))))))) -(f*s*dt)*P 
    dR <- (hG*dt)*C*(3 - R) - (hE*dt)*(1 - C)*(R - 1) - dommR
    #dY <- 0 - dM - dC 
    
    #on effectue les changements
    M <- pmax(M + dM,0.0001)
    C <- pmax(C + dC,0.0001)
    #Y <- Y + dY
    Y <- 1 - M - C
    R <- R + dR
    P <- pmax(P + dP, 0.0001)
    
    #on sauve les nouvelles valeurs dans leur vecteurs
    vecC <- c(vecC,C)
    vecM <- c(vecM,M)
    vecY <- c(vecY,Y)
    vecR <- c(vecR,R)
    vecP <- c(vecP,P)
    
    #on fait avancer le temps
    t <- t + dt
  }
  plot(seq(0,tmax+dt,dt),vecC, type="l", col="red", ylim=c(0,1), xlab= "temps", ylab="valeurs")
  axis(side = 4,at= seq(0,1,0.1),labels = seq(1,3,0.2), )
  lines(seq(0,tmax+dt,dt),vecM, col="blue")
  lines(seq(0,tmax+dt,dt),vecP, col="green")
  lines(seq(0,tmax+dt,dt),(vecR-1)/2, col="pink")
  #print(vecM)
}


corailstochastique(R0=2.4, f= 0.1, lambda = 1/9,impact = 2,dt=0.1)

