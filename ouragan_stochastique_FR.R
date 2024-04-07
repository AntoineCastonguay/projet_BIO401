#Ceci est un test pour un modèle avec les ouragans sans utiliser la fonction ode()

corailstochastique <- function (R0,f,lambda=1/30,impact=2,dt=1){
  #on instaure une étendue de temps
  tmax <- 100
  t <- 0
  #le temps au dernier ouragan
  to <- 0
  
  #on pose les valeurs de départ et les paramètres
  C <- 0.56
  M <- 0.01
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
  # plot(seq(0,tmax+dt,dt),vecC, type="l", col="red", ylim=c(0,1), xlab= "temps", ylab="valeurs")
  # axis(side = 4,at= seq(0,1,0.1),labels = seq(1,3,0.2), )
  # lines(seq(0,tmax+dt,dt),vecM, col="blue")
  # lines(seq(0,tmax+dt,dt),vecP, col="green")
  # lines(seq(0,tmax+dt,dt),(vecR-1)/2, col="pink")
  # legend("topright", legend = c("corail", "macroalgues","poissons","rugosité"), 
  #        col = c("red","blue","green","pink"), lty = c(1, 1), 
  #        lwd = 2, bg = "white",cex = 0.75)
  #print(vecM)
  return(tail(vecC,1))
}

corailstochastique(R0=2.4, f= 0.58, lambda = 1/30,impact = 2,dt=0.1)

matC <- matrix(0, nrow = 41,ncol = 51)
matsd <- matrix(0, nrow = 41,ncol = 51)

for (i in 1:41) {
  for (j in 1:51){
    res <- vector(length = 0)
    for (z in 1:10) {
      R1 <- i*0.05+0.9
      f0 <- j*0.02-0.02
      w <- corailstochastique(R0=R1,f=f0, lambda = 1/30,impact = 1,dt=1)
      res<-c(res,w)
    }
    matC[i,j]<-mean(res)
    matsd[i,j]<-sd(res)
  }
}

library(fields)

# mat<-matrix(0,nrow = 61,ncol = 51)
# mat[21:61,]<-matriceres
# colnames(mat)<-paste(seq(0,1,0.02))
# row.names(mat)<-paste(seq(0,3,0.05))

x_values <- seq(0, 1, length.out = ncol(matC))
y_values <- seq(1, 3, length.out = nrow(matC))

# Plot the matrix as an image with different colored pixels
image.plot(x = x_values, y = y_values, t(matC), col = viridis(20), xlab = "X", ylab = "Y", main="")



corailstochastique <- function (R0,f,lambda=1/30,impact=2,dt=1){
  #on instaure une étendue de temps
  tmax <- 100
  t <- 0
  #le temps au dernier ouragan
  to <- 0
  
  #on pose les valeurs de départ et les paramètres
  C <- 0.56
  M <- 0.01
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
  # plot(seq(0,tmax+dt,dt),vecC, type="l", col="red", ylim=c(0,1), xlab= "temps", ylab="valeurs")
  # axis(side = 4,at= seq(0,1,0.1),labels = seq(1,3,0.2), )
  # lines(seq(0,tmax+dt,dt),vecM, col="blue")
  # lines(seq(0,tmax+dt,dt),vecP, col="green")
  # lines(seq(0,tmax+dt,dt),(vecR-1)/2, col="pink")
  # legend("topright", legend = c("corail", "macroalgues","poissons","rugosité"), 
  #        col = c("red","blue","green","pink"), lty = c(1, 1), 
  #        lwd = 2, bg = "white",cex = 0.75)
  #print(vecM)
  return(tail(vecC,1))
}

corailstochastique(R0=2.4, f= 0.58, lambda = 1/30,impact = 2,dt=0.1)

matC <- matrix(0, nrow = 21,ncol = 26)
matsd <- matrix(0, nrow = 21,ncol = 26)

for (i in 1:21) {
  for (j in 1:26){
    res <- vector(length = 0)
    for (z in 1:10) {
      R1 <- i*0.1+0.9
      f0 <- j*0.04-0.04
      w <- corailstochastique(R0=R1,f=f0, lambda = 1/15,impact = 1,dt=1)
      res<-c(res,w)
    }
    matC[i,j]<-mean(res)
    matsd[i,j]<-sd(res)
  }
}

x_values <- seq(0, 1, length.out = ncol(matC))
y_values <- seq(1, 3, length.out = nrow(matC))
# Plot the matrix as an image with different colored pixels

image.plot(x = x_values, y = y_values, t(matC), col = viridis(20), xlab = "effort de pêche", ylab = "Rugosité", main="moyenne de C final")
image.plot(x = x_values, y = y_values, t(matsd), col = viridis(20), xlab = "effort de pêche", ylab = "Rugosité", main="écart-type de C final")
