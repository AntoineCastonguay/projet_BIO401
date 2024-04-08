corailstochastique <- function (R0,f,lambda=1/30,impact=2,dt=1,graph=FALSE){
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
  P <- 0.5
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
  
  
  while (t<tmax){
    
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
    dM <- ((a)*M*C - ((1)*P*M/(M+Y)) + (g)*M*Y)*dt - dommM
    dC <- ((r)*Y*C - (d)*C - (a)*M*C)*dt - dommC
    dP <- ((s)*P*(1 - (P/((((ca+cb*R)/kmax)*((delta*(M+Y))/(1+v*(M+Y))))))) -(f*s)*P)*dt 
    dR <- ((hG)*C*(3 - R) - (hE)*(1 - C)*(R - 1))*dt - dommR
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
  if (graph==TRUE){
    plot(seq(0,tmax,length.out = length(vecC)),vecC, type="l", col="red", ylim=c(0,1), xlab= "temps", ylab="valeurs")
    axis(side = 4,at= seq(0,1,0.1),labels = seq(1,3,0.2), )
    lines(seq(0,tmax,length.out = length(vecM)),vecM, col="blue")
    lines(seq(0,tmax,length.out = length(vecP)),vecP, col="green")
    lines(seq(0,tmax,length.out = length(vecR)),(vecR-1)/2, col="pink")
    legend("topright", legend = c("corail", "macroalgues","poissons","rugosité"),
           col = c("red","blue","green","pink"), lty = c(1, 1),
           lwd = 2, bg = "white",cex = 0.75)
  }
  return(c(tail(vecC,1),tail(vecR,1)))
}


#les fonctions pour faire les images

graph.imageC <- function(impact,freq,it=10,pas=0.1){
  
  matC <- matrix(0, nrow = 21,ncol = 26)
  matsd <- matrix(0, nrow = 21,ncol = 26)
  
  
  imp<-impact
  l<-freq
  for (i in 1:21) {
    for (j in 1:26){
      res <- vector(length = 0)
      for (z in 1:it) {
        R1 <- i*0.1+0.9
        #f0 <- j
        F0 <- j*0.04-0.04
        w <- corailstochastique(R0=R1,f=0.2, lambda = l,impact = imp,dt=pas, graph = F)[1]
        res<-c(res,w)
      }
      matC[i,j]<-mean(res)
      matsd[i,j]<-sd(res)
    }
  }
  
  
  #library(fields)
  
  
  x_values <- seq(0, 1, length.out = ncol(matC))
  y_values <- seq(1, 3, length.out = nrow(matC))
  # Plot the matrix as an image with different colored pixels
  
  par(mfrow=c(1,2))
  image.plot(x = x_values, y = y_values, t(matC), col = viridis(20), 
             xlab = "effort de pêche", ylab = "Rugosité initiale", 
             main=c("moyenne de C finale",paste(c("fréq:1/",1/l,", impact:",imp),collapse = "")))
  image.plot(x = x_values, y = y_values, t(matsd), col = viridis(20),
             xlab = "effort de pêche", ylab = "Rugosité initiale",
             main=c("moyenne de C finale",paste(c("fréq:1/",1/l,", impact:",imp),collapse = "")))
}
system.time(graph.imageC(1,1/9,5,0.25))
system.time(graph.imageC(2,1/9,5,0.1))
system.time(graph.imageC(3,1/9,5,0.1))
system.time(graph.imageC(1,1/15,5,0.1))
system.time(graph.imageC(2,1/15,5,0.1))
system.time(graph.imageC(3,1/15,5,0.1))
system.time(graph.imageC(1,1/30,5,0.1))
system.time(graph.imageC(2,1/30,5,0.1))
system.time(graph.imageC(3,1/30,5,0.1))


graph.imageR <- function(impact,freq,it=10,pas=0.1){
  
  matC <- matrix(0, nrow = 21,ncol = 26)
  matsd <- matrix(0, nrow = 21,ncol = 26)
  
  
  imp<-impact
  l<-freq
  for (i in 1:21) {
    for (j in 1:26){
      res <- vector(length = 0)
      for (z in 1:it) {
        R1 <- i*0.1+0.9
        f0 <- j*0.04-0.04
        w <- corailstochastique(R0=R1,f=f0, lambda = l,impact = imp,dt=pas, graph = F)
        res<-c(res,w)
      }
      matC[i,j]<-mean(res)
      matsd[i,j]<-sd(res)
    }
  }
  
  
  #library(fields)
  
  
  x_values <- seq(0, 1, length.out = ncol(matC))
  y_values <- seq(1, 3, length.out = nrow(matC))
  # Plot the matrix as an image with different colored pixels
  
  par(mfrow=c(1,2))
  image.plot(x = x_values, y = y_values, t(matC), col = viridis(20), 
             xlab = "effort de pêche", ylab = "Rugosité initiale", 
             main=c("moyenne de R finale",paste(c("fréq:1/",1/l,", impact:",imp),collapse = "")))
  image.plot(x = x_values, y = y_values, t(matsd), col = viridis(20),
             xlab = "effort de pêche", ylab = "Rugosité initiale",
             main=c("moyenne de R finale",paste(c("fréq:1/",1/l,", impact:",imp),collapse = "")))
}
system.time(graph.imageR(1,1/9,5,0.25))
system.time(graph.imageR(2,1/9,5,0.1))
system.time(graph.imageR(3,1/9,5,0.1))
system.time(graph.imageR(1,1/15,5,0.1))
system.time(graph.imageR(2,1/15,5,0.1))
system.time(graph.imageR(3,1/15,5,0.1))
system.time(graph.imageR(1,1/30,5,0.1))
system.time(graph.imageR(2,1/30,5,0.1))
system.time(graph.imageR(3,1/30,5,0.1))

