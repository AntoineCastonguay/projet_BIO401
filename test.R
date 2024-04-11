corailstochastiquev2 <- function (R0,f,lambda=1/15,impact=1,dt=1,graph=FALSE){
  #on instaure une étendue de temps
  tmax <- 100
  t <- 0
  #le temps au dernier ouragan, en tenant compte que la simulation ne commence pas nécessairement juste après un ouragan
  to <- runif(1,0,25.5)
  te <- t - to
  
  #on pose les valeurs de départ et les paramètres
  C <- 0.8
  M <- 0.2
  Y <- 1 - C - M
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
      F <- pmin(0.0512 + 0.0176*(t/200+1)*te, 0.20)
      G <- (R-1)*0.1
      
      #différents dommage selon valeur d'impact
      
      #on calcule les dommages d'un fort ouragan
      dommC <- C*(F + o*runif(1))*impact
      dommM <- M*mH
      dommR <- (R-1)*(G + o*runif(1))*impact
      
    }else if (0 < te & te < 8) {
      te <- te + 1
      dommC <- C*0.08*(impact/2)
      dommM <- 0
      dommR <- 0
    }else{
      dommC <- 0
      dommM <- 0
      dommR <- 0
    }
    
    #on calcule les changements
    dC <- ((r)*Y*C - (d)*C - (a)*M*C)*dt - dommC
    dM <- ((a)*M*C - ((1)*P*M/(M+Y)) + (g)*M*Y)*dt - dommM
    dP <- ((s)*P*(1 - (P/((((ca+cb*R)/kmax)*((delta*(M+Y))/(1+v*(M+Y))))))) -(f*s)*P)*dt 
    dR <- ((hG)*C*(3 - R) - (hE)*(1 - C)*(R - 1))*dt - dommR
    
    #on effectue les changements
    C <- pmax(C + dC,0.0001)
    M <- pmax(M + dM,0.0001)
    P <- pmax(P + dP, 0.0001)
    R <- R + dR
    Y <- 1 - M - C
    
    #on sauve les nouvelles valeurs dans leur vecteurs
    vecC <- c(vecC,C)
    vecM <- c(vecM,M)
    vecY <- c(vecY,Y)
    vecR <- c(vecR,R)
    vecP <- c(vecP,P)
    
    #on fait avancer le temps
    t <- t + dt
  }
  #si on veut afficher le graphique
  if (graph==TRUE){
    par(mfrow=c(1,1))
    plot(seq(0,tmax,length.out = length(vecC)),vecC, type="l", col="red", ylim=c(0,1), xlab= "temps", ylab="valeurs")
    axis(side = 4,at= seq(0,1,0.1),labels = seq(1,3,0.2), )
    lines(seq(0,tmax,length.out = length(vecM)),vecM, col="blue")
    lines(seq(0,tmax,length.out = length(vecP)),vecP, col="green")
    lines(seq(0,tmax,length.out = length(vecR)),(vecR-1)/2, col="pink")
    legend("topright", legend = c("corail", "macroalgues","poissons","rugosité"),
           col = c("red","blue","green","pink"), lty = c(1, 1),
           lwd = 2, bg = "white",cex = 0.75)
  }
  #retourne la valeur finale de C et de R
  return(c(tail(vecC,1),tail(vecR,1)))
}


corailstochastiquev2(R0=2.6,f=0.4,lambda=1/15,impact=1,dt=0.1,graph = T)

graph.imageCIf <- function(R0,freq,it=10,pas=0.1){
  
  matC <- matrix(0, nrow = 21,ncol = 26)
  matsd <- matrix(0, nrow = 21,ncol = 26)
  
  R1<-R0
  l<-freq
  for (i in 1:21) {
    for (j in 1:26){
      res <- vector(length = 0)
      for (z in 1:it) {
        f0 <- i*0.05-0.05
        imp <- j*0.08-0.08
        w <- corailstochastiquev2(R0=R1,f=f0, lambda = l,impact = imp,dt=pas, graph = F)[1]
        res<-c(res,w)
      }
      matC[i,j]<-mean(res)
      matsd[i,j]<-sd(res)
    }
  }
  
  x_values <- seq(0, 2, length.out = ncol(matC))
  y_values <- seq(0, 1, length.out = nrow(matC))
  # Plot the matrix as an image with different colored pixels
  
  par(mfrow=c(1,2))
  image.plot(x = x_values, y = y_values, t(matC), col = viridis(20), 
             xlab = "impact", ylab = "Effort de pêche", 
             main=c("moyenne de C finale",paste(c("fréq:1/",1/l,", rug. ini:",R0),collapse = "")))
  image.plot(x = x_values, y = y_values, t(matsd), col = viridis(20),
             xlab = "impact", ylab = "Effort de pêche",
             main=c("écart-type de C finale",paste(c("fréq:1/",1/l,", rug. ini:",R0),collapse = "")))
}
graph.imageCIf(R0=1.2,freq = 1/15,it=2,pas = 0.5)
graph.imageCIf(R0=1.4,freq = 1/15,it=2,pas = 0.5)
graph.imageCIf(R0=1.6,freq = 1/15,it=2,pas = 0.5)
graph.imageCIf(R0=1.8,freq = 1/15,it=2,pas = 0.5)
graph.imageCIf(R0=2.0,freq = 1/15,it=2,pas = 0.5)
graph.imageCIf(R0=2.2,freq = 1/15,it=5,pas = 0.1)
graph.imageCIf(R0=2.4,freq = 1/15,it=2,pas = 0.5)
graph.imageCIf(R0=2.6,freq = 1/15,it=2,pas = 0.5)
graph.imageCIf(R0=2.8,freq = 1/15,it=2,pas = 0.5)

#graphique de R finale en fonction de I et f
graph.imageRIf <- function(R0,freq,it=10,pas=0.1){
  
  matC <- matrix(0, nrow = 21,ncol = 26)
  matsd <- matrix(0, nrow = 21,ncol = 26)
  
  R1<-R0
  l<-freq
  for (i in 1:21) {
    for (j in 1:26){
      res <- vector(length = 0)
      for (z in 1:it) {
        f0 <- i*0.05-0.05
        imp <- j*0.04-0.04
        w <- corailstochastiquev2(R0=R1,f=f0, lambda = l,impact = imp,dt=pas, graph = F)[2]
        res<-c(res,w)
      }
      matC[i,j]<-mean(res)
      matsd[i,j]<-sd(res)
    }
  }
  
  x_values <- seq(0, 1, length.out = ncol(matC))
  y_values <- seq(0, 1, length.out = nrow(matC))
  # Plot the matrix as an image with different colored pixels
  
  par(mfrow=c(1,2))
  image.plot(x = x_values, y = y_values, t(matC), col = viridis(20), 
             xlab = "impact", ylab = "Effort de pêche", 
             main=c("moyenne de R finale",paste(c("fréq:1/",1/l,", rug. ini:",R0),collapse = "")))
  image.plot(x = x_values, y = y_values, t(matsd), col = viridis(20),
             xlab = "impact", ylab = "Effort de pêche",
             main=c("écart-type de R finale",paste(c("fréq:1/",1/l,", rug. ini:",R0),collapse = "")))
}
graph.imageRIf(R0=1.8,freq = 1/15,it=5,pas = 0.5)
graph.imageRIf(R0=2.1,freq = 1/15,it=5,pas = 0.5)
graph.imageRIf(R0=2.4,freq = 1/15,it=5,pas = 0.5)

graph.imageCIR <- function(fishing,freq,it=10,pas=0.1){
  
  matC <- matrix(0, nrow = 21,ncol = 26)
  matsd <- matrix(0, nrow = 21,ncol = 26)
  
  
  f0<-fishing
  l<-freq
  for (i in 1:21) {
    for (j in 1:26){
      res <- vector(length = 0)
      for (z in 1:it) {
        R1 <- i*0.1+0.9
        imp <- j*0.08-0.08
        w <- corailstochastiquev2(R0=R1,f=f0, lambda = l,impact = imp,dt=pas, graph = F)[1]
        res<-c(res,w)
      }
      matC[i,j]<-mean(res)
      matsd[i,j]<-sd(res)
    }
  }
  
  x_values <- seq(0, 2, length.out = ncol(matC))
  y_values <- seq(1, 3, length.out = nrow(matC))
  # Plot the matrix as an image with different colored pixels
  
  par(mfrow=c(1,2))
  image.plot(x = x_values, y = y_values, t(matC), col = viridis(20), 
             xlab = "impact", ylab = "Rugosité initiale", 
             main=c("moyenne de C finale",paste(c("fréq:1/",1/l,", fishing:",f0),collapse = "")))
  image.plot(x = x_values, y = y_values, t(matsd), col = viridis(20),
             xlab = "impact", ylab = "Rugosité initiale",
             main=c("écart-type de C finale",paste(c("fréq:1/",1/l,", fishing:",f0),collapse = "")))
}
graph.imageCIR(fishing = 0.1, freq = 1/15, it=2, pas = 0.5)
graph.imageCIR(fishing = 0.4, freq = 1/15, it=2, pas = 0.5)
graph.imageCIR(fishing = 0.7, freq = 1/15, it=2, pas = 0.5)

#graphique de R finale en fonction de I et R
graph.imageRIR <- function(fishing,freq,it=10,pas=0.1){
  
  matC <- matrix(0, nrow = 21,ncol = 26)
  matsd <- matrix(0, nrow = 21,ncol = 26)
  
  
  f0<-fishing
  l<-freq
  for (i in 1:21) {
    for (j in 1:26){
      res <- vector(length = 0)
      for (z in 1:it) {
        R1 <- i*0.1+0.9
        imp <- j*0.04-0.04
        w <- corailstochastiquev2(R0=R1,f=f0, lambda = l,impact = imp,dt=pas, graph = F)[2]
        res<-c(res,w)
      }
      matC[i,j]<-mean(res)
      matsd[i,j]<-sd(res)
    }
  }
  
  x_values <- seq(0, 1, length.out = ncol(matC))
  y_values <- seq(1, 3, length.out = nrow(matC))
  # Plot the matrix as an image with different colored pixels
  
  par(mfrow=c(1,2))
  image.plot(x = x_values, y = y_values, t(matC), col = viridis(20), 
             xlab = "impact", ylab = "Rugosité initiale", 
             main=c("moyenne de R finale",paste(c("fréq:1/",1/l,", fishing:",f0),collapse = "")))
  image.plot(x = x_values, y = y_values, t(matsd), col = viridis(20),
             xlab = "impact", ylab = "Rugosité initiale",
             main=c("écart-type de R finale",paste(c("fréq:1/",1/l,", fishing:",f0),collapse = "")))
}
graph.imageRIR(fishing = 0.1, freq = 1/15, it=2, pas = 0.5)
graph.imageRIR(fishing = 0.4, freq = 1/15, it=2, pas = 0.5)
graph.imageRIR(fishing = 0.7, freq = 1/15, it=2, pas = 0.5)