# projet BIO401

library(deSolve)

interaction_coral <- function(t, vars, parms){
  with(as.list(c(parms, vars)), {
    # Modèle interaction coral
    dY <- Y - M - C # dY/dt
    dM <- a*M*C - ((P/beta)*M/(M+Y)) + M*Y # dM/dt
    dC <- r*Y*C - d*C - a*M*C  # dC/dt
    dP <- s*P*(1 - (P/(beta*(((ca*R+cb)/kmax)*((delta*(M+Y)/(1+v*(M+Y)))))))) # dP/dt
    dR <- hG*C*(3 - R) - hE*(1 - C)*(R - 1)
    
    # Résultat
    res <- c(dY=dY, dM=dM, dC=dC, dP=dP, dR=dR)
    return(list(res))
  })
}

dessinSol <- function(ic=c(Y=1000,M=500,C=500,P=1000,R=1000), times=seq(1:1000),func=interaction_coral, 
                      parms=c(a=0.1,
                              beta=50,
                              d=0.44,
                              r=1,
                              s=0.49,
                              ca=-(3.21),
                              cb=3.65,
                              kmax=10000,
                              delta=4,557,
                              v=0.9877,
                              hG=0.03,
                              hE=0.01)) {
  soln <- ode(ic, times, func, parms)
  
  #plot(x = times, y = soln[,"M"], type = "l", col = "blue",
  # xlab = "Temps", ylab = "Données")
}

dessinSol()
