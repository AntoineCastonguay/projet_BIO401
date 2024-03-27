# projet BIO401

library(deSolve)

interaction_coral <- function(t, vars, parms){
  with(as.list(c(parms, vars)), {
    # Modèle interaction coral
    Y <- 1 - M - C # dY/dt
    dM <- a*M*C - ((P/beta)*M/(M+Y)) + gamma*M*Y # dM/dt
    dC <- r*Y*C - d*C - a*M*C  # dC/dt
    dP <- s*P*(1 - (P/(beta*kmax))) - f*P # dP/dt
    
    # Résultat
    res <- c(dM=dM, dC=dC, dP=dP)
    return(list(res))
  })
}

dessinSol <- function(ic=c(M=0.01,C=0.8,P=1000), times=seq(1:100),func=interaction_coral, 
                      parms=c(a=0.1,
                              beta=1000,
                              d=0.44,
                              r=1,
                              f=0.8,
                              s=0.49,
                              ca=-(3.21),
                              cb=3.65,
                              kmax=0.5,
                              gamma = 0.8,
                              delta=4,557,
                              v=0.9877,
                              hG=0.03,
                              hE=0.01)) {
  soln <- ode(ic, times, func, parms)
  
  plot(x = times, y = soln[,"M"], type = "l", col = "blue",
        xlab = "Temps", ylab = "Données")
  lines(x = times, y = soln[,"C"], type = "l", col = "red")
}

dessinSol()
