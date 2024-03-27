# projet BIO401

library(deSolve)

interaction_coral <- function(t, vars, parms){
  with(as.list(c(parms, vars)), {
    # Modèle interaction coral
    dT <- - M - C
    dM <- a*M*C - ((P/beta)*M/(M+T)) + M*T # dM/dt
    dC <- r*T*C - dC - a*M*C  # dC/dt
    dP <- s*P*(1 - (P/(beta*(((ca*R+cb)/kmax)*((gamma(M+T)/(1+v(M+T)))))))) # dP/dt
    dR <- hG*C*(3 - R) - hE(1 - C)(R - 1)
    
    # Résultat
    res <- c(dM=dM, dC=dC,dP=dP,dR=dR)
    return(list(res))
  })
}

dessinSol <- function(ic=c(R=25,C=20), times=1,func=interaction_coral, 
                      parms=c(a=0.1,
                              r=1,
                              s=0.49,
                              ca=0.75,
                              cb=1,
                              kmax=1,
                              gamma=1,
                              v=1,
                              hG=1,
                              hE=1)) {
  soln <- ode(ic, times, func, parms)
  
  plot(x = times, y = soln[,"R"], type = "l", col = "blue", 
       xlab = "Temps", ylab = "Données")
  
  # Ajouter la deuxième ligne
  lines(x = times, y = soln[,"C"], col = "red")
  
  
}