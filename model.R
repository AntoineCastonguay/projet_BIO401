# projet BIO401

library(deSolve)

interaction_coral <- function(t, vars, parms){
  with(as.list(c(parms, vars)), {
    # Modèle interaction coral
    Y <- 1 - M - C # dY/dt
    dM <- a*M*C - ((P/beta)*M/(M+Y)) + M*Y # dM/dt
    dC <- r*Y*C - d*C - a*M*C  # dC/dt
    dP <- s*P*(1 - (P/(beta*(((ca+cb*R)/kmax)*((delta*(M+Y))/(1+v*(M+Y))))))) -f*P # dP/dt
    dR <- hG*C*(3 - R) - hE*(1 - C)*(R - 1)
    print(paste("C",C))
    print(paste("R",R))
    print(paste("P",P))
    print(paste("k", (((ca+cb*R)/kmax)*((delta*(M+Y))/(1+v*(M+Y))))))
    
    # Résultat
    res <- c(dM=dM, dC=dC, dP=dP, dR=dR)
    return(list(res))
  })
}

dessinSol <- function(ic=c(M=0.5,C=0.5,P=1000,R=2), times=seq(1:10),func=interaction_coral, 
                      parms=c(a=0.1,
                              beta=10000,
                              d=0.44,
                              r=1,
                              f=0.1,
                              s=0.49,
                              ca=-(3.21),
                              cb=3.65,
                              kmax=17.745,
                              delta=4,557,
                              v=0.9877,
                              hG=0.03,
                              hE=0.01)) {
  soln <- ode(ic, times, func, parms)
  
  plot(x = times, y = soln[,"M"], type = "l", col = "blue",
       xlab = "Temps", ylab = "Données",ylim=c(0,1))
  lines(x = times, y = soln[,"C"], type = "l", col = "red")
  lines(x = times, y = soln[,"P"]/10000, type = "l", col = "green")
}

dessinSol()
