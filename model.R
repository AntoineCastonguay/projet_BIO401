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

f_values <- seq(0, 1, by = 0.05)
R0 <- 1

bifurcation_data <- lapply(f_values, function(f) {
  parms=c(a=0.1,
          beta=1,
          d=0.44,
          r=1,
          f=f,
          s=0.49,
          ca=-(3.21),
          cb=3.65,
          kmax=17.745,
          delta=4,557,
          v=0.9877,
          hG=0.03,
          hE=0.01)
  
  C.lim <- unname((parms["r"] - parms["d"])/parms["r"])
  M.lim <- 0
  P.lim <- unname((1-parms["f"])*(((parms["ca"]+parms["cb"]*1)/parms["kmax"])*((parms["delta"]*(1-C.lim))/(1+parms["v"]*(1-C.lim)))))
  ic=c(M=0,C=C.lim,P=P.lim,R=1)
  times = seq(1,100)
  
  out <- as.data.frame(ode(ic, times, interaction_coral, parms))
  out$R
})

plot(rep(f_values, each = length(seq(1,100))), unlist(bifurcation_data),
     xlab = "f", ylab = "R", main = "Diagramme de bifurcation - Équation logistique",
     col = "blue", pch = 20, cex = 0.5)

dessinSol <- function(ic=c(M=0.3,C=0.3,P=600,R=2.5), times=seq(1:100),func=interaction_coral, 
                      parms=c(a=0.1,
                              beta=1000,
                              d=0.44,
                              r=1,
                              f=0.2,
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
  lines(x = times, y = soln[,"P"]/parms[2], type = "l", col = "green")
  plot(x = times, y = soln[,"R"], type = "l", col = "pink",ylim=c(1,3))
  
}

dessinSol()
