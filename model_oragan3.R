interaction_coral <- function(t, vars, parms){
  with(as.list(c(parms, vars)), {
    # Modèle interaction coral
    
    #test <- sample(c(0.2,0.5,0.8),1,replace = T)
    Ra.C <- 0.5
    Ra.R <- 0.5
    print(t)
    
    ppois(1,lambda = 2)
    
    if(runif(1) < (1-30*dt)){
      dq <- 1
      print("oragan")
    }else{
      dq <- 0
    }
    
    Y <- 1 - M - C # dY/dt
    
    dC <- ((r/s)*Y*C - (d/s)*C - (a/s)*M*C) - C*(Ra.C + 0.05*runif(1))*dq # dC/dt    
    dM <- ((a/s)*M*C - ((1/s)*P*M/(M+Y)) + (gamma/s)*M*Y) - M*mh*dq# dM/dt
    dP <- P*(1 - (P/(((ca+cb*R)/kmax)*((delta*(M+Y))/(1+v*(M+Y)))))) - (f/s)*P # dP/dt
    dR <- hG*C*(3 - R) - hE*(1 - C)*(R - 1) - (R-1)*(Ra.R + 0.05*runif(1))*dq
    
    # Résultat
    res <- c(dC=dC,dM=dM,dP=dP,dR=dR)
    return(list(res))
  })
}
ic=c(C=0.56,M=0.001,P=0.95,R=3)
times=seq(1,100,by=1)
func=interaction_coral
#oragan.vec <- sample(c(0,1),length(times),replace = T, prob = c(0.9, 0.1))
parms=c(a=0.1,
        d=0.44,
        r=0.8,
        f=0.3,
        s=0.49,
        ca=-(3.21),
        cb=3.65,
        kmax=17.745,
        delta=4.557,
        v=0.9877,
        hG=0.03,
        hE=0.01,
        gamma = 0.8,
        mh=0.8)
soln <- ode(ic, times, func, parms)

plot(x = times, y = soln[,"M"], type = "l", col = "blue",
     xlab = "Temps", ylab = "Données",ylim=c(0,1))
lines(x = times, y = soln[,"C"], type = "l", col = "red")
lines(x = times, y = soln[,"P"], type = "l", col = "green")
