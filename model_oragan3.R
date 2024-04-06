interaction_coral <- function(t, vars, parms){
  with(as.list(c(parms, vars)), {
    # Modèle interaction coral
    
    if(runif(1)< 0.2){
      dq <- 1
    }else{
      dq <- 0
    }
    
    Y <- 1 - M - C # dY/dt
    
    dM <- ((a/s)*M*C - ((1/s)*P*M/(M+Y)) + (gamma/s)*M*Y) - M*mh*dq# dM/dt
    dC <- ((r/s)*Y*C - (d/s)*C - (a/s)*M*C) - C*(0.5)*dq # dC/dt
    dP <- P*(1 - (P/(((ca+cb*R)/kmax)*((delta*(M+Y))/(1+v*(M+Y)))))) - (f/s)*P # dP/dt
    dR <- hG*C*(3 - R) - hE*(1 - C)*(R - 1) - R*(0.5)*dq
    
    if(M < 0) dM <- 0
    if(M > 1) dM <- 1
    if(C < 0) dC <- 0
    if(C > 1) dC <- 1
    
    # Résultat
    res <- c(dM=dM, dC=dC, dP=dP, dR=dR)
    return(list(res))
  })
}

ic=c(M=0.001,C=0.75,P=0.95,R=3)
times=seq(1:100)
func=interaction_coral
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
