# figue 2
library(deSolve)
library(viridis)

interaction_coral4 <- function(t, vars, parms){
  with(as.list(c(parms, vars)), {
    
    Y <- 1 - M - C # dY/dt
    R <- rug
    
    dM <- (a*M*C - P*M/(M+Y) + gamma*M*Y)/s # dM/dt
    dC <- (r*Y*C - d*C - a*M*C)/s  # dC/dt
    dP <- P*(1 - (P/(((ca+cb*R)/kmax)*((delta*(M+Y))/(1+v*(M+Y)))))) - (f)*P # dP/dt
    
    # Résultat
    #res <- c(dM=dM, dC=dC, dP=dP, dR=dR)
    res <- c(dM=dM, dC=dC, dP=dP)
    return(list(res))
  })
}

# Fonction pour effectuer un simulation du modele
dessinSol4 <- function(ic=c(M=0.95,C=0.01,P=0.75), 
                       times=seq(1:30),func=interaction_coral4, 
                       parms=c(a=0.1,
                               d=0.44,
                               r=1,
                               f=0.1,
                               s=0.49,
                               ca=-(3.21),
                               cb=3.65,
                               kmax=17.745,
                               delta=4.557,
                               v=0.9877,
                               hG=0.03,
                               hE=0.01,
                               gamma = 0.8,
                               rug=2.2)) {
  soln <- ode(ic, times, func, parms)
  
  # plot(x = times, y = soln[,"M"], type = "l", col = "blue",
  #      xlab = "Temps", ylab = "Données",ylim=c(0,1))
  # lines(x = times, y = soln[,"C"], type = "l", col = "red")
  # lines(x = times, y = soln[,"P"], type = "l", col = "green")
  
  return(soln[30,"C"])
  
}

fish.vec <- data.frame()
for (i in 1:50) {
  fish <- (i - 1) * (1 / 50)  # Calcul de la valeur de fish
  res <- dessinSol4(ic = c(M = 0.2, C = 0.8, P = 0.75),
                    parms = c(a = 0.1, d = 0.44, r = 1, f = fish,
                              s = 0.49, ca = -3.21, cb = 3.65, kmax = 17.745,
                              delta = 4.557, v = 0.9877, hG = 0.03, hE = 0.01,
                              gamma=0.8,rug=1.6))
  res2 <- dessinSol4(ic = c(M = 0.8, C = 0.2, P = 0.75),
                     parms = c(a = 0.1, d = 0.44, r = 1, f = fish,
                               s = 0.49, ca = -3.21, cb = 3.65, kmax = 17.745,
                               delta = 4.557, v = 0.9877, hG = 0.03, hE = 0.01,
                               gamma=0.8,rug=1.6))
  fish.vec[i,1] <- fish
  fish.vec[i,2] <- res
  fish.vec[i,3] <- res2
}

plot(x = fish.vec[,1], y = fish.vec[,2],type = "l",ylim = c(0,0.7),xlab = "fishing effort", ylab = "Proportion de Corail apres 30 ans", main = "Graphique")
lines(x = fish.vec[,1], y = fish.vec[,3])
