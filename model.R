# projet BIO401

library(deSolve)
library(viridis)
# test
# Partie 1
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

dessinSol <- function(ic=c(M=0.25,C=0.25,P=600,R=2), times=seq(1:100),func=interaction_coral, 
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
  # plot(x = times, y = soln[,"R"], type = "l", col = "pink",ylim=c(1,3))
  
}

# Partie 2
dessinSol()
interaction_coral2 <- function(t, vars, parms){
  with(as.list(c(parms, vars)), {
    # Modèle interaction coral
    Y <- 1 - M - C # dY/dt
    R <- rug
    
    dM <- a*M*C - ((P/beta)*M/(M+Y)) + M*Y # dM/dt
    dC <- r*Y*C - d*C - a*M*C  # dC/dt
    dP <- s*P*(1 - (P/(beta*(((ca+cb*R)/kmax)*((delta*(M+Y))/(1+v*(M+Y))))))) -f*P # dP/dt
    
    # Résultat
    res <- c(dM=dM, dC=dC, dP=dP)
    return(list(res))
  })
}

dessinSol2 <- function(ic=c(M=0.1,C=0.8,P=750), times=seq(1:30),func=interaction_coral2, 
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
                               hE=0.01,
                               rug=1)) {
  soln <- ode(ic, times, func, parms)
  
  # plot(x = times, y = soln[,"M"], type = "l", col = "blue",
  #      xlab = "Temps", ylab = "Données",ylim=c(0,1))
  # lines(x = times, y = soln[,"C"], type = "l", col = "red")
  # lines(x = times, y = soln[,"P"]/parms[2], type = "l", col = "green")
  
  return(c(soln[30,"C"],soln[30,"M"]))
}

# Initialisation de la matrice coral
coral <- matrix(data = NA, nrow = 41, ncol = 41)
algue <- matrix(data = NA, nrow = 41, ncol = 41)

# Boucle sur rugo
for (i in 1:41) {
  fish <- (i - 1) * (1 / 40)  # Calcul de la valeur de fish

  # Boucle sur fish
  for (j in 1:41) {
    rugo <- 1 + (j - 1) * (2 / 40)  # Calcul de la valeur de rugo
    
    # Appel de la fonction dessinSol2 et mise à jour de la matrice coral
    res <- dessinSol2(ic = c(M = 0.005, C = 0.995, P = 95),
                      parms = c(a = 0.1, beta = 100, d = 0.44, r = 1, f = fish,
                                s = 0.49, ca = -3.21, cb = 3.65, kmax = 17.745,
                                delta = 4.557, v = 0.9877, hG = 0.03, hE = 0.01,
                                rug = rugo))
    coral[j, i] <- res[1]  # Met à jour la matrice coral
    algue[j,i] <- res[2]
  }
}

# Définition des noms de lignes et de colonnes
row_names <- paste(seq(1, 3, by = (2 / 40)))
col_names <- paste(seq(0, 1, by = (1 / 40)))

row.names(coral) <- row_names
colnames(coral) <- col_names

# Créer le graphique avec la matrice en tant que données
image(1:nrow(coral), 1:ncol(coral), t(coral), col = plasma(100),zlim = c(0,0.56),
      xlab = "X", ylab = "Y", main = "Graphique avec couleurs basées sur les valeurs de la matrice",axes=F);axis(side = 1, at = 1:ncol(coral), labels = col_names);axis(side = 2, at = 1:nrow(coral), labels = row_names)

image(1:nrow(algue), 1:ncol(algue), t(algue), col = plasma(100),
      xlab = "X", ylab = "Y", main = "Graphique avec couleurs basées sur les valeurs de la matrice",axes=F);axis(side = 1, at = 1:ncol(coral), labels = col_names);axis(side = 2, at = 1:nrow(coral), labels = row_names)
