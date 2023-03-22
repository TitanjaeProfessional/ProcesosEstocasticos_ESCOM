"Instituto Politecnico Nacional 
Escuela Superior de Computo
Procesos Estocasticos
Integrantes:
Angeles Lopez Nadya America
Garcia de Arcos Jose Angel Eduardo
Garcia Ochoa Christian Jaqueline
Guzman Toscano Sebastián
Padrón Hernández Andres
Ríos Rivera Jesús Rubén"


#instalacion de paquetes y importacion de librerias

install.packages("diagram")
install.packages("expm")
install.packages("markovchain")

library(diagram)
library(expm)
library(markovchain)

#Proceso de ejecucion del programa 
#----------------------------------------------------------------- Codigo 1
#Resolucion de Ejemplos 1.1-1.4 y 1.6
#specifying transition probability matrix
tm<- matrix(c(0.7, 0.1, 0.2, 0.0, 0.6, 0.4, 0.5, 0.2, 0.3),
            nrow=3, ncol=3, byrow=TRUE)
#transposing transition probability matrix
tm.tr<- t(tm)
#plotting diagram
plotmat(tm.tr, pos=c(1,2), arr.length=0.3, arr.width=0.1,
        box.col="light blue", box.lwd=1, box.prop=0.5, box.size=0.12,
        box.type="circle", cex.txt=0.8, lwd=1, self.cex=0.6,
        self.shiftx=0.17, self.shifty=-0.01)
#computing three-step transition probability matrix
tm3<- tm%^%3
print(tm3) 
#computing unconditional distribution after three steps
init.p<- c(1/3, 1/3, 1/3)
init.p%*%tm3
#creating Markov chain object
mc<- new("markovchain", transitionMatrix=tm, states=c("1","2", "3"))
#computing Markov chain characteristics
recurrentClasses(mc)
transientClasses(mc)
absorbingStates(mc)
period(mc)
round(steadyStates(mc), digits=4)
#---------------------------------------------------------------- COdigo 2
#Resolucion de ejemplos 1.5-1.7
#specifying transition probability matrix
tm<- matrix(c(0.3,0.7,0,0,0,0,1,0,0,0,0,0,0.5,0,0,0,0,0.5,
              0,0,0.6,0,0,0.4,0,0,0,0,0.1,0.9,0,0,0,0,0.7,0.3), nrow=6,
            ncol=6, byrow=TRUE)
tm.tr<- t(tm)
plotmat(tm.tr, arr.length=0.3, arr.width=0.1, box.col="lightblue", box.lwd=1, box.prop=0.5, box.size=0.09,
        box.type="circle", cex.txt=0.8, lwd=1, self.cex=0.3,
        self.arrpos=0.3, self.shiftx=0.09, self.shifty=-0.05)


mc<- new("markovchain", transitionMatrix=tm, states=c("1", "2", "3", "4", "5", "6"))
recurrentClasses(mc)
transientClasses(mc)
#finding periods of irreducible Markov chains
tm12.ir<- matrix(c(0.3,0.7,1,0), nrow=2, ncol=2, byrow=TRUE)
mc12.ir<- new("markovchain", transitionMatrix=tm12.ir,
              states=c("1","2"))
period(mc12.ir)
tm56.ir<- matrix(c(0.1,0.9,0.7,0.3), nrow=2, ncol=2,
                 byrow=TRUE)
mc56.ir<- new("markovchain", transitionMatrix=tm56.ir,
              states=c("5","6"))
period(mc56.ir)
#finding steady-state distribution
round(steadyStates(mc), digits=4)

#--------------------------------------------------------------------------------------------- Codigo 3 
#Resolucion Simulacion 1.1
#--------------------------------Simulation step 1


#specifying transition probability matrix
tm<- matrix(c(0.7, 0.1, 0.2, 0.0, 0.6, 0.4, 0.5, 0.2, 0.3),
            nrow=3, ncol=3, byrow=TRUE)
#creating Markov chain object
mc<- new("markovchain", transitionMatrix=tm, states=c("1",
                                                      "2", "3"))
#specifying total number of steps
nsteps<- 25
#specifying initial probability
p0<- c(1/3, 1/3, 1/3)
#specifying matrix containing states
MC.states<- matrix(NA, nrow=nsteps, ncol=2)
#specifying seed
set.seed(2443927)
#simulating trajectories
for (i in 1:2)
  state0<- sample(1:3, 1, prob=p0)
MC.states[,i]<- rmarkovchain(n=nsteps-1, object=mc,
                             t0=state0, include.t0=TRUE)
#plotting simulated trajectories
matplot(MC.states, type="l", lty=1, lwd=2, col=3:4,
        xaxt="n", yaxt="n", ylim=c(1,3), xlab="Step", ylab="State",
        panel.first=grid())
axis(side=1, at=c(1,5,10,15,20,25))
axis(side=2, at=1:3)
points(1:nsteps, MC.states[,1], pch=16, col=3)
points(1:nsteps, MC.states[,2], pch=16, col=4)

#--------------------------------Simulation step 2

# creating user-defined function
MC <- function(tm, p0, nsteps) {
  states <- numeric(nsteps)
  states[1] <- sample(1:3, 1, prob=p0)
  
  for(t in 2:nsteps) {
    p <- tm[states[t-1],]
    states[t] <- sample(1:3, 1, prob=p)
  }
  return(states)
}

# specifying seed
set.seed(2443927)

# simulating trajectories
MC.states2 <- matrix(NA, nrow=nsteps, ncol=2)
for (j in 1:2) {
  state0 <- sample(1:3, 1, prob=p0)
  MC.states2[,j] <- MC(tm, p0, nsteps)
}

# plotting simulated trajectories
matplot(MC.states2, type="l", lty=1, lwd=2, col=3:4,
        xaxt="n", yaxt="n", ylim=c(1,3), xlab="Step", ylab="State",
        panel.first= grid())
axis(side=1, at=c(1,5,10,15,20,25))
axis(side=2, at=c(1,2,3))
points(1:nsteps, MC.states2[,1], pch=16, col=3)
points(1:nsteps, MC.states2[,2], pch=16, col=4)

#--------------------------------Simulation step 3

#specifying matrix containing probabilities
probs<- matrix(NA, nrow=nsteps, ncol=3)
#specifying total number of steps
nsteps<- 15
#computing probabilities p_n
probs[1,] <- p0
for(n in 2:nsteps)
  probs[n,]<- probs[n-1,]%*%tm
#plotting probabilities against steps by state
matplot(probs, type="l", lty=1, lwd=2, col=2:4,
        ylim=c(0.2,0.5), xlab="Step ", ylab="Probability",
        panel.first=grid())
legend("right", c("State 1 ", "State 2", "State 3"), lty=1,
       lwd=2, col=2:4)



#--------------------------------------------------------------------------------------------- Codigo 4
# Resolucion de simulacion 1.2
#--------------------------------Simulation 1.2 step 1


tm<- matrix(c(0.3,0.7,0,0,0,0,1,0,0,0,0,0,0.5,0,0,0,0,0.5,
              0,0,0.6,0,0,0.4,0,0,0,0,0.1,0.9,0,0,0,0,0.7,0.3), nrow=6,
            ncol=6, byrow=TRUE)

mc<- new("markovchain", transitionMatrix=tm, states=c("1", "2", "3", "4", "5", "6"))

#specifying total number of steps
nsteps<- 20
#specifying initial probability
p0<- c(1/6, 1/6, 1/6, 1/6, 1/6, 1/6)
#specifying matrix containing states
MC.states<- matrix(NA, nrow=nsteps, ncol=3)
#specifying seed
set.seed(765881)
#simulating trajectories
for (i in 1:3) {
  state0<- sample(1:6, 1, prob=p0)
  MC.states[,i]<- rmarkovchain(n=nsteps-1, object=mc,
                               t0=state0,
                               include.t0=TRUE)
}
#plotting simulated trajectories
matplot(MC.states, type="l", lty=1, lwd=2, col=2:4, xaxt="n",
        ylim=c(1,6), xlab="Step", ylab="State", panel.first=grid())
axis(side=1, at=c(1,5,10,15,20))
points(1:nsteps, MC.states[,1], pch=16, col=2)
points(1:nsteps, MC.states[,2], pch=16, col=3)
points(1:nsteps, MC.states[,3], pch=16, col=4)

#---------------------------Simulation 1.2 step 2

#specifying initial state distribution (state 1)
#p0<- c(1,0,0,0,0,0)
#(state 2) 
#p0<- c(0,1,0,0,0,0)
#(state 3) 
#p0<-c(0,0,1,0,0,0)
#(state 4) 
#p0<- c(0,0,0,1,0,0)
#(state 5) 
#p0<-c(0,0,0,0,1,0)
#(state 6) 
p0<- c(0,0,0,0,0,1)
#specifying matrix containing probabilities
probs<- matrix(NA, nrow=nsteps, ncol=6)
#computing probabilities p_n
probs[1,] <- p0
for(n in 2:nsteps)
  probs[n,]<- probs[n-1,]%*%tm
#plotting probabilities vs. step by state
matplot(probs, main="Initial State 6", type="l", lty=1,
        lwd=2, col=1:6, ylim=c(-0.05, 1.1), xlab="Step",
        ylab="Probability", panel.first = grid())
legend("topright", c("State 1", "State 2", "State 3", "State
4", "State 5", "State 6"), lty=1, lwd=2, col=1:6)


#-----------------------------------------------------------------------------------------Codigo 5 
#Step 1
#specifying the transition probability matrix
tm<- matrix(c(0.1278, 0.8722, 0.6631, 0.3369), nrow=2,
            ncol=2,
            byrow=TRUE)
#creating Markov chain object

mc<- new("markovchain", transitionMatrix=tm, states=c("v",
                                                      "c"))
#computing limiting probabilities
steadyStates(mc)
#----------------------------------------------------------------------------------------------Codigo 6
#Aplication 1.3 step 1
#-------------------------
#specifying the transition probability matrix
tm<- matrix(c(1, 0, 0, 0.5, 0.5, 0, 0, 1, 0), nrow=3, ncol=3,
            byrow=TRUE)
#transposing the transition probability matrix
tm.tr<- t(tm)
plotmat(tm.tr, pos=c(1,2), name=c("AA", "Aa", "aa"),
        arr.length=0.3, arr.width=0.1, box.col="light blue",
        box.lwd=1, box.prop=0.5, box.size=0.12, box.type="circle",
        cex.txt=0.8, lwd=1, self.cex=0.6, self.shiftx=0.17,
        self.shifty=-0.01)

mc<- new("markovchain", transitionMatrix=tm,
         states=c("AA", "Aa", "aa"))
#computing stationary distribution
steadyStates(mc)

gen1<- c(0.99, 0, 0.01)
gen<- gen1%*%tm
for (n in 2:10) {
  print(n)
  print(round(gen, digits=10))
  gen<- gen%*%tm
}



