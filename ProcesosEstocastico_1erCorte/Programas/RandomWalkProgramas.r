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
install.packages("plot3D")
library(plot3D)

#Proceso de ejecucion del programa 
#----------------------------------------------------------------- Codigo 1
#specifying parameters
ntraj<- 3
p<- 0.6
nsteps<- 25
#specifying seed
set.seed(45568223)
#defining walk as matrix
walk<- matrix(NA, nrow=nsteps, ncol=ntraj)
#simulating trajectories
for (j in 1:ntraj) {
  walk[1,j]<- 0
  for (i in 2:nsteps)
    walk[i,j]<- ifelse(runif(1)<p, walk[i-1,j]+1, walk[i-1,j]-1)
}
#plotting trajectories
matplot(walk, type="l", lty=1, lwd=2, col=2:4,
        ylim=c(range(walk)), xlab="Step", ylab="Position",
        panel.first=grid())
points(1:nsteps, walk[,1], pch=16, col=2)
points(1:nsteps, walk[,2], pch=16, col=3)
points(1:nsteps, walk[,3], pch=16, col=4)

#-----------------------------------------------------------------------------Codigo 2
#specifying number of steps
nsteps<- 10000
#specifying seed
set.seed(607335)
#defining walk as matrix
walk<- matrix(NA, nrow=nsteps, ncol=2)
#setting starting point
walk[1,]<- c(0,0)
#definiting random steps
rstep<- matrix(c(1, 0, -1, 0, 0, 1, 0, -1), nrow=4, ncol=2,
               byrow=TRUE)

#simulating trajectories
for (i in 2:nsteps)
  walk[i,]<- walk[i-1,] + rstep[sample(1:4, size=1),]
#plotting trajectories
plot(x=walk[,1], y=walk[,2], type="l", col="blue",
     xlim=range(walk[,1]), ylim=range(walk[,2]), xlab="x",
     ylab="y", panel.first=grid())
#adding starting point
points(cbind(walk[1,1], walk[1,2]), pch=16, col="green",
       cex=2)
#adding ending point
points(cbind(walk[nsteps,1],walk[nsteps,2]), pch=16,
       col="red", cex=2)
#------------------------------------------------------------------------- Codigo 3
#specifying number of steps
nsteps<- 5000
#specifying seed
set.seed(830126)
#defining walk as matrix
walk<- matrix(NA, nrow=nsteps, ncol=3)
#setting starting point
walk[1,]<- c(0,0,0)
#defining random steps
rstep<- matrix(c(1,0,0,-1,0,0,0,1,0,0,-1,0,0,0,1,0,0,-1),
               nrow=6, ncol=3, byrow=TRUE)
#simulating trajectories
for (i in 2:nsteps)
  walk[i,]<- walk[i-1,]+rstep[sample(1:6, size=1),]
#plotting trajectories

lines3D(walk[,1], walk[,2], walk[,3], col="blue",
        xlim=range(walk[,1]), ylim=range(walk[,2]),
        zlim=range(walk[,3]), xlab="x", ylab="y", zlab="z", bty="b2",
        ticktype="detailed")
#adding starting point
points3D(x=walk[1,1], y=walk[1,2], z=walk[1,3], add=TRUE,
         pch=16, col="green", cex=2)
#adding ending point
points3D(walk[nsteps,1], walk[nsteps,2], walk[nsteps,3],
         add=TRUE, pch=16, col="red", cex=2)

#-------------------------------------------------------------------codigo 4
#---------------------------------step1

#specifying parameters
p<- 0.55
i<- 10
N<- 20
ntraj<- 100000
#defining walk as vector
walk<- c()
#setting counters
nNs<- 0
nzeros<- 0
ngames<- 0
#setting seed number
set.seed(30112443)
#simulating trajectories until hitting N or 0
for (j in 1:ntraj) {
  walk[1]<- i
  k<- 2
  repeat {
    walk[k]<- ifelse(runif(1)<p, walk[k-1]+1, walk[k-1]-1)
    ngames<- ngames + 1
    if (walk[k]==N) {
      nNs<- nNs+1
      break
    }
    else if(walk[k]==0) {
      nzeros<- nzeros+1
      break
    }
    k<- k+1
  }
}
print(prob.Ns<- nNs/ntraj)
print(prob.zeros <- nzeros/ntraj)
print(mean.ngames<- ngames/ntraj)

#-----------------step 2

p <- seq(0.35, 0.65, 0.001)
i <- 10
N <- 20
q <- 1 - p
p.ruin <- ifelse(p == 0.5, (N - i) / N, ((q / p)^i - (q / p)^N) / (1 - (q / p)^N))

#ploting the graphs
plot(p, p.ruin, type = "l", lwd = 2, col = "red", xlab = "p",
     ylab = "Probability", panel.first = grid())
lines(p, 1 - p.ruin, lwd = 2, col = "green")
legend("right", c("Probability of $0", "Probability of $N"),
       lty = 1, col = 2:3)

#-----------------Step 3


p<- seq(0.01,1,0.01)
i<- 10
N<- 20
q<- 1-p
E.ngames <- ifelse(p == 0.5, i * (N - i), (i - N * ((1 - (q / p)^i) / (1 - (q / p)^N))) / (1 - 2 * p))
plot(p, E.ngames, type="l", lwd=2, col="green", xlab="p",
     ylab="Expected number of games", panel.first=grid())

#----------------------------------------------------------------------------------------------Codigo 5

#specifying transition probability matrix
tm<- matrix(c(1,0,0,0,0,0,1/4,0,1/4,1/4,0,1/4,0,1/2,0,
              1/2,0,0,0,1/4,1/4,0,1/4,1/4,1/3,0,0,1/3,0,1/3,0,1/3,0,1/3,1/
                3,0), nrow=6, ncol=6, byrow=TRUE)
#setting counter
nsec<- 0
#estimating expected number of seconds
p<- matrix(NA, nrow=172, ncol=6)
p[1,]<- c(0, 0, 0, 1, 0, 0)
for (i in 2:172) {
  p[i,]<- p[i-1,]%*%tm
  nsec<- nsec+(i-1)*(p[i,1]-p[i-1,1])
}
print(nsec)


