library(hexbin)
######### Functions ##########
#Complex normal
Cnorm <- function(n) complex(real = rnorm(n),imaginary = rnorm(n))/sqrt(2)

#NxN GOE
GOE <- function(N){
  G <- matrix(rnorm(N*N),nrow = N)
  return((G+t(G))/sqrt(2))
}

#Simulation of M samples of NxN IID matrices with specified distribution
simIID <- function(N,M,dist = Cnorm,returnL = FALSE){
  L <- c()
  for (i in 1:M) {
    G <- matrix(dist(N*N),nrow = N)
    a <- eigen(G,only.values = TRUE)
    L <- append(L, a$values/sqrt(N))
  }
  plot(hexbin(L),
       xlab ="real",
       ylab ="imaginary",
       colramp=BTC)
  if(returnL){return(L)}
}

#Simulate M samples of the NxN GOE
simGOE <- function(N,M){
  L <- c()
  for (i in 1:M) {
    H <- GOE(N)
    a <- eigen(H,only.values = TRUE)
    L <- append(L, a$values/sqrt(N))
  }
  hist(L,prob=TRUE,
       main=sprintf("Eigenvalues from %d samples of GOE(%d)",M,N), 
       xlab = "Eigenvalues",
       border="blue", 
       col="grey",
       breaks=70)
}

#Create a Wigner matrix
Wigner <- function(N,dist=Cnorm){
  X <- matrix(dist(N*N),nrow = N)
  return(X+Conj(t(X)))
}

#Simulate M samples of NxN Wigner with specified distribution
simWigner <- function(N,M,dist=rnorm){
  L <- c()
  for (i in 1:M) {
    H <- wigner(N,dist)
    a <- eigen(H,only.values = TRUE)
    L <- append(L, a$values/sqrt(N))
  }
  hist(L,prob=TRUE,
       main=sprintf("Eigenvalues of %d samples of specified wigner matrices",M), 
       border="blue", 
       col="pink",
       breaks=100)
}

########## Random plot calls
#Zoom on real eigenvalues of real Ginibre
L1 <- simIID(500,100,rnorm,returnL = TRUE)
hist(Re(L1[Im(L1)==0]),breaks = 30,prob=T) 

#Added mean to a Ginibre and Wigner
L1 <- simIID(100,500,function(n){complex(real=rnorm(n,sqrt(2)/10,1),
                                         imaginary = rnorm(n,sqrt(2)/10,1))/sqrt(2)},returnL = T)
hexbinplot(Im(L1)~Re(L1),colramp=BTC)
simWigner(50,500,function(n){rexp(n)})
simWigner(100,500,function(n){rbinom(n,1,0.5)})

simGOE(100,100)
simGOE(1000,1)
simIID(100,100)
simIID(1000,1)

L<- simIID(10000,1,returnL = T) #This command takes time to execute (~25min)
hexbinplot(Im(L)~Re(L),colramp=BTC)
simIID(100,500,function(n){2*rbinom(n,1,0.5)-1})
simWigner(100,200,function(n){2*rbinom(n,1,0.5)-1})
simWigner(500,1,rcauchy)
