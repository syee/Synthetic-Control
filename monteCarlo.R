
setwd("/Users/stevenyee/Documents/UCSD/UCSDEconomics/Winter2021/ECON220E/monteCarlo")
rm(list = ls())
# install.packages("pracma")
# install.packages("matrixcalc")
# install.packages("limSolve")
library("pracma")
library("matrixcalc")
library("limSolve")

T  = 25
T0 = 15
J = 19
J1 = J + 1
rho = 0.99

sims = 100
lambda.vals = 5
# lambda.sequence = c(0, 0.000000000000000005, 0.00000000005, 1)
lambda.sequence = linspace(0, 2, lambda.vals+1)
p.value.mat = matrix(NA,sims,lambda.vals+1)

sc <- function(Y1, Y0){
    e = matrix(1, 1, J)
    f = 1
    g = diag(x=1, J, J)
    h = matrix(0, J, 1)
    w.hat = lsei(A=Y0, B=Y1,E=e,F=f,G=g,H=h, type=1)$X ###WHAT IS TYPE 1 VS TYPE 2
    u.hat <- Y1-Y0%*%w.hat
    return(list(u.hat=u.hat,w.hat=w.hat))
}

ri <- function(YT,T0,T){
    #T T01 <- T0+T1
    #YT Y1Y0 <- cbind(Y1,Y0)
    #J <- dim(Y0)[2]
    S.vec = rep(NA,J)
    gaps = matrix(NA,J+1,T-T0)
    pregaps = matrix(NA,J+1,T0)
    for (j in 1:(J+1)){
        Y1.temp <- YT[,j] # Making jth unit the "treated"
        Y0.temp <- YT[,-j] # Excludes the jth column
        w.hat <- sc(Y1.temp[1:T0],Y0.temp[1:T0,])$w.hat
        u.hat <- Y1.temp-Y0.temp%*%w.hat
        u.hat.pre <- u.hat[1:T0]
        u.hat.post <- u.hat[(T0+1):T]
        # test statistic is the ration of RMSPEs
        gaps[j,] <- u.hat.post
        pregaps[j,] <- u.hat.pre
        S.vec[j] <- sqrt(mean(u.hat.post^2))/sqrt(mean(u.hat.pre^2))
    }
    p.value <- mean(S.vec>=S.vec[1])
    return(list(p.value=p.value,gaps=gaps))
}


for (iter1 in 1:(lambda.vals+1)){
    lambda.test = lambda.sequence[iter1]
    for (iter2 in 1:sims){
        # Creating Tx(J+1) iid shocks
        eta = matrix(rnorm(J+1,0,1 - rho^2), T, J1)
        
        eps = matrix(0, T, J1)
        # Initializing t=0 epsilons
        eps[1,] = eta[1,]
        
        # Epsilons are AR(1) process
        for (i in 2:T){
            eps[i,] = rho * eps[i-1,] + eta[i,]
        }
        
        # Defining time fixed effects
        F = matrix(rnorm(T,0,1), T, 1)
        
        # Defining unit fixed effects
        alpha = matrix(0, 1, J1)
        for (i in 1:J1){
            alpha[,i] = rnorm(1,i/J1,1)
        }
        
        # Generate common factors
        F.mat = repmat(F, 1, J1)
        alpha.mat = repmat(alpha, T, 1)
        
        # Generate untreated outcomes
        YN = alpha.mat + F.mat + eps
        
        # Generating treatments
        # This is where we would make the change for spillover effects
        lambda = lambda.test # Generating fixed treatment effect
        lambda = lambda * matrix(1,T-T0, 1) # Creating vector of treatment effects in post treatment
        zeros = matrix(0, T0, 1) # Creating vector of zeros for treatment effects in pre treatment
        lambda.vec = rbind(zeros, lambda) # Stacking vectors to have treatment effect for all periods
        l.vec = linspace(-T0+1, T-T0,T) # Vector of time offsets from T0
        treatment.vec = hadamard.prod(lambda.vec, l.vec) # Creates a treatment vector that is time varying like Firpo
        
        # Generates the matrix of all outcomes with the treatment included
        # This is where we would make the change if we made treatment assignment non random
        YT = YN
        YT[,1] = YT[,1] + treatment.vec # This adds the treatment vector to the first column of untreated values. This makes the first column the treated unit
        
        # Y1=YT[1:T0,1] # Treated observation
        # Y0=YT[1:T0,2:J1] # Donor pool

        p.value = ri(YT, T0, T)$p.value
        p.value.mat[iter2, iter1] = p.value
    }
    
}

#Rejection rates
size = 0.10
rejection.mat = p.value.mat <= size
colMeans(rejection.mat)

p.value.avg = colMeans(p.value.mat)
p.value.avg

