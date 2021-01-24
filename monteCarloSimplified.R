# Authors Ali Uppal and Steven Yee
# ECON220E
# Synthetic Control Monte Carlo Simulation

rm(list = ls())
# install.packages("pracma")
# install.packages("matrixcalc")
# install.packages("limSolve")
library("pracma")
library("matrixcalc")
library("limSolve")

sc <- function(treated, control){
    J <- dim(control)[2] # Number of units in donor pool (Y0)
    e = matrix(1, 1, J)
    f = 1
    g = diag(x=1, J, J)
    h = matrix(0, J, 1)
    weights = lsei(A=control, B=treated,E=e,F=f,G=g,H=h, type=2)$X
    return(list(weights=weights))
}

ri <- function(observed,T0,T1){
    J1 <- dim(observed)[2]
    
    #Add other test statistics
    RMSPE = rep(NA,J1)

    for (j in 1:(J1)){
        treated <- observed[,j] # Making jth unit the "treated"
        control <- observed[,-j] # Excludes the jth column and uses everything else as control
        
        # Caluculate synthetic control and weights
        weights <- sc(treated[1:T0],control[1:T0,])$weights
        synth_control = control%*%weights
        
        # Calculate treatment effect
        effect <- treated - synth_control
        effect_pre <- effect[1:T0]
        effect_post <- effect[(T0+1):T1]
        # test statistic is the ration of RMSPEs
        RMSPE[j] <- sqrt(mean(effect_post^2))/sqrt(mean(effect_pre^2))
    }
    
    pvalue_RMPSE <- mean(RMSPE>=RMSPE[1])
    return(list(pvalue_RMSPE=pvalue_RMPSE))
}

gen_data<- function(case, J1, T1){
    data_mat = matrix(NA,T1,J1)
    # if(case != 1){
    #     # Creating Tx(J+1) iid shocks
    #     eta = matrix(rnorm(J+1,0,1 - rho^2), T1, J1)
    #     
    #     eps = matrix(0, T1, J1)
    #     # Initializing t=0 epsilons
    #     eps[1,] = eta[1,]
    #     
    #     # Epsilons are AR(1) process
    #     for (i in 2:T1){
    #         eps[i,] = rho * eps[i-1,] + eta[i,]
    #     }
    #     
    #     # Defining time fixed effects
    #     F = matrix(rnorm(T1,0,1), T1, 1)
    #     
    #     # Defining unit fixed effects
    #     alpha = matrix(0, 1, J1)
    #     for (i in 1:J1){
    #         alpha[,i] = rnorm(1,i/J1,1)
    #     }
    #     
    #     # Generate common factors
    #     F.mat = repmat(F, 1, J1)
    #     alpha.mat = repmat(alpha, T1, 1)
    #     
    #     # Generate untreated outcomes
    #     data_mat = alpha.mat + F.mat + eps
    #     
    # }
    if (case == 4){
        data_mat = matrix(rexp(J1*T1),T1,J1)
    }
    else if(case == 3){
        # Creating T1xJ1 iid shocks
        eps = matrix(rnorm(J1*T1,0,1), T1, J1)
        
        # Defining time fixed effects
        time_fe = matrix(rnorm(T1,0,1), T1, 1)
        
        # Defining unit fixed effects
        unit_fe = matrix(NA, 1, J1)
        for (i in 2:J1){
            unit_fe[,i] = rnorm(1,0,1)
        }
        unit_fe[,1] = -1.1 * abs(min(unit_fe[2:J1])) # We are setting the treated unit to be at the bottom of the distribution. We are taking the minimum of the control units and making sure teh treated unit_fe is 10% lower than the minimum
        
        # Generate common factors
        time_fe_mat = repmat(time_fe, 1, J1)
        unit_fe = repmat(unit_fe, T1, 1)
        
        # Generate untreated outcomes
        data_mat = unit_fe + time_fe_mat + eps
        
    }
    else if(case == 2){
        # Creating T1xJ1 iid shocks
        eps = matrix(rnorm(J1*T1,0,1), T1, J1)
        
        # Defining time fixed effects
        time_fe = matrix(rnorm(T1,0,1), T1, 1)
        
        # Defining unit fixed effects
        unit_fe = matrix(NA, 1, J1)
        for (i in 2:J1){
            unit_fe[,i] = rnorm(1,0,1)
        }
        unit_fe[,1] = mean(unit_fe[2:J1]) # We are setting the treated unit to be in the middle of the distribution
        
        # Generate common factors
        time_fe_mat = repmat(time_fe, 1, J1)
        unit_fe = repmat(unit_fe, T1, 1)
        
        # Generate untreated outcomes
        data_mat = unit_fe + time_fe_mat + eps
        
    }
    else if (case == 1){
        data_mat = matrix(rnorm(J1*T1),T1,J1)
    }
    return(data_mat)
}

# gen_treatment <- function(case, lambda, T0, T1){
#     lambda_vec = rep(0, T1)
#     lambda_vec[(T0+1):T1] = lambda
#     treatment_vec = rep(0, T1)
#     
#     if(case == 3){
#         treatment_vec = lambda_vec
#         # l.vec = linspace(-T0+1, T1-T0,T1) # Vector of time offsets from T0
#         # treatment.vec = hadamard.prod(lambda.vec, l.vec) # Creates a treatment vector that is time varying like Firpo
#     }
#     else if (case == 2){
#         treatment_vec = lambda_vec
#     }
#     else if (case == 1){
#         treatment_vec = lambda_vec
#     }
#     return(treatment_vec)
# }

gen_treatment <- function(case, lambda, T0, T1){
    treatment_vec = rep(0, T1)
    treatment_vec[(T0+1):T1] = lambda
    return(treatment_vec)
}

apply_treatment <- function(case, observed, treatment_vec){
    treated_data = observed
    treated_data[,1] = treated_data[,1] + treatment_vec # This adds the treatment vector to the first column of untreated values. This makes the first column the treated unit
    return(treated_data)
}

# apply_treatment <- function(case, observed, treatment_vec){
#     treated_data = observed
#     if(case == 3){
#         treated_data[,1] = treated_data[,1] + treatment_vec # This adds the treatment vector to the first column of untreated values. This makes the first column the treated unit
#         # l.vec = linspace(-T0+1, T1-T0,T1) # Vector of time offsets from T0
#         # treatment.vec = hadamard.prod(lambda.vec, l.vec) # Creates a treatment vector that is time varying like Firpo
#     }
#     else if (case == 2){
#         treated_data[,1] = treated_data[,1] + treatment_vec # This adds the treatment vector to the first column of untreated values. This makes the first column the treated unit
#     }
#     else if (case == 1) || (case == 2){
#         treated_data[,1] = treated_data[,1] + treatment_vec # This adds the treatment vector to the first column of untreated values. This makes the first column the treated unit
#     }
#     return(treated_data)
# }


T1  = 25
T0 = 15
J = 19
J1 = J + 1
# rho = 0.99

# Simulation Parameters
sims = 1000
lambda_vals = 0
lambda_seq = linspace(0, 0, lambda_vals+1)

pvalue_RMSPE_mat = matrix(NA,sims,lambda_vals+1)

case = 3

for (iter1 in 1:(lambda_vals+1)){
    lambda_test = lambda_seq[iter1]
    for (iter2 in 1:sims){
        
        data_mat = gen_data(case, J1, T1)

        treatment_vec = gen_treatment(case, lambda_test, T0, T1)

        observed = apply_treatment(case, data_mat, treatment_vec)

        pvalue_RMSPE = ri(observed, T0, T1)$pvalue_RMSPE
        pvalue_RMSPE_mat[iter2, iter1] = pvalue_RMSPE
    }
    
}

# Plot p-values conditional on size
space = 10
size = linspace(0,1,space+1)
pvalue_plot = matrix(NA,space+1,2)

for (i in 0:space+1){
    pvalue_plot[i,1]=size[i]
    pvalue_plot[i,2]=colMeans(pvalue_RMSPE_mat <= size[i])
}

plot(pvalue_plot[,1],pvalue_plot[,2])

# #Rejection rates
# size = 0.10
# rejection.mat = pvalue_RMSPE_mat <= size
# colMeans(rejection.mat)
# # colMeans(pvalue_RMSPE_mat <= 0.1)
# 
# colMeans(pvalue_RMSPE_mat <= 0.1)
# 
# p.value.avg = colMeans(pvalue_RMSPE_mat)
# p.value.avg

