# Authors Cole Dreier, Ali Uppal, and Steven Yee
# ECON220E
# Synthetic Control Monte Carlo Simulation

rm(list = ls())
# install.packages("pracma")
# install.packages("matrixcalc")
# install.packages("limSolve")
library("pracma")
library("matrixcalc")
library("limSolve")

# ==============================
# (1) Synthetic control function
# ==============================

# This function delivers the optimal weights for the synthetic control
sc <- function(treated, control){
    # Number of units in donor pool is J (i.e., control units)
    J0 = dim(control)[2] 
    
    # Row vector of 1s (number of columns = units in donor pool)
    e = matrix(1, 1, J0)
    
    # Constraints that weights must sum to 1 (ex=f)
    f = 1
    
    # Square matrix (JxJ) with with ones along diagonal
    g = diag(x=1, J0, J0)
    
    # Weights are non-negative (gx>=h) (removing this improves fit, but by extrapolating!)
    h = matrix(0, J0, 1)
    
    # Least Squares with Equalities and Inequalities to solve weights
    weights = lsei(A=control, B=treated,E=e,F=f,G=g,H=h, type=2)$X
    return(list(weights=weights))
}

# ====================================
# (2) Randomisation inference function
# ====================================

# This function generates a p-value using randomisation inference
ri <- function(observed,T0,T1){
    # Number of total units (treated + control)
    J1 = dim(observed)[2]
    
    # Test statistics under consideration
    RMSPE = rep(NA,J1)
    T_stat = rep(NA,J1)
    Post_treatment_avg = rep(NA,J1)
    actual_synth_control = rep(NA,T1)
    actual_treated = rep(NA,T1)
    
    # RI requires looping over all units (assigning each one the treatment status)
    for (j in 1:(J1)){
        # Making jth unit the "treated"
        treated = observed[,j] 
        
        # Excludes the jth column and uses everything else as control
        control = observed[,-j] 
        
        # Caluculate the synthetic control and weights (uses sc function!)
        weights = sc(treated[1:T0],control[1:T0,])$weights
        synth_control = control%*%weights
        
        # Calculate treatment effect (treatment takes place at T0+1)
        effect = treated - synth_control
        effect_pre = effect[1:T0]
        effect_post = effect[(T0+1):T1]
        
        # Test Statistics
        RMSPE[j] = sqrt(mean(effect_post^2))/sqrt(mean(effect_pre^2))
        
        average_treat = mean(effect[(T0+1):T1])
        std_treat = std(effect[(T0+1):T1])
        T_stat[j] = abs((average_treat/(T1-T0))/(std_treat/sqrt(T1-T0)))

        Post_treatment_avg[j] = mean(abs(effect[(T0+1):T1]))
        
        if(j==1){
            actual_synth_control = synth_control
            actual_treated = treated
        }
    }
    
    # Pvalue for each test stat is proportion of placebos with a greater test stat
    # than actual treated unit which is the first entry in the matrix
    pvalue_RMPSE = mean(RMSPE>=RMSPE[1])
    pvalue_tstat = mean(T_stat >= T_stat[1])
    pvalue_post = mean(Post_treatment_avg >= Post_treatment_avg[1])
    return(list(pvalue_RMSPE=pvalue_RMPSE, pvalue_tstat=pvalue_tstat, pvalue_post=pvalue_post, actual_treated=actual_treated, actual_synth_control=actual_synth_control))
}

gen_data<- function(case, J1, T1){
    data_mat = matrix(NA,T1,J1)
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

gen_treatment <- function(case, lambda, T0, T1){
    treatment_vec = rep(0, T1)
    treatment_vec[(T0+1):T1] = lambda
    return(treatment_vec)
}

gen_treatment_varying <- function(case, lambda, T0, T1){
    treatment_vec = rep(0, T1)
    treatment_vec[(T0+1):T1] = lambda
    time_offset = linspace(-T0+1, T1-T0,T1) # Vector of time offsets from T0
    treatment_vec = hadamard.prod(treatment_vec, time_offset) # Creates a treatment vector that is time varying like Firpo
    return(treatment_vec)
}

apply_treatment <- function(case, observed, treatment_vec){
    treated_data = observed
    treated_data[,1] = treated_data[,1] + treatment_vec # This adds the treatment vector to the first column of untreated values. This makes the first column the treated unit
    return(treated_data)
}

simulate <- function(case,sims,lambda_vals,lambda_start,lambda_end,varying, T0,T1,J0,J1){
    
    lambda_seq = linspace(lambda_start, lambda_end, lambda_vals+1)
    pvalue_RMSPE_mat = matrix(NA,sims,lambda_vals+1)
    pvalue_tstat_mat = matrix(NA,sims,lambda_vals+1)
    pvalue_post_mat = matrix(NA,sims,lambda_vals+1)
    synth_sum = matrix(0, lambda_vals+1, T1)
    treated_sum = matrix(0, lambda_vals+1, T1)
    
    for (iter1 in 1:(lambda_vals+1)){
        lambda_test = lambda_seq[iter1]
        for (iter2 in 1:sims){
            data_mat = gen_data(case, J1, T1)
            if (varying){
                treatment_vec = gen_treatment_varying(case, lambda_test, T0, T1)
            }
            else{
                treatment_vec = gen_treatment(case, lambda_test, T0, T1)
            }
            observed = apply_treatment(case, data_mat, treatment_vec)
            
            results = ri(observed, T0, T1)
            pvalue_RMSPE_mat[iter2, iter1] = results$pvalue_RMSPE
            pvalue_tstat_mat[iter2, iter1] = results$pvalue_tstat
            pvalue_post_mat[iter2, iter1] = results$pvalue_post
            
            synth_sum[iter1,] = synth_sum[iter1,] + results$actual_synth_control
            treated_sum[iter1,] = treated_sum[iter1,] + results$actual_treated
        }
    }
    
    synth_avg = synth_sum / sims
    treated_avg = treated_sum / sims
    return(list(pvalue_RMSPE_mat=pvalue_RMSPE_mat, pvalue_tstat_mat=pvalue_tstat_mat,pvalue_post_mat=pvalue_post_mat, synth_avg=synth_avg, treated_avg=treated_avg))
}

gen_power_curve <- function(case, varying, lambda_start,lambda_end, pvalue_RMSPE_mat, pvalue_tstat_mat, pvalue_post_mat, size){
    lambda_vals = dim(pvalue_RMSPE_mat)[2] - 1
    lambda_seq = linspace(lambda_start, lambda_end, lambda_vals+1)
    pvalue_plot = matrix(NA,lambda_vals+1,4)
    
    for (i in 0:lambda_vals+1){
        pvalue_plot[i,1]=lambda_seq[i]
        pvalue_plot[i,2]=mean(pvalue_RMSPE_mat[,i] <= size)
        pvalue_plot[i,3]=mean(pvalue_tstat_mat[,i] <= size)
        pvalue_plot[i,4]=mean(pvalue_post_mat[,i] <= size)
    }
    
    plot(pvalue_plot[,1],pvalue_plot[,2],type="l",col="red", ylim=c(0,1))
    lines(pvalue_plot[,1],pvalue_plot[,3],type="l",col="blue")
    lines(pvalue_plot[,1],pvalue_plot[,4],type="l",col="black")
    legend("topleft", legend=c("RMSPE", "T-stat","Post-Treatment"),
           col=c("red", "blue","black"),lty=1:1, cex=0.8)
}

check_size_control <-function(case, varying,lambda_vals,lambda_start,lambda_end, pvalue_RMSPE_mat, pvalue_tstat_mat, pvalue_post_mat,size_vals){
    size_seq = linspace(0,1,size_vals+1)
    pvalue_plot = matrix(NA,size_vals+1,4)
    
    for (i in 0:size_vals+1){
        pvalue_plot[i,1]=size_seq[i]
        pvalue_plot[i,2]=mean(pvalue_RMSPE_mat[,1] <= size_seq[i])
        pvalue_plot[i,3]=mean(pvalue_tstat_mat[,1] <= size_seq[i])
        pvalue_plot[i,4]=mean(pvalue_post_mat[,1] <= size_seq[i])
    }
    
    plot(pvalue_plot[,1],pvalue_plot[,2],type="l",col="red")
    lines(pvalue_plot[,1],pvalue_plot[,3],type="l",col="blue")
    lines(pvalue_plot[,1],pvalue_plot[,4],type="l",col="black")
    legend("bottomright", legend=c("RMSPE", "T-stat","Post-Treatment"),
           col=c("red", "blue","black"),lty=1:1, cex=0.8)
}


T1  = 25
T0 = 15
J0 = 19
J1 = J0 + 1
sims = 1000
lambda_vals = 0
size_vals = J1
lambda_start = 1
lambda_end = 1
size = 0.10

case = 3
varying = FALSE


simulation = simulate(case,sims,lambda_vals,lambda_start,lambda_end, varying, T0,T1,J0,J1)

plot(linspace(1,T1,T1), simulation$treated_avg[1,], type = "l", ylim=c(-5,5))
lines(linspace(1,T1,T1), simulation$synth_avg[1,], type="l")

# check_size_control(case, varying,lambda_vals,lambda_start,lambda_end, simulation$pvalue_RMSPE_mat, simulation$pvalue_tstat_mat, simulation$pvalue_post_mat,size_vals)

# gen_power_curve(case, varying, lambda_start,lambda_end, simulation$pvalue_RMSPE_mat, simulation$pvalue_tstat_mat, simulation$pvalue_post_mat, size)
