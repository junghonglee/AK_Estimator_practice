############################
# Created By Jung Hong Lee #
# Date : 2018-03-14        #
############################

# Initial Check up
memory <- memory.limit()

# Function Creation
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Importing Base Data
setwd("C:\\Users\\JHL\\Desktop\\thesis")  # Varies from cmputer!!!
unemployment <- read.csv("unemployment.csv", header=T)
data_length <- length(col(unemployment))

# initial setting
population <- matrix(NA,15000,data_length) # Sampled from true value
popul <- NULL                      # True value

##############################
# Number of simulation times #

simul <- 100
boot_simul <- 100

#######################################################################

# Creating Population
for (i in 1:192) {
  population [ ,i] <- rbinom(15000,1,unemployment[1,i]/100)
  popul[i] <- as.matrix(unemployment)[i]
}

# population data input
data <- population
popul <- as.matrix(unemployment)

##################
# Used Variables #
##################

# A and K values
a_i <- seq(from =0, to=1, by = 0.1)
k_i <- seq(from =0, to=1, by = 0.1)

# Timing values
time_calc <- NULL

# Local AK Estimator variables#
local_ak_sim <-NULL
local_booth <- matrix(0,boot_simul,data_length)

all_variance_local <- matrix(0,10,10)
all_mean_local <- matrix(0,10,10)
CV_local <- matrix(0,10,10)
MSE_local <- matrix(0,10,10)
VAR_local <- matrix(0,10,10)

ak_mean_local <- NULL
ak_var_local <- NULL

# Result variables #
local_var <- matrix(0,simul,data_length)
local_mean <- matrix(0,simul,data_length)

# SRS variables #
SRS_var <- NULL
SRS_mean <- NULL


#  #  #  #  #  #  #  #  #  #
#    S I M U L A T I O N   #
#  #  #  #  #  #  #  #  #  #

######################
# Time counter start #
e_time <- proc.time()

for (i in 1:simul){
  
  ##########################################
  # sampling id numbers for the population #
  ##########################################
  
  # initial setting
  no <- seq(1,15000)
  giri <- 20                              # Samping length
  sample_no <- matrix(0,1500,giri)
  
  
  
  for (t in 1:giri){
    if (t==1) {
      sample_no[1:100,t] <- sample(no,100, replace=F)
      
      for (j in 1:100){
        no <- no[no[]!=sample_no[j,t]]
      }
    }
    
    else {
      temp <- sample_no[1401:1500,t-1]
      sample_no[101:1500,t] <- sample_no[1:1400,t-1]
      sample_no[1:100,t] <- sample(no,100, replace=F)
      
      for (j in 1:100){
        no <- no[no[]!=sample_no[j,t]]
      }
      
      # adding fallout from the sample
      no <- c(no, temp)
      no <- no[no[]!=0]  # deleating 0 value
      no <- sort(no)
    }
  }
  
  #################################################
  # sampling random id numbers for the population #
  #################################################
  
  # initial setting
  no <- seq(1,15000)
  giri <- 20                                      # Samping length
  cuchul <- round(runif(giri, min=70, max=130))   # Random sampling through uniform dist
  cuchul[ cuchul < 0 ] <- 0
  
  sample_no <- matrix(0,1500,giri)
  
  for (t in 1:giri){
    if (t==1) {
      sample_no[1:cuchul[t],t] <- sample(no,cuchul[t], replace=F)
      
      for (j in 1:cuchul[t]){
        no <- no[no[]!=sample_no[j,t]]
      }
    }
    
    else {
      temp <- sample_no[(cuchul[t]+1):1500,t-1]
      sample_no[(cuchul[t]+1):1500,t] <- sample_no[1:(1500-cuchul[t]),t-1]
      sample_no[1:cuchul[t],t] <- sample(no,cuchul[t], replace=F)
      
      for (j in 1:cuchul[t]){
        no <- no[no[]!=sample_no[j,t]]
      }
      
      # adding fallout from the sample
      no <- c(no, temp)
      no <- no[no[]!=0]  # deleating 0 value
      no <- sort(no)
    }
  }
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # Getting the beset value for a and k by minimizing criterior #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  ###########################
  # Local A,K Optimization #
  ###########################
  
  # initial setting time
  t <- 20
  
  #* Start Timing *#
  e_time2 <- proc.time()
  
  #initial value
  local_ak_sim[t-4] <- mean(data[,t-4])
  delta <- NULL
  
  ############
  # 1st step #
  ############
  
  for (a in 1:10){
    for (k in 1:10){
      for (boo in 1:boot_simul) {
        
        new_input <- sample(sample_no[1:100,t-3], 100, replace=T)
        old_input <- sample(sample_no[101:1500,t-3], 1400, replace=T)
        delta[t-3] <- mean(data[old_input,t-3]-data[old_input,t-4])
        
        local_ak_sim[t-3] <- 1/1500*(  (1-k_i[k]+a_i[a])*sum(data[new_input,t-3]) + (1-k_i[k]-(100/1400)*a_i[a])*sum(data[old_input,t-3])  )  + k_i[k]*(local_ak_sim[t-4] + delta[t-3])
        local_booth[boo,t-3] <- local_ak_sim[t-3]
      }
      
      # Optimizing A and K value #
      all_variance_local[a,k] <- var(local_booth[,t-3])
      all_mean_local[a,k] <- mean(local_booth[,t-3])
      
      CV_local[a,k] <- sqrt(all_variance_local[a,k])/all_mean_local[a,k]
      #MSE_local[a,k] <- ((unemployment[[t-3]]/100) - all_mean_local[a,k])^2 + all_variance_local[a,k]
      #VAR_local[a,k] <- (all_variance_local[a,k])
      
      
    }
  }
  
  ak_max_local <- which(CV_local == min(CV_local), arr.ind = TRUE)
  #ak_max_local <- which(MSE_local == min(MSE_local), arr.ind = TRUE)
  #ak_max_local <- which(VAR_local == min(VAR_local), arr.ind = TRUE)
  
  # AK estimator value for the optimized A and K value #
  ak_mean_local[t-3] <- all_mean_local[ak_max_local[1],ak_max_local[2]]
  ak_var_local[t-3] <- all_variance_local[ak_max_local[1],ak_max_local[2]]
  
  ############
  # 2nd step #
  ############
  
  for (a in 1:10){
    for (k in 1:10){
      for (boo in 1:boot_simul) {
        
        new_input <- sample(sample_no[1:100,t-2], 100, replace=T)
        old_input <- sample(sample_no[101:1500,t-2], 1400, replace=T)
        delta[t-2] <- mean(data[old_input,t-2]-data[old_input,t-3])
        
        local_ak_sim[t-2] <- 1/1500*(  (1-k_i[k]+a_i[a])*sum(data[new_input,t-2]) + (1-k_i[k]-(100/1400)*a_i[a])*sum(data[old_input,t-2])  )  + k_i[k]*(ak_mean_local[t-3] + delta[t-2])
        local_booth[boo,t-2] <- local_ak_sim[t-2]
      }
      
      # Optimizing A and K value #
      all_variance_local[a,k] <- var(local_booth[,t-2])
      all_mean_local[a,k] <- mean(local_booth[,t-2])
      
      CV_local[a,k] <- sqrt(all_variance_local[a,k])/all_mean_local[a,k]
      #MSE_local[a,k] <- ((unemployment[[t-2]]/100) - all_mean_local[a,k])^2 + all_variance_local[a,k]
      #VAR_local[a,k] <- (all_variance_local[a,k])
      
      
    }
  }
  
  ak_max_local <- which(CV_local == min(CV_local), arr.ind = TRUE)
  #ak_max_local <- which(MSE_local == min(MSE_local), arr.ind = TRUE)
  #ak_max_local <- which(VAR_local == min(VAR_local), arr.ind = TRUE)
  
  # AK estimator value for the optimized A and K value #
  ak_mean_local[t-2] <- all_mean_local[ak_max_local[1],ak_max_local[2]]
  ak_var_local[t-2] <- all_variance_local[ak_max_local[1],ak_max_local[2]]
  
  
  ###########
  # 3 Step  #
  ###########
  
  for (a in 1:10){
    for (k in 1:10){
      for (boo in 1:boot_simul) {
        
        new_input <- sample(sample_no[1:100,t-1], 100, replace=T)
        old_input <- sample(sample_no[101:1500,t-1], 1400, replace=T)
        delta[t-1] <- mean(data[old_input,t]-data[old_input,t-2])
        
        local_ak_sim[t-1] <- 1/1500*(  (1-k_i[k]+a_i[a])*sum(data[new_input,t-1]) + (1-k_i[k]-(100/1400)*a_i[a])*sum(data[old_input,t-1])  )  + k_i[k]*(ak_mean_local[t-2] + delta[t-1])
        local_booth[boo,t-1] <- local_ak_sim[t-1]
        
      }
      
      # Optimizing A and K value #
      all_variance_local[a,k] <- var(local_booth[,t-1])
      all_mean_local[a,k] <- mean(local_booth[,t-1])
      
      CV_local[a,k] <- sqrt(all_variance_local[a,k])/all_mean_local[a,k]
      #MSE_local[a,k] <- ((unemployment[[t-1]]/100) - all_mean_local[a,k])^2 + all_variance_local[a,k]
      #VAR_local[a,k] <- (all_variance_local[a,k])
      
      
    }
  }
  
  ak_max_local <- which(CV_local == min(CV_local), arr.ind = TRUE)
  #ak_max_local <- which(MSE_local == min(MSE_local), arr.ind = TRUE)
  #ak_max_local <- which(VAR_local == min(VAR_local), arr.ind = TRUE)
  
  # AK estimator value for the optimized A and K value #
  ak_mean_local[t-1] <- all_mean_local[ak_max_local[1],ak_max_local[2]]
  ak_var_local[t-1] <- all_variance_local[ak_max_local[1],ak_max_local[2]]
  
  
  ###########
  # 4 Step  #
  ###########
  
  for (a in 1:10){
    for (k in 1:10){
      for (boo in 1:boot_simul) {
        
        new_input <- sample(sample_no[1:100,t], 100, replace=T)
        old_input <- sample(sample_no[101:1500,t], 1400, replace=T)
        delta[t] <- mean(data[old_input,t]-data[old_input,t-1])
        
        local_ak_sim[t] <- 1/1500*(  (1-k_i[k]+a_i[a])*sum(data[new_input,t]) + (1-k_i[k]-(100/1400)*a_i[a])*sum(data[old_input,t])  )  + k_i[k]*(ak_mean_local[t-1] + delta[t])
        local_booth[boo,t] <- local_ak_sim[t]
        
      }
      
      # Optimizing A and K value #
      all_variance_local[a,k] <- var(local_booth[,t])
      all_mean_local[a,k] <- mean(local_booth[,t])
      
      CV_local[a,k] <- sqrt(all_variance_local[a,k])/all_mean_local[a,k]
      #MSE_local[a,k] <- ((unemployment[[t]]/100) - all_mean_local[a,k])^2 + all_variance_local[a,k]
      #VAR_local[a,k] <- (all_variance_local[a,k])
      
      
    }
  }
  
  ak_max_local <- which(CV_local == min(CV_local), arr.ind = TRUE)
  #ak_max_local <- which(MSE_local == min(MSE_local), arr.ind = TRUE)
  #ak_max_local <- which(VAR_local == min(VAR_local), arr.ind = TRUE)
  
  # AK estimator value for the optimized A and K value #
  ak_mean_local[t] <- all_mean_local[ak_max_local[1],ak_max_local[2]]
  ak_var_local[t] <- all_variance_local[ak_max_local[1],ak_max_local[2]]
  
  # Saving Data
  local_var[i,t] <- ak_var_local[t]
  local_var[i,t-1] <- ak_var_local[t-1]
  local_var[i,t-2] <- ak_var_local[t-2]
  local_var[i,t-3] <- ak_var_local[t-3]
  local_mean[i,t] <- ak_mean_local[t]
  local_mean[i,t-1] <- ak_mean_local[t-1]
  local_mean[i,t-2] <- ak_mean_local[t-2]
  local_mean[i,t-3] <- ak_mean_local[t-3]
  
  #* end timing *#
  time_temp <- proc.time()-e_time2
  time_calc[i] <- time_temp[[3]]
  
  ##########################
  # Simple Random Sampling #
  ##########################
  
  local_boot <- sample(sample_no[,t], 1500, replace=T)
  
  SRS_mean[i] <- mean(data[local_boot])
  SRS_var[i] <- var(data[local_boot])
  
  #####################
  # Checking Progress #
  #####################
  
  cat( "[overall =", i/simul*100 ,"%] \t")
  cat( "[elapsed time =", (proc.time()-e_time)[[3]]/60, "minutes] \n")
}

time_elap <- proc.time()-e_time


# # # # # # # #
# R E S U L T #
# # # # # # # #

result_show <- function(simul, boot_simul, time_elap, SRS_var, SRS_mean, local_var, local_mean, time_calc){
  
  cat("\n")
  cat("# S I M U L A T I O N   S U M M A R Y ( L O C A L ) # \n\n")
  cat("simulation =", simul," ");cat("bootstrap =", boot_simul," "); cat("Time =", time_elap[[3]]/60,"min", "\n")
  cat("avg. time per AK calculation = ",mean(time_calc), "seconds")
  cat("\n\n")
  
  cat("# Variance Performance\n")
  cat("simple / ak2 =",mean(sqrt(SRS_var)) / mean(sqrt(local_var[,t-3])), "\n")
  cat("simple / ak3 =",mean(sqrt(SRS_var)) / mean(sqrt(local_var[,t-2])), "\n")
  cat("simple / ak4 =",mean(sqrt(SRS_var)) / mean(sqrt(local_var[,t-1])), "\n")
  cat("simple / ak5 =",mean(sqrt(SRS_var)) / mean(sqrt(local_var[,t])), "\n")
  cat("\n\n")
  
  cat("# Mean Performance \n")
  cat("Population / SRS =", popul[t]/100 / mean(SRS_mean), "\n")
  cat("Population / ak2 =", popul[t]/100 / mean(local_mean[,t-3]), "\n")
  cat("Population / ak3 =", popul[t]/100 / mean(local_mean[,t-2]), "\n")
  cat("Population / ak4 =", popul[t]/100 / mean(local_mean[,t-1]), "\n")
  cat("Population / ak5 =", popul[t]/100 / mean(local_mean[,t]), "\n")
  cat("\n\n")
  
  cat("#Bootstrap Test\n")
  cat("Bootstrap test ak2 =", mean(sqrt(local_var[,t-3])) / sd(local_mean[,t-3]), "\n")
  cat("Bootstrap test ak3 =", mean(sqrt(local_var[,t-3])) / sd(local_mean[,t-2]), "\n")
  cat("Bootstrap test ak4 =", mean(sqrt(local_var[,t-3])) / sd(local_mean[,t-1]), "\n")
  cat("Bootstrap test ak5 =", mean(sqrt(local_var[,t-3])) / sd(local_mean[,t]), "\n")
}
result_show(simul, boot_simul, time_elap, SRS_var, SRS_mean, local_var, local_mean, time_calc)