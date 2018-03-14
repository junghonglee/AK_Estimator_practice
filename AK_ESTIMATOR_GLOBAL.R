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
population <- matrix(NA,15000,192) # Sampled from true value
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

# Global AK Estimator variables#
global_ak_sim <-NULL
global_booth <- NULL

all_variance_global <- matrix(0,10,10)
all_mean_global <- matrix(0,10,10)
CV_global <- matrix(0,10,10)
MSE_global <- matrix(0,10,10)
VAR_global <- matrix(0,10,10)

ak_mean_global <- NULL
ak_var_global <- NULL

# Result variables #
global_var <- NULL
global_mean <- NULL

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
  # Global A,K Optimization #
  ###########################
  
  #initial setting time
  t <-20
  
  #* Start Timing *#
  e_time2 <- proc.time()
  
  global_ak_sim[16] <- mean(data[,16])
  delta <- NULL
  
  for (a in 1:10){
    for (k in 1:10){
      for (boo in 1:boot_simul) {
        
        # 1th step
        new_input <- sample(sample_no[1:100,t-3], 100, replace=T)
        old_input <- sample(sample_no[101:1500,t-3], 1400, replace=T)
        delta[t-3] <- mean(data[old_input,t-3]-data[old_input,t-4])
        
        global_ak_sim[t-3] <- 1/1500*(  (1-k_i[k]+a_i[a])*sum(data[new_input,t-3]) + (1-k_i[k]-(100/1400)*a_i[a])*sum(data[old_input,t-3])  )  + k_i[k]*(global_ak_sim[t-4] + delta[t-3])
        
        # 2th step
        new_input <- sample(sample_no[1:100,t-2], 100, replace=T)
        old_input <- sample(sample_no[101:1500,t-2], 1400, replace=T)
        delta[t-2] <- mean(data[old_input,t-2]-data[old_input,t-3])
        
        global_ak_sim[t-2] <- 1/1500*(  (1-k_i[k]+a_i[a])*sum(data[new_input,t-2]) + (1-k_i[k]-(100/1400)*a_i[a])*sum(data[old_input,t-2])  )  + k_i[k]*(global_ak_sim[t-3] + delta[t-2])
        
        
        # 3th step
        new_input <- sample(sample_no[1:100,t-1], 100, replace=T)
        old_input <- sample(sample_no[101:1500,t-1], 1400, replace=T)
        delta[t-1] <- mean(data[old_input,t-1]-data[old_input,t-2])
        
        global_ak_sim[t-1] <- 1/1500*(  (1-k_i[k]+a_i[a])*sum(data[new_input,t-1]) + (1-k_i[k]-(100/1400)*a_i[a])*sum(data[old_input,t-1])  )  + k_i[k]*(global_ak_sim[t-2] + delta[t-1])
        
        
        # 4th step
        new_input <- sample(sample_no[1:100,t], 100, replace=T)
        old_input <- sample(sample_no[101:1500,t], 1400, replace=T)
        delta[t] <- mean(data[old_input,t]-data[old_input,t-1])
        
        global_ak_sim[t] <- 1/1500*(  (1-k_i[k]+a_i[a])*sum(data[new_input,t]) + (1-k_i[k]-(100/1400)*a_i[a])*sum(data[old_input,t])  )  + k_i[k]*(global_ak_sim[t-1] + delta[t])
        
        global_booth[boo] <- global_ak_sim[t]
        
      }
      
      # Optimizing A and K value #
      all_variance_global[a,k] <- var(global_booth)
      all_mean_global[a,k] <- mean(global_booth)
      CV_global[a,k] <- sqrt(all_variance_global[a,k])/all_mean_global[a,k]
      #MSE_global[a,k] <- ((unemployment[[t]]/100) - all_mean_global[a,k])^2 + all_variance_global[a,k]
      #VAR_global[a,k] <- (all_variance_global[a,k])
    }
  }
  
  ak_max_global <- which(CV_global == min(CV_global), arr.ind = TRUE)     # Criterior CV
  #ak_max_global <- which(MSE_global == min(MSE_global), arr.ind = TRUE)  # Criterior MSE
  #ak_max_global <- which(VAR_global == min(VAR_global), arr.ind = TRUE)  # Criterior VAR
  
  a <- a_i[ak_max_global[1]]
  k <- k_i[ak_max_global[2]]
  
  # AK estimator value for the optimized A and K value #
  ak_mean_global[t] <- all_mean_global[ak_max_global[1],ak_max_global[2]]
  ak_var_global[t] <- all_variance_global[ak_max_global[1],ak_max_global[2]]
  
  global_var[i] <- ak_var_global[t]
  global_mean[i] <- ak_mean_global[t]
  
  
  #* end timing *#
  time_temp <- proc.time()-e_time2
  time_calc[i] <- time_temp[[3]]
  
  ##########################
  # Simple Random Sampling #
  ##########################
  
  global_boot <- sample(sample_no[,t], 1500, replace=T)
  
  SRS_mean[i] <- mean(data[global_boot])
  SRS_var[i] <- var(data[global_boot])
  
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

result_show <- function(simul, boot_simul, time_elap, SRS_var, SRS_mean, global_var, global_mean, time_calc){
  
  cat("\n")
  cat("# S I M U L A T I O N   S U M M A R Y ( G L O B A L ) # \n\n")
  cat("simulation =", simul," ");cat("bootstrap =", boot_simul," "); cat("Time =", time_elap[[3]]/60,"min \n")
  cat("avg. time per AK calculation = ",mean(time_calc), "seconds")
  cat("\n\n")
  
  cat("# Variance Performance\n")
  cat("simple / global", mean(sqrt(SRS_var)) / mean(sqrt(global_var)), "\n")
  cat("\n\n")
  
  cat("# Mean Performance \n")
  cat("Population / SRS", popul[t]/100 / mean(SRS_mean), "\n")
  cat("Population / Global", popul[t]/100 / mean(global_mean), "\n")
  cat("\n")
  
  cat("#Bootstrap Test\n")
  cat("Bootstrap test", mean(sqrt(global_var)) / sd(global_mean), "\n")
}
result_show(simul, boot_simul, time_elap, SRS_var, SRS_mean, global_var, global_mean, time_calc)


