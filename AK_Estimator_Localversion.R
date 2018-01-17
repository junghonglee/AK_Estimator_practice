# Introduction Added 
# AK estimator using local and global 

#simulation 
memory.limit()

#use 4 cores
#library(doSNOW)
#cl <- makeCluster(3, type="SOCK") # for 4 cores machine
#registerDoSNOW (cl)

#Function declaration
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#####################
#creating population#
#####################

# from 2000 Jan ~ 2015 Dec from 통계청
setwd("C:\\Users\\JHL\\Desktop\\thesis")
unemployment <- read.csv("unemployment.csv", header=T)
population <- matrix(NA,15000,192) # Sampled from true value
popul <- NULL                      # True value


#Creating Population data
for (i in 1:192) {
  population [ ,i] <- rbinom(15000,1,unemployment[1,i]/100)
  popul[i] <- as.matrix(unemployment)[i]
}

#population data input
data <- population
popul <- as.matrix(unemployment)

# Number of simulation times
simul <- 100
boot_simul <- 100

#Result Data
result <- matrix(0,simul,7)
colnames(result) <- (c("mean-ak2","mean-ak3","mean-ak4","popul-mean","popul-ak2","popul-ak3","popul-ak4"))


# Initial variables
ak_est <- matrix(0,simul,2)
trial_ak_estimator2 <- NULL
trial_ak_estimator3 <- NULL
trial_ak_estimator4 <- NULL

simple_var1 <- NULL
simple_var2 <- NULL
simple_var3 <- NULL
simple_mean1 <- NULL
simple_mean2 <- NULL
simple_mean3 <- NULL

boo_sum1 <- NULL
boo_sum2 <- NULL
boo_sum3 <- NULL
boo_mean1 <- 0
boo_mean2 <- 0
boo_mean3 <- 0
boo_var1 <- NULL
boo_var2 <- NULL
boo_var3 <- NULL

ak_sim <- NULL
ak_sim2 <-NULL
ak_sim3 <-NULL

global_var <- NULL
global_mean <- NULL 


e_time <- proc.time()
###########################################################
###########################################################

for (i in 1:simul){
  
  
  #########################################
  # sampling id numbers for the population#
  #########################################
  
  # initial setting
  no <- seq(1,15000)
  giri <- 19
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
  
  #######################################################
  # Getting the beset value for a and k for minimum mse #
  #######################################################
  
  all_variance <- matrix(0,10,10)
  all_variance2 <- matrix(0,10,10)
  all_variance3 <- matrix(0,10,10)
  all_mean <- matrix(0,10,10)
  all_mean2 <- matrix(0,10,10)
  all_mean3 <- matrix(0,10,10)
  all_mse <- matrix(0,10,10)
  all_mse2 <- matrix(0,10,10)
  all_mse3 <- matrix(0,10,10)
  
  a_i <- seq(from =0, to=1, by = 0.1)
  k_i <- seq(from =0, to=1, by = 0.1)
  
  ak_mean <- NULL
  ak_mean2 <- NULL
  ak_mean3 <- NULL
  
  ak_var <- NULL
  ak_var2 <- NULL
  ak_var3 <- NULL
  
  
  # Global AK Estimator variables#
  global_ak_sim <-NULL
  global_booth <- NULL
  all_variance_global <- matrix(0,10,10)
  all_mean_global <- matrix(0,10,10)
  all_mse_global <- matrix(0,10,10)
  ak_mean_global <- NULL
  ak_var_global <- NULL
  
  # SRS variables #
  sample_mean <- NULL
  
  # Initial value for AK simulation is simple mean of the data #
  
  ak_sim[16] <- mean(data[,16])
  
  ############
  # 1st step #
  ############
  
    for (a in 1:10){
      for (k in 1:10){
        # from 16th full set is made #
        t <- 17
        for (boo in 1:boot_simul) {
          

          
          delta <- NULL
          
          new_input <- sample(sample_no[1:100,t], 100, replace=T)
          old_input <- sample(sample_no[101:1500,t], 1400, replace=T)
          delta[t] <- mean(data[old_input,t]-data[old_input,t-1])
          
          ak_sim[t] <- 1/1500*(  (1-k_i[k]+a_i[a])*sum(data[new_input,t]) + (1-k_i[k]-(100/1400)*a_i[a])*sum(data[old_input,t])  )  + k_i[k]*(ak_sim[t-1] + delta[t])
          boo_sum1[boo] <- ak_sim[t]
        }
        
        # Optimizing A and K value #
        all_variance[a,k] <- var(boo_sum1)
        all_mean[a,k] <- mean(boo_sum1)
        all_mse[a,k] <- ((unemployment[[t]]/100) - all_mean[a,k])^2 + all_variance[a,k]
      }
    }

  ak_max <- which(all_mse == min(all_mse), arr.ind = TRUE)    
  
  a <- a_i[ak_max[1]]
  k <- k_i[ak_max[2]]
  
  # AK estimator value for the optimized A and K value #
  ak_mean[t] <- all_mean[ak_max[1],ak_max[2]]
  ak_var[t] <- all_variance[ak_max[1],ak_max[2]]
    
  ############
  # 2nd step #
  ############
  
  delta <- NULL
  ak_sim2 <-NULL
  
  for (a in 1:10){
    for (k in 1:10){
      for (boo in 1:boot_simul) {
        
        
        t <-18 
        new_input <- sample(sample_no[1:100,t], 100, replace=T)
        old_input <- sample(sample_no[101:1500,t], 1400, replace=T)
        delta[t] <- mean(data[old_input,t]-data[old_input,t-1])
        
        ak_sim2[t] <- 1/1500*(  (1-k_i[k]+a_i[a])*sum(data[new_input,t]) + (1-k_i[k]-(100/1400)*a_i[a])*sum(data[old_input,t])  )  + k_i[k]*(ak_mean[t-1] + delta[t])
        boo_sum2[boo] <- ak_sim2[t]
      }
      
      # Optimizing A and K value #
      all_variance2[a,k] <- var(boo_sum2)
      all_mean2[a,k] <- mean(boo_sum2)
      all_mse[a,k] <- ((unemployment[[t]]/100) - all_mean2[a,k])^2 + all_variance2[a,k]
    }
  }
      
  ak_max <- which(all_mse == min(all_mse), arr.ind = TRUE)    
  
  a <- a_i[ak_max[1]]
  k <- k_i[ak_max[2]]
  
  # AK estimator value for the optimized A and K value #
  ak_mean2[t] <- all_mean2[ak_max[1], ak_max[2]]
  ak_var2[t] <- all_variance2[ak_max[1],ak_max[2]]
      
  ###########
  # 3 Step  #
  ###########
  
  delta <- NULL
  ak_sim3 <-NULL
  
  for (a in 1:10){
    for (k in 1:10){
      for (boo in 1:boot_simul) {

        t <-19 
        
        new_input <- sample(sample_no[1:100,t], 100, replace=T)
        old_input <- sample(sample_no[101:1500,t], 1400, replace=T)
        delta[t] <- mean(data[old_input,t]-data[old_input,t-1])
        
        ak_sim3[t] <- 1/1500*(  (1-k_i[k]+a_i[a])*sum(data[new_input,t]) + (1-k_i[k]-(100/1400)*a_i[a])*sum(data[old_input,t])  )  + k_i[k]*(ak_mean2[t-1] + delta[t])
        boo_sum3[boo] <- ak_sim3[t]
        
      }
      
      # Optimizing A and K value #
      all_variance3[a,k] <- var(boo_sum3)
      all_mean3[a,k] <- mean(boo_sum3)
      all_mse[a,k] <- ((unemployment[[t]]/100) - all_mean3[a,k])^2 + all_variance3[a,k]
    }
  }
  
  ak_max <- which(all_mse == min(all_mse), arr.ind = TRUE)    
  
  a <- a_i[ak_max[1]]
  k <- k_i[ak_max[2]]
  
  # AK estimator value for the optimized A and K value #
  ak_mean3[t] <- all_mean3[ak_max[1],ak_max[2]]
  ak_var3[t] <- all_variance3[ak_max[1],ak_max[2]]
  
  ###########################
  # Global A,K Optimization #
  ###########################

  global_ak_sim[16] <- mean(data[,16])
  
  for (a in 1:10){
    for (k in 1:10){
      for (boo in 1:boot_simul) {
        t <-17 
        
        new_input <- sample(sample_no[1:100,t], 100, replace=T)
        old_input <- sample(sample_no[101:1500,t], 1400, replace=T)
        delta[t] <- mean(data[old_input,t]-data[old_input,t-1])
        
        global_ak_sim[t] <- 1/1500*(  (1-k_i[k]+a_i[a])*sum(data[new_input,t]) + (1-k_i[k]-(100/1400)*a_i[a])*sum(data[old_input,t])  )  + k_i[k]*(global_ak_sim[t-1] + delta[t])
        
        t <-18 
        new_input <- sample(sample_no[1:100,t], 100, replace=T)
        old_input <- sample(sample_no[101:1500,t], 1400, replace=T)
        delta[t] <- mean(data[old_input,t]-data[old_input,t-1])
        
        global_ak_sim[t] <- 1/1500*(  (1-k_i[k]+a_i[a])*sum(data[new_input,t]) + (1-k_i[k]-(100/1400)*a_i[a])*sum(data[old_input,t])  )  + k_i[k]*(global_ak_sim[t-1] + delta[t])
        
        
        t <-19
        new_input <- sample(sample_no[1:100,t], 100, replace=T)
        old_input <- sample(sample_no[101:1500,t], 1400, replace=T)
        delta[t] <- mean(data[old_input,t]-data[old_input,t-1])
        
        global_ak_sim[t] <- 1/1500*(  (1-k_i[k]+a_i[a])*sum(data[new_input,t]) + (1-k_i[k]-(100/1400)*a_i[a])*sum(data[old_input,t])  )  + k_i[k]*(global_ak_sim[t-1] + delta[t])
        
        
        global_booth[boo] <- global_ak_sim[t]
        
      }
      
      # Optimizing A and K value #
      all_variance_global[a,k] <- var(global_booth)
      all_mean_global[a,k] <- mean(global_booth)
      all_mse_global[a,k] <- ((unemployment[[t]]/100) - all_mean_global[a,k])^2 + all_variance_global[a,k]
    }
  }
  
  ak_max_global <- which(all_mse == min(all_mse), arr.ind = TRUE)    
  
  a <- a_i[ak_max_global[1]]
  k <- k_i[ak_max_global[2]]
  
  # AK estimator value for the optimized A and K value #
  ak_mean_global[t] <- all_mean_global[ak_max[1],ak_max[2]]
  ak_var_global[t] <- all_variance_global[ak_max[1],ak_max[2]]
  
  #########################
  #Simple Random Sampling #
  #########################
  
  sample_boot1 <- NULL
  sample_boot2 <- NULL
  sample_boot3 <- NULL
  sample_boot1 <- sample(sample_no[,17], 1500, replace=T)
  sample_boot2 <- sample(sample_no[,18], 1500, replace=T)
  sample_boot3 <- sample(sample_no[,19], 1500, replace=T)
  simple_var[i] <- var(data[sample_boot,t])
  simple_mean[i] <- mean(data[sample_boot,t])
  
  ######################
  # Result Saving Part #
  ######################
  
  ## MEAN and VAR of bootstrap
  boo_var1[i] <- ak_var[17]
  boo_var2[i] <- ak_var2[18]
  boo_var3[i] <- ak_var3[19]
  global_var[i] <- ak_var_global[19]
  simple_var1[i] <- var(data[sample_boot,17])
  simple_var2[i] <- var(data[sample_boot,18])
  simple_var3[i] <- var(data[sample_boot,19])
  
  
  boo_mean1[i] <- ak_mean[17]
  boo_mean2[i] <- ak_mean2[18]
  boo_mean3[i] <- ak_mean3[19]
  global_mean[i] <- ak_mean_global[19]
  simple_mean1[i] <- mean(data[sample_boot,17])
  simple_mean2[i] <- mean(data[sample_boot,18])
  simple_mean3[i] <- mean(data[sample_boot,19])
  
  ######################################
  
  ## Summarizing the performance result
  
  trial_ak_estimator2[i] <- ak_mean[17]
  trial_ak_estimator3[i] <- ak_mean2[18]
  trial_ak_estimator4[i] <- ak_mean3[19]
  
  result[i,1] <- simple_mean1[i]-boo_mean1[i]
  result[i,2] <- simple_mean2[i]-boo_mean2[i]
  result[i,3] <- simple_mean3[i]-boo_mean3[i]
  result[i,4] <- popul[19]/100 - simple_mean3[i]
  result[i,5] <- popul[17]/100 - boo_mean1[i]
  result[i,6] <- popul[18]/100 - boo_mean2[i]
  result[i,7] <- popul[19]/100 - boo_mean3[i]
  
  cat( "overall =", i/simul*100 ,"%\n")
}

time_elap <- proc.time()-e_time

cat( "Bootstraping Done, wish me a good luck!" );
cat("simulation =", simul," ");cat("bootstrap =", boot_simul," "); cat("Time =", time_elap[[3]]/60,"min")

##########################
#### SIMULATION  END #####
##########################


###########################################################

head(result)

#Variance differnce test
var_diff <- matrix(0,1,7)
var_diff[1,1] <- mean(sqrt(all_variance)) / mean(sqrt(boo_var1))
var_diff[1,2] <- mean(sqrt(all_variance2)) / mean(sqrt(boo_var2))
var_diff[1,3] <- mean(sqrt(all_variance3)) / mean(sqrt(boo_var3))
var_diff[1,4] <- mean(sqrt(simple_var1)) / mean(sqrt(boo_var1))
var_diff[1,5] <- mean(sqrt(simple_var2)) / mean(sqrt(boo_var2))
var_diff[1,6] <- mean(sqrt(simple_var3)) / mean(sqrt(boo_var3))
var_diff[1,7] <- mean(sqrt(global_var)) / mean(sqrt(boo_var3))
colnames(var_diff) <- (c("ak2/ak3","ak2/ak4","ak3/ak4","simple/ak2", "simple/ak3","simple/ak4", "Global/Local"))

#mean difference test
mean_diff <- colMeans(result)
print(mean_diff);print(var_diff)

#bootstrap check
boot_result <- matrix(0,1,3)
colnames(boot_result) <- (c("ak_3step","ak_2step","ak_1step"))

boot_result[1] <- mean(sqrt(boo_var3)) / sd(boo_mean3)
boot_result[2] <- mean(sqrt(boo_var2)) / sd(boo_mean2)
boot_result[3] <- mean(sqrt(boo_var1)) / sd(boo_mean1)


mean(boo_mean3); mean(sqrt(boo_var3))
mean(simple_mean); mean(sqrt(simple_var))