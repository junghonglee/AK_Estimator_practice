# Introduction Added 


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

# from 2000 Jan ~ 2015 Dec from ???계청 ?????률 ?????
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
simul <- 1000
boot_simul <- 1000

#Result Data
result <- matrix(0,simul,7)
colnames(result) <- (c("mean-ak2","mean-ak3","mean-ak4","popul-mean","popul-ak2","popul-ak3","popul-ak4"))


# Initial variables
ak_est <- matrix(0,simul,2)
trial_ak_estimator2 <- NULL
trial_ak_estimator3 <- NULL
trial_ak_estimator4 <- NULL
simple_var <- NULL
simple_mean <- NULL

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
  
  sample_mean <- NULL   # IS THIS NEEDED???????????????
  iter <- 0
  
  # Initial value for AK simulation is simple mean of the data #
  
  ak_sim[16] <- mean(data[,16]))


for (a in 1:10){
  for (k in 1:10){
    # from 16th full set is made #
    t <- 17
    for (boo in 1:boot_simul) {
      
      # 1st step #
      delta <- NULL
      ak_sim <- NULL
      
      
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
    ak_max <- which(all_mse == min(all_mse), arr.ind = TRUE)    
    
    a <- a_i[ak_max[1]]
    k <- k_i[ak_max[2]]
    
    # AK estimator value for the optimized A and K value #
    new_input <- sample_no[1:100,t]
    old_input <- sample_no[101:1500,t]
    delta[t] <- mean(data[old_input,t]-data[old_input,t-1])
    ak_sim[t] <- 1/1500*(  (1-k+a)*sum(data[new_input,t]) + (1-k-(100/1400)*a)*sum(data[old_input,t])  )  + k*(ak_sim[t-1] + delta[t])
    ak_mean[t] <- all_mean[a,k]
    ak_var[t] <- all_variance[a,k]
    
    # 2nd step #
    for (boo in 1:boot_simul) {
      delta <- NULL
      ak_sim2 <-NULL
      
      t <-18 
      new_input <- sample(sample_no[1:100,t], 100, replace=T)
      old_input <- sample(sample_no[101:1500,t], 1400, replace=T)
      delta[t] <- mean(data[old_input,t]-data[old_input,t-1])
      
      ak_sim2[t] <- 1/1500*(  (1-k_i[k]+a_i[a])*sum(data[new_input,t]) + (1-k_i[k]-(100/1400)*a_i[a])*sum(data[old_input,t])  )  + k_i[k]*(ak_sim[t-1] + delta[t])
      boo_sum2[boo] <- ak_sim2[t]
    }
    
    # Optimizing A and K value #
    all_variance2[a,k] <- var(boo_sum1)
    all_mean2[a,k] <- mean(boo_sum1)
    all_mse[a,k] <- ((unemployment[[t]]/100) - all_mean2[a,k])^2 + all_variance2[a,k]
    ak_max <- which(all_mse == min(all_mse), arr.ind = TRUE)    
    
    a <- a_i[ak_max[1]]
    k <- k_i[ak_max[2]]
    
    # AK estimator value for the optimized A and K value #
    new_input <- sample_no[1:100,t]
    old_input <- sample_no[101:1500,t]
    delta[t] <- mean(data[old_input,t]-data[old_input,t-1])
    ak_sim2[t] <- 1/1500*(  (1-k+a)*sum(data[new_input,t]) + (1-k-(100/1400)*a)*sum(data[old_input,t])  )  + k*(ak_sim[t-1] + delta[t])
    ak_mean2[t] <- all_mean2[a,k]
    ak_var2[t] <- all_variance2[a,k]
    
    ###########
    # 3 Step  #
    ###########
    
    for (boo in 1:boot_simul) {
      
      delta <- NULL
      ak_sim3 <-NULL
      
      t <-19 
      
      new_input <- sample(sample_no[1:100,t], 100, replace=T)
      old_input <- sample(sample_no[101:1500,t], 1400, replace=T)
      delta[t] <- mean(data[old_input,t]-data[old_input,t-1])
      
      ak_sim3[t] <- 1/1500*(  (1-k_i[k]+a_i[a])*sum(data[new_input,t]) + (1-k_i[k]-(100/1400)*a_i[a])*sum(data[old_input,t])  )  + k_i[k]*(ak_sim2[t-1] + delta[t])
      boo_sum3[boo] <- ak_sim3[t]
      
    }
    
    # Optimizing A and K value #
    all_variance3[a,k] <- var(boo_sum1)
    all_mean3[a,k] <- mean(boo_sum1)
    all_mse[a,k] <- ((unemployment[[t]]/100) - all_mean3[a,k])^2 + all_variance3[a,k]
    ak_max <- which(all_mse == min(all_mse), arr.ind = TRUE)    
    
    a <- a_i[ak_max[1]]
    k <- k_i[ak_max[2]]
    
    # AK estimator value for the optimized A and K value #
    new_input <- sample_no[1:100,t]
    old_input <- sample_no[101:1500,t]
    delta[t] <- mean(data[old_input,t]-data[old_input,t-1])
    ak_sim3[t] <- 1/1500*(  (1-k+a)*sum(data[new_input,t]) + (1-k-(100/1400)*a)*sum(data[old_input,t])  )  + k*(ak_sim2[t-1] + delta[t])
    ak_mean3[t] <- all_mean3[a,k]
    ak_var3[t] <- all_variance3[a,k]
  }
  
  iter = iter+1
  #cat( "bootsimul =", (iter)/(boot_simul) ,"% \t", "overall =", i/simul*100 ,"%\n")
  
}

#########################
#Simple Random Sampling #
#########################

sample_boot <- NULL
sample_boot <- sample(sample_no[,t], 1500, replace=T)
simple_var[i] <- var(data[sample_boot,t])
simple_mean[i] <- mean(data[sample_boot,t])


######################################
################# 요기부터 다시 시작 #


## MEAN and VAR of bootstrap
boo_var1[i] <- all_variance[ak_max[1],ak_max[2]]
boo_var2[i] <- all_variance2[ak_max[1],ak_max[2]]
boo_var3[i] <- all_variance3[ak_max[1],ak_max[2]]
simple_var[i] <- var(data[sample_boot,t])

boo_mean1[i] <- all_mean[ak_max[1],ak_max[2]]
boo_mean2[i] <- all_mean2[ak_max[1],ak_max[2]]
boo_mean3[i] <- all_mean3[ak_max[1],ak_max[2]]
simple_mean[i] <- mean(data[sample_boot,t])


######################################

## Summarizing the performance result
trial_ak_estimator2[i] <- ak_sim[t-2]
trial_ak_estimator3[i] <- ak_sim2[t-1]
trial_ak_estimator4[i] <- ak_sim3[t]


result[i,1] <- simple_mean[i]-boo_mean1[i]
result[i,2] <- simple_mean[i]-boo_mean2[i]
result[i,3] <- simple_mean[i]-boo_mean3[i]
result[i,4] <- popul[t]/100 - simple_mean[i]
result[i,5] <- popul[t]/100 - ak_sim[t-2]
result[i,6] <- popul[t]/100 - ak_sim2[t-1]
result[i,7] <- popul[t]/100 - ak_sim3[t]


cat( "overall =", i/simul*100 ,"%\n")
}

###########################################################
###########################################################

time_elap <- proc.time()-e_time

cat( "Bootstraping Done, wish me a good luck!" );
cat("simulation =", simul," ");cat("bootstrap =", boot_simul," "); cat("Time =", time_elap[[3]]/60,"min")


##########################
#### SIMULATION  END #####
##########################

head(result)
head(ak_est)

###########################################################
summary(ak_est)

Mode(ak_est[,1])
Mode(ak_est[,2])

mean(ak_est[,1])
mean(ak_est[,2])


#Variance differnce test
var_diff <- matrix(0,1,6)
var_diff[1,1] <- mean(sqrt(all_variance)) / mean(sqrt(boo_var1))
var_diff[1,2] <- mean(sqrt(all_variance)) / mean(sqrt(boo_var2))
var_diff[1,3] <- mean(sqrt(all_variance2)) / mean(sqrt(boo_var3))
var_diff[1,4] <- mean(sqrt(simple_var)) / mean(sqrt(boo_var1))
var_diff[1,5] <- mean(sqrt(simple_var)) / mean(sqrt(boo_var2))
var_diff[1,6] <- mean(sqrt(simple_var)) / mean(sqrt(boo_var3))
colnames(var_diff) <- (c("ak2/ak3","ak2/ak4","ak3/ak4","simple/ak2", "simple/ak3","simple/ak4"))

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