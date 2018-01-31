#collecting data from 2 simulations
a <- read.csv(file = "C:\\Users\\JHL\\Desktop\\thesis\\Local vs Global\\0130.csv", header=T)
b <- read.csv(file = "C:\\Users\\JHL\\Desktop\\thesis\\Local vs Global\\0131.csv", header=T)

#combining the data
simul_2000 <- rbind(a,b)

#Checking if data is well binded
ncol(simul_2000)
nrow(simul_2000)

#Testing the result of the bootstrap
boot_result <- matrix(0,1,4)
colnames(boot_result) <- (c("ak2","ak3","ak4","ak5"))

boot_result[1] <- mean(sqrt(simul_2000[[6]])) / sd(simul_2000[[2]])
boot_result[2] <- mean(sqrt(simul_2000[[7]])) / sd(simul_2000[[3]])
boot_result[3] <- mean(sqrt(simul_2000[[8]])) / sd(simul_2000[[4]])
boot_result[4] <- mean(sqrt(simul_2000[[9]])) / sd(simul_2000[[5]])
print(boot_result)

