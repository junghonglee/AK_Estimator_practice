
a <- read.csv(file = "C:\\Users\\JHL\\Desktop\\thesis\\Local vs Global\\0125.csv", header=T)
b <- read.csv(file = "C:\\Users\\JHL\\Desktop\\thesis\\Local vs Global\\0130.csv", header=T)
c <- read.csv(file = "C:\\Users\\JHL\\Desktop\\thesis\\Local vs Global\\0131.csv", header=T)

simul_3000 <- rbind( a,b,c)

nrow(simul_3000)

boot_result <- matrix(0,1,4)
colnames(boot_result) <- (c("ak2","ak3","ak4","ak5"))

boot_result[1] <- mean(sqrt(simul_3000[[5]])) / sd(simul_3000[[1]])
boot_result[2] <- mean(sqrt(simul_3000[[6]])) / sd(simul_3000[[2]])
boot_result[3] <- mean(sqrt(simul_3000[[7]])) / sd(simul_3000[[3]])
boot_result[4] <- mean(sqrt(simul_3000[[8]])) / sd(simul_3000[[4]])