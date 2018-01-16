
data <- read.csv("C:\\Inchon\\Q3_data.csv", header=T)
names(data)

data$숙박유형 <- trimws(data$숙박유형)   #중저가
data$숙박유형.1- trimws(data$숙박유형.1) #특2급이상

data$숙박일수       #중저가
data$숙박일수.1)    # 특2급이상

data$숙박일수[is.na(data$숙박일수)] <- 0
data$숙박일수.1[is.na(data$숙박일수.1)] <- 0

cbind(data$숙박유형, data$숙박숙박유형.1)

length(data$숙박일수)

place <- 0

#중저가 코딩 작업
for (i in 1:2099){
  if (data$숙박유형[i]=="") {
    place[i] <- 0}
  
  else if(data$숙박일수[i]=="0"){
    place[i] <- 0
  }
  
  else {
    place[i] <- 1
  }
}


place2 <- 0
#특2급이상 코딩 작업
for (i in 1:2099){
  if (data$숙박유형.1[i]=="") {
    place2[i] <- 0}
  
  else if(data$숙박일수.1[i]=="0"){
    place2[i] <- 0
  }
  
  else {
    place2[i] <- 2
  }
}

place_total <- place + place2 # 1: 중저가, 2:특2급 이상, 3: 둘다 


