setwd("/home/daniel/Dokumente/ueberlagerteSTA/")
library(gtools)
library(ggplot2)


##### GENERATE COMBINATIONS #########################################################

maxTrains <- 15

x <- 0:maxTrains
df <- data.frame(permutations(n=maxTrains+1,r=3,v=x,repeats.allowed=T))
df$SUM <- df$X1 + df$X2 + df$X3
df <- df[df$SUM <=15,]
df <- df[order(df$SUM),]

soll <- read.csv2(file = "./soll.csv", stringsAsFactors = F)

ist <- paste(df$X1, df$X2, df$X3, df$SUM, sep = "#")
soll <- paste(soll$ZCH1_SOLL, soll$ZCH2_SOLL, soll$ZCH3_SOLL, soll$GESAMT_SOLL, sep = "#")

comb <- expand.grid(soll, ist, stringsAsFactors = F)

write.csv2(comb, file = "./combination.csv", row.names = F)

#unlist(strsplit(comb$Var1, split = "#"))



###### CALCULATE ZFW ##########################################################

# soll * (1- ist / soll)²

combinations <- read.csv2(file = "./combination-sep.csv", stringsAsFactors = F, sep = ";")

y <- combinations$ZCH1_SOLL < combinations$ZCH1_IST | combinations$ZCH2_SOLL < combinations$ZCH2_IST | combinations$ZCH3_SOLL < combinations$ZCH3_IST
combinations <- combinations[!y,]


combinations$ZFW <- (combinations$ZCH1_SOLL)^0.75 * (1- 1.0*combinations$ZCH1_IST/combinations$ZCH1_SOLL)^2 +
  combinations$ZCH2_SOLL * (1- 1.0*combinations$ZCH2_IST/combinations$ZCH2_SOLL)^2 +
  combinations$ZCH3_SOLL * (1- 1.0*combinations$ZCH3_IST/combinations$ZCH3_SOLL)^2

k_alpha = 0.07
#k_alpha = 0.6
alpha <- log(2)/k_alpha
beta <- seq(0.1, 0.9, 0.2)
beta <- c(beta, 0.05, 0.03, 0.07)



combinations$ZFW_NEW <- exp(alpha/combinations$ZCH1_SOLL) * (1- exp(-alpha*combinations$ZCH1_IST/combinations$ZCH1_SOLL)) / (exp(alpha/combinations$ZCH1_SOLL) - 1) +
                        exp(alpha/combinations$ZCH2_SOLL) * (1- exp(-alpha*combinations$ZCH2_IST/combinations$ZCH2_SOLL)) / (exp(alpha/combinations$ZCH2_SOLL) - 1) +
                        exp(alpha/combinations$ZCH3_SOLL) * (1- exp(-alpha*combinations$ZCH3_IST/combinations$ZCH3_SOLL)) / (exp(alpha/combinations$ZCH3_SOLL) - 1)
  
combinations$ZFW_V061 <- exp(alpha/combinations$ZCH1_SOLL) * (exp(alpha*beta[1]/combinations$ZCH1_SOLL) - exp(alpha*(beta[1] - combinations$ZCH1_IST)/combinations$ZCH1_SOLL)) / (exp(alpha/combinations$ZCH1_SOLL) - 1) +
                         exp(alpha/combinations$ZCH2_SOLL) * (exp(alpha*beta[1]/combinations$ZCH2_SOLL) - exp(alpha*(beta[1] - combinations$ZCH2_IST)/combinations$ZCH2_SOLL)) / (exp(alpha/combinations$ZCH2_SOLL) - 1) +
                         exp(alpha/combinations$ZCH3_SOLL) * (exp(alpha*beta[1]/combinations$ZCH3_SOLL) - exp(alpha*(beta[1] - combinations$ZCH3_IST)/combinations$ZCH3_SOLL)) / (exp(alpha/combinations$ZCH3_SOLL) - 1)

combinations$ZFW_V063 <- exp(alpha/combinations$ZCH1_SOLL) * (exp(alpha*beta[2]/combinations$ZCH1_SOLL) - exp(alpha*(beta[2] - combinations$ZCH1_IST)/combinations$ZCH1_SOLL)) / (exp(alpha/combinations$ZCH1_SOLL) - 1) +
                         exp(alpha/combinations$ZCH2_SOLL) * (exp(alpha*beta[2]/combinations$ZCH2_SOLL) - exp(alpha*(beta[2] - combinations$ZCH2_IST)/combinations$ZCH2_SOLL)) / (exp(alpha/combinations$ZCH2_SOLL) - 1) +
                         exp(alpha/combinations$ZCH3_SOLL) * (exp(alpha*beta[2]/combinations$ZCH3_SOLL) - exp(alpha*(beta[2] - combinations$ZCH3_IST)/combinations$ZCH3_SOLL)) / (exp(alpha/combinations$ZCH3_SOLL) - 1)

combinations$ZFW_V065 <- exp(alpha/combinations$ZCH1_SOLL) * (exp(alpha*beta[3]/combinations$ZCH1_SOLL) - exp(alpha*(beta[3] - combinations$ZCH1_IST)/combinations$ZCH1_SOLL)) / (exp(alpha/combinations$ZCH1_SOLL) - 1) +
                         exp(alpha/combinations$ZCH2_SOLL) * (exp(alpha*beta[3]/combinations$ZCH2_SOLL) - exp(alpha*(beta[3] - combinations$ZCH2_IST)/combinations$ZCH2_SOLL)) / (exp(alpha/combinations$ZCH2_SOLL) - 1) +
                         exp(alpha/combinations$ZCH3_SOLL) * (exp(alpha*beta[3]/combinations$ZCH3_SOLL) - exp(alpha*(beta[3] - combinations$ZCH3_IST)/combinations$ZCH3_SOLL)) / (exp(alpha/combinations$ZCH3_SOLL) - 1)

combinations$ZFW_V067 <- exp(alpha/combinations$ZCH1_SOLL) * (exp(alpha*beta[4]/combinations$ZCH1_SOLL) - exp(alpha*(beta[4] - combinations$ZCH1_IST)/combinations$ZCH1_SOLL)) / (exp(alpha/combinations$ZCH1_SOLL) - 1) +
                         exp(alpha/combinations$ZCH2_SOLL) * (exp(alpha*beta[4]/combinations$ZCH2_SOLL) - exp(alpha*(beta[4] - combinations$ZCH2_IST)/combinations$ZCH2_SOLL)) / (exp(alpha/combinations$ZCH2_SOLL) - 1) +
                         exp(alpha/combinations$ZCH3_SOLL) * (exp(alpha*beta[4]/combinations$ZCH3_SOLL) - exp(alpha*(beta[4] - combinations$ZCH3_IST)/combinations$ZCH3_SOLL)) / (exp(alpha/combinations$ZCH3_SOLL) - 1)

combinations$ZFW_V069 <- exp(alpha/combinations$ZCH1_SOLL) * (exp(alpha*beta[5]/combinations$ZCH1_SOLL) - exp(alpha*(beta[5] - combinations$ZCH1_IST)/combinations$ZCH1_SOLL)) / (exp(alpha/combinations$ZCH1_SOLL) - 1) +
                         exp(alpha/combinations$ZCH2_SOLL) * (exp(alpha*beta[5]/combinations$ZCH2_SOLL) - exp(alpha*(beta[5] - combinations$ZCH2_IST)/combinations$ZCH2_SOLL)) / (exp(alpha/combinations$ZCH2_SOLL) - 1) +
                         exp(alpha/combinations$ZCH3_SOLL) * (exp(alpha*beta[5]/combinations$ZCH3_SOLL) - exp(alpha*(beta[5] - combinations$ZCH3_IST)/combinations$ZCH3_SOLL)) / (exp(alpha/combinations$ZCH3_SOLL) - 1)

combinations$ZFW_V0605 <- exp(alpha/combinations$ZCH1_SOLL) * (exp(alpha*beta[6]/combinations$ZCH1_SOLL) - exp(alpha*(beta[6] - combinations$ZCH1_IST)/combinations$ZCH1_SOLL)) / (exp(alpha/combinations$ZCH1_SOLL) - 1) +
                          exp(alpha/combinations$ZCH2_SOLL) * (exp(alpha*beta[6]/combinations$ZCH2_SOLL) - exp(alpha*(beta[6] - combinations$ZCH2_IST)/combinations$ZCH2_SOLL)) / (exp(alpha/combinations$ZCH2_SOLL) - 1) +
                          exp(alpha/combinations$ZCH3_SOLL) * (exp(alpha*beta[6]/combinations$ZCH3_SOLL) - exp(alpha*(beta[6] - combinations$ZCH3_IST)/combinations$ZCH3_SOLL)) / (exp(alpha/combinations$ZCH3_SOLL) - 1)

combinations$ZFW_V0603 <- exp(alpha/combinations$ZCH1_SOLL) * (exp(alpha*beta[7]/combinations$ZCH1_SOLL) - exp(alpha*(beta[7] - combinations$ZCH1_IST)/combinations$ZCH1_SOLL)) / (exp(alpha/combinations$ZCH1_SOLL) - 1) +
                          exp(alpha/combinations$ZCH2_SOLL) * (exp(alpha*beta[7]/combinations$ZCH2_SOLL) - exp(alpha*(beta[7] - combinations$ZCH2_IST)/combinations$ZCH2_SOLL)) / (exp(alpha/combinations$ZCH2_SOLL) - 1) +
                          exp(alpha/combinations$ZCH3_SOLL) * (exp(alpha*beta[7]/combinations$ZCH3_SOLL) - exp(alpha*(beta[7] - combinations$ZCH3_IST)/combinations$ZCH3_SOLL)) / (exp(alpha/combinations$ZCH3_SOLL) - 1)

combinations$ZFW_V0607 <- exp(alpha/combinations$ZCH1_SOLL) * (exp(alpha*beta[8]/combinations$ZCH1_SOLL) - exp(alpha*(beta[8] - combinations$ZCH1_IST)/combinations$ZCH1_SOLL)) / (exp(alpha/combinations$ZCH1_SOLL) - 1) +
                          exp(alpha/combinations$ZCH2_SOLL) * (exp(alpha*beta[8]/combinations$ZCH2_SOLL) - exp(alpha*(beta[8] - combinations$ZCH2_IST)/combinations$ZCH2_SOLL)) / (exp(alpha/combinations$ZCH2_SOLL) - 1) +
                          exp(alpha/combinations$ZCH3_SOLL) * (exp(alpha*beta[8]/combinations$ZCH3_SOLL) - exp(alpha*(beta[8] - combinations$ZCH3_IST)/combinations$ZCH3_SOLL)) / (exp(alpha/combinations$ZCH3_SOLL) - 1)



df <- data.frame(1.0*combinations$ZCH1_IST/combinations$ZCH1_SOLL, 
                 1.0*combinations$ZCH2_IST/combinations$ZCH2_SOLL,
                 1.0*combinations$ZCH3_IST/combinations$ZCH3_SOLL)

combinations$MIN_REL <- round(apply(df, 1, FUN=min), 3)
for(i in 1:length(combinations$ZCH1_SOLL)){
  val <- df[i,df[i,]>combinations$MIN_REL[i]+0.001]
  std_dev <- 0
  #if(length(val)==1 && combinations$MIN_REL[i] == 0){
  #  std_dev <- sqrt((val-combinations$MIN_REL[i])^2)
  #}
  if(length(val)>1){
    std_dev <- sd(val)
  }
  combinations$SD_RELgrMIN[i] <- round(std_dev, 3)
}

combinations$ZFW_V3 <- 100*combinations$ZFW_NEW - combinations$SD_RELgrMIN
combinations$ZFW_V4 <- combinations$ZFW_NEW * (2 - combinations$SD_RELgrMIN)
combinations$ZFW_V5 <- 100*combinations$ZFW_NEW - 1.2 * combinations$SD_RELgrMIN


beta2 <- -0.006*(combinations$SD_RELgrMIN)^2

combinations$ZFW_V07 <- exp(alpha/combinations$ZCH1_SOLL) * (exp(alpha*beta2/combinations$ZCH1_SOLL) - exp(alpha*(beta2 - combinations$ZCH1_IST)/combinations$ZCH1_SOLL)) / (exp(alpha/combinations$ZCH1_SOLL) - 1) +
                        exp(alpha/combinations$ZCH2_SOLL) * (exp(alpha*beta2/combinations$ZCH2_SOLL) - exp(alpha*(beta2 - combinations$ZCH2_IST)/combinations$ZCH2_SOLL)) / (exp(alpha/combinations$ZCH2_SOLL) - 1) +
                        exp(alpha/combinations$ZCH3_SOLL) * (exp(alpha*beta2/combinations$ZCH3_SOLL) - exp(alpha*(beta2 - combinations$ZCH3_IST)/combinations$ZCH3_SOLL)) / (exp(alpha/combinations$ZCH3_SOLL) - 1)


qplot(combinations$GESAMT_IST, combinations$ZFW_V5, color = combinations$SD_RELgrMIN) + theme_bw() +scale_color_gradient(low="blue", high="red")
qplot(combinations$MIN_REL, combinations$ZFW_V5, color = combinations$SD_RELgrMIN) + theme_bw() +scale_color_gradient(low="blue", high="red")
qplot(combinations$MIN_REL, combinations$ZFW_V4, color = combinations$SD_RELgrMIN) + theme_bw() +scale_color_gradient(low="blue", high="red")
qplot(combinations$MIN_REL, combinations$ZFW_NEW, color = combinations$SD_RELgrMIN) + theme_bw() +scale_color_gradient(low="blue", high="red")

qplot(combinations$GESAMT_IST, combinations$ZFW_V0605, color = combinations$SD_RELgrMIN) + theme_bw() +scale_color_gradient(low="blue", high="red")
qplot(combinations$MIN_REL, combinations$ZFW_V0605, color = combinations$SD_RELgrMIN) + theme_bw() +scale_color_gradient(low="blue", high="red")

qplot(combinations$GESAMT_IST, combinations$ZFW_V07, color = combinations$SD_RELgrMIN) + theme_bw() +scale_color_gradient(low="blue", high="red")
qplot(combinations$MIN_REL, combinations$ZFW_V07, color = combinations$SD_RELgrMIN) + theme_bw() +scale_color_gradient(low="blue", high="red")

combinations[c(1614,1997),]

#write.csv2(combinations, file = "./combinations-zfw-new-k0,07.csv", row.names = F)
combinations[combinations$GESAMT_SOLL == 15 & combinations$ZCH1_SOLL == 11 & combinations$ZCH2_SOLL == 2 & combinations$GESAMT_IST == 5 ,]

############ CHECK IF CORRECT FOR 1 STA #####################

resultFrame <- data.frame()

for(i_soll in 3:15){
  j_min <- min(combinations$ZCH1_SOLL[combinations$GESAMT_SOLL == i_soll])
  j_max <- max(combinations$ZCH1_SOLL[combinations$GESAMT_SOLL == i_soll])
  for(j in j_min:j_max){
    print(paste("soll", i_soll, "ZCH1", j))
    k_min <- min(combinations$ZCH2_SOLL[combinations$GESAMT_SOLL == i_soll & combinations$ZCH1_SOLL == j])
    k_max <- max(combinations$ZCH2_SOLL[combinations$GESAMT_SOLL == i_soll & combinations$ZCH1_SOLL == j])
    for(k in k_min:k_max){
      i_i_min <- min(combinations$GESAMT_IST[combinations$GESAMT_SOLL == i_soll & combinations$ZCH1_SOLL == j & combinations$ZCH2_SOLL == k])
      i_i_max <- max(combinations$GESAMT_IST[combinations$GESAMT_SOLL == i_soll & combinations$ZCH1_SOLL == j & combinations$ZCH2_SOLL == k])
      for(i_ist in i_i_min:i_i_max){
        tempDf <- combinations[combinations$GESAMT_SOLL == i_soll & 
                                 combinations$ZCH1_SOLL == j & 
                                 combinations$ZCH2_SOLL == k &
                                 combinations$GESAMT_IST == i_ist,]
        
        #id_zfw <- which.max(tempDf$ZFW_NEW)
        #id_zfw <- which.max(tempDf$ZFW_V3)
        #id_zfw <- which.max(tempDf$ZFW_V5)
        #id_zfw <- which.max(tempDf$ZFW_V061)
        #id_zfw <- which.max(tempDf$ZFW_V063)
        #id_zfw <- which.max(tempDf$ZFW_V065)
        #id_zfw <- which.max(tempDf$ZFW_V067)
        #id_zfw <- which.max(tempDf$ZFW_V069)
        #id_zfw <- which.max(tempDf$ZFW_V0605)
        #id_zfw <- which.max(tempDf$ZFW_V0603)
        #id_zfw <- which.max(tempDf$ZFW_V0607)
        id_zfw <- which.max(tempDf$ZFW_V07)
        
        idx11 <- integer(0)
        idx11_err <- integer(0)
        idx12 <- integer(0)
        idx12_err <- integer(0)
        idx21 <- integer(0)
        idx21_err <- integer(0)
        idx22 <- integer(0)
        idx22_err <- integer(0)
        
        # anz identic, min higher --> zfw higher
        for(n in 1:length(tempDf$ZCH1_SOLL)){
          x <- unlist(lapply(tempDf$MIN_REL, function(x) x > tempDf$MIN_REL[n]))
          #y <- unlist(lapply(tempDf$ZFW_NEW[x], function(x) x <= tempDf$ZFW_NEW[n]))
          #y <- unlist(lapply(tempDf$ZFW_V3[x], function(x) x <= tempDf$ZFW_V3[n]))
          #y <- unlist(lapply(tempDf$ZFW_V5[x], function(x) x <= tempDf$ZFW_V5[n]))
          #y <- unlist(lapply(tempDf$ZFW_V061[x], function(x) x <= tempDf$ZFW_V061[n]))
          #y <- unlist(lapply(tempDf$ZFW_V063[x], function(x) x <= tempDf$ZFW_V063[n]))
          #y <- unlist(lapply(tempDf$ZFW_V065[x], function(x) x <= tempDf$ZFW_V065[n]))
          #y <- unlist(lapply(tempDf$ZFW_V067[x], function(x) x <= tempDf$ZFW_V067[n]))
          #y <- unlist(lapply(tempDf$ZFW_V069[x], function(x) x <= tempDf$ZFW_V069[n]))
          #y <- unlist(lapply(tempDf$ZFW_V0605[x], function(x) x <= tempDf$ZFW_V0605[n]))
          #y <- unlist(lapply(tempDf$ZFW_V0603[x], function(x) x <= tempDf$ZFW_V0603[n]))
          #y <- unlist(lapply(tempDf$ZFW_V0607[x], function(x) x <= tempDf$ZFW_V0607[n]))
          y <- unlist(lapply(tempDf$ZFW_V07[x], function(x) x <= tempDf$ZFW_V07[n]))
          if(any(y)){
            z <- rownames(tempDf[x,][y,])
            for(z_t in z){
              idx11 <- c(idx11, which(rownames(combinations) == rownames(tempDf[n,])))
              idx11_err <- c(idx11_err, which(rownames(combinations) == z_t))
            }
            
          }
        }
        # anz identic, min identic, sd smaller --> zfw higher
        # no problem if MIN_REL == 0
        for(n in 1:length(tempDf$ZCH1_SOLL)){
          x1 <- unlist(lapply(tempDf$MIN_REL, function(x) x == tempDf$MIN_REL[n]))
          x2 <- unlist(lapply(tempDf$SD_RELgrMIN, function(x) x < tempDf$SD_RELgrMIN[n]))
          x3 <- unlist(lapply(tempDf$MIN_REL, function(x) x != 0))
          x <- x1 & x2 & x3
          #y <- unlist(lapply(tempDf$ZFW_NEW[x], function(x) x <= tempDf$ZFW_NEW[n]))
          #y <- unlist(lapply(tempDf$ZFW_V3[x], function(x) x <= tempDf$ZFW_V3[n]))
          #y <- unlist(lapply(tempDf$ZFW_V5[x], function(x) x <= tempDf$ZFW_V5[n]))
          #y <- unlist(lapply(tempDf$ZFW_V061[x], function(x) x <= tempDf$ZFW_V061[n]))
          #y <- unlist(lapply(tempDf$ZFW_V063[x], function(x) x <= tempDf$ZFW_V063[n]))
          #y <- unlist(lapply(tempDf$ZFW_V065[x], function(x) x <= tempDf$ZFW_V065[n]))
          #y <- unlist(lapply(tempDf$ZFW_V067[x], function(x) x <= tempDf$ZFW_V067[n]))
          #y <- unlist(lapply(tempDf$ZFW_V069[x], function(x) x <= tempDf$ZFW_V069[n]))
          #y <- unlist(lapply(tempDf$ZFW_V0605[x], function(x) x <= tempDf$ZFW_V0605[n]))
          #y <- unlist(lapply(tempDf$ZFW_V0603[x], function(x) x <= tempDf$ZFW_V0603[n]))
          #y <- unlist(lapply(tempDf$ZFW_V0607[x], function(x) x <= tempDf$ZFW_V0607[n]))
          y <- unlist(lapply(tempDf$ZFW_V07[x], function(x) x <= tempDf$ZFW_V07[n]))
          if(any(y)){
            z <- rownames(tempDf[x,][y,])
            for(z_t in z){
              idx12 <- c(idx12, which(rownames(combinations) == rownames(tempDf[n,])))
              idx12_err <- c(idx12_err, which(rownames(combinations) == z_t))
            }
          }
        }
        
        
        
        
        
        addDf <- combinations[combinations$GESAMT_SOLL == i_soll & 
                                 combinations$ZCH1_SOLL == j & 
                                 combinations$ZCH2_SOLL == k &
                                 combinations$GESAMT_IST == i_ist+1,]
        
        if(length(addDf$ZCH1_SOLL) > 0){
          # anz higher, min higher --> zfw higher
          for(n in 1:length(tempDf$ZCH1_SOLL)){
            x <- unlist(lapply(addDf$MIN_REL, function(x) x > tempDf$MIN_REL[n]))
            #y <- unlist(lapply(addDf$ZFW_NEW[x], function(x) x <= tempDf$ZFW_NEW[n]))
            #y <- unlist(lapply(addDf$ZFW_V3[x], function(x) x <= tempDf$ZFW_V3[n]))
            #y <- unlist(lapply(addDf$ZFW_V5[x], function(x) x <= tempDf$ZFW_V5[n]))
            #y <- unlist(lapply(addDf$ZFW_V061[x], function(x) x <= tempDf$ZFW_V061[n]))
            #y <- unlist(lapply(addDf$ZFW_V063[x], function(x) x <= tempDf$ZFW_V063[n]))
            #y <- unlist(lapply(addDf$ZFW_V065[x], function(x) x <= tempDf$ZFW_V065[n]))
            #y <- unlist(lapply(addDf$ZFW_V067[x], function(x) x <= tempDf$ZFW_V067[n]))
            #y <- unlist(lapply(addDf$ZFW_V069[x], function(x) x <= tempDf$ZFW_V069[n]))
            #y <- unlist(lapply(addDf$ZFW_V0605[x], function(x) x <= tempDf$ZFW_V0605[n]))
            #y <- unlist(lapply(addDf$ZFW_V0603[x], function(x) x <= tempDf$ZFW_V0603[n]))
            #y <- unlist(lapply(addDf$ZFW_V0607[x], function(x) x <= tempDf$ZFW_V0607[n]))
            y <- unlist(lapply(addDf$ZFW_V07[x], function(x) x <= tempDf$ZFW_V07[n]))
            if(any(y)){
              z <- rownames(addDf[x,][y,])
              for(z_t in z){
                idx21 <- c(idx21, which(rownames(combinations) == rownames(tempDf[n,])))
                idx21_err <- c(idx21_err, which(rownames(combinations) == z_t))
              }
            }
          }
          # anz higher, min identic, sd smaller --> zfw higher
          for(n in 1:length(tempDf$ZCH1_SOLL)){
            x1 <- unlist(lapply(addDf$MIN_REL, function(x) x == tempDf$MIN_REL[n]))
            x2 <- unlist(lapply(addDf$SD_RELgrMIN, function(x) x < tempDf$SD_RELgrMIN[n]))
            x3 <- unlist(lapply(addDf$MIN_REL, function(x) x != 0))
            x <- x1 & x2 & x3
            #y <- unlist(lapply(addDf$ZFW_NEW[x], function(x) x <= tempDf$ZFW_NEW[n]))
            #y <- unlist(lapply(addDf$ZFW_V3[x], function(x) x <= tempDf$ZFW_V3[n]))
            #y <- unlist(lapply(addDf$ZFW_V5[x], function(x) x <= tempDf$ZFW_V5[n]))
            #y <- unlist(lapply(addDf$ZFW_V061[x], function(x) x <= tempDf$ZFW_V061[n]))
            #y <- unlist(lapply(addDf$ZFW_V063[x], function(x) x <= tempDf$ZFW_V063[n]))
            #y <- unlist(lapply(addDf$ZFW_V065[x], function(x) x <= tempDf$ZFW_V065[n]))
            #y <- unlist(lapply(addDf$ZFW_V067[x], function(x) x <= tempDf$ZFW_V067[n]))
            #y <- unlist(lapply(addDf$ZFW_V069[x], function(x) x <= tempDf$ZFW_V069[n]))
            #y <- unlist(lapply(addDf$ZFW_V0605[x], function(x) x <= tempDf$ZFW_V0605[n]))
            #y <- unlist(lapply(addDf$ZFW_V0603[x], function(x) x <= tempDf$ZFW_V0603[n]))
            #y <- unlist(lapply(addDf$ZFW_V0607[x], function(x) x <= tempDf$ZFW_V0607[n]))
            y <- unlist(lapply(addDf$ZFW_V07[x], function(x) x <= tempDf$ZFW_V07[n]))
            if(any(y)){
              z <- rownames(addDf[x,][y,])
              for(z_t in z){
                idx22 <- c(idx22, which(rownames(combinations) == rownames(tempDf[n,])))
                idx22_err <- c(idx22_err, which(rownames(combinations) == z_t))
              }
            }
          }
        }
        
        
        if(length(idx11_err)==0){idx11 <- "x"; idx11_err <- "x"}
        if(length(idx12_err)==0){idx12 <- "x"; idx12_err <- "x"}
        if(length(idx21_err)==0){idx21 <- "x"; idx21_err <- "x"}
        if(length(idx22_err)==0){idx22 <- "x"; idx22_err <- "x"}
        
        
        resultFrame <- rbind(resultFrame, data.frame(ID = paste(i_soll, j, k, i_ist, sep = "#"),
                                                     CORRECT_ZFW_11 = length(idx11_err) == 1 && idx11_err=="x",
                                                     INDEX_ZFW_11 = paste(idx11, collapse = "#"),
                                                     ERR_ZFW_11 = paste(idx11_err, collapse= "#"),
                                                     CORRECT_ZFW_12 = length(idx12_err) == 1 && idx12_err=="x",
                                                     INDEX_ZFW_12 = paste(idx12, collapse = "#"),
                                                     ERR_ZFW_12 = paste(idx12_err, collapse= "#"),
                                                     CORRECT_T1_21 = length(idx21_err) == 1 && idx21_err=="x",
                                                     INDEX_T1_21 = paste(idx21, collapse = "#"),
                                                     ERR_T1_21 = paste(idx21_err, collapse = "#"),
                                                     CORRECT_T1_22 = length(idx22_err) == 1 && idx22_err=="x",
                                                     INDEX_T1_22 = paste(idx22, collapse = "#"),
                                                     ERR_T1_22 = paste(idx22_err, collapse = "#")
                                                     ))
      }
    }
  }
}

round(1.0 * sum(resultFrame$CORRECT_ZFW_11) / length(resultFrame$CORRECT_ZFW_11), 3)
round(1.0 * sum(resultFrame$CORRECT_ZFW_12) / length(resultFrame$CORRECT_ZFW_11), 3)
round(1.0 * sum(resultFrame$CORRECT_T1_21) / length(resultFrame$CORRECT_ZFW_11), 3)
round(1.0 * sum(resultFrame$CORRECT_T1_22) / length(resultFrame$CORRECT_ZFW_11), 3)

write.csv2(resultFrame, file = "./check1STA_k0,07.csv", row.names = F)

head(resultFrame[!resultFrame$CORRECT_ZFW_11,], 20)
combinations[c(2394,2490,2000),]

############ CHECK IF CORRECT FOR 2 STAs #####################

### generate 2 STA matrix

indG0 <- which(combinations$MIN_REL > 0)

for(j in 1:2){
  set.seed(j*123)
  ind1 <- sample(indG0,1000)
  ind2 <- sample(indG0,1000)
  
  ind_li <- expand.grid(ind1, ind2, stringsAsFactors = F)
  
  
  
  len <- length(ind_li$Var1)
  
  twoSTA <- data.frame(I = ind_li[, 1], J = ind_li[, 2], SOLL = integer(len), 
                       SUM_ZFW = integer(len), MIN_REL = integer(len), 
                       IST_GESAMT = integer(len), MEAN_SD = integer(len))
  
  
  
  for(i in 1:length(twoSTA$I)){
    if(i %% 1000 == 0){
      
      print(paste(j, i, "   Time:", Sys.time()))
      }
    t1 <- combinations[twoSTA$I[i],]
    t2 <- combinations[twoSTA$J[i],]
    twoSTA$SOLL[i] <- paste(t1$ZCH1_SOLL, t1$ZCH2_SOLL, t1$ZCH3_SOLL, t2$ZCH1_SOLL, t2$ZCH2_SOLL, t2$ZCH3_SOLL, sep = "#")
    twoSTA$SUM_ZFW[i] <- t1$ZFW_V0605 + t2$ZFW_V0605
    twoSTA$MIN_REL[i] <- min(t1$MIN_REL, t2$MIN_REL)
    twoSTA$IST_GESAMT[i] <- t1$GESAMT_IST + t2$GESAMT_IST
    twoSTA$MEAN_SD[i] <- mean(t1$SD_RELgrMIN, t2$SD_RELgrMIN)
  }
  
  write.csv2(twoSTA, file = paste0("./twoSTA_set0", j, "_V0605.csv"), row.names = F)
}



### check results

twoSTA <- read.csv2(file = "./twoSTA_set01_V0605.csv", stringsAsFactors = F)

head(twoSTA)

qplot(twoSTA$IST_GESAMT, twoSTA$SUM_ZFW, color = twoSTA$MEAN_SD) + theme_bw() +scale_color_gradient(low="blue", high="red")
qplot(twoSTA$MIN_REL, twoSTA$SUM_ZFW, color = twoSTA$MEAN_SD) + theme_bw() +scale_color_gradient(low="blue", high="red")

len <- length(twoSTA$I)

resultFrame <- data.frame(CORRECT_ZFW_11 = logical(len),
                          INDEX_ZFW_11 = integer(len),
                          ERR_ZFW_11 = integer(len),
                          ZFW_ERR_11 = integer(len),
                          CORRECT_ZFW_21 = logical(len),
                          INDEX_ZFW_21 = integer(len),
                          ERR_ZFW_21 = integer(len),
                          ZFW_ERR_21 = integer(len))

#write.csv2(resultFrame, file = "./resultFrame_set01_V0605.csv", row.names = F)
resultFrame <- read.csv2(file = "./resultFrame_set01_V0605.csv", stringsAsFactors = F)

for(i in 327037:len){
  if(i %% 200 == 0){print(paste(i, "   Time:", Sys.time()))}
  idx1 <- combinations$ZCH1_SOLL == combinations$ZCH1_SOLL[twoSTA$I[i]] & 
          combinations$ZCH2_SOLL == combinations$ZCH2_SOLL[twoSTA$I[i]] & 
          combinations$ZCH3_SOLL == combinations$ZCH3_SOLL[twoSTA$I[i]]
  idx2 <- combinations$ZCH1_SOLL == combinations$ZCH1_SOLL[twoSTA$J[i]] & 
          combinations$ZCH2_SOLL == combinations$ZCH2_SOLL[twoSTA$J[i]] & 
          combinations$ZCH3_SOLL == combinations$ZCH3_SOLL[twoSTA$J[i]]
  
  i1 <- (1:length(combinations$ZCH1_SOLL))[idx1]
  i2 <- (1:length(combinations$ZCH1_SOLL))[idx2]
  
  mergeSum_IST <- expand.grid(combinations$GESAMT_IST[idx1], combinations$GESAMT_IST[idx2], stringsAsFactors = F)
  idx3 <- apply(mergeSum_IST, 1, sum) == twoSTA$IST_GESAMT[i]
  mergeMin_REL <- apply(expand.grid(combinations$MIN_REL[idx1], combinations$MIN_REL[idx2], stringsAsFactors = F)[idx3,], 1, min)
  mergeZFW <- apply(expand.grid(combinations$ZFW_V5[idx1], combinations$ZFW_V5[idx2], stringsAsFactors = F)[idx3,], 1, sum)
  mergeIndex <- apply(expand.grid(i1, i2, stringsAsFactors = F)[idx3,], 1, paste, collapse = "§")
  
  idx11 <- integer(0)
  idx11_err <- integer(0)
  idx21 <- integer(0)
  idx21_err <- integer(0)
  z11 <- integer(0)
  z21 <- integer(0)
  
  if(sum(idx3) > 0){
    # anz identic min higher --> zfw higher
    for(j in 1:length(mergeZFW)){
      if(twoSTA$MIN_REL[i] < mergeMin_REL[j]){
        if(twoSTA$SUM_ZFW[i] > mergeZFW[j]){
          idx11 <- c(idx11, paste(twoSTA$I[i], twoSTA$J[i], sep = "§"))
          idx11_err <- c(idx11_err, mergeIndex[j])
          z11 <- c(z11, twoSTA$SUM_ZFW[i] / mergeZFW[j])
        }
      }
      if(twoSTA$MIN_REL[i] > mergeMin_REL[j]){
        if(twoSTA$SUM_ZFW[i] < mergeZFW[j]){
          idx11 <- c(idx11, paste(twoSTA$I[i], twoSTA$J[i], sep = "§"))
          idx11_err <- c(idx11_err, mergeIndex[j])
          z11 <- c(z11, mergeZFW[j] / twoSTA$SUM_ZFW[i])
        }
      }
    }
  }
  
  
  idx4 <- apply(mergeSum_IST, 1, sum) == (twoSTA$IST_GESAMT[i]+1)
  mergeMin_REL <- apply(expand.grid(combinations$MIN_REL[idx1], combinations$MIN_REL[idx2], stringsAsFactors = F)[idx4,], 1, min)
  mergeZFW <- apply(expand.grid(combinations$ZFW_V5[idx1], combinations$ZFW_V5[idx2], stringsAsFactors = F)[idx4,], 1, sum)
  mergeIndex <- apply(expand.grid(i1, i2, stringsAsFactors = F)[idx4,], 1, paste, collapse = "§")
  
  if(sum(idx4) > 0){
    # anz higher, min higher --> zfw higher
    for(j in 1:length(mergeZFW)){
      if(twoSTA$MIN_REL[i] < mergeMin_REL[j]){
        if(twoSTA$SUM_ZFW[i] > mergeZFW[j]){
          idx21 <- c(idx21, paste(twoSTA$I[i], twoSTA$J[i], sep = "§"))
          idx21_err <- c(idx21_err, mergeIndex[j])
          z21 <- c(z21, twoSTA$SUM_ZFW[i] / mergeZFW[j])
        }
      }
      if(twoSTA$MIN_REL[i] > mergeMin_REL[j]){
        if(twoSTA$SUM_ZFW[i] < mergeZFW[j]){
          idx21 <- c(idx21, paste(twoSTA$I[i], twoSTA$J[i], sep = "§"))
          idx21_err <- c(idx21_err, mergeIndex[j])
          z21 <- c(z21, mergeZFW[j] / twoSTA$SUM_ZFW[i])
        }
      }
    }
  }
  
  
  
  if(length(idx11_err)==0){idx11 <- "x"; idx11_err <- "x"}
  if(length(idx21_err)==0){idx21 <- "x"; idx21_err <- "x"}
  if(length(z11)==0){z11 <- "x"}
  if(length(z21)==0){z21 <- "x"}
  
  resultFrame$CORRECT_ZFW_11[i] <- length(idx11_err)==1 && idx11_err=="x"
  resultFrame$INDEX_ZFW_11[i] <- paste(idx11, collapse = "#")
  resultFrame$ERR_ZFW_11[i] <- paste(idx11_err, collapse= "#")
  resultFrame$ZFW_ERR_11[i] <- paste(z11, collapse= "#")
  resultFrame$CORRECT_ZFW_21[i] <- length(idx21_err)==1 && idx21_err=="x"
  resultFrame$INDEX_ZFW_21[i] <- paste(idx21, collapse = "#")
  resultFrame$ERR_ZFW_21[i] <- paste(idx21_err, collapse= "#")
  resultFrame$ZFW_ERR_21[i] <- paste(z21, collapse= "#")
}

tempFrame <- resultFrame[1:i,]

1.0 * sum(tempFrame$CORRECT_ZFW_11) / length(tempFrame$CORRECT_ZFW_11) # 0.00114054
1.0 * sum(tempFrame$CORRECT_ZFW_21) / length(tempFrame$CORRECT_ZFW_11) # 0.005773029

errZFW_11 <- tempFrame$ZFW_ERR_11[tempFrame$ZFW_ERR_11 != "x"]

errZFW_11 <- lapply(errZFW_11, strsplit, "#")
corr <- logical(length(errZFW_11))
for(j in 1:length(errZFW_11)){
  if(j %% 1000 == 0) {print(j)}
  corr[j] <- all(unlist(errZFW_11[[j]]) < 1.002)
}

1.0 * (sum(tempFrame$CORRECT_ZFW_11) + sum(corr)) / length(tempFrame$CORRECT_ZFW_11) # 0.7948234

errZFW_21 <- tempFrame$ZFW_ERR_21[tempFrame$ZFW_ERR_21 != "x"]

errZFW_21 <- lapply(errZFW_21, strsplit, "#")
corr <- logical(length(errZFW_21))
for(j in 1:length(errZFW_21)){
  if(j %% 1000 == 0) {print(j)}
  corr[j] <- all(unlist(errZFW_21[[j]]) < 1.002)
}

1.0 * (sum(tempFrame$CORRECT_ZFW_21) + sum(corr)) / length(tempFrame$CORRECT_ZFW_11) # 0.8153503

unlist(errZFW_11[[1]]) > 1.001

combinations[c(6056, 4119, 5000, 4845),]


############################################################
  

qplot(combinations$GESAMT_IST[combinations$GESAMT_SOLL == 15 & combinations$ZCH1_SOLL == 13], 
      combinations$ZFW_NEW[combinations$GESAMT_SOLL == 15 & combinations$ZCH1_SOLL == 13])


combinations[combinations$GESAMT_SOLL == 15 & combinations$ZCH1_SOLL == 10 & combinations$ZCH2_SOLL == 3 & combinations$GESAMT_IST == 9,]
