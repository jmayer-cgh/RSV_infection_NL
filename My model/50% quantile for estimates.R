setwd("LSHTM/Project/RSV_infection_NL/CSV files")
getwd()
trajsim2 <- read.csv("Daycare with addition/Updating contacts every 6 months/add_DEzs_trajsim_ll.csv")

Dquantiles <- plyr::ddply(.data=trajsim2, .variables="time", function(x) quantile(x[,"D"], prob = c(0.25, 0.5, 0.75), na.rm=T)) 
colnames(Dquantiles) <- c("agemid", "low50", "median", "up50")

Pquantiles <- plyr::ddply(.data=trajsim2, .variables="time", function(x) quantile(x[,"P"], prob = c(0.25, 0.5, 0.75), na.rm=T))
Pquantiles <- Pquantiles*30.41
Pquantiles <- Pquantiles + Mquantiles
colnames(Pquantiles) <- c("agemid", "low50", "median", "up50")

Mquantiles <- plyr::ddply(.data=trajsim2, .variables="time", function(x) quantile(x[,"M"], prob = c(0.25, 0.5, 0.75), na.rm=T)) 
Mquantiles <- Mquantiles*30.41
colnames(Mquantiles) <- c("agemid", "low50", "median", "up50")

Aquantiles <- plyr::ddply(.data=trajsim2, .variables="time", function(x) quantile(x[,"A"], prob = c(0.25, 0.5, 0.75), na.rm=T)) 
Aquantiles <- Aquantiles*30.41
Aquantiles <- Aquantiles + Mquantiles
colnames(Aquantiles) <- c("agemid", "low50", "median", "up50")

Wquantiles <- plyr::ddply(.data=trajsim2, .variables="time", function(x) quantile(x[,"W"], prob = c(0.25, 0.5, 0.75), na.rm=T)) 
Wquantiles <- Wquantiles*30.41
Wquantiles <- Wquantiles + Mquantiles
colnames(Wquantiles) <- c("agemid", "low50", "median", "up50")

Cquantiles <- plyr::ddply(.data=trajsim2, .variables="time", function(x) quantile(x[,"C"], prob = c(0.25, 0.5, 0.75), na.rm=T)) 
Cquantiles <- Cquantiles*30.41
colnames(Cquantiles) <- c("agemid", "low50", "median", "up50")

wquantiles <- plyr::ddply(.data=trajsim2, .variables="time", function(x) quantile(x[,"mu"], prob = c(0.25, 0.5, 0.75), na.rm=T)) 
wquantiles <- 1/wquantiles
colnames(wquantiles) <- c("agemid", "low50", "median", "up50")
