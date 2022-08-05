################################################################
# RSV seroconversion MSc project
# FOI ratios
# Author: Julia Mayer
# Last updated: 05.08.2022
################################################################
setwd("LSHTM/Project/RSV_infection_NL")
getwd()
trajsim <- read.csv("Daycare with DEzs/DEzs_trajsim.csv")
spring <- trajsim$P
winter <- trajsim$W
autumn <- trajsim$A
summer <- trajsim$M

trajsim$ratiowp <- winter/spring
ratio_winter_spring <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"ratiowp"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 

trajsim$ratiowm <- winter/summer
ratio_winter_summer <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"ratiowm"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 

trajsim$ratiowa <- winter/autumn
ratio_winter_autumn <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"ratiowa"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 

trajsim$ratiopm <- spring/summer
ratio_spring_summer <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"ratiopm"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 

trajsim$ratiopa <- spring/autumn
ratio_spring_autumn <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"ratiopa"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 

trajsim$ratioam <- autumn/summer
ratio_autumn_summer <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"ratioam"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 

write.csv(ratio_winter_autumn, "ratio_winter_autumn.csv")
write.csv(ratio_winter_spring, "ratio_winter_spring.csv")
write.csv(ratio_winter_summer, "ratio_winter_summer.csv")
write.csv(ratio_autumn_summer, "ratio_autumn_summer.csv")
write.csv(ratio_spring_autumn, "ratio_autumn_spring.csv")
write.csv(ratio_spring_summer, "ratio_spring_summer.csv")


winter_FOI <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"W"], prob = c(0.05, 0.5, 0.95), na.rm=T)) 
