################################################################
# RSV seroconversion MSc project
# FOI ratios
# Author: Julia Mayer
# Last updated: 05.08.2022
################################################################
setwd("LSHTM/Project/RSV_infection_NL/CSV files")
getwd()
trajsim2 <- read.csv("Daycare with addition/Updating contacts every 6 months/add_DEzs_trajsim_ll.csv")
spring <- trajsim2$P
winter <- trajsim2$W
autumn <- trajsim2$A
summer <- trajsim2$M

trajsim2$ratiowp <- winter/spring
ratio_winter_spring <- plyr::ddply(.data=trajsim2, .variables="time", function(x) quantile(x[,"ratiowp"], prob = c(0.25, 0.5, 0.75), na.rm=T)) 

trajsim2$ratiowm <- summer/winter
ratio_winter_summer <- plyr::ddply(.data=trajsim2, .variables="time", function(x) quantile(x[,"ratiowm"], prob = c(0.25, 0.5, 0.75), na.rm=T)) 

trajsim2$ratiowa <- winter/autumn
ratio_winter_autumn <- plyr::ddply(.data=trajsim2, .variables="time", function(x) quantile(x[,"ratiowa"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 

trajsim2$ratiopm <- summer/spring
ratio_spring_summer <- plyr::ddply(.data=trajsim2, .variables="time", function(x) quantile(x[,"ratiopm"], prob = c(0.25, 0.5, 0.75), na.rm=T)) 

trajsim2$ratiopa <- spring/autumn
ratio_spring_autumn <- plyr::ddply(.data=trajsim2, .variables="time", function(x) quantile(x[,"ratiopa"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 

trajsim2$ratioam <- summer/autumn
ratio_autumn_summer <- plyr::ddply(.data=trajsim2, .variables="time", function(x) quantile(x[,"ratioam"], prob = c(0.25, 0.5, 0.75), na.rm=T)) 

write.csv(ratio_winter_autumn, "add_daycare_ratio_winter_autumn.csv")
write.csv(ratio_winter_spring, "add__daycare_ratio_winter_spring.csv")
write.csv(ratio_winter_summer, "add_daycare_ratio_winter_summer.csv")
write.csv(ratio_autumn_summer, "add_ratio_autumn_summer.csv")
write.csv(ratio_spring_autumn, "new_ratio_autumn_spring_ll.csv")
write.csv(ratio_spring_summer, "new_ratio_spring_summer.csv")

