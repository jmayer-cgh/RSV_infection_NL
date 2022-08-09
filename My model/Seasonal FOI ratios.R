################################################################
# RSV seroconversion MSc project
# FOI ratios
# Author: Julia Mayer
# Last updated: 05.08.2022
################################################################
setwd("../")
getwd()
trajsim2 <- read.csv("Seasons with correct equation, return ll/S_trajsim_ll.csv")
spring <- trajsim2$P
winter <- trajsim2$W
autumn <- trajsim2$A
summer <- trajsim2$M

trajsim2$ratiowp <- winter/spring
ratio_winter_spring <- plyr::ddply(.data=trajsim2, .variables="time", function(x) quantile(x[,"ratiowp"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 

trajsim2$ratiowm <- winter/summer
ratio_winter_summer <- plyr::ddply(.data=trajsim2, .variables="time", function(x) quantile(x[,"ratiowm"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 

trajsim2$ratiowa <- winter/autumn
ratio_winter_autumn <- plyr::ddply(.data=trajsim2, .variables="time", function(x) quantile(x[,"ratiowa"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 

trajsim2$ratiopm <- spring/summer
ratio_spring_summer <- plyr::ddply(.data=trajsim2, .variables="time", function(x) quantile(x[,"ratiopm"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 

trajsim2$ratiopa <- spring/autumn
ratio_spring_autumn <- plyr::ddply(.data=trajsim2, .variables="time", function(x) quantile(x[,"ratiopa"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 

trajsim2$ratioam <- autumn/summer
ratio_autumn_summer <- plyr::ddply(.data=trajsim2, .variables="time", function(x) quantile(x[,"ratioam"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 

write.csv(ratio_winter_autumn, "new_ratio_winter_autumn_ll.csv")
write.csv(ratio_winter_spring, "new_ratio_winter_spring_ll.csv")
write.csv(ratio_winter_summer, "new_ratio_winter_summer_ll.csv")
write.csv(ratio_autumn_summer, "new_ratio_autumn_summer_ll.csv")
write.csv(ratio_spring_autumn, "new_ratio_autumn_spring_ll.csv")
write.csv(ratio_spring_summer, "new_ratio_spring_summer_ll.csv")

