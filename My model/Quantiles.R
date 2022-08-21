trajsim$FOI_summer <- trajsim$M * 30.41
trajsim$FOI_spring <- trajsim$P * 30.41 + FOI_summer
trajsim$FOI_autumn <- trajsim$A * 30.41 + FOI_summer
trajsim$FOI_winter <- trajsim$W * 30.41 + FOI_summer
trajsim$Contact <- trajsim$C * 30.41

spring_quantiles <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"FOI_spring"], prob = c(0.25, 0.5, 0.75), na.rm=T)) 
summer_quantiles <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"FOI_summer"], prob = c(0.25, 0.5, 0.75), na.rm=T)) 
autumn_quantiles <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"FOI_autumn"], prob = c(0.25, 0.5, 0.75), na.rm=T)) 
winter_quantiles <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"FOI_winter"], prob = c(0.25, 0.5, 0.75), na.rm=T)) 
Contcat_quantiles <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"Contact"], prob = c(0.25, 0.5, 0.75), na.rm=T))
Daycare_quantiles <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"D"], prob = c(0.06, 0.5, 0.94), na.rm=T)) 
immunity_quantiles <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"mu"], prob = c(0.25, 0.5, 0.75), na.rm=T)) 
immunity_quantiles <- 1/immunity_quantiles

trajsim$ratiowp <- trajsim$FOI_winter/trajsim$FOI_spring
ratio_winter_spring <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"ratiowp"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 

trajsim$ratiowm <- trajsim$FOI_summer/trajsim$FOI_winter
ratio_winter_summer <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"ratiowm"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 

trajsim$ratiowa <- trajsim$FOI_winter/trajsim$FOI_autumn
ratio_winter_autumn <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"ratiowa"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 

trajsim$ratiopm <- trajsim$FOI_summer/trajsim$FOI_spring
ratio_spring_summer <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"ratiopm"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 

trajsim$ratiopa <- trajsim$FOI_spring/trajsim$FOI_autumn
ratio_spring_autumn <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"ratiopa"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 

trajsim$ratioam <- trajsim$FOI_summer/trajsim$FOI_autumn
ratio_autumn_summer <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"ratioam"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 

immunity <- 1/wquantiles
head(immunity)

Cquantiles <- Cquantiles * 30.41
