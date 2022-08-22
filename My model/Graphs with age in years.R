setwd("LSHTM/Project/RSV_infection_NL/CSV files")
getwd()
trajsim2 <- read.csv("Daycare with addition/Updating contacts every 6 months/add_DEzs_trajsim_ll.csv")

trajquantiles <- plyr::ddply(.data=trajsim2, .variables="time", function(x) quantile(x[,"conv"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(trajquantiles) <- c("agemid", "low95", "median", "up95")

trajquantiles_sp <- plyr::ddply(.data=trajsim2, .variables="time", function(x) quantile(x[,"conv_spring"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(trajquantiles_sp) <- c("agemid", "low95", "median", "up95")

trajquantiles_sm <- plyr::ddply(.data=trajsim2, .variables="time", function(x) quantile(x[,"conv_summer"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(trajquantiles_sm) <- c("agemid", "low95", "median", "up95")

trajquantiles_au <- plyr::ddply(.data=trajsim2, .variables="time", function(x) quantile(x[,"conv_autumn"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(trajquantiles_au) <- c("agemid", "low95", "median", "up95")

trajquantiles_wt <- plyr::ddply(.data=trajsim2, .variables="time", function(x) quantile(x[,"conv_winter"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(trajquantiles_wt) <- c("agemid", "low95", "median", "up95")

fit <- ggplot() + theme_bw() + ggtitle("Overall model fit") +
  geom_point(data=data_no_season, aes(x=agemid, y=seroprev_mean)) +
  geom_linerange(data=data_no_season, aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95)) +
  geom_ribbon(data=trajquantiles, aes(x=agemid, ymin=low95, ymax=up95), fill="red", alpha=0.3) +
  geom_line(data=trajquantiles, aes(x=agemid, y=median), color="red") +
  xlab("age (years)") + ylab("proportion seroconverted") +
  scale_x_continuous(breaks=c(0,182.5, 365, 547.5, 730, 912.5, 1095, 1277.5, 1460, 1642.5, 1825), 
                     labels = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5 ,5))

fit

fit_sp2 <- ggplot() + theme_bw() + ggtitle("Model fit on spring birth cohort") +
  geom_point(data=spring_no_daycare, aes(x=agemid, y=seroprev_mean)) +
  geom_linerange(data=spring_no_daycare, aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95)) +
  geom_ribbon(data=spring_conv_quantiles, aes(x=agemid, ymin=low95, ymax=up95), fill="red", alpha=0.3) +
  geom_line(data=spring_conv_quantiles, aes(x=agemid, y=median), color="red") +
  xlab("age (years)") + ylab("proportion seroconverted") +
  scale_x_continuous(breaks=c(0,182.5, 365, 547.5, 730, 912.5, 1095, 1277.5, 1460, 1642.5, 1825), 
                     labels = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5 ,5))


fit_sp2

fit_sm2 <- ggplot() + theme_bw() + ggtitle("Model fit on summer birth cohort") +
  geom_point(data=summer_no_daycare, aes(x=agemid, y=seroprev_mean)) +
  geom_linerange(data=summer_no_daycare, aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95)) +
  geom_ribbon(data=trajquantiles_sm, aes(x=agemid, ymin=low95, ymax=up95), fill="red", alpha=0.3) +
  geom_line(data=trajquantiles_sm, aes(x=agemid, y=median), color="red") +
  xlab("age (years)") + ylab("proportion seroconverted") +
  scale_x_continuous(breaks=c(0,182.5, 365, 547.5, 730, 912.5, 1095, 1277.5, 1460, 1642.5, 1825), 
                     labels = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5 ,5))

fit_sm2

fit_au2 <- ggplot() + theme_bw() + ggtitle("Model fit on autumn birth cohort") +
  geom_point(data=autumn_no_daycare, aes(x=agemid, y=seroprev_mean)) +
  geom_linerange(data=autumn_no_daycare, aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95)) +
  geom_ribbon(data=trajquantiles_au, aes(x=agemid, ymin=low95, ymax=up95), fill="red", alpha=0.3) +
  geom_line(data=trajquantiles_au, aes(x=agemid, y=median), color="red") +
  xlab("age (years)") + ylab("proportion seroconverted") +
  scale_x_continuous(breaks=c(0,182.5, 365, 547.5, 730, 912.5, 1095, 1277.5, 1460, 1642.5, 1825), 
                     labels = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5 ,5))

fit_au2

fit_wt2 <- ggplot() + theme_bw() + ggtitle("Model fit on winter birth cohort") +
  geom_point(data=winter_no_daycare, aes(x=agemid, y=seroprev_mean)) +
  geom_linerange(data=winter_no_daycare, aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95)) +
  geom_ribbon(data=trajquantiles_wt, aes(x=agemid, ymin=low95, ymax=up95), fill="red", alpha=0.3) +
  geom_line(data=trajquantiles_wt, aes(x=agemid, y=median), color="red") +
  xlab("age (years)") + ylab("proportion seroconverted") +
  scale_x_continuous(breaks=c(0,182.5, 365, 547.5, 730, 912.5, 1095, 1277.5, 1460, 1642.5, 1825), 
                     labels = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5 ,5))

fit_wt2
