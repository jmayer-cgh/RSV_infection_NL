## Core equations for transitions between compartments:
# spring birth cohort
update(M1_sp) <- M1_sp - n_M1M2_sp
update(M2_sp) <- M2_sp + n_M1M2_sp - n_M2S_sp
update(S_sp) <- S_sp + n_M2S_sp - n_SR_sp
update(R_sp) <- R_sp + n_SR_sp

# summer birth cohort
update(M1_sm) <- M1_sm - n_M1M2_sm
update(M2_sm) <- M2_sm + n_M1M2_sm - n_M2S_sm
update(S_sm) <- S_sm + n_M2S_sm - n_SR_sm
update(R_sm) <- R_sm + n_SR_sm

# autumn birth cohort
update(M1_au) <- M1_au - n_M1M2_au
update(M2_au) <- M2_au + n_M1M2_au - n_M2S_au
update(S_au) <- S_au + n_M2S_au - n_SR_au
update(R_au) <- R_au + n_SR_au

# winter birth cohort
update(M1_wt) <- M1_wt - n_M1M2_wt
update(M2_wt) <- M2_wt + n_M1M2_wt - n_M2S_wt
update(S_wt) <- S_wt + n_M2S_wt - n_SR_wt
update(R_wt) <- R_wt + n_SR_wt

# Total seroprevalence
update(R_all) <- 0.26*R_sp + 
  0.29*R_sm + 
  0.24*R_au + 
  0.20*R_wt

## Individual probabilities of transition:
# spring birth cohort
n_M1M2_sp <- mu*2 * dt * M1_sp # M1 to M2
n_M2S_sp <- mu*2 * dt * M2_sp # M2 to S
n_SR_sp <- lambda_sp * dt * S_sp # S to R

# summer birth cohort
n_M1M2_sm <- mu*2 * dt * M1_sm # M1 to M2
n_M2S_sm <- mu*2 * dt * M2_sm # M2 to S
n_SR_sm <- lambda_sm * dt * S_sm # S to R

# autumn birth cohort
n_M1M2_au <- mu*2 * dt * M1_au # M1 to M2
n_M2S_au <- mu*2 * dt * M2_au # M2 to S
n_SR_au <- lambda_au * dt * S_au # S to R

# winter birth cohort
n_M1M2_wt <- mu*2 * dt * M1_wt# M1 to M2
n_M2S_wt <- mu*2 * dt * M2_wt # M2 to S
n_SR_wt <- lambda_wt * dt * S_wt# S to R

## Building the FOI
# Define booleans to know which season the cohort is in
spring_FOI_sp <- if ((time <= 30.41*3) ||  
                     ((time >= 365) && (time <= (365+30.41*3))) ||
                     ((time >= 2*365) && (time <= (2*365+30.41*3))) ||
                     ((time >= 3*365) && (time <= (3*365+30.41*3))) || 
                     ((time >= 4*365) && (time <= (4*365+30.41*3))) || 
                     ((time >= 5*365) && (time <= (5*365+30.41*3)))) 1 else 0
summer_FOI_sp <- if ((time > 30.41*3 && time <= 30.41*6 ) ||                      # FOI in summer for those born in spring
                     ((time > 365 + 30.41*3) && (time <= (365+30.41*6))) ||
                     ((time > 2*365 + 30.41*3) && (time <= (2*365+30.41*6))) ||
                     ((time > 3*365 + 30.41*3) && (time <= (3*365+30.41*6))) ||
                     ((time > 4*365 + 30.41*3) && (time <= (4*365+30.41*6))) ||
                     ((time > 5*365 + 30.41*3) && (time <= (5*365+30.41*6)))) 1 else 0 
autumn_FOI_sp <- if ( (time > 30.41*6 && time <= 30.41*9 ) ||                    
                      ((time > 365+30.41*6) && (time <= (365+30.41*9))) ||
                      ((time > 2*365 + 30.41*6) && (time <= (2*365+30.41*9))) ||
                      ((time > 3*365 + 30.41*6) && (time <= (3*365+30.41*9))) ||
                      ((time > 4*365 + 30.41*6) && (time <= (4*365+30.41*9))) ||
                      ((time > 5*365 + 30.41*6) && (time <= (5*365+30.41*9)))) 1 else 0
winter_FOI_sp <- if ( (time > 30.41*9 && time <= 30.41*12 ) ||                   
                      ((time > 365+30.41*9) && (time <= (365+30.41*12))) ||
                      ((time > 2*365 + 30.41*9) && (time <= (2*365+30.41*12))) ||
                      ((time > 3*365 + 30.41*9) && (time <= (3*365+30.41*12))) ||
                      ((time > 4*365 + 30.41*9) && (time <= (4*365+30.41*12))) ||
                      ((time > 5*365 + 30.41*9) && (time <= (5*365+30.41*12)))) 1 else 0
spring_FOI_sm <- if ( (time > 30.41*9 && time <= 30.41*12 )  ||                 
                      ((time > 365+30.41*9) && (time <= (365+30.41*12))) ||
                      ((time > 2*365 + 30.41*9) && (time <= (2*365+30.41*12))) ||
                      ((time > 3*365 + 30.41*9) && (time <= (3*365+30.41*12))) ||
                      ((time > 4*365 + 30.41*9) && (time <= (4*365+30.41*12))) ||
                      ((time > 5*365 + 30.41*9) && (time <= (5*365+30.41*12)))) 1 else 0
summer_FOI_sm <- if((time <= 30.41*3) ||                                     
                    ((time >= 365) && (time <= (365+30.41*3))) ||
                    ((time >= 2*365) && (time <= (2*365+30.41*3))) ||
                    ((time >= 3*365) && (time <= (3*365+30.41*3))) ||
                    ((time >= 4*365) && (time <= (4*365+30.41*3))) ||
                    ((time >= 5*365) && (time <= (5*365+30.41*3)))) 1 else 0
autumn_FOI_sm <- if ((time > 30.41*3 && time <= 30.41*6 ) ||                       # FOI in autumn for those born in summer
                     ((time > 365 + 30.41*3) && (time <= (365+30.41*6))) ||
                     ((time > 2*365 + 30.41*3) && (time <= (2*365+30.41*6))) || 
                     ((time > 3*365 + 30.41*3) && (time <= (3*365+30.41*6))) ||
                     ((time > 4*365 + 30.41*3) && (time <= (4*365+30.41*6))) ||
                     ((time > 5*365 + 30.41*3) && (time <= (5*365+30.41*6)))) 1 else 0
winter_FOI_sm <- if ( (time > 30.41*6 && time <= 30.41*9 ) ||                    
                      ((time > 365+30.41*6) && (time <= (365+30.41*9))) ||
                      ((time > 2*365 + 30.41*6) && (time <= (2*365+30.41*9))) ||
                      ((time > 3*365 + 30.41*6) && (time <= (3*365+30.41*9))) ||
                      ((time > 4*365 + 30.41*6) && (time <= (4*365+30.41*9))) ||
                      ((time > 5*365 + 30.41*6) && (time <= (5*365+30.41*9)))) 1 else 0
spring_FOI_au <- if ( (time > 30.41*6 && time <= 30.41*9 ) ||                    
                      ((time > 365+30.41*6) && (time <= (365+30.41*9))) ||
                      ((time > 2*365 + 30.41*6) && (time <= (2*365+30.41*9))) ||
                      ((time > 3*365 + 30.41*6) && (time <= (3*365+30.41*9))) ||
                      ((time > 4*365 + 30.41*6) && (time <= (4*365+30.41*9))) ||
                      ((time > 5*365 + 30.41*6) && (time <= (5*365+30.41*9)))) 1 else 0
summer_FOI_au <- if ( (time > 30.41*9 && time <= 30.41*12 ) ||                   
                      ((time > 365+30.41*9) && (time <= (365+30.41*12))) ||
                      ((time > 2*365 + 30.41*9) && (time <= (2*365+30.41*12))) ||
                      ((time > 3*365 + 30.41*9) && (time <= (3*365+30.41*12))) ||
                      ((time > 4*365 + 30.41*9) && (time <= (4*365+30.41*12))) ||
                      ((time > 5*365 + 30.41*9) && (time <= (5*365+30.41*12)))) 1 else 0
autumn_FOI_au <- if((time <= 30.41*3) ||                                      
                    ((time >= 365) && (time <= (365+30.41*3))) ||
                    ((time >= 2*365) && (time <= (2*365+30.41*3))) ||
                    ((time >= 3*365) && (time <= (3*365+30.41*3))) ||
                    ((time >= 4*365) && (time <= (4*365+30.41*3))) ||
                    ((time >= 5*365) && (time <= (5*365+30.41*3)))) 1 else 0
winter_FOI_au <- if ((time > 30.41*3 && time <= 30.41*6 ) ||                      # FOI in winter for those born in autumn
                     ((time > 365 + 30.41*3) && (time <= (365+30.41*6))) ||
                     ((time > 2*365 + 30.41*3) && (time <= (2*365+30.41*6))) || 
                     ((time > 3*365 + 30.41*3) && (time <= (3*365+30.41*6))) ||
                     ((time > 4*365 + 30.41*3) && (time <= (4*365+30.41*6))) ||
                     ((time > 5*365 + 30.41*3) && (time <= (5*365+30.41*6)))) 1 else 0
spring_FOI_wt <- if ((time > 30.41*3 && time <= 30.41*6 ) ||                      # FOI in spring for those born in winter
                     ((time > 365 + 30.41*3) && (time <= (365+30.41*6))) ||
                     ((time > 2*365 + 30.41*3) && (time <= (2*365+30.41*6))) ||
                     ((time > 3*365 + 30.41*3) && (time <= (3*365+30.41*6))) ||
                     ((time > 4*365 + 30.41*3) && (time <= (4*365+30.41*6))) ||
                     ((time > 5*365 + 30.41*3) && (time <= (5*365+30.41*6)))) 1 else 0
summer_FOI_wt <- if ( (time > 30.41*6 && time <= 30.41*9 ) ||                     
                      ((time > 365+30.41*6) && (time <= (365+30.41*9))) ||
                      ((time > 2*365 + 30.41*6) && (time <= (2*365+30.41*9))) ||
                      ((time > 3*365 + 30.41*6) && (time <= (3*365+30.41*9))) ||
                      ((time > 4*365 + 30.41*6) && (time <= (4*365+30.41*9))) ||
                      ((time > 5*365 + 30.41*6) && (time <= (5*365+30.41*9)))) 1 else 0
autumn_FOI_wt <- if ( (time > 30.41*9 && time <= 30.41*12 )     ||              
                      ((time > 365+30.41*9) && (time <= (365+30.41*12))) ||
                      ((time > 2*365 + 30.41*9) && (time <= (2*365+30.41*12))) ||
                      ((time > 3*365 + 30.41*9) && (time <= (3*365+30.41*12))) ||
                      ((time > 4*365 + 30.41*9) && (time <= (4*365+30.41*12))) ||
                      ((time > 5*365 + 30.41*9) && (time <= (5*365+30.41*12)))) 1 else 0
winter_FOI_wt <- if( (time <= 30.41*3) ||                                     
                     ((time >= 365) && (time <= (365+30.41*3))) ||
                     ((time >= 2*365) && (time <= (2*365+30.41*3))) ||
                     ((time >= 3*365) && (time <= (3*365+30.41*3))) ||
                     ((time >= 4*365) && (time <= (4*365+30.41*3))) ||
                     ((time >= 5*365) && (time <= (5*365+30.41*3)))) 1 else 0

# Putting it all together into four FOIs
lambda_sp = (summer_comp + spring_comp) * spring_FOI_sp + 
  summer_comp * summer_FOI_sp + 
  (summer_comp + autumn_comp) * autumn_FOI_sp +
  (summer_comp + winter_comp) * winter_FOI_sp 
lambda_sm = (summer_comp + spring_comp) * spring_FOI_sm + 
  summer_comp * summer_FOI_sm + 
  (summer_comp + autumn_comp) * autumn_FOI_sm + 
  (summer_comp + winter_comp) * winter_FOI_sm 
lambda_au = (summer_comp + spring_comp) * spring_FOI_au + 
  summer_comp * summer_FOI_au + 
  (summer_comp + autumn_comp) * autumn_FOI_au + 
  (summer_comp + winter_comp) * winter_FOI_au 
lambda_wt = (summer_comp + spring_comp) * spring_FOI_wt + 
  summer_comp * summer_FOI_wt +
  (summer_comp + autumn_comp) * autumn_FOI_wt + 
  (summer_comp + winter_comp) * winter_FOI_wt 

## Initial states:
# spring cohort
initial(M1_sp) <- prop * M1_sp_ini # only a proportion of children is born with maternal immunity
initial(M2_sp) <- M2_sp_ini
initial(S_sp) <- (1-prop) * M1_sp_ini
initial(R_sp) <- R_sp_ini

# summer cohort
initial(M1_sm) <- prop * M1_sm_ini
initial(M2_sm) <- M2_sm_ini
initial(S_sm) <- (1-prop) * M1_sm_ini
initial(R_sm) <- R_sm_ini

# autumn cohort
initial(M1_au) <- prop * M1_au_ini
initial(M2_au) <- M2_au_ini
initial(S_au) <- (1-prop) * M1_au_ini
initial(R_au) <- R_au_ini

# winter cohort
initial(M1_wt) <- prop * M1_wt_ini
initial(M2_wt) <- M2_wt_ini
initial(S_wt) <- (1-prop) * M1_wt_ini
initial(R_wt) <- R_wt_ini

# Total
initial(R_all) <- 0.26*R_sp_ini + 
  0.29*R_sm_ini + 
  0.24*R_au_ini + 
  0.20*R_wt_ini

## User defined parameters - default in parentheses:
M1_sp_ini <- parameter(1 - 2 * 1e-12)
M2_sp_ini <- parameter(1e-12)
R_sp_ini <- parameter(1e-12)

M1_sm_ini <- parameter(1 - 2 * 1e-12)
M2_sm_ini <- parameter(1e-12)
R_sm_ini <- parameter(1e-12)

M1_au_ini <- parameter(1 - 2 * 1e-12)
M2_au_ini <- parameter(1e-12)
R_au_ini <- parameter(1e-12)

M1_wt_ini <- parameter(1 - 2 * 1e-12)
M2_wt_ini <- parameter(1e-12)
R_wt_ini <- parameter(1e-12)


# transition parameters
mu_mean <- parameter(1/59.50121)
mu_sd <- parameter(1/59.50121)
spring_comp_mean <- parameter(1e-05)
spring_comp_sd <- parameter(1e-05)
summer_comp_mean <- parameter(0.02002)
summer_comp_sd <- parameter(0.02002)
autumn_comp_mean <- parameter(3e-05)
autumn_comp_sd <- parameter(3e-05)
winter_comp_mean <- parameter(4e-05)
winter_comp_sd <- parameter(4e-05)

mu <- Normal(mu_mean, mu_sd)
spring_comp <- Normal(spring_comp_mean, spring_comp_sd)
summer_comp <- Normal(summer_comp_mean, summer_comp_sd)
autumn_comp <- Normal(autumn_comp_mean, autumn_comp_sd)
winter_comp <- Normal(winter_comp_mean, winter_comp_sd)

# Proportion born with maternal immunity
prop_mean <- parameter(1)
prop_sd <- parameter(1)
prop <- Normal(prop_mean, prop_sd)