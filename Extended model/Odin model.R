## Definition of the time-step and output as "time"
dt <- user(1)
initial(time) <- 1
update(time) <- (step + 1) * dt


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
update(R_all) <- R_sp + R_sm + R_au + R_wt

## Individual probabilities of transition:
# spring birth cohort
n_M1M2_sp <- mu/2 * dt * M1_sp # M1 to M2
n_M2S_sp <- mu/2 * dt * M2_sp # M2 to S
n_SR_sp <- lambda_sp * dt * S_sp # S to R

# summer birth cohort
n_M1M2_sm <- mu/2 * dt * M1_sm # M1 to M2
n_M2S_sm <- mu/2 * dt * M2_sm # M2 to S
n_SR_sm <- lambda_sm * dt * S_sm # S to R

# autumn birth cohort
n_M1M2_au <- mu/2 * dt * M1_au # M1 to M2
n_M2S_au <- mu/2 * dt * M2_au # M2 to S
n_SR_au <- lambda_au * dt * S_au # S to R

# winter birth cohort
n_M1M2_wt <- mu/2 * dt * M1_wt# M1 to M2
n_M2S_wt <- mu/2 * dt * M2_wt # M2 to S
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


# time-dependent contact parameter for each birth cohort
# contacts_sp <- if((time > 30.41*6) &&  (time <= 30.41*12)) 5/4 else if 
#                          ((time > 30.41*12) && (time <= 30.41*18)) 5/4 else if
#                          ((time > 30.41*18) && (time <= 30.41*24)) 4.5/4 else if
#                          ((time > 30.41*24) && (time <= 30.41*30)) 6/4 else if
#                          ((time > 30.41*30) && (time <= 30.41*36)) 8/4 else if
#                          ((time > 30.41*36) && (time <= 30.41*42)) 20/4 else if
#                          ((time > 30.41*42) && (time <= 30.41*48)) 14/4 else if
#                          ((time > 30.41*48) && (time <= 30.41*54)) 29/4 else if
#                          ((time > 30.41*54)) 25/4 else  1
# contacts_sm <- if((time > 30.41*6) &&  (time <= 30.41*12) ~ 3.5/1,
#                          (time > 30.41*12) && (time <= 30.41*18) ~ 3.5/1,
#                          (time > 30.41*18) && (time <= 30.41*24) ~ 7/1,
#                          (time > 30.41*24) && (time <= 30.41*30) ~ 7.5/1,
#                          (time > 30.41*30) && (time <= 30.41*36) ~ 6/1,
#                          (time > 30.41*36) && (time <= 30.41*42) ~ 9/1,
#                          (time > 30.41*42) && (time <= 30.41*48) ~ 11/1,
#                          (time > 30.41*48) && (time <= 30.41*54) ~ 25/4,
#                          (time > 30.41*54) ~ 13/4,
#                          TRUE ~ 1)
# contacts_au <- if((time > 30.41*6) &&  (time <= 30.41*12) ~ 3/4,
#                          (time > 30.41*12) && (time <= 30.41*18) ~ 5/4,
#                          (time > 30.41*18) && (time <= 30.41*24) ~ 4/4,
#                          (time > 30.41*24) && (time <= 30.41*30) ~ 6/4,
#                          (time > 30.41*30) && (time <= 30.41*36) ~ 6/4,
#                          (time > 30.41*36) && (time <= 30.41*42) ~ 9/4,
#                          (time > 30.41*42) && (time <= 30.41*48) ~ 12/4,
#                          (time > 30.41*48) && (time <= 30.41*54) ~ 12/4,
#                          (time > 30.41*54) ~ 11.5/4,
#                          TRUE ~ 1)
# contacts_wt <- if((time > 30.41*6) &&  (time <= 30.41*12) ~ 4/4,
#                          (time > 30.41*12) && (time <= 30.41*18) ~ 4/4,
#                          (time > 30.41*18) && (time <= 30.41*24) ~ 4.5/4,
#                          (time > 30.41*24) && (time <= 30.41*30) ~ 6/4,
#                          (time > 30.41*30) && (time <= 30.41*36) ~ 6/4,
#                          (time > 30.41*36) && (time <= 30.41*42) ~ 7.5/4,
#                          (time > 30.41*42) && (time <= 30.41*48) ~ 11/4,
#                          (time > 30.41*48) && (time <= 30.41*54) ~ 21/4,
#                          (time > 30.41*54) ~ 26.5/4,
#                          TRUE ~ 1)

# Putting it all together into four FOIs
lambda_sp = (summer_comp + spring_comp) * spring_FOI_sp + 
  summer_comp * summer_FOI_sp + 
  (summer_comp + autumn_comp) * autumn_FOI_sp +
  (summer_comp + winter_comp) * winter_FOI_sp #+
 # contact_comp * contacts_sp
lambda_sm = (summer_comp + spring_comp) * spring_FOI_sm + 
  summer_comp * summer_FOI_sm + 
  (summer_comp + autumn_comp) * autumn_FOI_sm + 
  (summer_comp + winter_comp) * winter_FOI_sm #+
  #contact_comp * contacts_sm
lambda_au = (summer_comp + spring_comp) * spring_FOI_au + 
  summer_comp * summer_FOI_au + 
  (summer_comp + autumn_comp) * autumn_FOI_au + 
  (summer_comp + winter_comp) * winter_FOI_au #+
 # contact_comp * contacts_au
lambda_wt = (summer_comp + spring_comp) * spring_FOI_wt + 
  summer_comp * summer_FOI_wt +
  (summer_comp + autumn_comp) * autumn_FOI_wt + 
  (summer_comp + winter_comp) * winter_FOI_wt #+
  #contact_comp * contacts_wt

## Initial states:
# spring cohort
initial(M1_sp) <- M1_sp_ini
initial(M2_sp) <- M2_sp_ini
initial(S_sp) <- S_sp_ini
initial(R_sp) <- R_sp_ini

# summer cohort
initial(M1_sm) <- M1_sm_ini
initial(M2_sm) <- M2_sm_ini
initial(S_sm) <- S_sm_ini
initial(R_sm) <- R_sm_ini

# autumn cohort
initial(M1_au) <- M1_au_ini
initial(M2_au) <- M2_au_ini
initial(S_au) <- S_au_ini
initial(R_au) <- R_au_ini

# winter cohort
initial(M1_wt) <- M1_wt_ini
initial(M2_wt) <- M2_wt_ini
initial(S_wt) <- S_wt_ini
initial(R_wt) <- R_wt_ini

# Total
initial(R_all) <- R_sp_ini + R_sm_ini + R_au_ini + R_wt_ini

## User defined parameters - default in parentheses:
M1_sp_ini <- user (1191 * 0.26) # user(1-2*1e-12)
M2_sp_ini <- user(1e-12)
S_sp_ini <- user(1e-12)
R_sp_ini <- user(1e-12)

M1_sm_ini <- user (1911 * 0.29) # user(1-2*1e-12)
M2_sm_ini <- user(1e-12)
S_sm_ini <- user(1e-12)
R_sm_ini <- user(1e-12)

M1_au_ini <- user (1191 * 0.24) # user(1-2*1e-12)
M2_au_ini <- user(1e-12)
S_au_ini <- user(1e-12)
R_au_ini <- user(1e-12)

M1_wt_ini <- user (1191 * 0.20) # user(1-2*1e-12)
M2_wt_ini <- user(1e-12)
S_wt_ini <- user(1e-12)
R_wt_ini <- user(1e-12)


# transition parameters
mu <- user (1/182.5) # half-life is 6 months
spring_comp <- user(0.00001) # random values to be fitted
summer_comp <- user(0.02002)
autumn_comp <- user(0.00003)
winter_comp <- user(0.00004)
#contact_comp <- user(0.02)