## Definition of the time-step and output as "time"
dt <- user(1)
initial(time) <- 1
update(time) <- (step + 1) * dt


## Core equations for transitions between compartments:
# spring birth cohort
update(M_sp) <- M_sp - n_MS_sp
update(S_sp) <- S_sp + n_MS_sp - n_SR_sp
update(R_sp) <- R_sp + n_SR_sp

# summer birth cohort
update(M_sm) <- M_sm - n_MS_sm
update(S_sm) <- S_sm + n_MS_sm - n_SR_sm
update(R_sm) <- R_sm + n_SR_sm

# autumn birth cohort
update(M_au) <- M_au - n_MS_au
update(S_au) <- S_au + n_MS_au - n_SR_au
update(R_au) <- R_au + n_SR_au

# winter birth cohort
update(M_wt) <- M_wt - n_MS_wt
update(S_wt) <- S_wt + n_MS_wt - n_SR_wt
update(R_wt) <- R_wt + n_SR_wt

# Total seroprevalence
update(R_all) <- R_sp + R_sm + R_au + R_wt

# Get FOI in output
update(lambda_spring) <- lambda_sp # FOI for spring cohort
update(lambda_summer) <- lambda_sm
update(lambda_autumn) <- lambda_au
update(lambda_winter) <- lambda_wt

update(spring_component) <- spring_comp # FOI in the spring
update(summer_component) <- summer_comp
update(autumn_component) <- autumn_comp
update(winter_component) <- winter_comp

## Individual probabilities of transition:
# spring birth cohort
n_MS_sp <- mu * dt * M_sp # M to S
n_SR_sp <- lambda_sp * dt * S_sp # S to R

# summer birth cohort
n_MS_sm <- mu * dt * M_sm # M to S
n_SR_sm <- lambda_sm * dt * S_sm # S to R

# autumn birth cohort
n_MS_au <- mu * dt * M_au # M to S
n_SR_au <- lambda_au * dt * S_au # S to R

# winter birth cohort
n_MS_wt <- mu * dt * M_wt # M to S
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
initial(M_sp) <- M_sp_ini
initial(S_sp) <- S_sp_ini
initial(R_sp) <- R_sp_ini

# summer cohort
initial(M_sm) <- M_sm_ini
initial(S_sm) <- S_sm_ini
initial(R_sm) <- R_sm_ini

# autumn cohort
initial(M_au) <- M_au_ini
initial(S_au) <- S_au_ini
initial(R_au) <- R_au_ini

# winter cohort
initial(M_wt) <- M_wt_ini
initial(S_wt) <- S_wt_ini
initial(R_wt) <- R_wt_ini

# Total
initial(R_all) <- R_sp_ini + R_sm_ini + R_au_ini + R_wt_ini

# FOI for the output
initial(lambda_spring) <- 0.0001
initial(lambda_summer) <- 0.0001
initial(lambda_autumn) <- 0.0001
initial(lambda_winter) <- 0.0001

initial(spring_component) <- 0.00001
initial(summer_component) <- 0.00001
initial(autumn_component) <- 0.00001
initial(winter_component) <- 0.00001

## User defined parameters - default in parentheses:
M_sp_ini <- user(1-2*1e-12) # user (682 * 0.26) # 
S_sp_ini <- user(1e-12)
R_sp_ini <- user(1e-12)

M_sm_ini <- user(1-2*1e-12) # user (682 * 0.29) # 
S_sm_ini <- user(1e-12)
R_sm_ini <- user(1e-12)

M_au_ini <- user(1-2*1e-12) # user (682 * 0.24) # 
S_au_ini <- user(1e-12)
R_au_ini <- user(1e-12)

M_wt_ini <- user(1-2*1e-12) # user (682 * 0.20) # 
S_wt_ini <- user(1e-12)
R_wt_ini <- user(1e-12)


# transition parameters
mu <- user (1/182.5) # half-life is 6 months
spring_comp <- user(0.00001) # random values to be fitted
summer_comp <- user(0.02002)
autumn_comp <- user(0.00003)
winter_comp <- user(0.00004)
#contact_comp <- user(0.02)