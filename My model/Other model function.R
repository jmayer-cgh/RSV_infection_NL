model <- function(theta, age, inits, data) {
  
  catalytic <- function(age, state, param, data) {
    
    # FOI / seroconversion rate
    # Define booleans to know which season the cohort is in
    spring_FOI_sp = 0
    summer_FOI_sp = 0
    autumn_FOI_sp = 0
    winter_FOI_sp = 0
    spring_FOI_sm = 0
    summer_FOI_sm = 0
    autumn_FOI_sm = 0
    winter_FOI_sm = 0
    spring_FOI_au = 0
    summer_FOI_au = 0
    autumn_FOI_au = 0
    winter_FOI_au = 0
    spring_FOI_wt = 0
    summer_FOI_wt = 0
    autumn_FOI_wt = 0
    winter_FOI_wt = 0
    
    # Baseline contact parameter for each birth cohort
    contacts_sp = 1
    contacts_sm = 1
    contacts_au = 1
    contacts_wt = 1
    
    # Baseline percentage of daycare attendance 
    daycare = 0
    if (age >= 30.41*9){ # assume that all children who go to day-care start at 9 months
      daycare = 0.4075758
    }
    
    # Get seasonal FOI
    if ( (age <= 30.41*3)                                      # FOI of the season in which the children were born
         || ((age >= 365) & (age <= (365+30.41*3))) 
         || ((age >= 2*365) & (age <= (2*365+30.41*3))) 
         || ((age >= 3*365) & (age <= (3*365+30.41*3))) 
         || ((age >= 4*365) & (age <= (4*365+30.41*3))) 
         || ((age >= 5*365) & (age <= (5*365+30.41*3))) ){
      spring_FOI_sp = 1
      summer_FOI_sm = 1
      autumn_FOI_au = 1
      winter_FOI_wt = 1
    } 
    
    if ( (age>30.41*3 & age <=30.41*6 )                       # FOI of the season after the children were born
         || ((age > 365 + 30.41*3) & (age <= (365+30.41*6))) 
         || ((age > 2*365 + 30.41*3) & (age <= (2*365+30.41*6))) 
         || ((age > 3*365 + 30.41*3) & (age <= (3*365+30.41*6))) 
         || ((age > 4*365 + 30.41*3) & (age <= (4*365+30.41*6))) 
         || ((age > 5*365 + 30.41*3) & (age <= (5*365+30.41*6))) ){
      summer_FOI_sp = 1
      autumn_FOI_sm = 1
      winter_FOI_au = 1
      spring_FOI_wt = 1
    } 
    
    if ( (age > 30.41*6 & age <= 30.41*9 )                     # FOI 2 seasons after birth
         || ((age > 365+30.41*6) & (age <= (365+30.41*9))) 
         || ((age > 2*365 + 30.41*6) & (age <= (2*365+30.41*9))) 
         || ((age > 3*365 + 30.41*6) & (age <= (3*365+30.41*9))) 
         || ((age > 4*365 + 30.41*6) & (age <= (4*365+30.41*9))) 
         || ((age > 5*365 + 30.41*6) & (age <= (5*365+30.41*9))) ){
      autumn_FOI_sp = 1
      winter_FOI_sm = 1
      spring_FOI_au = 1
      summer_FOI_wt = 1
    } 
    
    if ( (age > 30.41*9 & age <= 30.41*12 )                     # FOI 3 seasons after birth
         || ((age > 365+30.41*9) & (age <= (365+30.41*12))) 
         || ((age > 2*365 + 30.41*9) & (age <= (2*365+30.41*12))) 
         || ((age > 3*365 + 30.41*9) & (age <= (3*365+30.41*12))) 
         || ((age > 4*365 + 30.41*9) & (age <= (4*365+30.41*12))) 
         || ((age > 5*365 + 30.41*9) & (age <= (5*365+30.41*12))) ){
      winter_FOI_sp = 1
      spring_FOI_sm = 1
      summer_FOI_au = 1
      autumn_FOI_wt = 1
    } 
    
    # Update contact and daycare parameters 
    '#the median number of contacts for children aged 9-12 months is x/y*(median 
    number of contacts of children aged 0-6 months)#'
    if ((age > 30.41*6) & (age <= 30.41*12)){ 
      contacts_sp = 5/4
      contacts_sm = 3.5/1
      contacts_au = 3/4
      contacts_wt = 4/4
    } 
    if ((age > 30.41*12) & (age <= 30.41*18)){
      contacts_sp = 5/4
      contacts_sm = 3.5/1
      contacts_au = 5/4
      contacts_wt = 4/4
    }
    if ((age > 30.41*18) & (age <= 30.41*24)){
      contacts_sp = 4.5/4
      contacts_sm = 7/1
      contacts_au = 4/4
      contacts_wt = 4.5/4
    }
    if ((age > 30.41*24) & (age <= 30.41*30)){
      contacts_sp = 6/4
      contacts_sm = 7.5/1
      contacts_au = 6/4
      contacts_wt = 6/4
    }
    if ((age > 30.41*30) & (age <= 30.41*36)){
      contacts_sp = 8/4
      contacts_sm = 6/1
      contacts_au = 6/4
      contacts_wt = 6/4
    }
    if ((age > 30.41*36) & (age <= 30.41*42)){
      contacts_sp = 20/4
      contacts_sm = 9/1
      contacts_au = 9/4
      contacts_wt = 7.5/4
    }
    if ((age > 30.41*42) & (age <= 30.41*48)){
      contacts_sp = 14/4
      contacts_sm = 11/1
      contacts_au = 12/4
      contacts_wt = 11/4
    }
    if ((age > 30.41*48) & (age <= 30.41*54)){
      contacts_sp = 29/4
      contacts_sm = 25/1
      contacts_au = 12/4  # NB: no children born in autumn are aged 48-54 months
      contacts_wt = 21/4
    }
    if (age > 30.41*54){
      contacts_sp = 25/4
      contacts_sm = 13/1
      contacts_au = 11.5/4
      contacts_wt = 26.5/4
    }
    
    
    # FOI for each birth cohort
    # participants not attending daycare
    # FOI for each birth cohort
    lambda_sp = (param[["M"]]+param[["P"]])*spring_FOI_sp + param[["M"]]*summer_FOI_sp + 
      (param[["M"]]+param[["A"]])*autumn_FOI_sp + (param[["M"]]+param[["W"]])*winter_FOI_sp +
      param[["C"]]*contacts_sp
    lambda_sm = (param[["M"]]+param[["P"]])*spring_FOI_sm + param[["M"]]*summer_FOI_sm + 
      (param[["M"]]+param[["A"]])*autumn_FOI_sm + (param[["M"]]+param[["W"]])*winter_FOI_sm +
      param[["C"]]*contacts_sm
    lambda_au = (param[["M"]]+param[["P"]])*spring_FOI_au + param[["M"]]*summer_FOI_au + 
      (param[["M"]]+param[["A"]])*autumn_FOI_au + (param[["M"]]+param[["W"]])*winter_FOI_au +
      param[["C"]]*contacts_au
    lambda_wt = (param[["M"]]+param[["P"]])*spring_FOI_wt + param[["M"]]*summer_FOI_wt +
      (param[["M"]]+param[["A"]])*autumn_FOI_wt + (param[["M"]]+param[["W"]])*winter_FOI_wt +
      param[["C"]]*contacts_wt
    
    # waning maternal immunity, same for all children
    mu = param[["B"]] 
    
    # states 
    # proportion with maternal immunity
    M_sp = exp(state[1]) # born in spring
    M_sm = exp(state[2]) # born in summer
    M_au = exp(state[3]) # born in autumn
    M_wt = exp(state[4]) # born in winter
    # Susceptible
    S_sp = exp(state[5]) # susceptible born in spring 
    S_sm = exp(state[6]) # susceptible born in summer
    S_au = exp(state[7]) # susceptible born in autumn
    S_wt = exp(state[8]) # susceptible born in winter 
    #Seroconverted
    Z_sp = exp(state[9]) # seroconverted after infection born in spring attending day-care
    Z_sm = exp(state[10]) # seroconverted after infection born in summer attending day-care
    Z_au = exp(state[11]) # seroconverted after infection born in autumn attending day-care
    Z_wt = exp(state[12]) # seroconverted after infection born in winter attending day-care
    
    
    # changes in states
    dM_sp = -mu*M_sp
    dM_sm = -mu*M_sm
    dM_au = -mu*M_au
    dM_wt = -mu*M_wt
    dS_sp = + mu*M_sp - (1-daycare)*lambda_sp*S_sp - daycare*param[["D"]]*lambda_sp*S_sp
    dS_sm = + mu*M_sm - (1-daycare)*lambda_sm*S_sm - daycare*param[["D"]]*lambda_sm*S_sm
    dS_au = + mu*M_au - (1-daycare)*lambda_au*S_au - daycare*param[["D"]]*lambda_au*S_au 
    dS_wt = + mu*M_wt - (1-daycare)*lambda_wt*S_wt - daycare*param[["D"]]*lambda_wt*S_wt
    dZ_sp = + (1-daycare)*lambda_sp*S_sp + daycare*param[["D"]]*lambda_sp*S_sp
    dZ_sm = + (1-daycare)*lambda_sm*S_sm + daycare*param[["D"]]*lambda_sm*S_sm
    dZ_au = + (1-daycare)*lambda_au*S_au + daycare*param[["D"]]*lambda_au*S_au
    dZ_wt = + (1-daycare)*lambda_wt*S_wt + daycare*param[["D"]]*lambda_wt*S_wt
    
    
    return(list(c(dM_sp/M_sp,dM_sm/M_sm, dM_au/M_au,dM_wt/M_wt,
                  dS_sp/S_sp, dS_sm/S_sm, dS_au/S_au,dS_wt/S_wt,
                  dZ_sp/Z_sp, dZ_sm/Z_sm, dZ_au/Z_au, dZ_wt/Z_wt),
                lambda_sp=lambda_sp, lambda_sm = lambda_sm, 
                lambda_au = lambda_au, lambda_wt = lambda_wt,
                mu=mu))
    
    
  }
  
  traj <- data.frame(ode(y=c(M_sp=log(inits[["M_sp"]]),
                             M_sm=log(inits[["M_sm"]]),
                             M_au=log(inits[["M_au"]]),
                             M_wt=log(inits[["M_wt"]]),
                             S_sp=log(inits[["S_sp"]]),
                             S_sm=log(inits[["S_sm"]]),
                             S_au=log(inits[["S_au"]]),
                             S_wt=log(inits[["S_wt"]]),
                             Z_sp=log(inits[["Z_sp"]]),
                             Z_sm=log(inits[["Z_sm"]]),
                             Z_au=log(inits[["Z_au"]]),
                             Z_wt=log(inits[["Z_wt"]])),
                         times=age, 
                         func=catalytic, 
                         parms=theta, 
                         data = data,
                         method="lsoda",
                         verbose=F))
  
  traj$conv_spring <- exp(traj$Z_sp) # cumulative seroconversion in spring cohort (=observed state)
  traj$inc_spring <- c(inits[["Z_sp"]], diff(exp(traj$Z_sp))) # incident seroconversion
  traj$conv_summer <- exp(traj$Z_sm) # cumulative seroconversion in summer cohort (=observed state)
  traj$inc_summer <- c(inits[["Z_sm"]], diff(exp(traj$Z_sm))) # incident seroconversion
  traj$conv_autumn <- exp(traj$Z_au) # cumulative seroconversion in autumn cohort (=observed state)
  traj$inc_autumn <- c(inits[["Z_au"]], diff(exp(traj$Z_au))) # incident seroconversion
  traj$conv_winter <- exp(traj$Z_wt) # cumulative seroconversion in winter cohort (=observed state)
  traj$inc_winter <- c(inits[["Z_wt"]], diff(exp(traj$Z_wt))) # incident seroconversion
  
  traj$Z_all = 0.26*traj$Z_sp + 0.29*traj$Z_sm + 0.24*traj$Z_au + 0.20*traj$Z_wt
  traj$conv <- exp(traj$Z_all)
  '#
  traj$conv_spring_d <- exp(traj$Z_sp_d) # cumulative seroconversion in spring cohort (=observed state)
  traj$inc_spring_d <- c(inits[["Z_sp_d"]], diff(exp(traj$Z_sp_d))) # incident seroconversion
  traj$conv_summer_d <- exp(traj$Z_sm_d) # cumulative seroconversion in summer cohort (=observed state)
  traj$inc_summer_d <- c(inits[["Z_sm_d"]], diff(exp(traj$Z_sm_d))) # incident seroconversion
  traj$conv_autumn_d <- exp(traj$Z_au_d) # cumulative seroconversion in autumn cohort (=observed state)
  traj$inc_autumn_d <- c(inits[["Z_au_d"]], diff(exp(traj$Z_au_d))) # incident seroconversion
  traj$conv_winter_d <- exp(traj$Z_wt_d) # cumulative seroconversion in winter cohort (=observed state)
  traj$inc_winter_d <- c(inits[["Z_wt_d"]], diff(exp(traj$Z_wt_d))) # incident seroconversion
  
  traj$Z_all_d = 0.26*traj$Z_sp_d + 0.29*traj$Z_sm_d + 0.24*traj$Z_au_d + 0.20*traj$Z_wt_d
  traj$conv_d <- exp(traj$Z_all_d)
  
  traj$conv_spring_n <- exp(traj$Z_sp_n) # cumulative seroconversion in spring cohort (=observed state)
  traj$inc_spring_n <- c(inits[["Z_sp_n"]], diff(exp(traj$Z_sp_n))) # incident seroconversion
  traj$conv_summer_n <- exp(traj$Z_sm_n) # cumulative seroconversion in summer cohort (=observed state)
  traj$inc_summer_n <- c(inits[["Z_sm_n"]], diff(exp(traj$Z_sm_n))) # incident seroconversion
  traj$conv_autumn_n <- exp(traj$Z_au_n) # cumulative seroconversion in autumn cohort (=observed state)
  traj$inc_autumn_n <- c(inits[["Z_au_n"]], diff(exp(traj$Z_au_n))) # incident seroconversion
  traj$conv_winter_n <- exp(traj$Z_wt_n) # cumulative seroconversion in winter cohort (=observed state)
  traj$inc_winter_n <- c(inits[["Z_wt_n"]], diff(exp(traj$Z_wt_n))) # incident seroconversion
  
  traj$Z_all_n = 0.26*traj$Z_sp_n + 0.29*traj$Z_sm_n + 0.24*traj$Z_au_n+ 0.20*traj$Z_wt_n
  traj$conv_n <- exp(traj$Z_all_n)
  
  daycare = 0.4075758
  
  traj$Z_spring_tot = traj$Z_sp_d*daycare + traj$Z_sp_n*(1-daycare)
  traj$Z_summer_tot = traj$Z_sm_d*daycare + traj$Z_sm_n*(1-daycare)
  traj$Z_autumn_tot = traj$Z_au_d*daycare + traj$Z_au_n*(1-daycare)
  traj$Z_winter_tot = traj$Z_wt_d*daycare + traj$Z_wt_n*(1-daycare)
  traj$conv_sp_tot = exp(traj$Z_spring_tot)
  traj$conv_sm_tot = exp(traj$Z_summer_tot)
  traj$conv_au_tot = exp(traj$Z_autumn_tot)
  traj$conv_wt_tot = exp(traj$Z_winter_tot)
  
  traj$Z_all_tot = 0.26*traj$Z_spring_tot + 0.29*traj$Z_summer_tot + 0.24*traj$Z_autumn_tot+ 0.20*traj$Z_winter_tot
  traj$conv_tot <- exp(traj$Z_all_tot)
  #'
  
  return(traj)
  
}

model <- function(theta, age, inits, data) {
  
  catalytic <- function(age, state, param, data) {
    
    # FOI / seroconversion rate
    # Define booleans to know which season the cohort is in
    spring_FOI_sp = 0
    summer_FOI_sp = 0
    autumn_FOI_sp = 0
    winter_FOI_sp = 0
    spring_FOI_sm = 0
    summer_FOI_sm = 0
    autumn_FOI_sm = 0
    winter_FOI_sm = 0
    spring_FOI_au = 0
    summer_FOI_au = 0
    autumn_FOI_au = 0
    winter_FOI_au = 0
    spring_FOI_wt = 0
    summer_FOI_wt = 0
    autumn_FOI_wt = 0
    winter_FOI_wt = 0
    
    # Baseline contact parameter for each birth cohort
    contacts_sp = 1
    contacts_sm = 1
    contacts_au = 1
    contacts_wt = 1
    
    # Baseline percentage of daycare attendance 
    daycare = 0
    if (age >= 30.41*9){ # assume that all children who go to day-care start at 9 months
      daycare = 0.4075758
    }
    
    # Get seasonal FOI
    if ( (age <= 30.41*3)                                      # FOI of the season in which the children were born
         || ((age >= 365) & (age <= (365+30.41*3))) 
         || ((age >= 2*365) & (age <= (2*365+30.41*3))) 
         || ((age >= 3*365) & (age <= (3*365+30.41*3))) 
         || ((age >= 4*365) & (age <= (4*365+30.41*3))) 
         || ((age >= 5*365) & (age <= (5*365+30.41*3))) ){
      spring_FOI_sp = 1
      summer_FOI_sm = 1
      autumn_FOI_au = 1
      winter_FOI_wt = 1
    } 
    
    if ( (age>30.41*3 & age <=30.41*6 )                       # FOI of the season after the children were born
         || ((age > 365 + 30.41*3) & (age <= (365+30.41*6))) 
         || ((age > 2*365 + 30.41*3) & (age <= (2*365+30.41*6))) 
         || ((age > 3*365 + 30.41*3) & (age <= (3*365+30.41*6))) 
         || ((age > 4*365 + 30.41*3) & (age <= (4*365+30.41*6))) 
         || ((age > 5*365 + 30.41*3) & (age <= (5*365+30.41*6))) ){
      summer_FOI_sp = 1
      autumn_FOI_sm = 1
      winter_FOI_au = 1
      spring_FOI_wt = 1
    } 
    
    if ( (age > 30.41*6 & age <= 30.41*9 )                     # FOI 2 seasons after birth
         || ((age > 365+30.41*6) & (age <= (365+30.41*9))) 
         || ((age > 2*365 + 30.41*6) & (age <= (2*365+30.41*9))) 
         || ((age > 3*365 + 30.41*6) & (age <= (3*365+30.41*9))) 
         || ((age > 4*365 + 30.41*6) & (age <= (4*365+30.41*9))) 
         || ((age > 5*365 + 30.41*6) & (age <= (5*365+30.41*9))) ){
      autumn_FOI_sp = 1
      winter_FOI_sm = 1
      spring_FOI_au = 1
      summer_FOI_wt = 1
    } 
    
    if ( (age > 30.41*9 & age <= 30.41*12 )                     # FOI 3 seasons after birth
         || ((age > 365+30.41*9) & (age <= (365+30.41*12))) 
         || ((age > 2*365 + 30.41*9) & (age <= (2*365+30.41*12))) 
         || ((age > 3*365 + 30.41*9) & (age <= (3*365+30.41*12))) 
         || ((age > 4*365 + 30.41*9) & (age <= (4*365+30.41*12))) 
         || ((age > 5*365 + 30.41*9) & (age <= (5*365+30.41*12))) ){
      winter_FOI_sp = 1
      spring_FOI_sm = 1
      summer_FOI_au = 1
      autumn_FOI_wt = 1
    } 
    
    # Update contact and daycare parameters 
    '#the median number of contacts for children aged 9-12 months is x/y*(median 
    number of contacts of children aged 0-6 months)#'
    if ((age > 30.41*6) & (age <= 30.41*12)){ 
      contacts_sp = 5/4
      contacts_sm = 3.5/1
      contacts_au = 3/4
      contacts_wt = 4/4
    } 
    if ((age > 30.41*12) & (age <= 30.41*18)){
      contacts_sp = 5/4
      contacts_sm = 3.5/1
      contacts_au = 5/4
      contacts_wt = 4/4
    }
    if ((age > 30.41*18) & (age <= 30.41*24)){
      contacts_sp = 4.5/4
      contacts_sm = 7/1
      contacts_au = 4/4
      contacts_wt = 4.5/4
    }
    if ((age > 30.41*24) & (age <= 30.41*30)){
      contacts_sp = 6/4
      contacts_sm = 7.5/1
      contacts_au = 6/4
      contacts_wt = 6/4
    }
    if ((age > 30.41*30) & (age <= 30.41*36)){
      contacts_sp = 8/4
      contacts_sm = 6/1
      contacts_au = 6/4
      contacts_wt = 6/4
    }
    if ((age > 30.41*36) & (age <= 30.41*42)){
      contacts_sp = 20/4
      contacts_sm = 9/1
      contacts_au = 9/4
      contacts_wt = 7.5/4
    }
    if ((age > 30.41*42) & (age <= 30.41*48)){
      contacts_sp = 14/4
      contacts_sm = 11/1
      contacts_au = 12/4
      contacts_wt = 11/4
    }
    if ((age > 30.41*48) & (age <= 30.41*54)){
      contacts_sp = 29/4
      contacts_sm = 25/1
      contacts_au = 12/4  # NB: no children born in autumn are aged 48-54 months
      contacts_wt = 21/4
    }
    if (age > 30.41*54){
      contacts_sp = 25/4
      contacts_sm = 13/1
      contacts_au = 11.5/4
      contacts_wt = 26.5/4
    }
    
    
    # FOI for each birth cohort
    # participants not attending daycare
    # FOI for each birth cohort
    lambda_sp = (param[["M"]]+param[["P"]])*spring_FOI_sp + param[["M"]]*summer_FOI_sp + 
      (param[["M"]]+param[["A"]])*autumn_FOI_sp + (param[["M"]]+param[["W"]])*winter_FOI_sp +
      param[["C"]]*contacts_sp
    lambda_sm = (param[["M"]]+param[["P"]])*spring_FOI_sm + param[["M"]]*summer_FOI_sm + 
      (param[["M"]]+param[["A"]])*autumn_FOI_sm + (param[["M"]]+param[["W"]])*winter_FOI_sm +
      param[["C"]]*contacts_sm
    lambda_au = (param[["M"]]+param[["P"]])*spring_FOI_au + param[["M"]]*summer_FOI_au + 
      (param[["M"]]+param[["A"]])*autumn_FOI_au + (param[["M"]]+param[["W"]])*winter_FOI_au +
      param[["C"]]*contacts_au
    lambda_wt = (param[["M"]]+param[["P"]])*spring_FOI_wt + param[["M"]]*summer_FOI_wt +
      (param[["M"]]+param[["A"]])*autumn_FOI_wt + (param[["M"]]+param[["W"]])*winter_FOI_wt +
      param[["C"]]*contacts_wt
    
    # waning maternal immunity, same for all children
    mu = param[["B"]] 
    
    # states 
    # proportion with maternal immunity
    M_sp = exp(state[1]) # born in spring
    M_sm = exp(state[2]) # born in summer
    M_au = exp(state[3]) # born in autumn
    M_wt = exp(state[4]) # born in winter
    # Susceptible
    S_sp = exp(state[5]) # susceptible born in spring 
    S_sm = exp(state[6]) # susceptible born in summer
    S_au = exp(state[7]) # susceptible born in autumn
    S_wt = exp(state[8]) # susceptible born in winter 
    #Seroconverted
    Z_sp_d = exp(state[9]) # seroconverted after infection born in spring attending day-care
    Z_sm_d = exp(state[10]) # seroconverted after infection born in summer attending day-care
    Z_au_d = exp(state[11]) # seroconverted after infection born in autumn attending day-care
    Z_wt_d = exp(state[12]) # seroconverted after infection born in winter attending day-care
    Z_sp_n = exp(state[13]) # seroconverted after infection born in spring not attending day-care
    Z_sm_n = exp(state[14]) # seroconverted after infection born in summer not attending day-care
    Z_au_n = exp(state[15]) # seroconverted after infection born in autumn not attending day-care
    Z_wt_n = exp(state[16]) # seroconverted after infection born in winter not attending day-care
    
    # changes in states
    dM_sp = -mu*M_sp
    dM_sm = -mu*M_sm
    dM_au = -mu*M_au
    dM_wt = -mu*M_wt
    dS_sp = + mu*M_sp - (1-daycare)*lambda_sp*S_sp - daycare*param[["D"]]*lambda_sp*S_sp
    dS_sm = + mu*M_sm - (1-daycare)*lambda_sm*S_sm - daycare*param[["D"]]*lambda_sm*S_sm
    dS_au = + mu*M_au - (1-daycare)*lambda_au*S_au - daycare*param[["D"]]*lambda_au*S_au 
    dS_wt = + mu*M_wt - (1-daycare)*lambda_wt*S_wt - daycare*param[["D"]]*lambda_wt*S_wt
    dZ_sp_d = + daycare*param[["D"]]*lambda_sp*S_sp
    dZ_sm_d = + daycare*param[["D"]]*lambda_sm*S_sm
    dZ_au_d = + daycare*param[["D"]]*lambda_au*S_au
    dZ_wt_d = + daycare*param[["D"]]*lambda_wt*S_wt
    dZ_sp_n = + (1-daycare)*lambda_sp*S_sp
    dZ_sm_n = + (1-daycare)*lambda_sm*S_sm
    dZ_au_n = + (1-daycare)*lambda_au*S_au
    dZ_wt_n = + (1-daycare)*lambda_wt*S_wt 
    
    
    return(list(c(dM_sp/M_sp,dM_sm/M_sm, dM_au/M_au,dM_wt/M_wt,
                  dS_sp/S_sp, dS_sm/S_sm, dS_au/S_au,dS_wt/S_wt,
                  dZ_sp_d/Z_sp_d, dZ_sm_d/Z_sm_d, dZ_au_d/Z_au_d, dZ_wt_d/Z_wt_d,
                  dZ_sp_n/Z_sp_n, dZ_sm_n/Z_sm_n, dZ_au_n/Z_au_n, dZ_wt_n/Z_wt_n),
                lambda_sp=lambda_sp, lambda_sm = lambda_sm, 
                lambda_au = lambda_au, lambda_wt = lambda_wt,
                mu=mu))
    
    
  }
  
  traj <- data.frame(ode(y=c(M_sp=log(inits[["M_sp"]]),
                             M_sm=log(inits[["M_sm"]]),
                             M_au=log(inits[["M_au"]]),
                             M_wt=log(inits[["M_wt"]]),
                             S_sp=log(inits[["S_sp"]]),
                             S_sm=log(inits[["S_sm"]]),
                             S_au=log(inits[["S_au"]]),
                             S_wt=log(inits[["S_wt"]]),
                             Z_sp_d=log(inits[["Z_sp_d"]]),
                             Z_sm_d=log(inits[["Z_sm_d"]]),
                             Z_au_d=log(inits[["Z_au_d"]]),
                             Z_wt_d=log(inits[["Z_wt_d"]]),
                             Z_sp_n=log(inits[["Z_sp_n"]]),
                             Z_sm_n=log(inits[["Z_sm_n"]]),
                             Z_au_n=log(inits[["Z_au_n"]]),
                             Z_wt_n=log(inits[["Z_wt_n"]])),
                         times=age, 
                         func=catalytic, 
                         parms=theta, 
                         data = data,
                         method="lsoda",
                         verbose=F))
  
  traj$conv_spring_d <- exp(traj$Z_sp_d) # cumulative seroconversion in spring cohort (=observed state)
  traj$inc_spring_d <- c(inits[["Z_sp_d"]], diff(exp(traj$Z_sp_d))) # incident seroconversion
  traj$conv_summer_d <- exp(traj$Z_sm_d) # cumulative seroconversion in summer cohort (=observed state)
  traj$inc_summer_d <- c(inits[["Z_sm_d"]], diff(exp(traj$Z_sm_d))) # incident seroconversion
  traj$conv_autumn_d <- exp(traj$Z_au_d) # cumulative seroconversion in autumn cohort (=observed state)
  traj$inc_autumn_d <- c(inits[["Z_au_d"]], diff(exp(traj$Z_au_d))) # incident seroconversion
  traj$conv_winter_d <- exp(traj$Z_wt_d) # cumulative seroconversion in winter cohort (=observed state)
  traj$inc_winter_d <- c(inits[["Z_wt_d"]], diff(exp(traj$Z_wt_d))) # incident seroconversion
  
  traj$Z_all_d = 0.26*traj$Z_sp_d + 0.29*traj$Z_sm_d + 0.24*traj$Z_au_d + 0.20*traj$Z_wt_d
  traj$conv_d <- exp(traj$Z_all_d)
  
  traj$conv_spring_n <- exp(traj$Z_sp_n) # cumulative seroconversion in spring cohort (=observed state)
  traj$inc_spring_n <- c(inits[["Z_sp_n"]], diff(exp(traj$Z_sp_n))) # incident seroconversion
  traj$conv_summer_n <- exp(traj$Z_sm_n) # cumulative seroconversion in summer cohort (=observed state)
  traj$inc_summer_n <- c(inits[["Z_sm_n"]], diff(exp(traj$Z_sm_n))) # incident seroconversion
  traj$conv_autumn_n <- exp(traj$Z_au_n) # cumulative seroconversion in autumn cohort (=observed state)
  traj$inc_autumn_n <- c(inits[["Z_au_n"]], diff(exp(traj$Z_au_n))) # incident seroconversion
  traj$conv_winter_n <- exp(traj$Z_wt_n) # cumulative seroconversion in winter cohort (=observed state)
  traj$inc_winter_n <- c(inits[["Z_wt_n"]], diff(exp(traj$Z_wt_n))) # incident seroconversion
  
  traj$Z_all_n = 0.26*traj$Z_sp_n + 0.29*traj$Z_sm_n + 0.24*traj$Z_au_n+ 0.20*traj$Z_wt_n
  traj$conv_n <- exp(traj$Z_all_n)
  
  daycare = 0.4075758
  '#for (i in seq_along(traj$conv_spring_d)){
    if (traj$conv_spring_d[i] < 2e-12){
      daycare = 0.4075758
    }
 }#'
  traj$Z_spring_tot = traj$Z_sp_d*daycare + traj$Z_sp_n*(1-daycare)
  traj$Z_summer_tot = traj$Z_sm_d*daycare + traj$Z_sm_n*(1-daycare)
  traj$Z_autumn_tot = traj$Z_au_d*daycare + traj$Z_au_n*(1-daycare)
  traj$Z_winter_tot = traj$Z_wt_d*daycare + traj$Z_wt_n*(1-daycare)
  
  traj$conv_sp_tot = exp(traj$Z_spring_tot)
  traj$conv_sm_tot = exp(traj$Z_summer_tot)
  traj$conv_au_tot = exp(traj$Z_autumn_tot)
  traj$conv_wt_tot = exp(traj$Z_winter_tot)
  
  traj$Z_all_tot = 0.26*traj$Z_spring_tot + 0.29*traj$Z_summer_tot + 0.24*traj$Z_autumn_tot+ 0.20*traj$Z_winter_tot
  traj$conv_tot <- exp(traj$Z_all_tot)
  traj$daycare = daycare
  
  return(traj)
  
}