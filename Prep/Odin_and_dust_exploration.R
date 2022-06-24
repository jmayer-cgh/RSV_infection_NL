#odin package
lorenz <- odin::odin({
  ## Derivatives
  deriv(y1) <- sigma * (y2 - y1)
  deriv(y2) <- R * y1 - y2 - y1 * y3
  deriv(y3) <- -b * y3 + y1 * y2
  
  ## Initial conditions
  initial(y1) <- 10.0
  initial(y2) <- 1.0
  initial(y3) <- 1.0
  
  ## parameters
  sigma <- 10.0
  R     <- 28.0
  b     <-  8.0 / 3.0
})

mod <- lorenz()
t <- seq(0, 100, length.out = 50000)
y <- mod$run(t)

a <- if (2<10) 5 else 7
b[] <- runif (10,20)

#dust package-----------------------------------------------------------------------------------------------------------------------------------------------------
#a simple example
walk <- dust::dust_example("walk")
model <- walk$new(list(sd = 1), 0, 20)
#at each time step we move our position with a draw from a normal distribution with mean 0 and standard deviation 1, step is 0 and 20 particles
model
model$state() #initial model
model$run(100) #run the model on 100 steps

model <- walk$new(list(sd = 1), 0, 20000) #now 20000 particles instead of 20, all are run one after the other
invisible(model$run(100))
hist(model$state(), freq = FALSE, las = 1, col = "steelblue2", main = "",
     ylim = c(0., 0.04), xlab = "State")
curve(dnorm(x, 0, 10), col = "orange", add = TRUE, lwd = 2)

model <- walk$new(list(sd = 1), 0, 20, n_threads = 2) #n_threads means we're running the particles in parallel on 2 threads
model$run(100)

#a more complicated example
sir <- dust::dust_example("sir") #SIR model (Susceptible - Infected - Recovered)
coef(sir)
model <- sir$new(list(), 0, 20)
model$info()
model$state() #five states per particle, but WHY?
model$set_index(2L) #will only return the state of the I variable when run, WHY? is 1L = S, 2L= I, 3L = R?
model$run(10)
model$set_index(c(I = 2L)) #same as above only rows will be labelled
model$run(30)
model$state() #still gives you all states

model$state(c(S = 1L, R = 3L)) #S and R variables
model$run(20)

model <- sir$new(list(), 0, 200)
model$set_index(2L)
steps <- seq(0, 600, by = 5) #run over steps 0 to 600 and record all the states
state <- model$simulate(steps) #The output here is a 1 x 200 x 121 matrix (n state x n particles x n steps)
traces <- t(drop(state)) #drop first dimension and transpose for ease of plotting
time <- steps / 4  #Plotting this over time (with 4 steps per day)
matplot(time, traces, type = "l", lty = 1, col = "#00000022",
        xlab = "Time", ylab = "Number infected (I)")
lines(time, rowMeans(traces), col = "red", lwd = 2)

#OTHER METHODS
#Reodering particles
model <- walk$new(list(sd = 1), 0, 20) #create a model
model$run(1)
index <- order(model$state()) #particles are now in increasing state
index
model$reorder(index) #particles will now be in the new order
model$state()
#Set particle state, three mutable things; pars, state and step
model <- sir$new(list(), 0, 20)
model$update_state(state = c(1000, 1, 0, 0, 0)) #reset step
model$state()
steps <- seq(0, 600, by = 5) #plot the epidemic
state <- model$simulate(steps)
time <- steps / 4
matplot(time, t(state[2, , ]), type = "l", lty = 1, col = "#00000022",
        xlab = "Time", ylab = "Number infected (I)")

step0 <- sample(0:30, 20, replace = TRUE) #10 individuals in the infected state but we do so with a number of possible starting points sampled from between step 0 and step 30:
model <- sir$new(list(), 0, 20)
model$update_state(state = c(1000, 10, 0, 0, 0), step = step0)
model$step() #model should be at step 29
model$state() #initial states are a mix of states, despite having been seeded with the same 10 infected individuals, as different particles have run for different numbers of steps to reach this common time point
steps <- seq(max(step0), 600, by = 5)
state <- model$simulate(steps)
time <- steps / 4
matplot(time, t(state[2, , ]), type = "l", lty = 1, col = "#00000022",
        xlab = "Time", ylab = "Number infected (I)")

I0 <- rpois(20, 10) #et the initial number of infections to be Poisson distributed with a mean of 10
state0 <- rbind(1010 - I0, I0, 0, 0, 0, deparse.level = 0)
model$update_state(state = state0, step = 0L)
model$step()
model$state()

#reset the model
model <- walk$new(list(sd = 1), 0, 10, seed = 1L) #initial model with sd = 1
y1 <- model$run(100) #run for 100 steps
model$update_state(pars = list(sd = 2), step = 0) #set new parameters into the model and set the time back to zero
y2 <- model$run(100) #run again