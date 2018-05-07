#'@title Gillespie Stochastic Simulation Algorithm ordinary differential equation simulation
#'
#'@description
#'\code{ode} Run an ODE simulation (population version of the reaction rate equations). Homodimerization (e.g. A+A type) propensities are converted to Reaction Rate Equation rates (c*A*A instead of c*A*(A-1)/2) for mass action reactions. Propensities of customized propensities are not converted. Uses the 5th order Dormand-Prince Runge-Kutta method.
#'
#'@param modelFile Character string with path to StochKit2 .xml model file
#'@param time Simulation time
#'@param intervals Number of output intervals. Default 0 outputs at end time only. 1=keep data at start and end time, 2=keep data at start, middle, and end times, etc. Note data is stored at (intervals+1) equally spaced time points.
#'@return Data frame containing time and species population
#'@examples
#'\dontrun{
#'#example using included dimer_decay.xml file
#'out <- ode(system.file("dimer_decay.xml",package="StochKit2R"),time=10,intervals=20)
#'}
ode <- function(modelFile,time,intervals) {

  #checks on modelFile  
  if (!file.exists(modelFile)) {
    stop("ERROR: Model file does not exist")
  }

  # process the xml model file, put in R Lists
  StochKit2Rmodel <- buildStochKit2Rmodel(modelFile)

  # everything is set up, ready to run the simulation
  odecpp(StochKit2Rmodel,time,intervals)
}