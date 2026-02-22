#' @docType package
#' @name StochKit2R
#'
#' @title StochKit2R
#' @description
#' StochKit2R: a package for efficient discrete stochastic simulation
#' using the Gillespie Algorithm (SSA) and tau-leaping. The package
#' provides an R interface to stochastic simulation functions and
#' provides functions for generating plots of the simulation data.
#'
#' @useDynLib StochKit2R, .registration=TRUE
#'
#' @import reshape
#' @import ggplot2
#' @import XML
#' @importFrom Rcpp evalCpp
#' @importFrom RcppXPtrUtils cppXPtr
#' @importFrom R.utils mkdirs
#' @importFrom stats runif time
#' @importFrom utils read.table write.table
#' @export ssa
#' @export ssaSingle
#' @export tauLeaping
#' @export plotStats
#' @export plotTrajectories
#' @export plotHistogram
#' @export histogramDistance
#' @export ode
#' @export writeEnsembleData
#' @export readEnsembleData
"_PACKAGE"
