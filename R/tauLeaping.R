#'@title Explicit Adaptive Tau-leaping simulation
#'
#'@description
#'\code{tauLeaping} Runs an ensemble of Tau-leaping \code{realizations} (trajectories). Will revert to SSA if Tau-leaping is not beneficial.
#'
#'@param modelFile Character string with path to StochKit2 .xml model file
#'@param time Simulation time of each realization
#'@param realizations Number of realizations
#'@param intervals Number of output intervals. Default 0 outputs at end time only. 1=keep data at start and end time, 2=keep data at start, middle, and end times, etc. Note data is stored at (intervals+1) equally spaced time points.
#'@param noStats Do not keep means and variances data. Default \code{FALSE} creates stats output.
#'@param keepTrajectories Keep trajectory data.
#'@param keepHistograms Keep histogram data.
#'@param bins Number of histogram bins
#'@param seed Seed the random number generator. By default the seed is determined by the R random number generator, so the seed can also be set by calling \code{set.seed} in R immediately before calling \code{tauLeaping}
#'@param p Override default and specify the number of processes (threads) to use. By default (=0), the number of processes will be determined automatically (recommended). Ignored on systems without OpenMP support.
#'@param epsilon Set the tolerance (applicable to tauLeaping only), default is 0.03. Valid values: must be greater than 0.0 and less than 1.0
#'@param threshold Set the threshold (minimum number of reactions per leap before switching to ssa) for tauLeaping
#'@return List of statistics, trajectory and histogram output data.
#'@examples
#'\dontrun{
#'#example using included dimer_decay.xml file
#'#store path to StochKit2R dimer_decay.xml file in a variable
#'model <- system.file("dimer_decay.xml", package="StochKit2R")
#'out <- tauLeaping(modelFile=model, time=10, realizations=100, intervals=20)
#'
#'#more typical example where model file is stored elsewhere
#'#(must be valid path to existing .xml StochKit2 model file)
#'#also, keep trajectory data
#'out <- tauLeaping("Desktop/dimer_decay.xml", 10, 100, 20, keepTrajectories=TRUE)
#'}
tauLeaping <- function(modelFile,time,realizations,intervals=0,noStats=FALSE,keepTrajectories=FALSE,keepHistograms=FALSE,bins=32,seed=NULL,p=0,epsilon=0.03,threshold=10) {
  # can set seed in R with set.seed()
  
  #checks on modelFile  
  if (!file.exists(modelFile)) {
    stop("ERROR: Model file does not exist")
  }
  #simple checks
  if (time<=0) {
    stop("ERROR: simulation time must be positive")
  }
  if (realizations<=0) {
    stop("ERROR: realizations must be positive")
  }
  if (intervals<0) {
    stop("ERROR: intervals must not be negative")
  }
  if (bins<=0) {
    stop("ERROR: number of bins must be positive")
  }
  if (p<0) {
    stop("ERROR: number of processes must not be negative (when p=0, processes will be determined automatically)")
  }
  
  # if(!is.null(outputDir)){
  #   #expand path
  #   outputDir <- tryCatch(suppressWarnings(normalizePath(outputDir)), error = function(e) {stop("Invalid or missing outputDir output directory path, terminating StochKit2R")}, finally = NULL)
  #   #remove tailing slashes or backslashes
  #   #because file.exists returns false if directory ends in slash or backslash
  #   outputDir <- gsub("//*$","",outputDir)
  #   outputDir <- gsub("\\\\*$","",outputDir)
  #   if (file.exists(outputDir)) {
  #     if (!force) {
  #       stop("ERROR: Output directory already exists. Delete existing directory, specify a unique directory name, or run with force=TRUE to overwrite")
  #     }
  #     else {
  #       #delete stats, trajectories, histograms directories, if they exist
  #       unlink(file.path(outputDir,"stats"),recursive=TRUE)
  #       unlink(file.path(outputDir,"trajectories"),recursive=TRUE)
  #       unlink(file.path(outputDir,"histograms"),recursive=TRUE)
  #     }
  #   }
  #   else {
  #     dir.create(outputDir,recursive=TRUE)
  #   }
    
    #must keep some output
    if (noStats & !keepTrajectories & !keepHistograms) {
      stop("ERROR: must output at least one of stats, trajectories, or histograms.")
    }
    
    #create stats, trajectories, histograms output directories, if needed
  #   if (!noStats) {
  #     dir.create(file.path(outputDir,"stats"))
  #   }
  #   if (keepTrajectories) {
  #     dir.create(file.path(outputDir,"trajectories"))
  #   }
  #   if (keepHistograms) {
  #     dir.create(file.path(outputDir,"histograms"))
  #   }
  # }
  
  if (is.null(seed)) {
      seed=floor(runif(1,-.Machine$integer.max,.Machine$integer.max))
    }

  # if(is.null(outputDir)){
  #   outputDir=""
  # }

  # process the xml model file, put in R Lists
  StochKit2Rmodel <- buildStochKit2Rmodel(modelFile)
  
  #tauLeaping specific checks
  if (epsilon>=1.0 | epsilon <= 0.0) {
    stop("epsilon must be greater than 0.0 and less than 1.0")
  }
  
  if (threshold<0) {
    stop("threshold must be a positive integer")
  }
  
  if (as.integer(threshold)!=threshold) {
    warning("Converting threshold to integer")
    threshold=as.integer(threshold)
  }
  # everything is set up, ready to run the simulation
  tauLeapingStochKit2RInterface(StochKit2Rmodel,time,realizations,intervals,!noStats,keepTrajectories,keepHistograms,bins,"",seed,p,epsilon,threshold)
}