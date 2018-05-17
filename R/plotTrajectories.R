#'@title Plot StochKit2R simulation trajectory output data
#'
#'@description
#'\code{plotTrajectories} Plots trajectories \code{trajectoryIndexes} of populations specified in \code{species} in the StochKit2R trajectories element of \code{data}
#'
#'@param data ensemble output from ssa or tauLeaping (trajectories object must exist, i.e. ensemble must have been run with keepTrajectories=TRUE).
#'@param trajectoryIndexes Integer vector containing indexes for which simulations will be plotted. By default (=NULL) all trajectories are plotted.
#'@param species A character vector of species names or a numeric vector of species indexes of the species that will be plotted. For numeric indexes, the first species is index 1. By default =NULL, all species are plotted.
#'@return The ggplot object
#'@examples
#'\dontrun{
#'#example using included dimer_decay.xml file
#'model <- system.file("dimer_decay.xml",package="StochKit2R")
#'#output written to ex_out directory (created in current working directory)
#'out <- ssa(model,time=10,realizations=100,intervals=20,keepTrajectories=TRUE)
#'#plot the data for species S2 and S3 for trajectories 2,3,4 and 5
#'plotTrajectories(out,trajectoryIndexes=2:5,species=c("S2","S3"))
#'}
plotTrajectories <- function(data,trajectoryIndexes=NULL,species=NULL) {
  
  # if (is.null(data$trajectories)) {
  #   stop("data does not contain trajectories element. Run ensemble with keepTrajectories=TRUE.")
  # }

  if (length(data$trajectories)==0) {
    stop("trajectories element of input data has length 0. Run ensemble with keepTrajectories=TRUE.")
  }  
  
  # checks on trajectoryIndexes
  if (!is.null(trajectoryIndexes)) {
    if (! (class(trajectoryIndexes)=="numeric" || class(trajectoryIndexes)=="integer")) {
      stop("trajectoryIndexes must be a character vector of species names or a numeric vector")
    }
    
    # Check to make sure that the indices are integers
    if (sum((round(trajectoryIndexes)==trajectoryIndexes)) != length(trajectoryIndexes)) {
      stop('trajectoryIndexes indexes must be integers')
    }
    
    # Check to make sure that the indices are greater than 0
    if (min(trajectoryIndexes)<=0) {
      stop('trajectoryIndexes must be positive integers')
    }
    
    # Check to eliminate duplicate species
    trajectoryIndexes2 = unique(trajectoryIndexes);
    if (length(trajectoryIndexes2)!=length(trajectoryIndexes)) {
      warning(paste("removing ",length(trajectoryIndexes) - length(trajectoryIndexes2), " duplicate trajectoryIndexes"))
      trajectoryIndexes=trajectoryIndexes2;
    }
  } else {
    trajectoryIndexes=1:length(data$trajectories)
  }
  
  # checks on species
  if (!is.null(species)) {
    if (! (class(species)=="character" || class(species)=="numeric" || class(species)=="integer")) {
      stop("Species must be a character vector of species names or a numeric vector")
    }
    
    if (class(species)=="numeric" || class(species)=="integer") {
      # Check to make sure that the indices are integers
      if (sum((round(species)==species)) != length(species)) {
        stop('Species indexes must be integers')
      }
      
      # Check to make sure that the indices are greater than 0
      if (min(species)<=0) {
        stop('Species indexes must be positive integers')
      }
    }
    
    # Check to eliminate duplicate species
    species2 = unique(species);
    if (length(species2)!=length(species)) {
      warning(paste("removing ",length(species) - length(species2), " duplicate species."))
      species=species2;
    }
    
    #if numeric, add 1 (if character, handle it later)
    if (class(species)=="numeric" || class(species)=="integer") {
      indices = species+1; # add 1 to account for the time column
      indices = sort(indices);
      len_indices =length(indices);
    }
  }
  
  #convert species names to indices
  if (class(species)=="character") {
    indices = sapply(species,function(s) which(s==names(data$trajectories[[1]])))
    # if one or more names is not found
    if (class(indices)=="list") {
      bad_names = paste(species[sapply(indices,function(i) length(i)==0)],collapse=",")
      stop(paste("ERROR: invalid species name(s) (",bad_names,") entered.",sep=""))
    }
  }
  
  if (is.null(species)) {
    indices=1:(ncol(data$trajectories[[1]])-1)+1
  }
  
  #plot data...
  
  trajectoriesData=data$trajectories
  trajData <- trajectoriesData[[1]][,c(1,indices)]
  
  #put data into plottable form
  value=NULL#only to appease R CMD CHECK for CRAN
  variable=NULL#only to appease R CMD CHECK for CRAN
  trajDataMelted <- melt(trajData,id="time")
  
  p <- qplot(time,value,data=trajDataMelted,colour=variable,geom="line") + ylab("population") + ggtitle("trajectory plot")
  
  # plot the rest of the trajectories!
  if (length(trajectoryIndexes)>1) {
    
    for (i in trajectoryIndexes[2:length(trajectoryIndexes)]) {
      trajData <- trajectoriesData[[i]][,c(1,indices)]

      trajDataMelted2 <- melt(trajData,id="time")
      p <- p + geom_line(data=trajDataMelted2)
    }
  }  
  
  return(p)
}
