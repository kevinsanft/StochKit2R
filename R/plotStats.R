#'@title Plot StochKit2R simulation statistics output data
#'
#'@description
#'\code{plotStats} Plots means and mean +/- one standard deviation of populations specified in \code{species}
#'
#'@param data ensemble output from ssa or tauLeaping (stats object must exist, i.e. ensemble must NOT have been run with noStats=TRUE).
#'@param species A character vector of species names or a numeric vector of species indexes of the species that will be plotted. For numeric indexes, the first species is index 1. By default =NULL, all species are plotted.
#'@return The ggplot object
#'@examples
#'\dontrun{
#'#example using included dimer_decay.xml file
#'model <- system.file("dimer_decay.xml",package="StochKit2R")
#'#output written to ex_out directory (created in current working directory)
#'out <- ssa(model,time=10,realizations=100,intervals=20)
#'#plot the data for all species
#'plotStats(out)
#'#plot the data for species S2 and S3
#'plotStats(out,species=c("S2","S3"))
#'}
plotStats <- function(data,species=NULL) {

  # if (is.null(data$stats)) {
  #   stop("data does not contain stats element.")
  # }
  if (nrow(data$stats$means)==0) {
    stop("data is missing stats information. Run ensemble with noStats=FALSE.")
  }
    
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
    indices = sapply(species,function(s) which(s==names(data$stats$means)))
    # if one or more names is not found
    if (class(indices)=="list") {
      bad_names = paste(species[sapply(indices,function(i) length(i)==0)],collapse=",")
      stop(paste("ERROR: invalid species name(s) (",bad_names,") entered.",sep=""))
    }
  }
  
  if (is.null(species)) {
    indices=1:(ncol(data$stats$means)-1)+1
  }
  
  if (max(indices)>ncol(data$stats$means)) {
    stop("ERROR: invalid species index")
  }
  #get the means data
  meansData <- data$stats$means[,c(1,indices)]
  #give (slightly) more meaningful labels if none
  if (!(names(data$stats$means[1])=="time")) {
    names(meansData) <- c("time",names(meansData)[1:(length(meansData)-1)])
  }
  #get the variances data
  variancesData <- data$stats$variances[,c(1,indices)]
  #give (slightly) more meaningful labels if none
  if (!(names(data$stats$variances[1])=="time")) {
    names(variancesData) <- names(meansData)
  }
  
  #put data into plottable form
  value=NULL#only to appease R CMD CHECK for CRAN
  variable=NULL#only to appease R CMD CHECK for CRAN
  meansDataMelted <- melt(meansData,id="time")
  p <- qplot(time,value,data=meansDataMelted,colour=variable) + geom_line() + ylab("population") + ggtitle("stats plot")
  # create mean+stdev data
  plusStdevData <- meansData
  plusStdevData[,2:ncol(plusStdevData)] <- plusStdevData[,2:ncol(plusStdevData)] + sqrt(variancesData[,2:ncol(variancesData)])
  plusStdevDataMelted <- melt(plusStdevData,id="time")
  p2 <- p + geom_line(data=plusStdevDataMelted,linetype="dashed")
  # create mean-stdev data
  minusStdevData <- meansData
  minusStdevData[,2:ncol(minusStdevData)] <- minusStdevData[,2:ncol(minusStdevData)] - sqrt(variancesData[,2:ncol(variancesData)])
  minusStdevDataMelted <- melt(minusStdevData,id="time")
  p3 <- p2 + geom_line(data=minusStdevDataMelted,linetype="dashed")
  return(p3)
}