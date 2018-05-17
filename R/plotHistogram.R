#'@title Plot StochKit2R simulation histogram output data
#'
#'@description
#'\code{plotHistogram} Plots histogram of data stored in StochKit2R output object or histogram output file.
#'
#'@param data ensemble output from ssa or tauLeaping (histograms object must exist, i.e. ensemble must have been run with keepHistograms=TRUE).
#'@param species The species name or index. For numeric indexes, the first species is index 1.
#'@param timeIndex The time index (the initial condition is index 1, the index of the end time is equal to the number of output intervals+1). By default =NULL the last index is used.
#'@return The ggplot object
#'@examples
#'\dontrun{
#'#example using included dimer_decay.xml file
#'model <- system.file("dimer_decay.xml",package="StochKit2R")
#'#output written to ex_out directory (created in current working directory)
#'out <- ssa(model,time=10,realizations=100,intervals=20,keepHistograms=TRUE)
#'#plot the histogram for species "S2" at time point 5 (t=2.0)
#'plotHistogram(out,species="S2",timeIndex=5)
#'}
plotHistogram <- function(data,species,timeIndex=NULL) {
  
  # if (is.null(data$histograms)) {
  #   stop("data does not contain histograms element. Run ensemble with keepHistograms=TRUE.")
  # }
  
  if (length(data$histograms)==0) {
    stop("data does not contain histogram data. Run ensemble with keepHistograms=TRUE.")
  }
  
  if (!(class(species)=="integer" || class(species)=="numeric" || class(species)=="character")) {
    stop("species must be a species name or index")
  }
  if (length(species)!=1) {
    stop("species must be a single name or index.")
  }
  
  if (class(species)=="integer" || class(species)=="numeric") {
    # convert to species name
    if (!round(species)==species) {
      stop('Species index must be an integer')
    }
    if (species<1 || species>length(names(data$histograms))) {
      stop('Invalid species index.')
    }
    species = names(data$histograms)[species]
  } else {
    # verify species name
    if (!any(species==names(data$histograms))) {
      stop("Invalid species name.")
    }
  }
  speciesIndex = which(names(data$histograms)==species)
  speciesHistogramData = data$histograms[[speciesIndex]]
  
  # verify timeIndex
  if (is.null(timeIndex)) {
    timeIndex=length(speciesHistogramData)
  } else {
    if (class(timeIndex)!="integer" && class(timeIndex)!="numeric") {
      stop("timeIndex must be an integer.")
    }
    if (!round(timeIndex)==timeIndex) {
      stop('timeIndex must be an integer')
    }
    if (timeIndex<1 || timeIndex>length(speciesHistogramData)) {
      stop("Invalid timeIndex")
    }
  }
  
  histogramData = speciesHistogramData[[timeIndex]]
  
  lines <- strsplit(histogramData,split="\t")

  sID <- lines[[1]][1] # species ID
  time = as.numeric(lines[[1]][2]) # time
  sInd = as.numeric(lines[[1]][2]) # species index
  tInd = as.numeric(lines[[1]][3]) # time index

  lb = as.numeric(lines[[2]][1]) # lower bound
  ub = as.numeric(lines[[2]][2]) # upper bound
  width = as.numeric(lines[[2]][3]) # width
  bsize = as.numeric(lines[[2]][4]) # number of bins

  # read in the third line (bin counts)
  dat = as.numeric(lines[[3]])
  # error checking
  if (bsize!=length(dat)) {
    stop('retrieved size of bins are not equal to the number of bins in the data line')
  }

  dataRange <- seq(lb,ub,by=width)
  bincenters <- (dataRange[-1] + dataRange[-length(dataRange)])/2
  centers=NULL#only to appease R CMD CHECK for CRAN
  counts=NULL#only to appease R CMD CHECK for CRAN
  df <- data.frame(centers=bincenters,counts=dat)
  
  ggplot(data=df,aes(x=centers,y=counts)) + geom_bar(stat="identity",width=width,colour="dark blue",fill="dark blue") + ylab("Bin Counts") + xlab("Bin Centers") + ggtitle(paste("Species:",sID," Time:",time))
}
