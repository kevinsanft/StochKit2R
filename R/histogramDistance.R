#'@title Plot two StochKit2R simulation histograms and display distance
#'
#'@description
#'\code{histogramDistance} Plots histograms of data stored in two StochKit2R ensemble (ssa or tauLeaping) output objects and calculates histogram distance.
#'
#'@param data1 ensemble output from ssa or tauLeaping (histograms object must exist, i.e. ensemble must have been run with keepHistograms=TRUE).
#'@param species1 The species name or index in data1 to use for the first histogram. For numeric indexes, the first species is index 1.
#'@param timeIndex1 The time index in data1 to use for the second histogram (the initial condition is index 1, the index of the end time is equal to the number of output intervals+1). By default =NULL the last index is used.
#'@param data2 ensemble output (could be the same as data1) from ssa or tauLeaping (histograms object must exist, i.e. ensemble must have been run with keepHistograms=TRUE).
#'@param species2 The species name or index in data2 to use for the second histogram. For numeric indexes, the first species is index 1.
#'@param timeIndex2 The time index in data2 to use for the second histogram (the initial condition is index 1, the index of the end time is equal to the number of output intervals+1). By default =NULL the last index is used.
#'@return The ggplot object
#'@examples
#'\dontrun{
#'#example using included dimer_decay.xml file
#'model <- system.file("dimer_decay.xml",package="StochKit2R")
#'#output written to ex_out directory (created in current working directory)
#'result1 = ssa(model,10,100,20,keepHistograms=TRUE)
#'#another ensemble
#'result2 = ssa(model,10,100,20,keepHistograms=TRUE)
#'#plot the histograms for species 2 ("S2") at time point 5 (t=2.0) for the two runs above
#'histogramDistance(result1,"S2",5,result2,"S2",5)
#'#compare species from the same ensemble at different time points
#'histogramDistance(result1,"S2",5,result1,"S2",6) #warning about time mismatch
#'}
histogramDistance <- function(data1,species1,timeIndex1=NULL,data2,species2,timeIndex2=NULL) {
  
  # if (is.null(data1$histograms)) {
  #   stop("data1 does not contain histograms element. Run ensemble with keepHistograms=TRUE.")
  # }
  
  if (length(data1$histograms)==0) {
    stop("data1 does not contain histogram data. Run ensemble with keepHistograms=TRUE.")
  }
  
    
  if (!(class(species1)=="integer" || class(species1)=="numeric" || class(species1)=="character")) {
    stop("species1 must be a species name or index")
  }
  if (length(species1)!=1) {
    stop("species1 must be a single name or index.")
  }
  
  if (class(species1)=="integer" || class(species1)=="numeric") {
    # convert to species name
    if (!round(species1)==species1) {
      stop('species1 index must be an integer')
    }
    if (species1<1 || species1>length(names(data1$histograms))) {
      stop('Invalid species1 index.')
    }
    species1 = names(data1$histograms)[species1]
  } else {
    # verify species name
    if (!any(species1==names(data1$histograms))) {
      stop("Invalid species1 name.")
    }
  }
  species1Index = which(names(data1$histograms)==species1)
  speciesHistogramData1 = data1$histograms[[species1Index]]
  
  # verify timeIndex
  if (is.null(timeIndex1)) {
    timeIndex1=length(speciesHistogramData1)
  } else {
    if (class(timeIndex1)!="integer" && class(timeIndex1)!="numeric") {
      stop("timeIndex1 must be an integer.")
    }
    if (!round(timeIndex1)==timeIndex1) {
      stop('timeIndex1 must be an integer')
    }
    if (timeIndex1<1 || timeIndex1>length(speciesHistogramData1)) {
      stop("Invalid timeIndex1")
    }
  }
  
  histogramData1 = speciesHistogramData1[[timeIndex1]]
  
  lines1 <- strsplit(histogramData1,split="\t")

  sID1 <- lines1[[1]][1] # species ID
  time1 = as.numeric(lines1[[1]][2]) # time
  sInd1 = as.numeric(lines1[[1]][2]) # species index
  tInd1 = as.numeric(lines1[[1]][3]) # time index

  lb1 = as.numeric(lines1[[2]][1]) # lower bound
  ub1 = as.numeric(lines1[[2]][2]) # upper bound
  width1 = as.numeric(lines1[[2]][3]) # width
  bsize1 = as.numeric(lines1[[2]][4]) # number of bins

  # read in the third line (bin counts)
  data1 = as.numeric(lines1[[3]])
  # error checking
  if (bsize1!=length(data1)) {
    stop('retrieved size of bins are not equal to the number of bins in the data line')
  }

  dataRange1 <- seq(lb1,ub1,by=width1)
  bincenters1 <- (dataRange1[-1] + dataRange1[-length(dataRange1)])/2
  df1 <- data.frame(centers=bincenters1,counts=data1)

  # second histogram
  # if (is.null(data2$histograms)) {
  #   stop("data2 does not contain histograms element. Run ensemble with keepHistograms=TRUE.")
  # }
  
  if (length(data2$histograms)==0) {
    stop("data2 does not contain histogram data. Run ensemble with keepHistograms=TRUE.")
  }
  
  if (!(class(species2)=="integer" || class(species2)=="numeric" || class(species2)=="character")) {
    stop("species2 must be a species name or index")
  }
  if (length(species2)!=1) {
    stop("species2 must be a single name or index.")
  }
  
  if (class(species2)=="integer" || class(species2)=="numeric") {
    # convert to species name
    if (!round(species2)==species2) {
      stop('species2 index must be an integer')
    }
    if (species2<1 || species2>length(names(data2$histograms))) {
      stop('Invalid species2 index.')
    }
    species2 = names(data2$histograms)[species2]
  } else {
    # verify species name
    if (!any(species2==names(data2$histograms))) {
      stop("Invalid species2 name.")
    }
  }
  species2Index = which(names(data2$histograms)==species2)
  speciesHistogramData2 = data2$histograms[[species2Index]]
  
  # verify timeIndex
  if (is.null(timeIndex2)) {
    timeIndex2=length(speciesHistogramData2)
  } else {
    if (class(timeIndex2)!="integer" && class(timeIndex2)!="numeric") {
      stop("timeIndex2 must be an integer.")
    }
    if (!round(timeIndex2)==timeIndex2) {
      stop('timeIndex2 must be an integer')
    }
    if (timeIndex2<1 || timeIndex2>length(speciesHistogramData2)) {
      stop("Invalid timeIndex2")
    }
  }
  
  histogramData2 = speciesHistogramData2[[timeIndex2]]
  
  lines2 <- strsplit(histogramData2,split="\t")

  sID2 <- lines2[[1]][1] # species ID
  time2 = as.numeric(lines2[[1]][2]) # time
  sInd2 = as.numeric(lines2[[1]][2]) # species index
  tInd2 = as.numeric(lines2[[1]][3]) # time index
  
  lb2 = as.numeric(lines2[[2]][1]) # lower bound
  ub2 = as.numeric(lines2[[2]][2]) # upper bound
  width2 = as.numeric(lines2[[2]][3]) # width
  bsize2 = as.numeric(lines2[[2]][4]) # number of bins
  
  # read in the third line (bin counts)
  data2 = as.numeric(lines2[[3]])
  # error checking
  if (bsize2!=length(data2)) {
    stop('retrieved size of bins are not equal to the number of bins in the data line')
  }
  
  dataRange2 <- seq(lb2,ub2,by=width2)
  bincenters2 <- (dataRange2[-1] + dataRange2[-length(dataRange2)])/2
  centers=NULL#only to appease R CMD CHECK for CRAN
  counts=NULL#only to appease R CMD CHECK for CRAN
  df2 <- data.frame(centers=bincenters2,counts=data2)
 
  if (sID1!=sID2) {
    warning(paste('Species ID mismatch, ID1:',sID1,", ID2:",sID2))
  }
  if (time1!=time2) {
    warning(paste('Time point mismatch, Time1:',time1,", Time2:",time2))
  }
  if (bsize1!=bsize2) {
    stop('two histograms must have the same number of bins')
  }
  
  # Syncronize the histograms if necessary
  llb1 = lb1 + (which(data1>0)[1]-1)*width1
  llb2 = lb2 + (which(data2>0)[1]-1)*width2
  glb = min(llb1,llb2)
  lub1 = lb1 + (which(data1>0)[length(which(data1>0))])*width1
  lub2 = lb2 + (which(data2>0)[length(which(data2>0))])*width2
  gub = max(lub1,lub2)
  gwidth = max(width1,width2)

  if (glb<0 | glb>gub) {
    stop("error in bounds")
  }
  nlb = floor((glb/gwidth)*gwidth)
  nub = nlb+length(data1)*gwidth
  while (gub > nub) {
    gwidth = gwidth*2
    nlb = floor(glb / gwidth) * gwidth
    nub = nlb + length(data1) * gwidth
  }
  gIwidth = 1/gwidth
  
  his1 = rep(0,length(data1))
  his2 = rep(0,length(data2))

  #data1
  for (i in 0:(length(data1)-1)) {
    if (data1[i+1] != 0) {
      event = lb1+i*width1
      if (nlb>event || event>=nub) {
        stop('invalid new upper and/or lower bound')
      }
      his1[1+floor((event-nlb)*gIwidth)] = his1[1+floor((event-nlb)*gIwidth)]+data1[i+1];
    }
  }
  
  #data2
  for (i in 0:(length(data2)-1)) {
    if (data2[i+1] != 0) {
      event2 = lb2+i*width2
      if (nlb>event2 ||event2>=nub) {
        stop('invalid new upper and/or lower bound')
      }
      his2[1+floor((event2-nlb)*gIwidth)] = his2[1+floor((event2-nlb)*gIwidth)]+data2[i+1]
    }
  }
 
  # plot
  dataRange <- seq(nlb,nub,by=gwidth)
  bincenters <- (dataRange[-1] + dataRange[-length(dataRange)])/2
  df1 <- data.frame(centers=bincenters,counts=his1)
  
  p <- ggplot(data=df1,aes(x=centers,y=counts)) + geom_bar(stat="identity",width=gwidth,colour="red",fill="red",alpha=.65) + ylab("Bin Counts") + xlab("Bin Centers") + ggtitle(paste("Species:",sID1," Time:",time1))

  df2 <- data.frame(centers=bincenters,counts=his2)
  p2 <- p + geom_bar(data=df2,stat="identity",width=gwidth,colour="dark blue",fill="dark blue",alpha=.65)
  
  maxCount <- max(his1,his2)
  
  # compute histogram distance
  his1= his1/sum(his1);
  his2= his2/sum(his2);
  
  euclidean_distance = sqrt(sum(abs(his1-his2)^2))
  manhattan_distance = sum(abs(his1-his2))
  p3 <- p2 + annotate("text", size=4, x = (dataRange[length(dataRange)]+dataRange[1])/2, y = maxCount*1.1, label = sprintf("Euclidean d=%2.3f, Manhattan d=%2.3f",euclidean_distance,manhattan_distance))
  return(p3)
}
