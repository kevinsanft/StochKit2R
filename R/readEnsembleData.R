#'@title Read StochKit2R ensemble output files (reads data from files written by writeEnsembleData)
#'
#'@description
#'\code{writeEnsembleOutput} Reads StochKit2R ensemble output data from directory. If data is from StochKit2R C++ starting indexes will be converted from 0 to 1.
#'
#'@param outputDir Character string with path to output directory.
#'@examples
#'\dontrun{
#'#'#example using included dimer_decay.xml file
#'#run 100 simulations for 10 time units, keeping output at 20 time intervals
#'#store model file name in a variable first
#'model <- system.file("dimer_decay.xml",package="StochKit2R")
#'out <- ssa(modelFile=model,time=10,realizations=100,intervals=20)
#'#creates ex_out directory and writes output
#'writeEnsembleData(data=out,outputDir="ex_out")
#'out2=readEnsembleData("ex_out") #out and out2 should match
#'}
readEnsembleData <- function(outputDir) {

  #expand path
  outputDir <- tryCatch(suppressWarnings(normalizePath(outputDir)), error = function(e) {stop("Invalid or missing outputDir output directory path, terminating.")}, finally = NULL)
  #remove tailing slashes or backslashes
  #because file.exists returns false if directory ends in slash or backslash
  outputDir <- gsub("//*$","",outputDir)
  outputDir <- gsub("\\\\*$","",outputDir)

  if (!file.exists(outputDir)) {
    stop("ERROR: outputDir does not exist")
  }
  out=list(stats=list(means=data.frame(),variances=data.frame()),trajectories=list(),histograms=list())
    
  # see what subdirectories we have
  subdirs = dir(outputDir)

  for (dir in subdirs) {
    if (dir=="stats") {
      # verify means.txt and variances.txt files exist
      if (! (file.exists(file.path(outputDir,dir,"means.txt")) && file.exists(file.path(outputDir,dir,"variances.txt")))) {
        warning("stats directory exists in outputDir but means.txt or variances.txt does not exist. Statistics data not returned.")
      } else {
        # check for header row
        firstRow = readLines(file.path(outputDir,dir,"means.txt"),n=1)
        if (substr(firstRow,1,4)=="time") {
          hasHeader=T
        } else {
          hasHeader=F
          warning("stats data does not have header row. A generic one (X1, X2, ..., XN) will be created.")
        }
        ##
        out$stats$means = read.table(file.path(outputDir,dir,"means.txt"),header=hasHeader)
        out$stats$variances = read.table(file.path(outputDir,dir,"variances.txt"),header=hasHeader)
        if (!hasHeader) {
          speciesNames = paste("X",1:(ncol(out$stats$means)-1),sep="")
          names(out$stats$means) = c("time",speciesNames)
          names(out$stats$variances) = c("time",speciesNames)          
        }
      }     
    } else if (dir=="trajectories") {

      trajectoryFiles = dir(file.path(outputDir,dir))
      # verify all files start with "trajectory"
      if (length(grep("trajectory[0-9]+\\.txt",trajectoryFiles))!=length(trajectoryFiles)) {
        stop("Unexpected file in trajectories directory.")
      }
      # get indexes
      trajectoryIndexes = sort(as.integer(substring(trajectoryFiles,11,(regexpr("\\.",trajectoryFiles)-1))))
      minIndex = min(trajectoryIndexes)
      if (minIndex!=1) {
        message(paste("MESSAGE: Minimum trajectory index ",minIndex," will be index 1 in returned trajectories list.",sep=""))
      }
      if (!identical(minIndex:max(trajectoryIndexes),trajectoryIndexes)) {
        warning("Trajectory file indexes are not consecutive.")
      }
      # iterate over the files
      # check for header row
      firstRow = readLines(file.path(outputDir,dir,trajectoryFiles[1]),n=1)
      if (substr(firstRow,1,4)=="time") {
        hasHeader=T
      } else {
        hasHeader=F
        warning("trajectory data does not have header row. A generic one (X1, X2, ..., XN) will be created.")
      }
      
      for (i in 1:length(trajectoryIndexes)) {
        filename = paste("trajectory",trajectoryIndexes[i],".txt",sep="")
        out$trajectories[[i]] = read.table(file.path(outputDir,dir,filename),header=hasHeader)
        if (!hasHeader) {
          speciesNames = paste("X",1:(ncol(out$trajectories[[i]])-1),sep="")
          names(out$trajectories[[i]]) = c("time",speciesNames)
        }
      }
      
    } else if (dir=="histograms") {
      histogramFiles = dir(file.path(outputDir,dir))
      # verify all files start with "hist"
      if (length(grep("hist_[0-9]+_[0-9]+\\.dat",histogramFiles))!=length(histogramFiles)) {
        stop("Unexpected file in histogram directory.")
      }
      # get indexes
      fileSpeciesAndTimeIndexes = strsplit(substring(histogramFiles,6,(regexpr("\\.",histogramFiles)-1)),split="_")
      speciesIndexes = sort(as.integer(unique(sapply(fileSpeciesAndTimeIndexes,function(s) s[1]))))
      timeIndexes = sort(as.integer(unique(sapply(fileSpeciesAndTimeIndexes,function(s) s[2]))))
      
      minSpeciesIndexes = min(speciesIndexes)
      if (minSpeciesIndexes!=1) {
        message(paste("MESSAGE: Minimum histogram species index ",minSpeciesIndexes," will be index 1 in returned histogram data (renumbering).",sep=""))
      }
      if (!identical(minSpeciesIndexes:max(speciesIndexes),speciesIndexes)) {
        warning("Histogram file species indexes are not consecutive, renumbering!")
      }
      minTimeIndexes = min(timeIndexes)
      if (minTimeIndexes!=1) {
        message(paste("MESSAGE: Minimum histogram time index ",minTimeIndexes," will be index 1 in returned histogram data (renumbering).",sep=""))
      }
      if (!identical(minTimeIndexes:max(timeIndexes),timeIndexes)) {
        warning("Histogram file time indexes are not consecutive, renumbering!")
      }
      
      # iterate over species indexes
      for (i in 1:length(speciesIndexes)) {
        # get species name from first file
        firstfilename = paste("hist_",speciesIndexes[i],"_",timeIndexes[1],".dat",sep="")
        speciesName = strsplit(readLines(file.path(outputDir,dir,firstfilename),n=1),split="\t")[[1]][1]
        out$histograms[[i]] = list()
        names(out$histograms)[i]=speciesName
        for (j in 1:length(timeIndexes)) {
          filename = paste("hist_",speciesIndexes[i],"_",timeIndexes[j],".dat",sep="")
          out$histograms[[i]][[j]] = readLines(file.path(outputDir,dir,filename))
          # update species and time indexes in the data
          oldData = strsplit(out$histograms[[i]][[j]][1],split="\t")[[1]]
          newSpeciesIndex = i
          newTimeIndex = j
          out$histograms[[i]][[j]][1] = paste(c(oldData[1:2],as.character(newSpeciesIndex),as.character(newTimeIndex)),collapse="\t")
          
        }
      }
    } else {
      warning(paste("Unexpected directory in outputDir (",dir,").",sep=""))
    }
  }
  return(out)
}