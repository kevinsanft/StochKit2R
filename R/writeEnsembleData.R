#'@title Write StochKit2R data to files
#'
#'@description
#'\code{writeEnsembleData} Writes StochKit2R ensemble output data to files.
#'
#'@param data StochKit2R ensemble (from ssa or tauLeaping) output data
#'@param outputDir Character string with path to output directory. If specified output directory does not exist, it will be created. If output directory or file already exists, use \code{force=TRUE} to overwrite.
#'@param force Force overwriting of existing data. Warning: force=TRUE will delete all existing content in the output directory.
#'@param indexStart Starting index. Default =1 for R data. Set to 0 for compatibility with StochKit2 C++ version.
#'@examples
#'\dontrun{
#'#'#example using included dimer_decay.xml file
#'#run 100 simulations for 10 time units, keeping output at 20 time intervals
#'#store model file name in a variable first
#'model <- system.file("dimer_decay.xml",package="StochKit2R")
#'out <- ssa(modelFile=model,time=10,realizations=100,intervals=20)
#'#creates ex_out directory and writes output
#'writeEnsembleData(data=out,outputDir="ex_out")
#'}
writeEnsembleData <- function(data,outputDir,force=FALSE,indexStart=1) {
  if (!(indexStart==1 || indexStart==0)) {
    warning("indexStart is not 0 or 1")
  }
  
  #expand path
  outputDir <- tryCatch(suppressWarnings(normalizePath(outputDir)), error = function(e) {stop("Invalid or missing outputDir output directory path, terminating.")}, finally = NULL)
  #remove tailing slashes or backslashes
  #because file.exists returns false if directory ends in slash or backslash
  outputDir <- gsub("//*$","",outputDir)
  outputDir <- gsub("\\\\*$","",outputDir)
  
  # see what subdirectories we will have
  subdirs = character(0)
  if (nrow(data$stats$means)>0) subdirs = c(subdirs,"stats")
  if (length(data$trajectories)>0) subdirs = c(subdirs,"trajectories")
  if (length(data$histograms)>0) subdirs = c(subdirs,"histograms")
  
  #message(paste("creating outputDir",outputDir,", subdirs: ",paste(subdirs,collapse=",")))
  createOutputDirs(outputDir,subdirs,force)
  
  for (dir in subdirs) {
    if (dir=="stats") {
      write.table(data$stats$means,file=file.path(outputDir,dir,"means.txt"),sep="\t",quote=F,row.names=F)
      write.table(data$stats$variances,file=file.path(outputDir,dir,"variances.txt"),sep="\t",quote=F,row.names=F)
    }
    if (dir=="trajectories") {
      realizations = length(data$trajectories)
      indexOffset = indexStart-1
      cat("indexOffset=",indexOffset)
      for (i in 1:realizations) {
        index=i+indexOffset
        cat("index=",index)
        filename = paste("trajectory",index,".txt",sep="")
        write.table(data$trajectories[[i]],file=file.path(outputDir,dir,filename),sep="\t",quote=F,row.names=F)
      }
    }
    if (dir=="histograms") {
      # indexStart is used in histogram file names and in first line of data in each file
      # compute the start index in the data (should be 1 for StochKit2R data)
      firstLine=data$histograms[[1]][[1]][1]
      firstLineSplit=strsplit(firstLine,split="\t")[[1]]
      firstSpecies=firstLineSplit[1]
      firstSpeciesIndex=as.integer(firstLineSplit[3])
      firstTimeIndex=as.integer(firstLineSplit[4])
      if (!(firstTimeIndex==0 || firstTimeIndex==1)) {
        warning(paste("Time index in first histogram data object is not 0 or 1 (detected it is",firstTimeIndex,").",sep=""))
      }
      # firstSpecies should match name of first element in histogram list
      if (firstSpecies!=names(data$histograms)[1]) {
        stop("First species in histogram list does not match species name in data.")
      }
      if (!(firstSpeciesIndex==0 || firstSpeciesIndex==1)) {
        warning(paste("Species index in first histogram data object is not 0 or 1 (detected it is ",firstSpeciesIndex,").",sep=""))
      }
      if (firstSpeciesIndex!=firstTimeIndex) {
        stop("First species index and first time index do not match in first histogram data object.")
      }
      
      # difference between indexStart and detected firstIndex
      indexDiff = indexStart-firstTimeIndex
      
      # if firstIndex is not the same as indexStart, adjust the DATA
      if (indexDiff!=0) {
        message(paste("Adjusting histogram data for indexStart=",indexStart,sep=""))
        # for each species
        for (i in 1:length(data$histograms)) {
          # for each timepoint
          for (j in 1:length(data$histograms[[i]])) {
            oldData = strsplit(data$histograms[[i]][[j]][1],split="\t")[[1]]
            newSpeciesIndex = as.integer(oldData[3])+indexDiff
            newTimeIndex = as.integer(oldData[4])+indexDiff
            data$histograms[[i]][[j]][1] = paste(c(oldData[1:2],as.character(newSpeciesIndex),as.character(newTimeIndex)),collapse="\t")
          }
        }
        
      }
      
      # write all data to files
      # for each species
      for (i in 1:length(data$histograms)) {
        # for each timepoint
        for (j in 1:length(data$histograms[[i]])) {
          line1split = strsplit(data$histograms[[i]][[j]][1],split="\t")[[1]]
          speciesIndex = line1split[3]
          timeIndex = line1split[4]
          filename = paste("hist_",speciesIndex,"_",timeIndex,".dat",sep="")
          full_filename = file.path(outputDir,dir,filename)
          writeLines(data$histograms[[i]][[j]],full_filename)
        }
      }
    }
  }
}