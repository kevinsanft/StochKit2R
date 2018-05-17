#ifndef _HISTOGRAM_SINGLE_H_
#define _HISTOGRAM_SINGLE_H_
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <vector>
#include "boost/foreach.hpp"
#include "boost/tokenizer.hpp"
#include <limits>

#define to_string_macro( x ) static_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

//export LD_LIBRARY_PATH=libs/boost_1_41_0/stage/lib/

namespace STOCHKIT
{
 template<typename _populationValueType>
 class HistogramSingle
 {
public:

 protected:
  //! The closed lower bound is a multiple of the width.
  double _lowerBound;
  //! The open upper bound.
  /*! _upperBound == _lowerBound + _width * _bins.size() */
  double _upperBound;
  //! The width of a bin is a power of 2.
  double _width;
  //! The inverse of the bin width.
  double _inverseWidth;
  //! The number of bins.
  std::size_t _size;
  //! Time index of the histogram
  std::size_t _timeIndex;
  //! Species index of the histogram
  std::size_t _speciesIndex;
  //! The data container
  std::vector<int> _data;
  void setHistogramData(std::vector<int> data, double lb, double ub, double wd, double iwd){
    _lowerBound=lb;
    _upperBound=ub;
    _width=wd;
    _inverseWidth=iwd;
#ifndef DEBUG_StochKit 
    if (data.size()!=_size){
      Rcpp::Rcout<<"StochKit ERROR (HistogramSingle::setHistogramData): requires bin sizes to be the same\n";
      Rcpp::stop("Fatal error encountered, terminating StochKit2R");
    }
#endif
    for(std::size_t i=0;i<_size;i++){
      _data[i]=data[i];
    }
    // Rcpp::Rcout<<"lb: "<<_lowerBound<<" ub:"<<_upperBound<<" width:"<<_width<<" iwidth:"<<_inverseWidth<<'\n';
  }

 
 public:

 HistogramSingle() :
  _lowerBound(0),
    _upperBound(0),
    _width(0),
    _inverseWidth(0),
    _size(0) 
      { 
      }
    
  //! Construct from the number of bins 
 HistogramSingle(const std::size_t size,const std::size_t tIndex, const std::size_t sIndex) :
  _lowerBound(0),
    _upperBound(size),
    _width(1),
    _inverseWidth(1),
    _size(size),
    _timeIndex(tIndex),
    _speciesIndex(sIndex),
    _data(size,0){
#ifndef DEBUG_StochKit 
    if(size<1){
      Rcpp::Rcout<<"StochKit ERROR (HistogramSingle::HistogramSingle): requires bin size greater than 0\n";
      Rcpp::stop("Fatal error encountered, terminating StochKit2R");
    }
#endif
  }

  //! Copy constructor. Used in merge to avoid error by modifying constant (referenced) histogram object
 HistogramSingle(const HistogramSingle& x) :
  _lowerBound(x._lowerBound),
    _upperBound(x._upperBound),
    _width(x._width),
    _inverseWidth(x._inverseWidth),
    _size(x._size),
    _timeIndex(x._timeIndex),
    _speciesIndex(x._speciesIndex),
    _data(x._data) {
  }

  //! Destructor
  ~HistogramSingle() {}

  //! Initialize with the number of bins and a pointer to the bin data.
  
  void initialize(const std::size_t size, std::size_t tIndex, std::size_t sIndex) {
#ifndef DEBUG_StochKit
    if (size<1){
      Rcpp::Rcout<<"StochKit ERROR (HistogramSingle::initialize): requires bin size greater than 0\n";
      Rcpp::stop("Fatal error encountered, terminating StochKit2R");
    }
    if (_size>0){ // bin size is already determined. cannot be overwriiten unless syncro function is used
      Rcpp::Rcout<<"StochKit ERROR (HistogramSingle::initialize): bin size is already determined\n";
      Rcpp::stop("Fatal error encountered, terminating StochKit2R");
    }
#endif
    _lowerBound = 0;
    _width = 1;
    _inverseWidth = 1;
    _upperBound = size;
    _size = size;
    _timeIndex = tIndex;
    _speciesIndex = sIndex;
    std::size_t i;
    for(i=0; i!=size;i++)
      _data.push_back(0);
  }

  //! Return the number of bins.
  std::size_t size() const {
    return _size;
  }

  //! Return the closed lower bound.
  double getLowerBound() const {
    return _lowerBound;
  }

  //! Return the bin width.
  double getWidth() const {
    return _width;
  }

  //! Return the open upper bound.
  double getUpperBound() const {
    return _upperBound;
  }

  //! Return the inverse width, used in a copy constructor.
  double getInverseWidth() const {
    return _inverseWidth;
  }

  //! Return the time index
  std::size_t _getTimeIndex() const {
    return _timeIndex;
  }

  //! Return the species index
  std::size_t _getSpeciesIndex() const {
    return _speciesIndex;
  }
 
  //! Clear the histogram. Empty the bins but do not alter the lower bound or width.
  void clear() {
    std::size_t i;
    for(i=0; i!=_data.size();i++){
      _data[i]=0;
    }
  }

  //! Reset the histogram. Empty the bins. Reset the lower bound and width.
  void reset() {
    _lowerBound = 0;
    _width = 1;
    _inverseWidth = 1;
    _upperBound = size();
    clear();
  }

  //! Compute the sum of counts for all bins
  int computeSum() const {
    int tSum=0;
    for(std::size_t i=0;i<_size; i++)
      tSum += _data[i];
    return tSum;
  }

  //! Compute the minimum non-zero element
  double computeMinimumNonzero() const {
    for (std::size_t i = 0; i != _size; i++) {
      if (_data[i] != 0) {
	return _lowerBound + i * _width;
      }
    }
#ifndef DEBUG_StochKit
    Rcpp::Rcout<<"StochKit ERROR (HistogramSingle::computeMinimumNonzero)  \n";
    return 0.;
#endif
  }
 
  //! Compute the minimum non-zero element
  double  computeMaximumNonzero() const {
    for (std::size_t i = _size; i != 0; i--) {
      if (_data[i-1] != 0) {
	return _lowerBound + i * _width;
      }
    }
#ifndef DEBUG_StochKit
    Rcpp::Rcout<<"StochKit ERROR (HistogramSingle::computeMaximumNonzero)  \n";
    return 0.;
#endif
  }
  //! Compute the number of data in the histogram
  int getCounts(std::size_t index) const{
#ifndef DEBUG_StochKit
    if (index<0 || index>=_size){
      Rcpp::Rcout<<"StochKit ERROR (HistogramSingle::getCounts) requires 0 <= index < numberOfBins  \n";
      return 1;
    }
#endif
    return(_data[index]);
  }

  //! Add an event
  void accumulate(const _populationValueType event) {

    // The bin width is a power of 2.
    // The lower bound is a multiple of the bin width.
    // If an event lies outside of the range of the current histogram, the bin width is doubled (and the lower bound adjusted) until all events lie within the histogram's rang

#ifndef DEBUG_StochKit
    if (event<0){
      Rcpp::Rcout<<"StochKit ERROR (HistogramSingle::accumulate) requires population to be non-negative \n";
      Rcpp::stop("Fatal error encountered, terminating StochKit2R");
    }
#endif
    rebuild(event);
    _data[std::size_t((event-_lowerBound)*_inverseWidth)]+=1;
    return;
  }
   
  //! If necessary, rebuild the histogram so it can contain the specified event.
  void rebuild(const _populationValueType event) {
    // Do nothing if the event will be placed in the current histogram.
    if (_lowerBound <= event && event < _upperBound) {
      return;
    }

    // Determine the new lower bound.
    _populationValueType lower = event;
    if (_lowerBound < lower) {
      std::size_t i;
      for (i = 0; i != size(); i++) {
	if (_data[i] != 0) {
	  break;
	}
      }
      if (i != size() && _lowerBound + i * _width < lower) {
	lower = _lowerBound + i * _width;
      }
    }
    // Determine the new open upper bound.
    // Add one half to get an open upper bound.
    double upper = event + 0.5;
    if (_upperBound > upper) {
      int i;
      for (i = size() - 1; i >= 0; i--) {
	if (_data[i] != 0) {
	  break;
	}
      }
      if (i != -1 && (double)_lowerBound + ((double)i + 1.0) *(double) _width > upper) {
	upper =(double)_lowerBound + ((double)i + 1.0) * (double)_width;
      }
    }
    // Rebuild with the new lower and upper bounds.
    rebuild(lower, upper, _width);
  }

  //! Rebuild the histogram so it covers the specified range and has at least the specified minimum width.
  void rebuild( _populationValueType low,  double high, double newWidth) {
#ifndef DEBUG_StochKit
    if(low < 0 || low >= high){
      Rcpp::Rcout<<"StochKit ERROR (HistogramSingle::rebuild): invalid population counts\n";
      Rcpp::stop("Fatal error encountered, terminating StochKit2R");
    }
#endif
    // Determine the new bounds and a bin width.
    // Note that the width is only allowed to grow.
    _populationValueType newLowerBound = std::floor(low / newWidth) * newWidth;
    double newUpperBound = newLowerBound + size() * newWidth;
    while (high > newUpperBound) {
      newWidth *= 2;
      newLowerBound = std::floor(low / newWidth) * newWidth;
      newUpperBound = newLowerBound + size() * newWidth;
    }
    // Rebuild the histogram.
    double newInverseWidth = 1. / newWidth;
    std::vector<int> newData(size(), 0);
    for (std::size_t i = 0; i != size(); i++) {
      if (_data[i] != 0) {
	_populationValueType event = _lowerBound + i * _width;

#ifndef DEBUG_StochKit
	if(newLowerBound>event || event>= newUpperBound){
	  Rcpp::Rcout<<"StochKit ERROR (HistogramSingle::rebuild): invalid new upper and/or lower bound\n";
	  Rcpp::stop("Fatal error encountered, terminating StochKit2R");
	}
#endif	 
	newData[std::size_t((event - newLowerBound) * newInverseWidth)] += _data[i];
      }
    }
    for (std::size_t i = 0; i != size(); i++) 
      _data[i] = newData[i];

    // New bounds and width.
    _lowerBound = newLowerBound;
    _width = newWidth;
    _inverseWidth = newInverseWidth;
    _upperBound = newUpperBound;
  }


  void mergeHistogram( HistogramSingle<_populationValueType> x){
#ifndef DEBUG_StochKit
    if (size() != x.size()){
      Rcpp::Rcout<<"StochKit ERROR (HistogramSingle::mergeHistogram): bin size of two histograms must be equal\n";
      Rcpp::stop("Fatal error encountered, terminating StochKit2R");
    }
#endif   
    // check to make sure time and species index are identical
    if (x._getSpeciesIndex()!=_speciesIndex){
      Rcpp::Rcout<<"StochKit ERROR (HistogramSingle::mergeHistogram): species indices must match\n";
      Rcpp::stop("Fatal error encountered, terminating StochKit2R");
    }
    if (x._getTimeIndex()!= _timeIndex){
     Rcpp::Rcout<<"StochKit ERROR (HistogramSingle::mergeHistogram): time indices must match\n";
      Rcpp::stop("Fatal error encountered, terminating StochKit2R");
    }

 
    double lower, upper; 
    if (computeSum()==0){
      
      lower = x.getLowerBound();
      upper = x.getUpperBound();
      
    }
    else if (x.computeSum()==0){
      lower = _lowerBound;
      upper = _upperBound;
    }
    else {
      lower = std::min(computeMinimumNonzero(), x.computeMinimumNonzero());
      upper = std::max(computeMaximumNonzero(), x.computeMaximumNonzero());
    }
    
    double width = std::max(getWidth(), x.getWidth());
    rebuild(lower, upper, width);
    x.rebuild(lower, upper, width);
    if(x.getLowerBound()!=_lowerBound){
      Rcpp::Rcout << "StochKit ERROR (HistogramSingle::mergeHistogram): lower bounds do not match after rebuild.\n";
      Rcpp::stop("Fatal error encountered, terminating StochKit2R");
    }
    if(x.getWidth()!=_width){
      Rcpp::Rcout << "StochKit ERROR (HistogramSingle::mergeHistogram): lower bounds do not match after rebuild.\n";
      Rcpp::stop("Fatal error encountered, terminating StochKit2R");
    }
    // merge data
    for(std::size_t i=0; i<_size;i++)
      _data[i] += x.getCounts(i);
  }   

  void writeToFile(std::string filename, std::vector<double>& outputTimes) {
    std::ofstream outfile;
    //   std::size_t numberOfProcesses=commandLine.getProcesses();
    //   append the processor at the end for parallel simulation
    /*
    std::ostringstream buffer;
    buffer << _timeIndex;
    std::ostringstream buffer2;
    buffer2 << _speciesIndex;
    */

    //std::string fname = filename+"-"+buffer.str()+"-"+buffer2.str();

    outfile.open(filename.c_str());
    if (!outfile) {
      Rcpp::Rcout << "StochKit ERROR (HistogramSingle::writeToFile): Unable to open output file \"" << filename << "\".\n";
      Rcpp::stop("Fatal error encountered, terminating StochKit2R");
    }
    try {
      // time index and species index 
      outfile << _speciesIndex << "\t"<<outputTimes[_timeIndex] << "\t"<< _speciesIndex <<"\t"<<_timeIndex<<"\n";
      // other necessary information
      outfile << _lowerBound<<"\t"<<_upperBound<<"\t"<<_width<<"\t"<<_size<<"\t"<<_inverseWidth<<"\n";

      for (std::size_t i=0; i<_data.size(); i++) {
	  outfile << _data[i] << "\t";
      }
      outfile << "\n";
      outfile.close();
    }
    catch (...) {
      Rcpp::Rcout << "StochKit ERROR (HistogramSingle::writeToFile): error writing data to output file.\n";
      Rcpp::stop("Fatal error encountered, terminating StochKit2R");
    } 
  }
  
  void writeToFile(std::string filename, std::vector<double>& outputTimes, std::string speciesID) {
    std::ofstream outfile;

    outfile.open(filename.c_str());
    if (!outfile) {
      Rcpp::Rcout << "StochKit ERROR (HistogramSingle::writeToFile): Unable to open output file \"" << filename << "\".\n";
      Rcpp::stop("Fatal error encountered, terminating StochKit2R");
    }
    try {
      // time index and species index 
      outfile << speciesID <<"\t"<< outputTimes[_timeIndex] << "\t"<< _speciesIndex <<"\t"<<_timeIndex<<"\n";
      // other necessary information
      outfile << _lowerBound<<"\t"<<_upperBound<<"\t"<<_width<<"\t"<<_size<<"\t"<<_inverseWidth<<"\n";

      for (std::size_t i=0; i<_data.size(); i++) {
	  outfile << _data[i] << "\t";
      }
      outfile << "\n";
      outfile.close();
    }
    catch (...) {
      Rcpp::Rcout << "StochKit ERROR (HistogramSingle::writeToFile): error writing data to output file.\n";
      Rcpp::stop("Fatal error encountered, terminating StochKit2R");
    } 
  }

     //StochKit2R: use index start of 1 instead of 0!
	 std::vector<std::string> fileDataAsString(std::vector<double> outputTimes, std::string speciesID) const {
		 std::vector<std::string> outdat;//(1);
			 // time index and species index
		 std::string str = speciesID +"\t"+ to_string_macro(outputTimes[_timeIndex]) + "\t"+ to_string_macro(_speciesIndex+1) +"\t"+to_string_macro(_timeIndex+1);
		 outdat.push_back(str);
			 // other necessary information
		 str =  to_string_macro(_lowerBound)+"\t"+to_string_macro(_upperBound)+"\t"+to_string_macro(_width)+"\t"+to_string_macro(_size)+"\t"+to_string_macro(_inverseWidth);
		 outdat.push_back(str);
		 str = "";
		 for (std::size_t i=0; i<_data.size(); i++) {
			 str+= to_string_macro(_data[i])+"\t";
		 }
		 outdat.push_back(str);
		 return outdat;
	 }


 };
}

#endif
