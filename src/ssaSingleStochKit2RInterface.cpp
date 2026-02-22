// [[Rcpp::depends(BH)]]

#include <iostream>
#include <fstream>
#include <string>
#include "boost/numeric/ublas/vector_sparse.hpp"
#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/matrix_sparse.hpp"
#include "boost/numeric/ublas/matrix_proxy.hpp"
#include "CustomPropensitySet.h"
#include "SSA_Direct.h"
#include "StandardDriverUtilities.h"
#include <boost/random/mersenne_twister.hpp>
#include "MassActionModel.h"
#include "solver_helper_functions.h"

#include <Rcpp.h>

//'@title C++ Interface to Gillespie Stochastic Simulation Algorithm single trajectory
//'
//'@description
//'\code{ssa} Called by StochKit2R ssaSingle function, do not call this C++ interface directly
//'
//'@param StochKit2Rmodel R list (Rcpp List built from buildStochKit2Rmodel output)
//'@param startTime Simulation start time
//'@param endTime Simulation end time
//'@param outputFileNameString Character string with path to output file.
//'@param seed Seed for random number generator
//'@return Dataframe containing the time and population sizes
//'@keywords internal
// [[Rcpp::export]]
RcppExport SEXP ssaSingleStochKit2RInterface(Rcpp::List StochKit2Rmodel,
                                             double startTime, double endTime, std::string outputFileNameString, unsigned int seed)
{
  //create StochKit2R mass action model object
  //first, pull out pieces from list object
  Rcpp::List rParameterList=StochKit2Rmodel[0];
  Rcpp::List rSpeciesList=StochKit2Rmodel[1];
  Rcpp::List rReactionList=StochKit2Rmodel[2];
  Rcpp::List rCustomPropensityList=StochKit2Rmodel[3];
	
  //get species labels...
  std::vector<std::string> modelSpeciesList = getSpeciesList(rSpeciesList);

  STOCHKIT::MassActionModel<STOCHKIT::StandardDriverTypes::populationType,
                            STOCHKIT::StandardDriverTypes::denseStoichiometryType,
                            STOCHKIT::StandardDriverTypes::propensitiesType,
                            STOCHKIT::StandardDriverTypes::graphType> model(rParameterList,
                                                                            rReactionList,
                                                                            rSpeciesList,rCustomPropensityList);

  //create the output object
  std::vector< std::pair<double, STOCHKIT::StandardDriverTypes::populationType> > output;

  //use StochKit2 seeding (and parallel RNG) strategy
  std::vector<unsigned int> seeds;
  boost::mt19937 seedGenerator;
  seedGenerator.seed(seed);
  seeds.push_back(seedGenerator());//thread 0

  STOCHKIT::SSA_Direct<STOCHKIT::StandardDriverTypes::populationType,
                       STOCHKIT::StandardDriverTypes::denseStoichiometryType,
                       STOCHKIT::StandardDriverTypes::propensitiesType,
                       STOCHKIT::StandardDriverTypes::graphType> ssa(model.writeInitialPopulation(),
                                                                     model.writeStoichiometry(),
                                                                     model.writePropensities(),
                                                                     model.writeDependencyGraph(),
                                                                     seeds[0]);

  ssa.simulateSingle(startTime, endTime, output);

  // construct vector of proper row size
  std::vector<double> row(rSpeciesList.size()+1); 
  // this 2D double array will get converted to a dataframe
  std::vector< std::vector<double> > buffer; 
  
  //if user entered an output directory then write file and fill data buffer
  if(outputFileNameString != ""){
    //try to write labels to file
    STOCHKIT::IntervalOutput<STOCHKIT::StandardDriverTypes::populationType>::writeLabelsToFile(outputFileNameString, modelSpeciesList);
    // write output
    std::ofstream outfile;

    //open for appending
    outfile.open(outputFileNameString.c_str(), std::ios::out | std::ios::app);

    if (!outfile) {
      Rcpp::Rcout << "StochKit ERROR (ssaSingleStochKit2RInterface): Unable to open output file.\n";
      Rcpp::stop("Fatal error encountered, terminating StochKit2R");
    }

    try {
      int i = 0;
      for (std::size_t step=0; step!=output.size(); ++step) {
        i = 0;
        row[i++] = output[step].first;

        //write time
        outfile << output[step].first << "\t";

        for (size_t index=0; index!=output[step].second.size(); ++index) {
          row[i++] = output[step].second[index];
          outfile << std::setprecision(8) << output[step].second[index] << "\t";
        }
        outfile << "\n";
        buffer.push_back(row);
      }
      outfile.close();
    }
    catch (...) {
      Rcpp::Rcout << "StochKit ERROR (IntervalOutput::writeDataToFile): error writing data to output file.\n";
      Rcpp::stop("Fatal error encountered, terminating StochKit2R");
    }
  }

  // to fill data buffer without writing file
  else{
    int i=0;
    for(std::size_t step=0; step!=output.size(); ++step){
      i=0;
      row[i++] = output[step].first;
      for (size_t index=0; index!=output[step].second.size(); ++index) {
            row[i++] = output[step].second[index];
          }
        buffer.push_back(row);
    }
  }

  // write buffer to return data 
  int nr = buffer.size();
  long long unsigned int nc = buffer[0].size();
  Rcpp::NumericMatrix m( nr, nc );
  for (int i=0; i<nr; i++) {
    std::vector<double>& result_i = buffer[i];

    if (result_i.size() != nc) {
      //the length of data in a row is inconsistent, will not be able to convert to a data frame.
      Rcpp::stop("Fatal error encountered, terminating StochKit2R");
    }

    for (long long unsigned int j=0; j<nc; j++)
      m(i, j) = result_i[j];
  }


  // name the columns accordingly
  Rcpp::CharacterVector col_names;
  col_names.push_back("time");

  for (long long unsigned int i = 0; i < modelSpeciesList.size(); i++)
    col_names.push_back(modelSpeciesList[i]);

  Rcpp::DataFrame df = Rcpp::DataFrame(m);
  df.attr("names")= col_names;

  

  return df;
}
