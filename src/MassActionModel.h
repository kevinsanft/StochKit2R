/*!
 \brief MassActionModel model class (Note: no longer restricted to mass-action!)
*/

#ifndef STOCHKIT_MASS_ACTION_MODEL_H
#define STOCHKIT_MASS_ACTION_MODEL_H

#define MAXPARAMETERNUM 200 // maximum parameter number there could be
#define BADRESULT -32678 // indicate bad result such as divisor to be 0 or something in calculation

//#include <iostream>
//#include <sstream>
//#include <algorithm>
//#include <stdio.h>
//#include <string.h>
//#include <stdlib.h>
//#include <libxml/xmlmemory.h>
//#include <libxml/parser.h>
//#include "boost/shared_ptr.hpp"
//#include <vector>
//#include <limits>
#include "StandardDriverTypes.h"
#include "Parameter.h"
#include "StringCalculator.h"
//#include "VectorManipulation.h"
#include "CustomPropensity.h"
//#include "CustomSimplePropensity.h"
#include "CustomPropensitySet.h"
#include <Rcpp.h>

namespace STOCHKIT
{
 template<typename _populationVectorType,
	typename _stoichiometryType,
	typename _propensitiesFunctorType,
	typename _dependencyGraphType>
 class MassActionModel
 {
 protected:
	typedef typename _populationVectorType::value_type _populationValueType;

	int NumberOfReactions, NumberOfSpecies;

	int NumberOfParameters;

	//! class to store parameters and related information
	ListOfParameters ParametersList;
	
	//! class to store species and related information
	class Species{
		public:
			std::string Id; // Id of the species
			std::string InitialPopulation;  // initial population, only valid in SpeciesList
			std::vector<int> AffectReactions; // only valid in SpeciesList, record what reactions does one Species affect
	};
	std::vector<Species> SpeciesList;

	class SpeciesReference{
		public:
			std::string Id;
			int Index;
			int Stoichiometry;
	};
	
	//! class to store reactions and related information
	class Reaction{
		public:
			std::string Id;  // Id of the reaction
			int Type; // type of the reaction, 0 = mass-action, 1 = michaelis-menten, 2 = customized
			std::string Rate; // for mass-action
			std::string Vmax;  // for michealis-menten
			std::string Km;  // for michealis-menten
			//std::string Customized; // for customized propensity
		SEXP Customized;
			std::vector<SpeciesReference> Reactants;  // reactant list
			std::vector<SpeciesReference> Products;  // product list
	};
	std::vector<Reaction> ReactionsList;

	//! class to handle calculation of simple math expression strings
	StringCalculator simpleCalculator;

 public:
	 
	 MassActionModel(Rcpp::List rParametersList, Rcpp::List rReactionList, Rcpp::List rSpeciesList, Rcpp::List rCustomPropensityList=NULL) : simpleCalculator()
	{
		NumberOfReactions = 0;
		NumberOfSpecies = 0;
		NumberOfParameters = 0;
		SpeciesList.clear();
		ReactionsList.clear();
        NumberOfSpecies=rSpeciesList.size();
        NumberOfReactions=rReactionList.size();
        recordParametersList(rParametersList);
        recordReactionsList(rReactionList,rCustomPropensityList);
        recordSpeciesList(rSpeciesList);
        
        if(!ParametersList.linkParameters()) {
            //unreachable
            Rcpp::stop("Unexpected error in linkParameters function. Please report bug");
		}

		if(!ParametersList.calculateParameters()){
            Rcpp::stop("Fatal error encountered, terminating StochKit2R");
		}

		if(!linkSpeciesAndReactions()){
			Rcpp::stop("Fatal error encountered, terminating StochKit2R");
		}

		if(!checkUniqueID()){
			Rcpp::stop("Fatal error encountered, terminating StochKit2R");
		}
	}

 protected:
     bool recordParametersList(Rcpp::List rParameterList);
     bool recordReactionsList(Rcpp::List rReactionList, Rcpp::List rCustomPropensityList);
     bool recordSpeciesList(Rcpp::List rSpeciesList);
     bool linkSpeciesAndReactions();

	_populationValueType populationCalculation(std::string equation);

	//! check to see if all IDs are unique, return true if they are, return false if there are duplicate IDs
	bool checkUniqueID();
     
     double rateCalculation(std::string equation);
     
public:
     _populationVectorType writeInitialPopulation();
     _stoichiometryType writeStoichiometry();
     _stoichiometryType writeMatrixStoichiometry();
     _propensitiesFunctorType writePropensities();
     _dependencyGraphType writeDependencyGraph();
     _dependencyGraphType writeDependencyGraphMatrixStoichiometry();
	 
 }; // end of MassActionModel class definition

//include function definitions
template<typename _populationVectorType,
typename _stoichiometryType,
typename _propensitiesFunctorType,
typename _dependencyGraphType>
bool
MassActionModel<_populationVectorType,
_stoichiometryType,
_propensitiesFunctorType,
_dependencyGraphType>::
    recordParametersList(Rcpp::List rParameterList)
{
    for (int i=0; i!=rParameterList.size(); ++i) {
        Rcpp::CharacterVector thisParam = rParameterList[i];
        std::string pId = Rcpp::as<std::string>(thisParam[0]);
        std::string pExpression = Rcpp::as<std::string>(thisParam[1]);
        //Rcout << "parameter " << pId << " has expression " << pExpression << "\n";
        
        ParametersList.push_back(Parameter());
        ParametersList.back().Id=pId;
        ParametersList.back().Expression=pExpression;
        ParametersList.back().Type=0;
        ParametersList.back().CalculateFlag=-1;
    }
    
    return true;
}

template<typename _populationVectorType,
typename _stoichiometryType,
typename _propensitiesFunctorType,
typename _dependencyGraphType>
bool
MassActionModel<_populationVectorType,
_stoichiometryType,
_propensitiesFunctorType,
_dependencyGraphType>::
    recordReactionsList(Rcpp::List rReactionList, Rcpp::List rCustomPropensityList)
{

	//Rcpp::Rcout << "rReactionList.size()="<<rReactionList.size()<<"\n";


    for (int i=0; i!=rReactionList.size(); ++i) {
        //contents: List of "reactions"
		//each reaction contains Id, Type, Rate, PropensityFunction,Reactants, Products
        //where Reactants and Products are lists containing...
        Rcpp::List thisReaction = rReactionList[i];
		//Rcpp::Rcout << "thisReaction.size()="<<thisReaction.size()<<"\n";
		
		
        std::string rId = Rcpp::as<std::string>(thisReaction[0]);
        std::string rRate = Rcpp::as<std::string>(thisReaction[2]);
		//Rcpp::Rcout << "reaction " << i << " has Id " << rId << ", Rate=" << rRate << std::endl;
//        Rcpp::List reactants = thisReaction[2];
//        Rcpp::List products = thisReaction[3];
		Rcpp::List reactants = thisReaction[4];
		Rcpp::List products = thisReaction[5];

		
        ReactionsList.push_back(Reaction());
        ReactionsList.back().Id=rId;
		
		std::string rxnTypeString = thisReaction[1];
		if (rxnTypeString.compare("mass-action")==0) {
			ReactionsList.back().Type=0;
			ReactionsList.back().Rate=rRate;
		}
		else if (rxnTypeString.compare("customized")==0) {
			ReactionsList.back().Type=2;
			ReactionsList.back().Customized=(SEXP)rCustomPropensityList[i];
		}
		else {
			Rcpp::stop("Detected invalid model reaction Type. Terminating.");
		}
		
        //process Reactants
        int reactionOrder=0;//counter to check for >3rd order reaction
        for (int j=0; j!=reactants.size(); ++j) {
            Rcpp::List reactant=reactants[j];
            std::string id=reactant[0];
            int stoich=reactant[1];
            //Rcpp::Rcout << "  reactant " << j << " id="<<id<<", stoichiometry="<<stoich<<"\n";
            ReactionsList.back().Reactants.push_back(SpeciesReference());
            ReactionsList.back().Reactants.back().Id=id;
            ReactionsList.back().Reactants.back().Stoichiometry=-stoich;//negative for Reactants!
            reactionOrder+=stoich;
        }
        if (reactionOrder>3) {
            Rcpp::Rcout << "ERROR: reaction order >3 for mass-action model\n";
            Rcpp::stop("Fatal error encountered, terminating StochKit2R");
        }
        
        //process Products
        for (int j=0; j!=products.size(); ++j) {
            Rcpp::List product=products[j];
            std::string id=product[0];
            int stoich=product[1];
            //Rcpp::Rcout << "  product " << j << " id="<<id<<", stoichiometry="<<stoich<<"\n";
            ReactionsList.back().Products.push_back(SpeciesReference());
            ReactionsList.back().Products.back().Id=id;
            ReactionsList.back().Products.back().Stoichiometry=stoich;
        }
    }

    return true;
}
    
template<typename _populationVectorType,
typename _stoichiometryType,
typename _propensitiesFunctorType,
typename _dependencyGraphType>
bool
MassActionModel<_populationVectorType,
_stoichiometryType,
_propensitiesFunctorType,
_dependencyGraphType>::
    recordSpeciesList(Rcpp::List rSpeciesList)
{
    for (int i=0; i!=rSpeciesList.size(); ++i) {
        Rcpp::CharacterVector thisSpec = rSpeciesList[i];
        std::string sId = Rcpp::as<std::string>(thisSpec[0]);
        std::string sIP = Rcpp::as<std::string>(thisSpec[1]);
        //Rcpp::Rcout << "species " << sId << " has InitialPopulation = " << sIP << "\n";
        SpeciesList.push_back(Species());
        SpeciesList.back().Id=sId;
        SpeciesList.back().InitialPopulation=sIP;
    }
    
    return true;
}

template<typename _populationVectorType,
typename _stoichiometryType,
typename _propensitiesFunctorType,
typename _dependencyGraphType>
bool
MassActionModel<_populationVectorType,
_stoichiometryType,
_propensitiesFunctorType,
_dependencyGraphType>::linkSpeciesAndReactions()
{
    class Reaction *cur_reaction;
    class SpeciesReference *cur_reactant, *cur_product;
    int flag = 0;
    
    for( unsigned int i=0; i<ReactionsList.size(); ++i){
        cur_reaction = &ReactionsList[i];
        
        if(cur_reaction->Reactants.size() != 0){
            for( unsigned int j=0; j<cur_reaction->Reactants.size(); ++j){
                cur_reactant = &cur_reaction->Reactants[j];
                flag = 0;
                
                unsigned int k=0;
                while( (k<SpeciesList.size()) && (flag==0) ){
                    if(cur_reactant->Id.compare(SpeciesList[k].Id) == 0){
                        flag = 1;
                        //if(cur_reaction->Type == 0 || cur_reaction->Type == 1){
                            SpeciesList[k].AffectReactions.push_back(i);
                        //}
                        cur_reactant->Index = k;
                    }
                    ++k;
                }
                
                if( flag == 0 ){
                    Rcpp::Rcout << "StochKit ERROR (Input::linkSpeciesAndReactions): reactant " << cur_reactant->Id << " in reaction " << cur_reaction->Id << " does not exist in SpeciesList\n";
                    return false;
                }
            }
        }

        
        if(cur_reaction->Products.size() != 0){
            for( unsigned int j=0; j<cur_reaction->Products.size(); ++j){
                cur_product = &cur_reaction->Products[j];
                flag = 0;
                
                unsigned int k=0;
                while( (k<SpeciesList.size()) && (flag==0) ){
                    if(cur_product->Id.compare(SpeciesList[k].Id) == 0){
                        flag = 1;
                        cur_product->Index = k;
                    }
                    ++k;
                }
                
                if( flag == 0 ){
                    Rcpp::Rcout << "StochKit ERROR (Input::linkSpeciesAndReactions): product " << cur_product->Id << " in reaction " << cur_reaction->Id << " does not exist in SpeciesList\n";
                    return false;
                }
            }
        }
    }
    
    return true;
}

template<typename _populationVectorType,
typename _stoichiometryType,
typename _propensitiesFunctorType,
typename _dependencyGraphType>
bool
MassActionModel<_populationVectorType,
_stoichiometryType,
_propensitiesFunctorType,
_dependencyGraphType>::checkUniqueID()
{
    for( int i=0; i< NumberOfSpecies; ++i ){
        for( int j=i+1; j<NumberOfSpecies; ++j){
            if( SpeciesList[i].Id.compare(SpeciesList[j].Id) == 0 ){
                Rcpp::Rcout << "StochKit ERROR (Input::checkUniqueID): there are two species having the same ID \"" << SpeciesList[i].Id << "\"\n";
                return false;
            }
        }
        
        for( int j=0; j<NumberOfParameters; ++j){
            if( SpeciesList[i].Id.compare(ParametersList[j].Id) == 0 ){
                Rcpp::Rcout << "StochKit ERROR (Input::checkUniqueID): there are one species and one parameter having the same ID \"" << SpeciesList[i].Id << "\"\n";
                return false;
            }
        }
    }
    
    for( int i=0; i< NumberOfParameters; ++i ){
        for( int j=i+1; j<NumberOfParameters; ++j){
            if( ParametersList[i].Id.compare(ParametersList[j].Id) == 0 ){
                Rcpp::Rcout << "StochKit ERROR (Input::checkUniqueID): there are two parameters having the same ID \"" << ParametersList[i].Id << "\"\n";
                return false;
            }
        }
    }
    
    return true;
}

// calculate rate based on the value stored in parameterslist
template<typename _populationVectorType,
typename _stoichiometryType,
typename _propensitiesFunctorType,
typename _dependencyGraphType>
typename
MassActionModel<_populationVectorType,
_stoichiometryType,
_propensitiesFunctorType,
_dependencyGraphType>::
_populationValueType
MassActionModel<_populationVectorType,
_stoichiometryType,
_propensitiesFunctorType,
_dependencyGraphType>::
populationCalculation(std::string equation)
{
    std::vector<unsigned int> ParametersAffectRate;
    std::vector<unsigned int>::iterator para_it; // iterator of parameters in link graph
    
    ParametersAffectRate = ParametersList.analyzeParameterExpression(equation);
    
    bool calculationStatus = false;
    
    for( para_it = ParametersAffectRate.begin(); para_it < ParametersAffectRate.end(); ++para_it ){
        if( ParametersList[*para_it].CalculateFlag == -1 ){
            calculationStatus = ParametersList.calculateParameter(*para_it);
            if(!calculationStatus){
                Rcpp::Rcout << "StochKit ERROR (Input::populationCalculation): while calculating rate " << equation << std::endl;
                return BADRESULT;
            }
        }
    }
    
    std::string substitutedEquation = ParametersList.parameterSubstitution(equation);
    if( substitutedEquation.empty() ){
        Rcpp::Rcout << "StochKit ERROR (Input::populationCalculation): while calculating initial population " << equation << std::endl;
        return BADRESULT;
    }
    
    double calculatedPopulation = simpleCalculator.calculateString(substitutedEquation);
    double roundedPopulation = floor(calculatedPopulation+0.5);
    if( fabs(roundedPopulation-calculatedPopulation) > 1e-7 ){
        Rcpp::Rcout.precision(10);
        Rcpp::Rcout << "StochKit WARNING (Input::populationCalculation): population was rounded from " << calculatedPopulation << " to " << roundedPopulation << std::endl;
    }
    
    return (_populationValueType)roundedPopulation;
}
    
template<typename _populationVectorType,
typename _stoichiometryType,
typename _propensitiesFunctorType,
typename _dependencyGraphType>
_populationVectorType
MassActionModel<_populationVectorType,
_stoichiometryType,
_propensitiesFunctorType,
_dependencyGraphType>::
writeInitialPopulation()
{
    _populationValueType cur_population;
    _populationVectorType X(NumberOfSpecies);
    
    for(int i=0; i<NumberOfSpecies; ++i){
        cur_population = populationCalculation(SpeciesList[i].InitialPopulation);
        if( cur_population == BADRESULT ){
            Rcpp::Rcout << "StochKit ERROR (Input::writeInitialPopulation): while calculating initial population of " << SpeciesList[i].Id << std::endl;
            Rcpp::stop("Fatal error encountered. Terminating StochKit2R");
        }
        
        X[i] = cur_population;
    }
    
    return X;
}

template<typename _populationVectorType,
typename _stoichiometryType,
typename _propensitiesFunctorType,
typename _dependencyGraphType>
_stoichiometryType
MassActionModel<_populationVectorType,
_stoichiometryType,
_propensitiesFunctorType,
_dependencyGraphType>::
writeStoichiometry()
{
    typedef typename _stoichiometryType::value_type vectorType;
    _stoichiometryType nu(NumberOfReactions, vectorType(NumberOfSpecies));
    // initialize values to 0
    for(unsigned int i=0; i<nu.size(); ++i) {
        for(unsigned int j=0; j<nu[i].size(); ++j) {
            if( nu[i][j] != 0 ){
                nu[i][j] = 0;
            }
        }
    }
    
    //calculate Stoichiometry and write
    class SpeciesReference *cur_species = NULL;
    
    for(int i=0; i<NumberOfReactions; ++i){ // the i-th reaction
        // write reactants stoichiometry
        for(unsigned int j=0; j<ReactionsList[i].Reactants.size(); ++j){
            cur_species = &ReactionsList[i].Reactants[j];
            nu[i][cur_species->Index] += cur_species->Stoichiometry;
        }
        
        // write products stoichiometry
        for(unsigned int j=0; j<ReactionsList[i].Products.size(); ++j){
            cur_species = &ReactionsList[i].Products[j];
            nu[i][cur_species->Index] += cur_species->Stoichiometry;
        }
    }
	
//	Rcpp::Rcout << "writeStoichiometry:"<<std::endl;
//	for(unsigned int i=0; i<nu.size(); ++i) {
//		for(unsigned int j=0; j<nu[i].size(); ++j) {
//			Rcpp::Rcout << nu[i][j] <<"\t";
//		}
//		Rcpp::Rcout <<"\n";
//	}

    return nu;
}

//specifically for matrix stoichiometry type
template<typename _populationVectorType,
typename _stoichiometryType,
typename _propensitiesFunctorType,
typename _dependencyGraphType>
_stoichiometryType
MassActionModel<_populationVectorType,
_stoichiometryType,
_propensitiesFunctorType,
_dependencyGraphType>::
writeMatrixStoichiometry()
{
    _stoichiometryType nu(NumberOfReactions, NumberOfSpecies);
    
    //calculate Stoichiometry and write
    class SpeciesReference *cur_species = NULL;
    
    for(int i=0; i<NumberOfReactions; ++i){ // the i-th reaction
        // write reactants stoichiometry
        for(unsigned int j=0; j<ReactionsList[i].Reactants.size(); ++j){
            cur_species = &ReactionsList[i].Reactants[j];
            nu(i,cur_species->Index) += cur_species->Stoichiometry;
        }
        
        // write products stoichiometry
        for(unsigned int j=0; j<ReactionsList[i].Products.size(); ++j){
            cur_species = &ReactionsList[i].Products[j];
            nu(i,cur_species->Index) += cur_species->Stoichiometry;
        }
    }
    
    return nu;
}

template<typename _populationVectorType,
typename _stoichiometryType,
typename _propensitiesFunctorType,
typename _dependencyGraphType>
_propensitiesFunctorType
MassActionModel<_populationVectorType,
_stoichiometryType,
_propensitiesFunctorType,
_dependencyGraphType>::writePropensities()
{
    _propensitiesFunctorType propensitiesList;
    double rate;
    
    typename MassActionModel<_populationVectorType, _stoichiometryType, _propensitiesFunctorType, _dependencyGraphType>::Reaction *cur_reaction;
    std::vector<int> reactantsList;
    
    for(int i=0; i<this->NumberOfReactions; ++i){
        cur_reaction = &this->ReactionsList[i];
        
        if(cur_reaction->Type == 0){
            rate = rateCalculation(cur_reaction->Rate);
            if( rate == BADRESULT ){
                Rcpp::Rcout << "StochKit ERROR (Input_mass_action::writePropensities): while calculating rate of reaction " << cur_reaction->Id << "\n";
                Rcpp::stop("Fatal error encountered, terminating StochKit2R");
            }
            switch ( cur_reaction->Reactants.size() ){
                case 0:
                    propensitiesList.pushSimplePropensity(rate);
                    break;
                case 1:
                    if( cur_reaction->Reactants[0].Stoichiometry == -1 ) {
                        propensitiesList.pushSimplePropensity(rate, cur_reaction->Reactants[0].Index);
                    }
                    else if( cur_reaction->Reactants[0].Stoichiometry == -2 )
                        propensitiesList.pushSimplePropensity(rate, cur_reaction->Reactants[0].Index, cur_reaction->Reactants[0].Index);
                    else if( cur_reaction->Reactants[0].Stoichiometry == -3 )
                        propensitiesList.pushSimplePropensity(rate, cur_reaction->Reactants[0].Index, cur_reaction->Reactants[0].Index, cur_reaction->Reactants[0].Index);
                    else{
                        Rcpp::Rcout << "StochKit ERROR (Input_mass_action::writePropensities): currently the highest order mass-action reaction supported is tri-molecular reaction while Reaction " << cur_reaction->Id << " is not\n";
                        Rcpp::stop("Fatal error encountered, terminating StochKit2R");
                    }
                    break;
                case 2:
                    if( cur_reaction->Reactants[0].Stoichiometry == -1 ){
                        if (cur_reaction->Reactants[1].Stoichiometry == -1) {
                            propensitiesList.pushSimplePropensity(rate, cur_reaction->Reactants[0].Index, cur_reaction->Reactants[1].Index);
                        }
                        else if (cur_reaction->Reactants[1].Stoichiometry == -2)
                            propensitiesList.pushSimplePropensity(rate, cur_reaction->Reactants[0].Index, cur_reaction->Reactants[1].Index, cur_reaction->Reactants[1].Index);
                        else{
                            Rcpp::Rcout << "StochKit ERROR (Input_mass_action::writePropensities): currently the highest order mass-action reaction supported is tri-molecular reaction while Reaction " << cur_reaction->Id << " is not\n";
                            Rcpp::stop("Fatal error encountered, terminating StochKit2R");
                        }
                    }
                    else if(cur_reaction->Reactants[0].Stoichiometry == -2){
                        if (cur_reaction->Reactants[1].Stoichiometry == -1)
                            propensitiesList.pushSimplePropensity(rate, cur_reaction->Reactants[0].Index, cur_reaction->Reactants[0].Index, cur_reaction->Reactants[1].Index);
                        else{
                            Rcpp::Rcout << "StochKit ERROR (Input_mass_action::writePropensities): currently the highest order mass-action reaction supported is tri-molecular reaction while Reaction " << cur_reaction->Id << " is not\n";
                            Rcpp::stop("Fatal error encountered, terminating StochKit2R");
                        }
                    }
                    else{
                        Rcpp::Rcout << "StochKit ERROR (Input_mass_action::writePropensities): currently the highest order mass-action reaction supported is tri-molecular reaction while Reaction " << cur_reaction->Id << " is not\n";
                        Rcpp::stop("Fatal error encountered, terminating StochKit2R");
                    }
                    break;
                case 3:
                    if( cur_reaction->Reactants[0].Stoichiometry != -1 || cur_reaction->Reactants[1].Stoichiometry != -1 || cur_reaction->Reactants[2].Stoichiometry != -1 ){
                        Rcpp::Rcout << "StochKit ERROR (Input_mass_action::writePropensities): currently the highest order mass-action reaction supported is tri-molecular reaction while Reaction " << cur_reaction->Id << " is not\n";
                        Rcpp::stop("Fatal error encountered, terminating StochKit2R");
                    }
                    else
                        propensitiesList.pushSimplePropensity(rate, cur_reaction->Reactants[0].Index, cur_reaction->Reactants[1].Index, cur_reaction->Reactants[2].Index);
                    
                    break;
                default:
                    Rcpp::Rcout << "StochKit ERROR (Input_mass_action::writePropensities): more than 3 reactants in mass-action reaction " << cur_reaction->Id << "\n";
                    Rcpp::stop("Fatal error encountered, terminating StochKit2R");
            }
        }
        else if(cur_reaction->Type == 2){
			//
			propensitiesList.pushCustomPropensity(*Rcpp::XPtr<STOCHKIT::CustomPropensity<STOCHKIT::StandardDriverTypes::populationType>::customPropensityFunction>(cur_reaction->Customized));
		}
		else {
            Rcpp::Rcout << "StochKit ERROR (Input_mass_action::writePropensities): reaction " << cur_reaction->Id << " has invalid Type.\n";
            Rcpp::stop("Fatal error encountered, terminating StochKit2R");
        }
    }
    
    return propensitiesList;
}

template<typename _populationVectorType,
typename _stoichiometryType,
typename _propensitiesFunctorType,
typename _dependencyGraphType>
double MassActionModel<_populationVectorType,
    _stoichiometryType,
    _propensitiesFunctorType,
    _dependencyGraphType>::rateCalculation(std::string equation)
{
    std::vector<unsigned int> ParametersAffectRate;
    std::vector<unsigned int>::iterator para_it; // iterator of parameters in link graph
    
    ParametersAffectRate = this->ParametersList.analyzeParameterExpression(equation);
    
    bool calculationStatus = false;
    
    for( para_it = ParametersAffectRate.begin(); para_it < ParametersAffectRate.end(); ++para_it ){
        if( this->ParametersList[*para_it].CalculateFlag == -1 ){
            calculationStatus = this->ParametersList.calculateParameter(*para_it);
            if(!calculationStatus){
                Rcpp::Rcout << "StochKit ERROR (Input_mass_action::rateCalculation): while calculating rate " << equation << std::endl;
                return BADRESULT;
            }
        }
    }
    
    std::string substitutedEquation = this->ParametersList.parameterSubstitution(equation);
    if( substitutedEquation.empty() ){
        Rcpp::Rcout << "StochKit ERROR (Input_mass_action::rateCalculation): while calculating rate " << equation << std::endl;
        return BADRESULT;
    }
    
    return this->simpleCalculator.calculateString(substitutedEquation);
}

template<typename _populationVectorType,
typename _stoichiometryType,
typename _propensitiesFunctorType,
typename _dependencyGraphType>
_dependencyGraphType
MassActionModel<_populationVectorType,
_stoichiometryType,
_propensitiesFunctorType,
_dependencyGraphType>::
writeDependencyGraph()
{
    typedef typename _dependencyGraphType::value_type _vectorType;
    typedef typename _vectorType::value_type _memberType;
    
    _dependencyGraphType dg(NumberOfReactions);
    
    _stoichiometryType nu = writeStoichiometry();
    
    class SpeciesReference *cur_species = NULL;
    
    for(int i=0; i<NumberOfReactions; ++i){
        for(int j=0; j<NumberOfSpecies; ++j){
                if(nu[i][j] != 0)
                    mergeToSortedArray<std::vector<int>, std::vector<_memberType> >(SpeciesList[j].AffectReactions, dg[i]);
        }
        
        // add reactions i to species s's affectreactions list if species s is a reactant of reaction i, even if the stiochiometry is 0
        for(unsigned int j=0; j<ReactionsList[i].Reactants.size(); ++j){
            cur_species = &ReactionsList[i].Reactants[j];
            mergeToSortedArray<std::vector<int>, std::vector<_memberType> >(SpeciesList[cur_species->Index].AffectReactions, dg[i]);
        }
        
    }
	
//	Rcpp::Rcout << "writeDependencyGraph...\n";
//	for (int i=0; i<dg.size(); i++) {
//		for (int j=0; j<dg[i].size(); j++) {
//			Rcpp::Rcout << dg[i][j] << "\t";
//		}
//		Rcpp::Rcout << "\n";
//	}
    return dg;
}

template<typename _populationVectorType,
typename _stoichiometryType,
typename _propensitiesFunctorType,
typename _dependencyGraphType>
_dependencyGraphType
MassActionModel<_populationVectorType,
_stoichiometryType,
_propensitiesFunctorType,
_dependencyGraphType>::
writeDependencyGraphMatrixStoichiometry()
{
    typedef typename _dependencyGraphType::value_type _vectorType;
    typedef typename _vectorType::value_type _memberType;
    
    _dependencyGraphType dg(NumberOfReactions);
    
    _stoichiometryType nu = writeStoichiometry();
    
    class SpeciesReference *cur_species = NULL;
    
    for(int i=0; i<NumberOfReactions; ++i){
        for(int j=0; j<NumberOfSpecies; ++j){
            if(nu(i,j) != 0)
                    mergeToSortedArray<std::vector<int>, std::vector<_memberType> >(SpeciesList[j].AffectReactions, dg[i]);
        }
        
        // add reactions i to species s's affectreactions list if species s is a reactant of reaction i, even if the stiochiometry is 0
        for(unsigned int j=0; j<ReactionsList[i].Reactants.size(); ++j){
            cur_species = &ReactionsList[i].Reactants[j];
            mergeToSortedArray<std::vector<int>, std::vector<_memberType> >(SpeciesList[cur_species->Index].AffectReactions, dg[i]);
        }
        
    }
	
	
	
    return dg;
}

}// end namespace STOCHKIT
#endif

