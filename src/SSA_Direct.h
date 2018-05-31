/*!
\brief the direct method of the Stochastic Simulation Algorithm (SSA)
*/

#ifndef _SSA_DIRECT_H_
#define _SSA_DIRECT_H_

#include <iostream>
#include <vector>
#include <limits>
#include <ctime>
#include "Random.h"
#include "StandardDriverTypes.h"
#include <utility>

/*! 
\file SSA_Direct.h
\brief Gillespie's Stochastic Simulation Algorithm (SSA): direct method.
\param _populationVectorType the population vector type, should be dense
\param _stoichiometryType
\param _propensitiesFunctorType functor takes reaction index and _populationVectorType population
and returns the propensity for that reaction
\param _dependencyGraphType [expand]
*/
namespace STOCHKIT
{
	template<typename _populationVectorType, 
		typename _stoichiometryType,
		typename _propensitiesFunctorType,
		typename _dependencyGraphType>
	class SSA_Direct
	{
	public:	
		typedef _populationVectorType populationVectorType;
		typedef _stoichiometryType stoichiometryType;
		typedef _propensitiesFunctorType propensitiesType;
		typedef _dependencyGraphType dependencyGraphType;

#ifdef MATRIX_STOICHIOMETRY
		typedef StandardDriverTypes::stoichiometryRow matrixrow;
#endif

	protected:
		//! the class that implements all random number generator functions
		/*! change the RandomGenerator class to swap generators */
        Random randomGenerator;

		//! the initial population
		/*! currentPopulation should be set to initialPopulation at the beginning of each
		realization as in initialize()
		\see initialize()
		*/
		_populationVectorType initialPopulation;

		//! the stoichiometric matrix
		/*!
		should actually be a dense vector of (usually sparse) vectors of dimension NumberOfReactions x NumberOfSpecies
		so that currentPopulation+=stoichiometry[reactionNumber] modifies the population based on
		the stoichiometry of the given reaction number
		*/
		_stoichiometryType stoichiometry;

		//! the propensities functor
		/*! propensities(rxn, pop) returns the propensity of reaction number rxn based on population pop
		after a simulation step, currentPropensities[rxn] should equal propensities(rxn, currentPopulation)
		but since propensities() is a function call, accessing current propensities should be done with
		currentPropensities
		\see currentPropensities
		*/
		_propensitiesFunctorType propensities;
		//! the dependency graph which describes the propensities that are affected by each reaction
		/*!
		should be a dense vector (size=NumberOfReactions) of variable length vectors
		where dependencyGraph[rxn] returns the variable length vector of reaction indices that are affected by reaction rxn
		e.g. if dependencyGraph[4][0]=2 and dependencyGraph[4][1]=4 and dependencyGraph[4][2]=5 then dependencyGraph[4].size()=3 and
		reaction 4 affects reactions 2, 4, and 5
		\see fireReaction
		*/
		_dependencyGraphType dependencyGraph;


		//! number of species in the system
		std::size_t NumberOfSpecies;

		//! number of reactions in the system
		std::size_t NumberOfReactions;

		//! current time of the simulation, should be incremented at each simulation time step
		double currentTime;

		//! current population of the simulation
		_populationVectorType currentPopulation;

		//! current propensities of the simulation
		//std::vector<double> currentPropensities;
		_populationVectorType currentPropensities;

		//! current sum of propensities of the simulation, used to determine time step
		double propensitySum;

		//! index of the last reaction that fired
		/*!
		default and error value is -1
		*/
		int previousReactionIndex;

		//! counter for simulation steps taken since the last time calculateAllPropensities was called
		/*!
		is used to ensure that roundoff errors in propensities do not accumulate
		selectReaction() uses a simple strategy to recalculate all propensities using maxStepsCalculateAllPropensities
		for conservative strategy, see S. Mauch, M. Stalzer "Efficient formulations for
		exact stochastic simulation of chemical systems" IEEE/ACM Trans. on Comp. Bio. and Bioinformatics, 30 April 2009
		\see calculateAllPropensities
		\see defaultMaxStepsCalculateAllPropensities
		\see maxStepsCalculateAllPropensities
		\see selectReaction()
		*/
		std::size_t stepsSinceCalculateAllPropensities;

		//! default value for maxStepsCalculateAllPropensities
		static const std::size_t defaultMaxStepsCalculateAllPropensities=10000;

		//! maximum number of steps allowed before calling calculateAllPropensities
		/*!
		\see stepsSinceCalculateAllPropensities
		*/
		std::size_t maxStepsCalculateAllPropensities;

	private:
		//! default constructor not implemented
		SSA_Direct();

	public:

		//! Constructor
		SSA_Direct(const _populationVectorType& initialPop,
			const _stoichiometryType& stoich,
			const _propensitiesFunctorType& propensitiesFunctor,
			const _dependencyGraphType& depGraph,
			int seed=time(NULL));

		//! compiler-generated copy constructor OK
		//! compiler-generated assignment operator OK

		//! destructor
		virtual ~SSA_Direct() {
		}

		/*!
		\brief seed the random number generator
		*/
		void seed(int seed);

		/*!
		\brief run an ensemble simulation with output recorded at fixed time intervals

		output must have a conforming initialize, getOutputTimes, and record method
		outputTimes should be set in output prior to calling simulate
		if doValidate=true (the default) calls validate before ensemble
		calls initialize before each realization

		\param realizations number of simulations in the ensemble
		\param startTime the initial value of currentTime for each realization
		\param endTime the end time of each realization
		\param Output the class that handles storing the output for the simulation
		*/
		template<typename IntervalOutputType>
		void simulate(std::size_t realizations, double startTime, double endTime, IntervalOutputType& output, bool doValidate=true);

        //run a single trajectory, keep output for every reaction event
        //output is vector of pairs, first is simulation time, second is species population vector
        void simulateSingle(double startTime, double endTime, std::vector<std::pair<double, _populationVectorType> >& output, bool doValidate=true);

		/*!
		\brief initialize the state for a new simulation realization, this should be called before each realization
		*/
		virtual void initialize(double startTime=0.0);

		/*!
		\brief consistency checks to validate that the class is set up properly for a simulation, should be called before an ensemble as in simulate(...,doValidate=true)
		*/
		bool validate(double startTime, double endTime);

		/*!
		\brief update all the propensities

		updates currentPropensities by calling the propensities functor
		for each reaction using currentPopulation
		updates propensitySum
		resets stepsSinceCalculateAllPropensities to 0
		*/
		void calculateAllPropensities();

		/*!
		\brief selects the step size based on the propensitySum

		returns infinity if propensitySum is less than or equal to 0
		issues a warning if propensitySum is less than 0
		*/
		double selectStepSize();

		/*!
		\brief selects the index of the next reaction to fire based on currentPropensities

		returns -1 if there is an error
		calls calculateAllPropensities if stepsSinceCalculateAllPropensities is greater than maxStepsCalculateAllPropensities
		*/
		int selectReaction();

		/*!
		\brief fire a reaction

		updates currentPopulation
		updates currentPropensities for all affected reactions (determined by dependencyGraph[reactionIndex])
		updates propensitiesSum
		increments stepsSinceCalculateAllPropensities

		\param reactionIndex the index of the reaction to fire (-1 is an error value)
		*/
		bool fireReaction(int reactionIndex);

		/*!
		\brief take one step (select step size, increment time, fire a reaction)		
		*/
		bool step();

		double getCurrentTime();

		bool setCurrentTime(double newCurrentTime);

		_populationVectorType getCurrentPopulation();

		bool detectedVerySmallPropensity;

	};//end SSA_Direct class
}

#define _SSA_DIRECT_IPP_	
/*!
	\brief the direct method of the Stochastic Simulation Algorithm (SSA)
 */

#if !defined(_SSA_DIRECT_IPP_)
#error This file is the implementation of SSA_Direct
#endif

namespace STOCHKIT
{
    
    template<typename _populationVectorType,
    typename _stoichiometryType,
    typename _propensitiesFunctorType,
    typename _dependencyGraphType>
    SSA_Direct<_populationVectorType,
    _stoichiometryType,
    _propensitiesFunctorType,
    _dependencyGraphType>::
    SSA_Direct(const _populationVectorType& initialPop,
               const _stoichiometryType& stoich,
               const _propensitiesFunctorType& propensitiesFunctor,
               const _dependencyGraphType& depGraph,
               int seed) :
    initialPopulation(initialPop),
    stoichiometry(stoich),
    propensities(propensitiesFunctor),
    dependencyGraph(depGraph),
    NumberOfSpecies(initialPop.size()),
#ifdef MATRIX_STOICHIOMETRY
    NumberOfReactions(stoich.size1()),
    currentPropensities(stoichiometry.size1()),
#else
    NumberOfReactions(stoich.size()),
    currentPropensities(stoichiometry.size()),
#endif
    previousReactionIndex(-1),
    maxStepsCalculateAllPropensities(defaultMaxStepsCalculateAllPropensities),
    detectedVerySmallPropensity(false)
    {
        randomGenerator.seed(seed);
    }
    
    template<typename _populationVectorType,
    typename _stoichiometryType,
    typename _propensitiesFunctorType,
    typename _dependencyGraphType>
    void
    SSA_Direct<_populationVectorType,
    _stoichiometryType,
    _propensitiesFunctorType,
    _dependencyGraphType>::
    seed(int seed) {
        randomGenerator.seed(seed);
    }
    
    template<typename _populationVectorType,
    typename _stoichiometryType,
    typename _propensitiesFunctorType,
    typename _dependencyGraphType>
    void
    SSA_Direct<_populationVectorType,
    _stoichiometryType,
    _propensitiesFunctorType,
    _dependencyGraphType>::
    initialize(double startTime) {
        previousReactionIndex=-1;
        currentTime=startTime;
        currentPopulation=initialPopulation;
        calculateAllPropensities();
    }
    
    template<typename _populationVectorType,
    typename _stoichiometryType,
    typename _propensitiesFunctorType,
    typename _dependencyGraphType>
    bool
    SSA_Direct<_populationVectorType,
    _stoichiometryType,
    _propensitiesFunctorType,
    _dependencyGraphType>::
    validate(double startTime, double endTime) {
        if (startTime>=endTime) {
            Rcpp::Rcout << "StochKit ERROR (SSA_Direct::validate): startTime not before endTime\n";
            return false;
        }
        
        std::size_t N=initialPopulation.size();
#ifdef MATRIX_STOICHIOMETRY
        std::size_t M=stoichiometry.size1();
#else
        std::size_t M=stoichiometry.size();
#endif
        if (N==0) {
            Rcpp::Rcout << "StochKit ERROR (SSA_Direct::validate): initial population size=0\n";
            return false;
        }
        if (N!=NumberOfSpecies) {
            Rcpp::Rcout << "StochKit ERROR (SSA_Direct::validate): Number of species does not equal initial population size\n";
            return false;
        }
        if (M!=NumberOfReactions) {
            Rcpp::Rcout << "StochKit ERROR (SSA_Direct::validate): Number of reactions does not equal stoichiometry size\n";
            return false;
        }
        if (M!=propensities.size()) {
            Rcpp::Rcout << "StochKit ERROR (SSA_Direct::validate): Number of reactions does not equal propensities size\n";
            return false;
        }
        
        //check initial populations are all non-negative
        for (std::size_t i=0; i!=NumberOfSpecies; ++i) {
            if (initialPopulation[i]<0) {
                Rcpp::Rcout << "StochKit ERROR (SSA_Direct::validate): negative value detected in initial population\n";
                return false;
            }
        }
        
        //check that propensities, evaluated with initial population, are all non-negative
        for (std::size_t i=0; i!=NumberOfReactions; ++i) {
            if (propensities(i,initialPopulation)<0.0) {
                Rcpp::Rcout << "StochKit ERROR (SSA_Direct::validate): negative propensity detected based on initial population\n";
                return false;
            }
        }
        
        return true;
    }
    
    template<typename _populationVectorType,
    typename _stoichiometryType,
    typename _propensitiesFunctorType,
    typename _dependencyGraphType>
    void
    SSA_Direct<_populationVectorType,
    _stoichiometryType,
    _propensitiesFunctorType,
    _dependencyGraphType>::
    calculateAllPropensities() {
        double smallestNonzeroPropensity=std::numeric_limits<double>::max();
        
        propensitySum=0.0;
        for (std::size_t i=0; i!=NumberOfReactions; ++i) {
            currentPropensities[i]=propensities(i,currentPopulation);
			//Rcpp::Rcout << "propensities["<<i<<"]="<<currentPropensities[i]<<std::endl;
#ifdef DEBUG
            if (currentPropensities[i]!=currentPropensities[i]) {
                Rcpp::Rcout << "StochKit DEBUG (SSA_Direct::calculateAllPropensities): detected 'NaN (not a number)' propensity for reaction index "<<i<<".\n";
                Rcpp::Rcout << "currentPopulation was:\n";
                for (std::size_t j=0; j!=currentPopulation.size(); ++j) {
                    Rcpp::Rcout << "currentPopulation["<<j<<"]="<<currentPopulation[j]<<"\n";
                }
                Rcpp::Rcout << "Terminating.\n";
                Rcpp::stop("Fatal error encountered, terminating StochKit2R");
            }
            if (currentPropensities[i]==std::numeric_limits<double>::infinity()) {
                Rcpp::Rcout << "StochKit DEBUG (SSA_Direct::calculateAllPropensities): detected 'infinity' propensity for reaction index "<<i<<".\n";
                Rcpp::Rcout << "currentPopulation was:\n";
                for (std::size_t j=0; j!=currentPopulation.size(); ++j) {
                    Rcpp::Rcout << "currentPopulation["<<j<<"]="<<currentPopulation[j]<<"\n";
                }
                Rcpp::Rcout << "Terminating.\n";
                Rcpp::stop("Fatal error encountered, terminating StochKit2R");
            }
            if (currentPropensities[i]<0.0) {
                Rcpp::Rcout << "StochKit DEBUG (SSA_Direct::calculateAllPropensities): detected negative propensity ("<<currentPropensities[i]<<") for reaction index "<<i<<"\n";
                Rcpp::Rcout << "currentPopulation was:\n";
                for (std::size_t j=0; j!=currentPopulation.size(); ++j) {
                    Rcpp::Rcout << "currentPopulation["<<j<<"]="<<currentPopulation[j]<<"\n";
                }
                Rcpp::Rcout << "Terminating.\n";
                Rcpp::stop("Fatal error encountered, terminating StochKit2R");
            }
#endif
            propensitySum+=currentPropensities[i];
            if (currentPropensities[i]>0.0 && currentPropensities[i]<smallestNonzeroPropensity) {
                smallestNonzeroPropensity=currentPropensities[i];
            }
        }
		
        stepsSinceCalculateAllPropensities=0;
        
        if (propensitySum>0.0 && smallestNonzeroPropensity/propensitySum<2E-10) { //per S.Mauch, M.Stalzer. (2009) "Efficient Formulations for Exact..."
            if (detectedVerySmallPropensity==false) {
                detectedVerySmallPropensity=true;
                Rcpp::Rcout << "StochKit WARNING (SSA_Direct::calculateAllPropensities): detected very small propensity value, biased sampling of small propensity reactions may occur\n";
            }
        }
#ifdef DEBUG
        //a reasonable place to check for possible time step inaccuracy
        if (propensitySum>0.0 && currentTime>2E21/propensitySum) { //per Mauch, Stalzer. (2009) "Efficient Formulations..."
            Rcpp::Rcout << "StochKit DEBUG (SSA_Direct::calculateAllPropensities): ratio of average time step size to simulation currentTime is very small, may lead to step size inaccuracies\n";
        }
#endif
    }
    
    template<typename _populationVectorType,
    typename _stoichiometryType,
    typename _propensitiesFunctorType,
    typename _dependencyGraphType>
    inline
    double
    SSA_Direct<_populationVectorType,
    _stoichiometryType,
    _propensitiesFunctorType,
    _dependencyGraphType>::
    selectStepSize() {
        
        if (propensitySum<0.0) {
#ifdef DEBUG
            Rcpp::Rcout << "StochKit DEBUG (SSA_Direct::selectStepSize): detected negative propensitySum, recalculating\n";
#endif
            //if propensitySum negative, recalculate all propensities
            calculateAllPropensities();
            //if still negative, give warning and return infinity
            if (propensitySum<0.0) {
                Rcpp::Rcout << "StochKit WARNING (SSA_Direct::selectStepSize): propensitySum<0, returning step size=infinity\n";
                return std::numeric_limits<double>::infinity();
            }
            if (propensitySum==0) {
                return std::numeric_limits<double>::infinity();
            }
        }
        
        return randomGenerator.exponential(propensitySum);
    }
    
    template<typename _populationVectorType,
    typename _stoichiometryType,
    typename _propensitiesFunctorType,
    typename _dependencyGraphType>
    inline
    int
    SSA_Direct<_populationVectorType,
    _stoichiometryType,
    _propensitiesFunctorType,
    _dependencyGraphType>::
    selectReaction() {
		
//		Rcpp::Rcout << "selectReaction...\n";
//		for (int i = 0; i<currentPropensities.size(); i++) {
//			Rcpp::Rcout << "currentPropensities["<<i<<"]="<<currentPropensities[i]<<std::endl;
//		}
		
        previousReactionIndex=-1;
        if (stepsSinceCalculateAllPropensities>maxStepsCalculateAllPropensities) {
            calculateAllPropensities();
        }
        
        //check for 0 propensity sum
        if (propensitySum<=0) {
            //if negative, may be due to numerical error
            if (propensitySum<0) {
                calculateAllPropensities();
            }
            if (propensitySum<=0) {
                return -1;
            }
        }
        
        //generate a uniform random number between (0,propensitySum)
        double r=0;
        while (r==0) {
            r=randomGenerator.continuousUniform(0,1)*propensitySum;
        }
        double jsum=0;
        while (jsum < r) {
            ++previousReactionIndex;
            //test that we don't run off end of array
            if (previousReactionIndex==(int)NumberOfReactions) {
#ifdef DEBUG
                Rcpp::Rcout << "StochKit DEBUG (SSA_Direct::selectReaction): detected numerical error in propensities, recalculating\n";
#endif
                calculateAllPropensities();
                return selectReaction();
            }
            else {
#ifdef DEBUG
                if (currentPropensities[previousReactionIndex]<0.0) {
                    Rcpp::Rcout << "StochKit DEBUG (SSA_Direct::selectReaction): detected negative propensity ("<<currentPropensities[previousReactionIndex]<<") for reaction index "<<previousReactionIndex<<"\n";
                    Rcpp::Rcout << "currentPopulation was:\n";
                    for (std::size_t i=0; i!=currentPopulation.size(); ++i) {
                        Rcpp::Rcout << "currentPopulation["<<i<<"]="<<currentPopulation[i]<<"\n";
                    }
                    Rcpp::Rcout << "Terminating.\n";
                    Rcpp::stop("Fatal error encountered, terminating StochKit2R");
                }
#endif
                jsum+=currentPropensities[previousReactionIndex];
            }
        }
		
		//Rcpp::Rcout << "firing "<<previousReactionIndex<<std::endl;
        return previousReactionIndex;
    }
    
    template<typename _populationVectorType,
    typename _stoichiometryType,
    typename _propensitiesFunctorType,
    typename _dependencyGraphType>
    inline
    bool
    SSA_Direct<_populationVectorType,
    _stoichiometryType,
    _propensitiesFunctorType,
    _dependencyGraphType>::
    fireReaction(int reactionIndex) {
        if (reactionIndex==-1) {
#ifdef DEBUG
            Rcpp::Rcout << "StochKit DEBUG (SSA_Direct::fireReaction): attempt to fire reaction index = -1. Terminating.\n";
            Rcpp::stop("Fatal error encountered, terminating StochKit2R");
#endif
            return false;
        }
        else {
            //if not -1, assumes valid reactionIndex
#ifdef MATRIX_STOICHIOMETRY
            matrixrow dX(stoichiometry, reactionIndex);
            typename matrixrow::iterator it;
            for(it=dX.begin();it!=dX.end();it++) {
                currentPopulation[it.index()]+=*it;
            }
#else
            currentPopulation+=stoichiometry[reactionIndex];
#endif
			
			//Rcpp::Rcout << "currentPopulation[0]="<<currentPopulation[0]<<std::endl;
			
            int affectedReactionIndex;
            double oldPropensity;
            for (std::size_t i=0; i!=dependencyGraph[reactionIndex].size(); ++i) {
                affectedReactionIndex=dependencyGraph[reactionIndex][i];
                oldPropensity=currentPropensities[affectedReactionIndex];
                currentPropensities[affectedReactionIndex]=propensities(affectedReactionIndex,currentPopulation);
#ifdef DEBUG
                if (currentPropensities[affectedReactionIndex]!=currentPropensities[affectedReactionIndex]) {
                    Rcpp::Rcout << "StochKit DEBUG (SSA_Direct::fireReaction): detected 'NaN (not a number)' propensity for reaction index "<<affectedReactionIndex<<" while updating after firing reaction index "<<reactionIndex<<".\n";
                }
                if (currentPropensities[affectedReactionIndex]==std::numeric_limits<double>::infinity()) {
                    Rcpp::Rcout << "StochKit DEBUG (SSA_Direct::fireReaction): detected 'infinity' propensity for reaction index "<<affectedReactionIndex<<" while updating after firing reaction index "<<reactionIndex<<".\n";
                }
                if (currentPropensities[affectedReactionIndex]<0.0) {
                    Rcpp::Rcout << "StochKit DEBUG (SSA_Direct::fireReaction): detected negative propensity for reaction index "<<affectedReactionIndex<<" while updating after firing reaction index "<<reactionIndex<<".\n";
                    Rcpp::Rcout << "updated currentPopulation is:\n";
                    for (std::size_t i=0; i!=currentPopulation.size(); ++i) {
                        Rcpp::Rcout << "currentPopulation["<<i<<"]="<<currentPopulation[i]<<"\n";
                    }
                    Rcpp::Rcout << "population before firing reaction was:\n";
                    _populationVectorType oldPop=currentPopulation-=stoichiometry[reactionIndex];
                    for (std::size_t i=0; i!=oldPop.size(); ++i) {
                        Rcpp::Rcout << "oldPopulation["<<i<<"]="<<oldPop[i]<<"\n";
                    }
                    Rcpp::Rcout << "Terminating.\n";
                    Rcpp::stop("Fatal error encountered, terminating StochKit2R");
                }
#endif
                propensitySum+=currentPropensities[affectedReactionIndex]-oldPropensity;
            }
            stepsSinceCalculateAllPropensities++;
            return true;
        }
    }
    
    template<typename _populationVectorType,
    typename _stoichiometryType,
    typename _propensitiesFunctorType,
    typename _dependencyGraphType>
    inline
    bool
    SSA_Direct<_populationVectorType,
    _stoichiometryType,
    _propensitiesFunctorType,
    _dependencyGraphType>::
    step() {
        currentTime+=selectStepSize();
        return fireReaction(selectReaction());
    }
    
    template<typename _populationVectorType,
    typename _stoichiometryType,
    typename _propensitiesFunctorType,
    typename _dependencyGraphType>
    template<typename IntervalOutputType>
    void
    SSA_Direct<_populationVectorType,
    _stoichiometryType,
    _propensitiesFunctorType,
    _dependencyGraphType>::
    simulate(std::size_t realizations, double startTime, double endTime, IntervalOutputType& output, bool doValidate) {
        
        if (doValidate) {
            if (!validate(startTime,endTime)) {
                Rcpp::Rcout << "StochKit ERROR (SSA_Direct::simulate): validate() failed, simulation aborted\n";
                Rcpp::stop("Fatal error encountered, terminating StochKit2R");
            }		
        }
        
        if (!output.initialize(realizations,startTime,endTime,initialPopulation)) {
            Rcpp::Rcout << "StochKit ERROR (SSA_Direct::simulate): initialization of output object failed, simulation aborted\n";
            Rcpp::stop("Fatal error encountered, terminating StochKit2R");
        }
        
        std::vector<double> outputTimes = output.getOutputTimes();
        std::size_t totalIntervals=outputTimes.size();
        
        std::size_t currentInterval;
        
        for (std::size_t currentRealization=0; currentRealization!=realizations; ++currentRealization) {
            
            initialize(startTime);
            currentInterval=0;
            
            
            currentTime+=selectStepSize();
            
            
            while (currentTime<endTime) {
                
                while (currentInterval<totalIntervals && currentTime >=outputTimes[currentInterval]){
                    output.record(currentRealization,currentInterval,currentPopulation);
                    currentInterval++;
                }
                
                fireReaction(selectReaction());
                currentTime+=selectStepSize();
            }
            while (currentInterval<totalIntervals && currentTime>=outputTimes[currentInterval]){
                output.record(currentRealization,currentInterval,currentPopulation);
                currentInterval++;
            }	
        }
        
    }//end simulate

    template<typename _populationVectorType,
    typename _stoichiometryType,
    typename _propensitiesFunctorType,
    typename _dependencyGraphType>
    void
    SSA_Direct<_populationVectorType,
    _stoichiometryType,
    _propensitiesFunctorType,
    _dependencyGraphType>::
    simulateSingle(double startTime, double endTime, std::vector<std::pair<double, _populationVectorType> >& output, bool doValidate) {
        
        if (doValidate) {
            if (!validate(startTime,endTime)) {
                Rcpp::Rcout << "StochKit ERROR (SSA_Direct::simulateSingle): validate() failed, simulation aborted\n";
                Rcpp::stop("Fatal error encountered, terminating StochKit2R");
            }
        }
            
        initialize(startTime);
        
        output.push_back(std::make_pair(startTime,initialPopulation));
        
        currentTime+=selectStepSize();
        
        while (currentTime<endTime) {
            fireReaction(selectReaction());
            output.push_back(std::make_pair(currentTime,currentPopulation));
            currentTime+=selectStepSize();
        }
        
        output.push_back(std::make_pair(endTime,currentPopulation));

        
    }//end simulateSingle

    template<typename _populationVectorType,
    typename _stoichiometryType,
    typename _propensitiesFunctorType,
    typename _dependencyGraphType>
    inline
    double
    SSA_Direct<_populationVectorType, 
    _stoichiometryType,
    _propensitiesFunctorType,
    _dependencyGraphType>::
    getCurrentTime(){
        return currentTime;
    }
    
    
    template<typename _populationVectorType, 
    typename _stoichiometryType,
    typename _propensitiesFunctorType,
    typename _dependencyGraphType>
    inline
    bool
    SSA_Direct<_populationVectorType, 
    _stoichiometryType,
    _propensitiesFunctorType,
    _dependencyGraphType>::
    setCurrentTime(double newCurrentTime){
        currentTime = newCurrentTime;
        return true;
    }
    
    template<typename _populationVectorType, 
    typename _stoichiometryType,
    typename _propensitiesFunctorType,
    typename _dependencyGraphType>
    inline
    _populationVectorType
    SSA_Direct<_populationVectorType, 
    _stoichiometryType,
    _propensitiesFunctorType,
    _dependencyGraphType>::
    getCurrentPopulation(){
        return currentPopulation;
    }
    
}


#undef _SSA_DIRECT_IPP_

#endif
