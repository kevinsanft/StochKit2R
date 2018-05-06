// [[Rcpp::depends(BH)]]
#include "CustomPropensity.h"
#include "MassActionModel.h"
#include "StandardDriverTypes.h"
#include "solver_helper_functions.h"
#include <Rcpp.h>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

using namespace std;

const double sigmac = 10.0;
const double Rc = 28.0;
const double bc = 8.0 / 3.0;

//typedef boost::array< double , 3 > state_type;
typedef boost::numeric::ublas::vector< double > state_type;
void lorenz(const state_type &x , state_type &dxdt , double t )
{
	dxdt(0) = sigmac * ( x(1) - x(0) );
	dxdt(1) = Rc * x(0) - x(1) - x(0) * x(2);
	dxdt(2) = -bc * x(2) + x(0) * x(1);
}

void write_lorenz( const state_type &x , const double t )
{
	Rcpp::Rcout << t << '\t' << x(0) << '\t' << x(1) << '\t' << x(2) << endl;
}

class RHS {

	STOCHKIT::StandardDriverTypes::propensitiesType props;
	STOCHKIT::StandardDriverTypes::denseStoichiometryType nu;
	int N;
public:
	RHS(const STOCHKIT::StandardDriverTypes::propensitiesType& props, const STOCHKIT::StandardDriverTypes::denseStoichiometryType& nu) : props(props),nu(nu)
	{
		N=nu[0].size();
	}
	
	//propensities/rate vector same type as population
	void operator()( const state_type &x, state_type &dxdt, const double t ) {

		std::fill(dxdt.begin(),dxdt.end(),0.0);
		for (int i=0; i!=props.size(); i++) {
			dxdt += props(i,x)*nu[i];
		}
	}
};

class Recorder {
public:
	std::vector<double>& times;
	std::vector<state_type>& state;
	
	
	Recorder(std::vector<double>& t ,std::vector<state_type>& s) : times(t), state(s)
	{}
	
public:
	void operator()(const state_type &x , const double t ) {
//		Rcpp::Rcout << "recording: " << t << '\t';
//		for (int i=0; i<x.size(); i++) {
//			Rcpp::Rcout << x(i) << '\t';
//		}
//		Rcpp::Rcout << endl;
		times.push_back(t);
		//Rcpp::Rcout << "times.size()="<<times.size()<<endl;
		state.push_back(x);
	}
};

//'@title C++ Interface to Gillespie Stochastic Simulation Algorithm single trajectory
//'
//'@description
//'\code{ode} Called by StochKit2R ode function, do not call this C++ interface directly
//'
//'@param StochKit2Rmodel R list (Rcpp List built from buildStochKit2Rmodel output)
//'@param endTime Simulation end time
//'@return Dataframe containing the time and population sizes
//'@keywords internal
// [[Rcpp::export]]
RcppExport SEXP odecpp(Rcpp::List StochKit2Rmodel, double endTime, int intervals)
{
	//create StochKit2R mass action model object
	//first, pull out pieces from list object
	Rcpp::List rParameterList=StochKit2Rmodel[0];
	Rcpp::List rSpeciesList=StochKit2Rmodel[1];
	Rcpp::List rReactionList=StochKit2Rmodel[2];
	Rcpp::List rCustomPropensityList=StochKit2Rmodel[3];
	
	//get species labels...
	std::vector<std::string> modelSpeciesList = getSpeciesList(rSpeciesList);

	//create stochastic model
	STOCHKIT::MassActionModel<STOCHKIT::StandardDriverTypes::populationType,
	STOCHKIT::StandardDriverTypes::denseStoichiometryType,
	STOCHKIT::StandardDriverTypes::propensitiesType,
	STOCHKIT::StandardDriverTypes::graphType> model(rParameterList,
													rReactionList,
													rSpeciesList,rCustomPropensityList);

	//get stochastic propensities
	STOCHKIT::StandardDriverTypes::propensitiesType propensityFunctions=model.writePropensities();
	
	//convert to ODE model by iterating over reactions
	//convert A+A type propensities to ODE rates
	propensityFunctions.convertToODE();

	//get stoichiometry matrix
	STOCHKIT::StandardDriverTypes::denseStoichiometryType stoichiometry=model.writeStoichiometry();
	
	RHS rhs(model.writePropensities(),model.writeStoichiometry());
	
	state_type x = model.writeInitialPopulation();
	
	int outputrows=0;
	std::vector<double> times(intervals+1);
	if (intervals==0) {
		times[0]=0.0;
		times.push_back(endTime);
		outputrows=1;
	}
	else {
		times[0] = 0.0;
		for (int i=1; i<intervals; ++i) {
			times[i] = (double)i*endTime/(double)intervals;
		}
		times[intervals] = endTime;
		outputrows=times.size();
	}
	
	//typedef boost::numeric::odeint::runge_kutta_cash_karp54< state_type > error_stepper_type;
	//typedef runge_kutta_fehlberg78< state_type > error_stepper_type;
	//typedef boost::numeric::odeint::controlled_runge_kutta< error_stepper_type > controlled_stepper_type;
	//double abs_err = 1.0e-14 , rel_err = 1.0e-10 , a_x = 1.0 , a_dxdt = 1.0;
	//	double abs_err = 1.0e-8 , rel_err = 1.0e-7 , a_x = 1.0 , a_dxdt = 1.0;
	//controlled_stepper_type controlled_stepper(boost::numeric::odeint::default_error_checker< double ,  boost::numeric::odeint::range_algebra ,  boost::numeric::odeint::default_operations >( abs_err , rel_err , a_x , a_dxdt ) );
	typedef boost::numeric::odeint::runge_kutta_dopri5< state_type > stepper_type;

	std::vector<double> times2;
	std::vector<state_type> state;
	Recorder recorder(times2,state);
	
	Rcpp::NumericMatrix df(outputrows,modelSpeciesList.size()+1);
	
	//Rcpp::Rcout << "modelSpeciesList.size()=" << modelSpeciesList.size() << std::endl;

	double initialDeltaT = 1e-8;
	if (intervals==0) {
		//only output end time
		integrate_times( boost::numeric::odeint::make_controlled( 1E-8 , 1E-7 , stepper_type() ), rhs , x , times.begin() , times.end() , initialDeltaT, recorder);
		df(0,0) = recorder.times[recorder.times.size()-1];
		for (int i=0; i<modelSpeciesList.size();i++) {
			df(0,i+1)=recorder.state[recorder.state.size()-1](i);
		}
	}
	else {

		
		integrate_times( boost::numeric::odeint::make_controlled( 1E-8 , 1E-7 , stepper_type() ), rhs , x , times.begin() , times.end() , initialDeltaT, recorder);

		for (int i=0; i<recorder.times.size(); i++) {
			df(i,0)=recorder.times[i];
			for (int j=0; j<modelSpeciesList.size();j++) {
				df(i,j+1)=recorder.state[i](j);
			}
		}
	}
	//boost::numeric::odeint::integrate_times(boost::numeric::odeint::rosenbrock4_controller<state_type>(), lorenz , x , 0.0 , 25.0 , 25.0/20.0 , write_lorenz );
	
	
	// name the columns accordingly
	Rcpp::CharacterVector col_names;
	col_names.push_back("time");
	
	for (int i = 0; i < modelSpeciesList.size(); i++)
		col_names.push_back(modelSpeciesList[i]);
	Rcpp::colnames(df) = col_names;
	return df;
}
