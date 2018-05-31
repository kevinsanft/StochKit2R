#include "Random.h"

//#include <random>//for random_device, for seed only

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <ctime>

Random::Random() {
    seed();
}

void Random::seed()
{
//    std::random_device rd;
    generator.seed((boost::uint32_t) time(NULL));
}

void Random::seed(unsigned int seed)
{
    generator.seed(seed);
}

// continuous uniform random number in (a,b)
double Random::continuousUniform(double a, double b)
{
    boost::uniform_01<> continuousZeroOneOpen;
    boost::variate_generator<boost::mt19937&, boost::uniform_01<> > generatorContinuousZeroOneOpen(generator, continuousZeroOneOpen);
    return (generatorContinuousZeroOneOpen()*(b-a) + a);
}

// mean = 1/rate
double Random::exponential(double rate)
{
    if (rate<=0.0) {
        return std::numeric_limits<double>::infinity();
    }
    else {
        boost::exponential_distribution<> exponentialDistribution(rate);
        boost::variate_generator<boost::mt19937&, boost::exponential_distribution<> > generatorExponential(generator, exponentialDistribution);
        return generatorExponential();
    }
}

double Random::poisson(double mean)
{
    boost::poisson_distribution<> poissonDistribution(mean);
    boost::variate_generator<boost::mt19937&, boost::poisson_distribution<> > generatorPoisson(generator, poissonDistribution);
    return generatorPoisson();
}
