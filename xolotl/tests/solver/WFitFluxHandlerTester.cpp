#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include "WFitFluxHandler.h"

using namespace std;
using namespace xolotlSolver;

/**
 * The test suite is responsible for testing the WFitFluxHandler.
 */
BOOST_AUTO_TEST_SUITE (FitFluxHandlerTester_testSuite)

BOOST_AUTO_TEST_CASE(checkgetIncidentFlux) {

	// Specify the number of grid points that will be used
	int nGridpts = 8;
	// Specify the step size between grid points
	double step = 1.142857142857143;

    auto testFitFlux = std::make_shared<xolotlSolver::WFitFluxHandler>();
    // Initialize the flux handler
    testFitFlux->initializeFluxHandler(nGridpts, step);

	// Create a composition vector
	std::vector<int> compVec = {1, 0, 0};
	// Create a time
	double currTime = 1.0;

	// Create a vector representing the position of the cluster
	std::vector<double> x = {1.142857142857143, 0.0, 0.0};
	// x is a gridpoint (position) in PetscSolver, RHSFunction
	//x=1.142857142857143=8/7 where x = xi * hx with xi=1, hx=1.142857142857143
	double fitFunction = 6.0 * x[0] * x[0] * x[0] - 87.0 * x[0] * x[0] + 300.0 * x[0];
	//fitFunction = 0.02381807580174927

	auto testFlux = testFitFlux->getIncidentFlux(compVec, x, 1);

	BOOST_TEST_MESSAGE( "\nWFitFluxHandlerTester Message: \n"
						<< "incidentFlux = " << testFlux << " with composition "
						<< "(" << x[0] << "," << x[1] << "," << x[2] << ") "
						<< "at time = " << currTime << "\n"
						<< "fitFunction = " << fitFunction << "\n");
	BOOST_REQUIRE_EQUAL(testFlux, fitFunction);

}

BOOST_AUTO_TEST_CASE(checkheFluence) {

	// Specify the number of grid points that will be used
	int nGridpts = 8;
	// Specify the step size between grid points
	double step = 1.142857142857143;

    auto testFitFlux = std::make_shared<xolotlSolver::WFitFluxHandler>();
    // Initialize the flux handler
    testFitFlux->initializeFluxHandler(nGridpts, step);

	// Create a composition vector
	std::vector<int> compVec = {1, 0, 0};
	// Create a time
	double currTime = 1.0;

	// Create a vector representing the position of the cluster
	std::vector<double> x = {1.142857142857143, 0.0, 0.0};
	// x is a gridpoint (position) in PetscSolver, RHSFunction
	//x=1.142857142857143=8/7 where x = xi * hx with xi=1, hx=1.142857142857143
	double fitFunction = 6.0 * x[0] * x[0] * x[0] - 87.0 * x[0] * x[0] + 300.0 * x[0];
	//fitFunction = 238.1807580174927

	auto testFlux = testFitFlux->getIncidentFlux(compVec, x, 1);
	BOOST_TEST_MESSAGE( "\nWFitFluxHandlerTester Message: \n"
						<< "incidentFlux = " << testFlux << " with composition "
						<< "(" << x[0] << "," << x[1] << "," << x[2] << ") "
						<< "at time = " << currTime << "\n"
						<< "fitFunction = " << fitFunction << "\n");
	BOOST_REQUIRE_EQUAL(testFlux, fitFunction);

	// Set the maximum helium fluence value
	testFitFlux->setMaxHeFluence(9e-06);
	BOOST_REQUIRE_EQUAL(testFitFlux->getMaxHeFluence(), 9e-06);

	// Increment the helium fluence
	testFitFlux->incrementHeFluence(1.0e-8, 4/3);
	BOOST_REQUIRE_EQUAL(testFitFlux->getHeFluence(), 9.1521865889212782e-06);
	// Increment the helium fluence again, although the helium fluence value should not change
	testFitFlux->incrementHeFluence(1.0e-7, 4/3);
	// Check to make sure the helium fluence value did not increment a second time
	BOOST_REQUIRE_EQUAL(testFitFlux->getHeFluence(), 9.1521865889212782e-06);
	// Since the maximum helium fluence value has been exceeded, an incident flux of 0 should be returned
	BOOST_REQUIRE_EQUAL(testFitFlux->getIncidentFlux(compVec, x, 1), 0.0);

}


BOOST_AUTO_TEST_SUITE_END()
