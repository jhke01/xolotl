/*
 * PSIClusterTester.cpp
 *
 *  Created on: May 6, 2013
 *      Author: Jay Jay Billings
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <PSICluster.h>
#include "SimpleReactionNetwork.h"
#include <HeVCluster.h>
#include <memory>
#include <typeinfo>
#include <limits>
#include <algorithm>

using namespace std;
using namespace xolotlCore;
using namespace testUtils;


/**
 * This suite is responsible for testing the HeVCluster.
 */
BOOST_AUTO_TEST_SUITE(HeVCluster_testSuite)


BOOST_AUTO_TEST_CASE(getSpeciesSize) {
	HeVCluster cluster(4, 5);
	BOOST_REQUIRE_EQUAL(cluster.getSpeciesSize("He"), 4);
	BOOST_REQUIRE_EQUAL(cluster.getSpeciesSize("V"), 5);
	BOOST_REQUIRE_EQUAL(cluster.getSpeciesSize("I"), 0);
}

/**
 * This operation checks the ability of the HeVCluster to describe
 * its connectivity to other clusters.
 */
BOOST_AUTO_TEST_CASE(checkConnectivity) {

	// Local Declarations
	shared_ptr<ReactionNetwork> network = getSimpleReactionNetwork();
	int maxClusterSize = 10;
	int numClusters = maxClusterSize;
	
	shared_ptr<vector<shared_ptr<Reactant>>> reactants = network->reactants;
	shared_ptr<map<string, string>> props = network->properties;

	// Write the cluster information to stdout
	BOOST_TEST_MESSAGE("Sizes of clusters in network:");
	
	for (auto reactantIt = reactants->begin(); reactantIt != reactants->end(); reactantIt++) {
		// Write the size of the psi cluster to stdout
		shared_ptr<PSICluster> psiCluster = static_pointer_cast<PSICluster>(*reactantIt);
		BOOST_TEST_MESSAGE(psiCluster->getSize());
	}
	
	BOOST_TEST_MESSAGE("Maximum He Cluster Size = " << (*props)["maxHeClusterSize"]);
	BOOST_TEST_MESSAGE("Maximum V Cluster Size = " << (*props)["maxVClusterSize"]);
	BOOST_TEST_MESSAGE("Number of He clusters = " << (*props)["numHeClusters"]);
	BOOST_TEST_MESSAGE("Number of V clusters = " << (*props)["numVClusters"]);
	BOOST_TEST_MESSAGE("Number of mixed clusters = " << (*props)["numMixedClusters"]);

	// Get the connectivity of the HeV cluster with 5 He and 2 V
	std::map<std::string, int> speciesMap;
	speciesMap["He"] = 5;
	speciesMap["V"] = 2;
	int clusterIndex = network->toClusterIndex(speciesMap);
	
	shared_ptr<Reactant> reactant = reactants->at(clusterIndex);
	shared_ptr<HeVCluster> cluster = std::dynamic_pointer_cast<HeVCluster>(reactant);
	
	// Make sure we're retreiving the correct reactant
	
	// The HeVCluster should have 5 He and 2 V, by the ordering of
	// the SimpleReactionNetwork.
	
	BOOST_REQUIRE_EQUAL(cluster->getSpeciesSize("He"), 5);
	BOOST_REQUIRE_EQUAL(cluster->getSpeciesSize("V"), 2);
	
	vector<int> connectivityArray = cluster->getConnectivity();
	
	// The connectivity array should be the same size as the reactants array
	BOOST_TEST_MESSAGE("Connectivity Array Size = " << connectivityArray.size());
	BOOST_REQUIRE_EQUAL(connectivityArray.size(), reactants->size());
	
	// The following are reactions of the HeVCluster
	
	// xHe * yV + zHe --> (x + z)He + yV
	// xHe * yV + V   --> xHe + (y + 1)V
	// xHe * yV + zI  --> xHe + (y - z)V
	
	int connectivityExpected[] = {
		// He
		1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
		
		// V
		1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		
		// I
		1, 0, 0, 0, 0, 0, 0, 0, 0, 0
	};
	
	// Check the He, V, and I reaction
	
	for (int i = 0; i < 30; i++) {
		BOOST_REQUIRE_EQUAL(connectivityArray.at(i), connectivityExpected[i]);
	}
	
	// The rest should be zero
	
	for (int i = 30; i < reactants->size(); i++) {
		BOOST_REQUIRE_EQUAL(connectivityArray.at(i), 0);
	}
}

BOOST_AUTO_TEST_CASE(checkGetFlux) {
	/*shared_ptr<ReactionNetwork> network = getSimpleReactionNetwork();
	int maxClusterSize = 10;
	int numClusters = maxClusterSize;

	shared_ptr<vector<shared_ptr<Reactant>>> reactants = network->reactants;
	shared_ptr<map<string, string>> props = network->properties;

	for (int i = 0; i < reactants->size(); i++) {
		reactants->at(i)->setConcentration(2.0*i);
		(dynamic_pointer_cast<PSICluster>(reactants->at(i)))->setMigrationEnergy(1.0);
		(dynamic_pointer_cast<PSICluster>(reactants->at(i)))->setDiffusionFactor(1.0);
	}
	// Get the connectivity of the 20th HeV cluster (index 59).
	shared_ptr<Reactant> reactant = reactants->at(3 * numClusters + 20 - 1);
	shared_ptr<HeVCluster> cluster = dynamic_pointer_cast<HeVCluster>(reactant);

	double flux = cluster->getProductionFlux(1.0);

	std::cout << "Flux is " << flux << "\n";*/
}

/**
 * This operation checks the reaction radius for HeVCluster.
 */
BOOST_AUTO_TEST_CASE(checkReactionRadius) {

	std::vector<std::shared_ptr<HeVCluster>> clusters;
	std::shared_ptr<HeVCluster> cluster;
	double expectedRadii[] = { 0.4330127019, 0.5609906819, 0.6507642333,
			0.7222328328, 0.7825853415 };

	for (int i = 1; i <= 5; i++) {
		cluster = std::shared_ptr<HeVCluster>(new HeVCluster(1, i));
		BOOST_CHECK_CLOSE(expectedRadii[i - 1], cluster->getReactionRadius(),
				.000001);
	}
}

BOOST_AUTO_TEST_SUITE_END()
