#ifndef INTERSTITIALCLUSTER_H
#define INTERSTITIALCLUSTER_H

// Includes
#include <sstream>
#include "PSICluster.h"

namespace xolotlCore {

/**
 * This class represents a cluster composed entirely of interstitial defects.
 */
class InterstitialCluster: public PSICluster {

    static
    std::string buildName(IReactant::SizeType nI) {
        // Set the reactant name appropriately
        std::stringstream nameStream;
        nameStream << "I_" << nI;
        return nameStream.str();
    }

public:

    /**
     * Default constructor, deleted because we require info to construct.
     */
    InterstitialCluster() = delete;

	/**
	 * The constructor. All InterstitialClusters must be initialized with
	 * a size.
	 *
	 * @param nI The number of interstitial defect in this cluster
	 * @param registry The performance handler registry
	 */
	InterstitialCluster(int nI,
            IReactionNetwork& _network,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

    /**
     * Copy constructor, deleted to prevent use.
     */
    InterstitialCluster(const InterstitialCluster& other) = delete;

	/**
	 * The Destructor
	 */
	~InterstitialCluster() {
	}

	/**
	 * This operation returns the total change in this cluster due its
	 * own dissociation.
	 *
	 * @return The flux due to its dissociation
	 */
	double getEmissionFlux() const;

	/**
	 * This operation computes the partial derivatives due to emission
	 * reactions.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 */
	void getEmissionPartialDerivatives(std::vector<double> & partials) const;

};
//end class InterstitialCluster

} /* end namespace xolotlCore */
#endif
