#ifndef FEHECLUSTER_H
#define FEHECLUSTER_H

// Includes
#include <sstream>
#include "FeCluster.h"
#include <xolotlPerf.h>
#include "FeClusterReactionNetwork.h"

namespace xolotlCore {

/**
 * This class represents a cluster composed entirely of helium.
 */
class FeHeCluster: public FeCluster {

private:
	static std::string buildName(IReactant::SizeType size) {
		// Set the reactant name appropriately
		std::stringstream nameStream;
		nameStream << "He_" << size;
		return nameStream.str();
	}

public:

	/**
	 * Default constructor, deleted because we require info to construct.
	 */
	FeHeCluster() = delete;

	/**
	 * The constructor. All HeClusters must be initialized with a size.
	 *
	 * @param nHe the number of helium atoms in the cluster
	 * @param registry The performance handler registry
	 */
	FeHeCluster(int nHe, IReactionNetwork& _network,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
			FeCluster(_network, registry, buildName(nHe)) {

		// Set the size
		size = nHe;
		// Update the composition map
		composition[toCompIdx(Species::He)] = size;

		// Set the typename appropriately
		type = ReactantType::He;

		// Compute the reaction radius
		constexpr auto FourPi = 4 * xolotlCore::pi;
		constexpr auto aCubed = ipow<3>(xolotlCore::ironLatticeConstant);
		double termOne = std::cbrt((3 / FourPi) * (1.0 / 10) * aCubed * size);
		double termTwo = std::cbrt((3 / FourPi) * (1.0 / 10) * aCubed);
		reactionRadius = 0.3 + termOne - termTwo;

		// Bounds on He and V
		heBounds = IntegerRange<IReactant::SizeType>(
				static_cast<IReactant::SizeType>(size),
				static_cast<IReactant::SizeType>(size + 1));
		vBounds = IntegerRange<IReactant::SizeType>(
				static_cast<IReactant::SizeType>(0),
				static_cast<IReactant::SizeType>(1));

		return;
	}

	/**
	 * Copy constructor, deleted to prevent use.
	 */
	FeHeCluster(const FeHeCluster& other) = delete;

	/**
	 * Destructor
	 */
	~FeHeCluster() {
	}

};
//end class FeHeCluster

} /* namespace xolotlCore */
#endif
