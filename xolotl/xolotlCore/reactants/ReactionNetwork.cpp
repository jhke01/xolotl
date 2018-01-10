#include "ReactionNetwork.h"
#include <xolotlPerf.h>
#include <iostream>
#include <cassert>

namespace xolotlCore {

ReactionNetwork::ReactionNetwork(
		const std::set<ReactantType>& _knownReactantTypes,
		ReactantType _superClusterType,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> _registry) :
		knownReactantTypes(_knownReactantTypes), superClusterType(
				_superClusterType), handlerRegistry(_registry), temperature(
				0.0), dissociationsEnabled(true) {

	// Ensure our per-type cluster map can store Reactants of the types
	// we support.
	for (auto const& currType : knownReactantTypes) {
		clusterTypeMap.emplace(
				std::make_pair(currType, IReactionNetwork::ReactantMap()));
	}

	// Ensure we have a baseline for determining max cluster size for
	// the types we support.
	for (auto const& currType : knownReactantTypes) {
		maxClusterSizeMap.insert( { currType, 0 });
	}
	return;
}

void ReactionNetwork::add(std::unique_ptr<IReactant> reactant) {

	// Ensure we have a valid object to work with.
	assert(reactant);

	// Get the composition
	auto& composition = reactant->getComposition();

	// Check if we already know about this reactant.
	auto& currTypeMap = clusterTypeMap[reactant->getType()];
	auto iter = currTypeMap.find(composition);
	if (iter == currTypeMap.end()) {

		// Set the id for this cluster
		// (It is networkSize+1 because we haven't added
		// it to the network yet.)
		reactant->setId(size() + 1);

		// Update the max cluster size for the new cluster's type.
		maxClusterSizeMap[reactant->getType()] = std::max(reactant->getSize(),
				maxClusterSizeMap[reactant->getType()]);

		// Note the reactant in our flat list of all reactants.
		allReactants.emplace_back(*reactant);

		// Give reactant to the appropriate per-type map.
		currTypeMap.emplace(composition, std::move(reactant));

	} else {
		std::stringstream errStream;
		errStream << "ReactionNetwork: not adding duplicate "
				<< toString(reactant->getType()) << ':' << composition;

		throw errStream.str();
	}

	return;
}

// TODO have this return an iterator to the appropriate map.  Probably
// will have to have it return a pair with iter and bool, to ease
// detection of whether we found it or not.  (Figuring out which
// map to test the end against is probable less easy.)
// Alternative is to remove 'get()' entirely and rely on callers
// to getAll() for particular type, then do their find/test to see if
// we found it.  That would be (slightly?) faster since now caller
// has to test return for nullptr anyway, and we're already testing
// against end() of the appropriate map.
IReactant * ReactionNetwork::get(ReactantType type,
		const IReactant::Composition& comp) const {

	IReactant* ret = nullptr;

	// Check if the reactant is in the map
	auto const& currTypeMap = clusterTypeMap.at(type);
	auto iter = currTypeMap.find(comp);
	if (iter != currTypeMap.end()) {
		ret = iter->second.get();
	}

	return ret;
}

double ReactionNetwork::calculateReactionRateConstant(
		const ProductionReaction& reaction) const {

	// Get the reaction radii
	double r_first = reaction.first.getReactionRadius();
	double r_second = reaction.second.getReactionRadius();

	// Get the diffusion coefficients
	double firstDiffusion = reaction.first.getDiffusionCoefficient();
	double secondDiffusion = reaction.second.getDiffusionCoefficient();

	// Calculate and return
	double k_plus = 4.0 * xolotlCore::pi
			* (r_first + r_second + xolotlCore::reactionRadius)
			* (firstDiffusion + secondDiffusion);
	return k_plus;
}

void ReactionNetwork::fillConcentrationsArray(double * concentrations) {

	// Fill the array
	std::for_each(allReactants.begin(), allReactants.end(),
			[&concentrations](const IReactant& currReactant) {
				auto id = currReactant.getId() - 1;
				concentrations[id] = currReactant.getConcentration();
			});

	return;
}

void ReactionNetwork::updateConcentrationsFromArray(double * concentrations) {

	std::for_each(allReactants.begin(), allReactants.end(),
			[&concentrations](IReactant& currReactant) {
				auto id = currReactant.getId() - 1;
				currReactant.setConcentration(concentrations[id]);
			});

	return;
}

void ReactionNetwork::setTemperature(double temp) {
	// Set the temperature
	temperature = temp;

	// Update the temperature for all of the clusters
	std::for_each(allReactants.begin(), allReactants.end(),
			[&temp](IReactant& currReactant) {

				// This part will set the temperature in each reactant
				// and recompute the diffusion coefficient
				currReactant.setTemperature(temp);
			});

	return;
}

#if READY
static std::shared_ptr<xolotlPerf::IEventCounter> addPRCounter_calls;
static std::shared_ptr<xolotlPerf::IEventCounter> addPRCounter_adds;
static std::shared_ptr<xolotlPerf::IEventCounter> addDRCounter_calls;
static std::shared_ptr<xolotlPerf::IEventCounter> addDRCounter_adds;
#endif // READY

ProductionReaction& ReactionNetwork::add(
		std::unique_ptr<ProductionReaction> reaction) {

#if READY
	if(not addPRCounter_calls) {
		addPRCounter_calls = handlerRegistry->getEventCounter("addPR_calls");
	}
	addPRCounter_calls->increment();
#endif // READY

	// Ensure we know about the reaction.
	// Map's emplace() returns a pair (iter, bool) where
	// iter points to the item in the map and the bool indicates
	// whether it was added by this emplace() call.
	auto key = reaction->descriptiveKey();
	auto eret = productionReactionMap.emplace(key, std::move(reaction));

#if READY
	if(eret.second == true) {
		if(not addPRCounter_adds) {
			addPRCounter_adds = handlerRegistry->getEventCounter("addPR_adds");
		}
		addPRCounter_adds->increment();
	}
#endif // READY

	// Regardless of whether we added it in this emplace() call or not,
	// the iter within eret refers to the desired reaction in the map.
	return *(eret.first->second);
}

DissociationReaction& ReactionNetwork::add(
		std::unique_ptr<DissociationReaction> reaction) {

#if READY
	if(not addDRCounter_calls) {
		addDRCounter_calls = handlerRegistry->getEventCounter("addDR_calls");
	}
	addDRCounter_calls->increment();
#endif // READY

	// Unlike addProductionReaction, we use a check-add approach
	// instead of emplace() because if the item isn't already in
	// the network, we will need to contruct its reverse reaction.

	// Check if we already know about this reaction.
	auto key = reaction->descriptiveKey();
	auto iter = dissociationReactionMap.find(key);
	if (iter != dissociationReactionMap.end()) {
		// We already knew about the reaction.
		// Return the existing one.
		return *(iter->second);
	}

	// Add the dissociation reaction to our set of known reactions.
	auto eret = dissociationReactionMap.emplace(key, std::move(reaction));

#if READY
	if(eret.second == true) {
		if(not addDRCounter_adds) {
			addDRCounter_adds = handlerRegistry->getEventCounter("addDR_adds");
		}
		addDRCounter_adds->increment();
	}
#endif // READY

	// Since we checked earlier and the reaction wasn't in the map,
	// our emplace() call should have added it.
	assert(eret.second);

	// Return the newly-added dissociation reaction.
	return *(eret.first->second);
}

void ReactionNetwork::removeReactants(
		const IReactionNetwork::ReactantVector& doomedReactants) {

	// Build a ReactantMatcher functor for the doomed reactants.
	// Doing this here allows us to construct the canonical composition
	// strings for the doomed reactants once and reuse them.
	// If we used an anonymous functor object in the std::remove_if
	// calls we would build these strings several times in this function.
	ReactionNetwork::ReactantMatcher doomedReactantMatcher(doomedReactants);

	// Remove the doomed reactants from our collection of all known reactants.
	auto ariter = std::remove_if(allReactants.begin(), allReactants.end(),
			doomedReactantMatcher);
	allReactants.erase(ariter, allReactants.end());

	// Remove the doomed reactants from the type-specific cluster vectors.
	// First, determine all cluster types used by clusters in the collection
	// of doomed reactants...
	std::set<ReactantType> typesUsed;
	for (IReactant const& reactant : doomedReactants) {
		typesUsed.insert(reactant.getType());
	}

	// ...Next, examine each type's collection of clusters and remove the
	// doomed reactants.
	for (auto currType : typesUsed) {
		auto& clusters = clusterTypeMap[currType];

		for (IReactant const& currDoomedReactant : doomedReactants) {
			auto iter = clusters.find(currDoomedReactant.getComposition());
			assert(iter != clusters.end());
			clusters.erase(iter);
		}
	}

	return;
}

void ReactionNetwork::dumpTo(std::ostream& os) const {

	// Dump flat view of reactants.
	os << size() << " reactants:\n";
	for (auto const& currReactant : allReactants) {
		os << currReactant << '\n';
	}

#if READY
	// Dump reactants of each type.
	// TODO what does this give us that the flat view doesn't?
	os << "per-type reactant map:\n";
	auto const& knownTypes = getKnownReactantTypes();
	for (auto const& currType : knownTypes) {
		auto const& currTypeReactantMap = clusterTypeMap.at(currType);
		os << currTypeReactantMap.size() << " " << toString(currType) << " reactants:\n";
		for (auto const& currMapItem : currTypeReactantMap) {
			os << *(currMapItem.second) << '\n';
		}
	}
#endif // READY

	// Dump ProductionReactions.
	os << productionReactionMap.size() << " production reactions:\n";
	for (auto const& currMapItem : productionReactionMap) {
		os << *(currMapItem.second) << '\n';
	}

	// Dump DissociationReactions.
	os << dissociationReactionMap.size() << " dissociation reactions:\n";
	for (auto const& currMapItem : dissociationReactionMap) {
		os << *(currMapItem.second) << '\n';
	}

	// For each reactant, dump coefficients it uses for reactions it
	// participates in.
	os << size() << " reactant coefficients:\n";
	for (IReactant const& currReactant : allReactants) {
		currReactant.outputCoefficientsTo(os);
	}
}

std::map<std::string, int> ReactionNetwork::getIDMap() const {
	// Create the map
	std::map<std::string, int> idMap;

	// Fill it
	for (IReactant const& currReactant : allReactants) {
		idMap[currReactant.getName()] = currReactant.getId() - 1;
	}

	return idMap;
}

std::pair<int, int> ReactionNetwork::getHighestClusterCoordinates() const {
	// Create the pair
	std::pair<int, int> coor = std::make_pair(0, 0);
	// Initialize the highest concentration
	double highestConc = 0.0;

	// Loop on the reactants
	for (IReactant const& currReactant : allReactants) {
		// Get its concentration
		double conc = currReactant.getConcentration();
		// Compare it to the highest conc
		if (conc > highestConc) {
			highestConc = conc;
			// Get its composition
			auto comp = currReactant.getComposition();
			coor.first = comp[toCompIdx(Species::He)];
			coor.second = comp[toCompIdx(Species::V)];
		}
	}

	return coor;
}

} // xolotlCore
