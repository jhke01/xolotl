#include "PSIClusterReactionNetwork.h"
#include "PSICluster.h"
#include "PSISuperCluster.h"
#include <xolotlPerf.h>
#include <Constants.h>
#include <MathUtils.h>

using namespace xolotlCore;

void PSIClusterReactionNetwork::setDefaultPropsAndNames() {

	// Initialize default properties
	dissociationsEnabled = true;
	numHeClusters = 0;
	numVClusters = 0;
	numIClusters = 0;
	numHeVClusters = 0;
	numHeIClusters = 0;
	numSuperClusters = 0;
	maxHeClusterSize = 0;
	maxVClusterSize = 0;
	maxIClusterSize = 0;
	maxHeVClusterSize = 0;
	maxHeIClusterSize = 0;

	// Initialize the current and last size to 0
	networkSize = 0;
	// Set the reactant names
	names.push_back(heType);
	names.push_back(vType);
	names.push_back(iType);
	// Set the compound reactant names
	compoundNames.push_back(heVType);
	compoundNames.push_back(heIType);
	compoundNames.push_back(PSISuperType);

	// Specify cluster types we know about.
    for (auto& currName : names) {
        clusterTypeMap.insert( { currName, ReactionNetwork::ReactantVector() } );
    }
    for (auto& currName : compoundNames) {
        clusterTypeMap.insert( { currName, ReactionNetwork::ReactantVector() } );
    }

	return;
}

PSIClusterReactionNetwork::PSIClusterReactionNetwork() :
		ReactionNetwork() {
	// Setup the properties map and the name lists
	setDefaultPropsAndNames();

	return;
}

PSIClusterReactionNetwork::PSIClusterReactionNetwork(
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		ReactionNetwork(registry) {
	// Setup the properties map and the name lists
	setDefaultPropsAndNames();

	return;
}

PSIClusterReactionNetwork::PSIClusterReactionNetwork(
		const PSIClusterReactionNetwork &other) :
		ReactionNetwork(other) {
	// The size and ids do not need to be copied. They will be fixed when the
	// reactants are added.

	// Reset the properties table so that it can be properly updated when the
	// network is filled.
	setDefaultPropsAndNames();
	// Get all of the reactants from the other network and add them to this one
	// Load the single-species clusters. Calling getAll() will not work because
	// it is not const.
	std::vector<std::shared_ptr<IReactant> > reactants;
	for (auto it = other.singleSpeciesMap.begin();
			it != other.singleSpeciesMap.end(); ++it) {
		reactants.push_back(it->second);
	}
	// Load the mixed-species clusters
	for (auto it = other.mixedSpeciesMap.begin();
			it != other.mixedSpeciesMap.end(); ++it) {
		reactants.push_back(it->second);
	}
	// Load the super-species clusters
	for (auto it = other.superSpeciesMap.begin();
			it != other.superSpeciesMap.end(); ++it) {
		reactants.push_back(it->second);
	}
	for (unsigned int i = 0; i < reactants.size(); i++) {
		add(reactants[i]->clone());
	}

	return;
}

double PSIClusterReactionNetwork::calculateDissociationConstant(
		DissociationReaction * reaction) const {
	// If the dissociations are not allowed
	if (!dissociationsEnabled)
		return 0.0;

	// The atomic volume is computed by considering the BCC structure of the
	// tungsten. In a given lattice cell in tungsten there are tungsten atoms
	// at each corner and a tungsten atom in the center. The tungsten atoms at
	// the corners are shared across a total of eight cells. The fraction of
	// the volume of the lattice cell that is filled with tungsten atoms is the
	// atomic volume and is a_0^3/(8*1/8 + 1) = 0.5*a_0^3.
	double atomicVolume = 0.5 * xolotlCore::tungstenLatticeConstant
			* xolotlCore::tungstenLatticeConstant
			* xolotlCore::tungstenLatticeConstant;

	// Get the rate constant from the reverse reaction
	double kPlus = reaction->reverseReaction->kConstant;

	// Calculate and return
	double bindingEnergy = computeBindingEnergy(reaction);
	double k_minus_exp = exp(
			-1.0 * bindingEnergy / (xolotlCore::kBoltzmann * temperature));
	double k_minus = (1.0 / atomicVolume) * kPlus * k_minus_exp;

	return k_minus;
}

void PSIClusterReactionNetwork::createReactionConnectivity() {
	// Initial declarations
	int firstSize = 0, secondSize = 0, productSize = 0;

	// Single species clustering (He, V, I)
	// We know here that only Xe_1 can cluster so we simplify the search
	// X_(a-i) + X_i --> X_a
	// Make a vector of types
	std::vector<string> typeVec = { heType, vType, iType };
	// Loop on it
	for (auto tvIter = typeVec.begin(); tvIter != typeVec.end(); ++tvIter) {
		string typeName = *tvIter;

		// Get all the reactants of this type
		auto allTypeReactants = getAll(typeName);
		// Loop on them
		for (auto firstIt = allTypeReactants.begin();
				firstIt != allTypeReactants.end(); firstIt++) {
			// Get its size
			firstSize = (*firstIt)->getSize();
			// Loop on the second cluster starting at the same pointer to avoid double counting
			for (auto secondIt = firstIt; secondIt != allTypeReactants.end();
					secondIt++) {
				// Get its size
				secondSize = (*secondIt)->getSize();
				productSize = firstSize + secondSize;
				// Get the product
				auto product = get(typeName, productSize);
				// Check that the reaction can occur
				if (product
						&& ((*firstIt)->getDiffusionFactor() > 0.0
								|| (*secondIt)->getDiffusionFactor() > 0.0)) {
					// Create a production reaction
					auto reaction = std::make_shared<ProductionReaction>(
							(*firstIt), (*secondIt));
					// Tell the reactants that they are in this reaction
					(*firstIt)->createCombination(reaction);
					(*secondIt)->createCombination(reaction);
					product->createProduction(reaction);

					// Check if the reverse reaction is allowed
					checkDissociationConnectivity(product, reaction);
				}
			}
		}
	}

	// Helium absorption by HeV clusters
	// He_(a) + (He_b)(V_c) --> [He_(a+b)](V_c)
	// Get all the He and HeV clusters
	auto allHeReactants = getAll(heType);
	auto allHeVReactants = getAll(heVType);
	auto allSuperReactants = getAll(PSISuperType);
	// Loop on the He clusters
	for (auto firstIt = allHeReactants.begin(); firstIt != allHeReactants.end();
			firstIt++) {
		// Skip if it can't diffuse
		if (xolotlCore::equal((*firstIt)->getDiffusionFactor(), 0.0))
			continue;
		// Get its size
		firstSize = (*firstIt)->getSize();
		// Loop on the HeV clusters
		for (auto secondIt = allHeVReactants.begin();
				secondIt != allHeVReactants.end(); secondIt++) {
			// Get its composition
			auto comp = (*secondIt)->getComposition();
			// Create the composition of the potential product
			std::vector<int> compositionVec = { comp[heType] + firstSize,
					comp[vType], 0 };
			// Get the product
			auto product = getCompound(heVType, compositionVec);

			// Check if the product can be a super cluster
			if (!product) {
				// Check if it is a super cluster from the map
				product = getSuperFromComp(compositionVec[0], compositionVec[1]);
			}

			// Check that the reaction can occur
			if (product
					&& ((*firstIt)->getDiffusionFactor() > 0.0
							|| (*secondIt)->getDiffusionFactor() > 0.0)) {
				// Create a production reaction
				auto reaction = std::make_shared<ProductionReaction>((*firstIt),
						(*secondIt));
				// Tell the reactants that they are in this reaction
				(*firstIt)->createCombination(reaction);
				(*secondIt)->createCombination(reaction);
				product->createProduction(reaction, compositionVec[0],
						compositionVec[1]);

				// Check if the reverse reaction is allowed
				checkDissociationConnectivity(product, reaction,
						compositionVec[0], compositionVec[1]);
			}
		}

		// Loop on the super clusters
		for (auto secondIt = allSuperReactants.begin();
				secondIt != allSuperReactants.end(); secondIt++) {
			auto superCluster = (PSISuperCluster *) *secondIt;
			IReactant * product = nullptr;
			// Get its boundaries
			auto boundaries = superCluster->getBoundaries();
			// Loop on them
			for (int i = boundaries[0]; i <= boundaries[1]; i++) {
				for (int j = boundaries[2]; j <= boundaries[3]; j++) {
					// Assume the product can only be a super cluster here
					product = getSuperFromComp(i + firstSize, j);

					// Check that the reaction can occur
					if (product
							&& ((*firstIt)->getDiffusionFactor() > 0.0
									|| (*secondIt)->getDiffusionFactor() > 0.0)) {
						// Create a production reaction
						auto reaction = std::make_shared<ProductionReaction>(
								(*firstIt), (*secondIt));
						// Tell the reactants that they are in this reaction
						(*firstIt)->createCombination(reaction, i, j);
						(*secondIt)->createCombination(reaction, i, j);
						product->createProduction(reaction, i + firstSize, j, i,
								j);

						// Check if the reverse reaction is allowed
						checkDissociationConnectivity(product, reaction,
								i + firstSize, j, i, j);
					}
				}
			}
		}
	}

	// Vacancy absorption by HeV clusters
	// (He_a)(V_b) + V_c --> (He_a)[V_(b+c)]
	// Get all the V clusters
	auto allVReactants = getAll(vType);
	// Loop on the V clusters
	// Loop on the HeV clusters
	for (auto firstIt = allVReactants.begin(); firstIt != allVReactants.end();
			firstIt++) {
		// Skip if it can't diffuse
		if (xolotlCore::equal((*firstIt)->getDiffusionFactor(), 0.0))
			continue;
		// Get the V size
		firstSize = (*firstIt)->getSize();
		// Loop on the HeV clusters
		for (auto secondIt = allHeVReactants.begin();
				secondIt != allHeVReactants.end(); secondIt++) {
			// Get its composition
			auto comp = (*secondIt)->getComposition();
			// Create the composition of the potential product
			std::vector<int> compositionVec = { comp[heType], comp[vType]
					+ firstSize, 0 };
			// Get the product
			auto product = getCompound(heVType, compositionVec);

			// Check if the product can be a super cluster
			if (!product) {
				product = getSuperFromComp(compositionVec[0], compositionVec[1]);
			}

			// Check that the reaction can occur
			if (product
					&& ((*firstIt)->getDiffusionFactor() > 0.0
							|| (*secondIt)->getDiffusionFactor() > 0.0)) {
				// Create a production reaction
				auto reaction = std::make_shared<ProductionReaction>((*firstIt),
						(*secondIt));
				// Tell the reactants that they are in this reaction
				(*firstIt)->createCombination(reaction);
				(*secondIt)->createCombination(reaction);
				product->createProduction(reaction, compositionVec[0],
						compositionVec[1]);

				// Check if the reverse reaction is allowed
				checkDissociationConnectivity(product, reaction,
						compositionVec[0], compositionVec[1]);
			}
		}

		// Loop on the super clusters
		for (auto secondIt = allSuperReactants.begin();
				secondIt != allSuperReactants.end(); secondIt++) {
			auto superCluster = (PSISuperCluster *) *secondIt;
			IReactant * product = nullptr;
			// Get its boundaries
			auto boundaries = superCluster->getBoundaries();
			// Loop on them
			for (int i = boundaries[0]; i <= boundaries[1]; i++) {
				for (int j = boundaries[2]; j <= boundaries[3]; j++) {
					// Assume the product can only be a super cluster here
					product = getSuperFromComp(i, j + firstSize);

					// Check that the reaction can occur
					if (product
							&& ((*firstIt)->getDiffusionFactor() > 0.0
									|| (*secondIt)->getDiffusionFactor() > 0.0)) {
						// Create a production reaction
						auto reaction = std::make_shared<ProductionReaction>(
								(*firstIt), (*secondIt));
						// Tell the reactants that they are in this reaction
						(*firstIt)->createCombination(reaction, i, j);
						(*secondIt)->createCombination(reaction, i, j);
						product->createProduction(reaction, i, j + firstSize, i,
								j);

						// Check if the reverse reaction is allowed
						checkDissociationConnectivity(product, reaction, i,
								j + firstSize, i, j);
					}
				}
			}
		}
	}

	// Helium-Vacancy clustering
	// He_a + V_b --> (He_a)(V_b)
	// Loop on the He clusters
	for (auto firstIt = allHeReactants.begin(); firstIt != allHeReactants.end();
			firstIt++) {
		// Get its size
		firstSize = (*firstIt)->getSize();
		// Loop on the HeV clusters
		for (auto secondIt = allVReactants.begin();
				secondIt != allVReactants.end(); secondIt++) {
			// Get its size
			secondSize = (*secondIt)->getSize();
			// Create the composition of the potential product
			std::vector<int> compositionVec = { firstSize, secondSize, 0 };
			// Get the product
			auto product = getCompound(heVType, compositionVec);

			// Check if the product can be a super cluster
			if (!product) {
				product = getSuperFromComp(compositionVec[0], compositionVec[1]);
			}

			// Check that the reaction can occur
			if (product
					&& ((*firstIt)->getDiffusionFactor() > 0.0
							|| (*secondIt)->getDiffusionFactor() > 0.0)) {
				// Create a production reaction
				auto reaction = std::make_shared<ProductionReaction>((*firstIt),
						(*secondIt));
				// Tell the reactants that they are in this reaction
				(*firstIt)->createCombination(reaction);
				(*secondIt)->createCombination(reaction);
				product->createProduction(reaction, compositionVec[0],
						compositionVec[1]);

				// Check if the reverse reaction is allowed
				checkDissociationConnectivity(product, reaction,
						compositionVec[0], compositionVec[1]);
			}
		}
	}

	// Vacancy reduction by Interstitial absorption in HeV clusters
	// (He_a)(V_b) + (I_c) --> (He_a)[V_(b-c)]
	// Get all the I clusters
	auto allIReactants = getAll(iType);
	// Loop on them
	for (auto firstIt = allIReactants.begin(); firstIt != allIReactants.end();
			firstIt++) {
		// Get its size
		firstSize = (*firstIt)->getSize();
		// Loop on the HeV clusters
		for (auto secondIt = allHeVReactants.begin();
				secondIt != allHeVReactants.end(); secondIt++) {
			// Get its composition
			auto comp = (*secondIt)->getComposition();
			// The product can be He or HeV
			IReactant * product = nullptr;
			if (comp[vType] == firstSize) {
				// The product is He
				product = get(heType, comp[heType]);
			} else {
				// The product is HeV
				// Create the composition of the potential product
				std::vector<int> compositionVec = { comp[heType], comp[vType]
						- firstSize, 0 };
				// Get the product
				product = getCompound(heVType, compositionVec);
			}
			// Check that the reaction can occur
			if (product
					&& ((*firstIt)->getDiffusionFactor() > 0.0
							|| (*secondIt)->getDiffusionFactor() > 0.0)) {
				// Create a production reaction
				auto reaction = std::make_shared<ProductionReaction>((*firstIt),
						(*secondIt));
				// Tell the reactants that they are in this reaction
				(*firstIt)->createCombination(reaction);
				(*secondIt)->createCombination(reaction);
				product->createProduction(reaction);

				// Check if the reverse reaction is allowed
				checkDissociationConnectivity(product, reaction);
			}
		}

		// Loop on the super clusters
		for (auto secondIt = allSuperReactants.begin();
				secondIt != allSuperReactants.end(); secondIt++) {
			auto superCluster = (PSISuperCluster *) *secondIt;
			IReactant * product = nullptr;
			// Get its boundaries
			auto boundaries = superCluster->getBoundaries();
			// Loop on them
			for (int i = boundaries[0]; i <= boundaries[1]; i++) {
				for (int j = boundaries[2]; j <= boundaries[3]; j++) {
					// The product might be HeV or He
					// Get the product
					if (j == firstSize) {
						// The product is He
						product = get(heType, i);
					}
					else {
						// Create the composition of the potential product
						std::vector<int> compositionVec = { i, j - firstSize, 0 };
						product = getCompound(heVType, compositionVec);

						// If the product doesn't exist check for super clusters
						if (!product) {
							product = getSuperFromComp(i, j - firstSize);
						}
					}

					// Check that the reaction can occur
					if (product
							&& ((*firstIt)->getDiffusionFactor() > 0.0
									|| (*secondIt)->getDiffusionFactor() > 0.0)) {
						// Create a production reaction
						auto reaction = std::make_shared<ProductionReaction>(
								(*firstIt), (*secondIt));
						// Tell the reactants that they are in this reaction
						(*firstIt)->createCombination(reaction, i, j);
						(*secondIt)->createCombination(reaction, i, j);
						product->createProduction(reaction, i, j - firstSize, i,
								j);

						// Check if the reverse reaction is allowed
						checkDissociationConnectivity(product, reaction, i,
								j - firstSize, i, j);
					}
				}
			}
		}
	}

//	// Helium clustering leading to trap mutation
//	// He_a + He_b --> [He_(a+b)](V_c) + I_c
//	// Loop on the He clusters
//	for (auto firstIt = allHeReactants.begin(); firstIt != allHeReactants.end();
//			firstIt++) {
//		// Get its size
//		firstSize = (*firstIt)->getSize();
//		// Loop on the second He cluster starting at the same pointer to avoid double counting
//		for (auto secondIt = firstIt; secondIt != allHeReactants.end();
//				secondIt++) {
//			// Get its size
//			secondSize = (*secondIt)->getSize();
//			// Get the simple product
//			productSize = firstSize + secondSize;
//			auto product = get(heType, productSize);
//			// Doesn't do anything if the product exist
//			if (product)
//				continue;
//
//			// Trap mutation is happening
//			// Loop on the possible I starting by the smallest
//			for (auto it = allIReactants.begin(); it != allIReactants.end();
//					it++) {
//				// Get the size of the I cluster
//				int iSize = (*it)->getSize();
//				// Create the composition of the potential product
//				std::vector<int> compositionVec = { firstSize + secondSize,
//						iSize, 0 };
//				product = getCompound(heVType, compositionVec);
//				// Check that the reaction can occur
//				if (product
//						&& ((*firstIt)->getDiffusionFactor() > 0.0
//								|| (*secondIt)->getDiffusionFactor() > 0.0)) {
//					// Create a production reaction
//					auto reaction = std::make_shared<ProductionReaction>(
//							(*firstIt), (*secondIt));
//					// Tell the reactants that they are in this reaction
//					(*firstIt)->createCombination(reaction);
//					(*secondIt)->createCombination(reaction);
//					product->createProduction(reaction);
//					(*it)->createProduction(reaction);
//
//					// Stop the loop on I clusters here
//					break;
//				}
//			}
//		}
//	}

//	// Helium absorption by HeV leading to trap mutation
//	// (He_a)(V_b) + He_c --> [He_(a+c)][V_(b+d)] + I_d
//	// Loop on the He clusters
//	for (auto firstIt = allHeReactants.begin(); firstIt != allHeReactants.end();
//			firstIt++) {
//		// Get its size
//		firstSize = (*firstIt)->getSize();
//		// Loop on the HeV clusters
//		for (auto secondIt = allHeVReactants.begin();
//				secondIt != allHeVReactants.end(); secondIt++) {
//			// Get its composition
//			auto comp = (*secondIt)->getComposition();
//			// Get the simple product
//			std::vector<int> compositionVec = { firstSize + comp[heType],
//					comp[vType], 0 };
//			auto product = getCompound(heVType, compositionVec);
//			// Doesn't do anything if the product exist
//			if (product)
//				continue;
//
//			// Trap mutation is happening
//			// Loop on the possible I starting by the smallest
//			for (auto it = allIReactants.begin(); it != allIReactants.end();
//					it++) {
//				// Get the size of the I cluster
//				int iSize = (*it)->getSize();
//				// Create the composition of the potential product
//				compositionVec[1] = comp[vType] + iSize;
//				product = getCompound(heVType, compositionVec);
//				// Check that the reaction can occur
//				if (product
//						&& ((*firstIt)->getDiffusionFactor() > 0.0
//								|| (*secondIt)->getDiffusionFactor() > 0.0)) {
//					// Create a production reaction
//					auto reaction = std::make_shared<ProductionReaction>(
//							(*firstIt), (*secondIt));
//					// Tell the reactants that they are in this reaction
//					(*firstIt)->createCombination(reaction);
//					(*secondIt)->createCombination(reaction);
//					product->createProduction(reaction);
//					(*it)->createProduction(reaction);
//
//					// Stop the loop on I clusters here
//					break;
//				}
//			}
//		}
//	}

	// Vacancy-Interstitial annihilation
	// I_a + V_b
	//        --> I_(a-b), if a > b
	//        --> V_(b-a), if a < b
	//        --> 0, if a = b
	// Loop on the I clusters
	for (auto firstIt = allIReactants.begin(); firstIt != allIReactants.end();
			firstIt++) {
		// Get its size
		firstSize = (*firstIt)->getSize();
		// Loop on the V clusters
		for (auto secondIt = allVReactants.begin();
				secondIt != allVReactants.end(); secondIt++) {
			// Get its size
			secondSize = (*secondIt)->getSize();
			// Check the possibilities
			if (firstSize > secondSize) {
				// Get the product
				productSize = firstSize - secondSize;
				auto product = get(iType, productSize);
				// Check that the reaction can occur
				if (product
						&& ((*firstIt)->getDiffusionFactor() > 0.0
								|| (*secondIt)->getDiffusionFactor() > 0.0)) {
					// Create a production reaction
					auto reaction = std::make_shared<ProductionReaction>(
							(*firstIt), (*secondIt));
					// Tell the reactants that they are in this reaction
					(*firstIt)->createCombination(reaction);
					(*secondIt)->createCombination(reaction);
					product->createProduction(reaction);
				}
			} else if (firstSize < secondSize) {
				// Get the product
				productSize = secondSize - firstSize;
				auto product = get(vType, productSize);
				// Check that the reaction can occur
				if (product
						&& ((*firstIt)->getDiffusionFactor() > 0.0
								|| (*secondIt)->getDiffusionFactor() > 0.0)) {
					// Create a production reaction
					auto reaction = std::make_shared<ProductionReaction>(
							(*firstIt), (*secondIt));
					// Tell the reactants that they are in this reaction
					(*firstIt)->createCombination(reaction);
					(*secondIt)->createCombination(reaction);
					product->createProduction(reaction);
				}

			} else {
				// Annihilation
				// Check that the reaction can occur
				if (((*firstIt)->getDiffusionFactor() > 0.0
						|| (*secondIt)->getDiffusionFactor() > 0.0)) {
					// Create a production reaction
					auto reaction = std::make_shared<ProductionReaction>(
							(*firstIt), (*secondIt));
					// Tell the reactants that they are in this reaction
					(*firstIt)->createCombination(reaction);
					(*secondIt)->createCombination(reaction);
				}
			}
		}
	}

//	// Helium absorption by HeI clusters
//	// He_(a) + (He_b)(I_c) --> [He_(a+b)](I_c)
//	// Get all the HeI clusters
//	auto allHeIReactants = getAll(heIType);
//	// Loop on the He clusters
//	for (auto firstIt = allHeReactants.begin(); firstIt != allHeReactants.end();
//			firstIt++) {
//		// Get its size
//		firstSize = (*firstIt)->getSize();
//		// Loop on the HeV clusters
//		for (auto secondIt = allHeIReactants.begin();
//				secondIt != allHeIReactants.end(); secondIt++) {
//			// Get its composition
//			auto comp = (*secondIt)->getComposition();
//			// Create the composition of the potential product
//			std::vector<int> compositionVec = { comp[heType] + firstSize, 0,
//					comp[iType] };
//			// Get the product
//			auto product = getCompound(heIType, compositionVec);
//			// Check that the reaction can occur
//			if (product
//					&& ((*firstIt)->getDiffusionFactor() > 0.0
//							|| (*secondIt)->getDiffusionFactor() > 0.0)) {
//				// Create a production reaction
//				auto reaction = std::make_shared<ProductionReaction>((*firstIt),
//						(*secondIt));
//				// Tell the reactants that they are in this reaction
//				(*firstIt)->createCombination(reaction);
//				(*secondIt)->createCombination(reaction);
//				product->createProduction(reaction);
//
//				// Check if the reverse reaction is allowed
//				checkDissociationConnectivity(product, reaction);
//			}
//		}
//	}
//
//	// Single Interstitial absorption by HeI clusters
//	// (He_a)(I_b) + I --> (He_a)[I_(b+1)]
//	// Get the single interstitial cluster
//	auto singleInterstitialCluster = get(iType, 1);
//	// Loop on the HeI clusters
//	for (auto secondIt = allHeIReactants.begin();
//			secondIt != allHeIReactants.end(); secondIt++) {
//		// Get its composition
//		auto comp = (*secondIt)->getComposition();
//		// Create the composition of the potential product
//		std::vector<int> compositionVec = { comp[heType], 0, comp[iType] + 1 };
//		// Get the product
//		auto product = getCompound(heIType, compositionVec);
//		// Check that the reaction can occur
//		if (product
//				&& (singleInterstitialCluster->getDiffusionFactor() > 0.0
//						|| (*secondIt)->getDiffusionFactor() > 0.0)) {
//			// Create a production reaction
//			auto reaction = std::make_shared<ProductionReaction>(
//					singleInterstitialCluster, (*secondIt));
//			// Tell the reactants that they are in this reaction
//			singleInterstitialCluster->createCombination(reaction);
//			(*secondIt)->createCombination(reaction);
//			product->createProduction(reaction);
//
//			// Check if the reverse reaction is allowed
//			checkDissociationConnectivity(product, reaction);
//		}
//	}
//
//	// Helium-Interstitial clustering
//	// He_a + I_b --> (He_a)(I_b)
//	// Loop on the He clusters
//	for (auto firstIt = allHeReactants.begin(); firstIt != allHeReactants.end();
//			firstIt++) {
//		// Get its size
//		firstSize = (*firstIt)->getSize();
//		// Loop on the I clusters
//		for (auto secondIt = allIReactants.begin();
//				secondIt != allIReactants.end(); secondIt++) {
//			// Get its size
//			secondSize = (*secondIt)->getSize();
//			// Create the composition of the potential product
//			std::vector<int> compositionVec = { firstSize, 0, secondSize };
//			// Get the product
//			auto product = getCompound(heIType, compositionVec);
//			// Check that the reaction can occur
//			if (product
//					&& ((*firstIt)->getDiffusionFactor() > 0.0
//							|| (*secondIt)->getDiffusionFactor() > 0.0)) {
//				// Create a production reaction
//				auto reaction = std::make_shared<ProductionReaction>((*firstIt),
//						(*secondIt));
//				// Tell the reactants that they are in this reaction
//				(*firstIt)->createCombination(reaction);
//				(*secondIt)->createCombination(reaction);
//				product->createProduction(reaction);
//
//				// Check if the reverse reaction is allowed
//				checkDissociationConnectivity(product, reaction);
//			}
//		}
//	}
//
//	// Interstitial reduction by Vacancy absorption in HeI clusters
//	// (He_a)(I_b) + (V_c) --> (He_a)[I_(b-c)]
//	// Loop on V clusters
//	for (auto firstIt = allVReactants.begin(); firstIt != allVReactants.end();
//			firstIt++) {
//		// Get its size
//		firstSize = (*firstIt)->getSize();
//		// Loop on the HeI clusters
//		for (auto secondIt = allHeIReactants.begin();
//				secondIt != allHeIReactants.end(); secondIt++) {
//			// Get its composition
//			auto comp = (*secondIt)->getComposition();
//			// The product can be He or HeI
//			IReactant * product = nullptr;
//			if (comp[iType] == firstSize) {
//				// The product is He
//				product = get(heType, comp[heType]);
//			} else {
//				// The product is HeI
//				// Create the composition of the potential product
//				std::vector<int> compositionVec = { comp[heType], 0, comp[iType]
//						- firstSize };
//				// Get the product
//				product = getCompound(heIType, compositionVec);
//			}
//			// Check that the reaction can occur
//			if (product
//					&& ((*firstIt)->getDiffusionFactor() > 0.0
//							|| (*secondIt)->getDiffusionFactor() > 0.0)) {
//				// Create a production reaction
//				auto reaction = std::make_shared<ProductionReaction>((*firstIt),
//						(*secondIt));
//				// Tell the reactants that they are in this reaction
//				(*firstIt)->createCombination(reaction);
//				(*secondIt)->createCombination(reaction);
//				product->createProduction(reaction);
//
//				// Check if the reverse reaction is allowed
//				checkDissociationConnectivity(product, reaction);
//			}
//		}
//	}

	return;
}

void PSIClusterReactionNetwork::checkDissociationConnectivity(
		IReactant * emittingReactant,
		std::shared_ptr<ProductionReaction> reaction, int a, int b, int c,
		int d) {
	// Check if at least one of the potentially emitted cluster is size one
	if (reaction->first->getSize() != 1 && reaction->second->getSize() != 1) {
		// Don't add the reverse reaction
		return;
	}
	// remove He+He
	if (reaction->first->getSize() == 1 && reaction->second->getSize() == 1
			&& reaction->first->getType() == heType
			&& reaction->second->getType() == heType) {
		// Don't add the reverse reaction
		return;
	}

//	// Check for trap mutations (with XOR)
//	if ((reaction->first->getType() == iType)
//			== !(reaction->second->getType() == iType)) {
//		// Don't add the reverse reaction
//		return;
//	}

	// The reaction can occur, create the dissociation
	// Create a dissociation reaction
	auto dissociationReaction = std::make_shared<DissociationReaction>(
			emittingReactant, reaction->first, reaction->second);
	// Set the reverse reaction
	dissociationReaction->reverseReaction = reaction.get();
	// Tell the reactants that their are in this reaction
	reaction->first->createDissociation(dissociationReaction, a, b, c, d);
	reaction->second->createDissociation(dissociationReaction, a, b, c, d);
	emittingReactant->createEmission(dissociationReaction, a, b, c, d);

	return;
}

void PSIClusterReactionNetwork::setTemperature(double temp) {
	ReactionNetwork::setTemperature(temp);

	computeRateConstants();

	return;
}

double PSIClusterReactionNetwork::getTemperature() const {
	return temperature;
}

IReactant * PSIClusterReactionNetwork::get(const std::string& type,
		const int size) const {
	// Local Declarations
	static std::map<std::string, int> composition = { { heType, 0 },
			{ vType, 0 }, { iType, 0 }, { xeType, 0 } };
	std::shared_ptr<IReactant> retReactant;

	// Initialize the values because it's static
	composition[heType] = 0;
	composition[vType] = 0;
	composition[iType] = 0;

	// Only pull the reactant if the name and size are valid
	if ((type == heType || type == vType || type == iType) && size >= 1) {
		composition[type] = size;
		//std::string encodedName = PSICluster::encodeCompositionAsName(composition);
		// Make sure the reactant is in the map
		std::string compStr = Reactant::toCanonicalString(type, composition);
		if (singleSpeciesMap.count(compStr)) {
			retReactant = singleSpeciesMap.at(compStr);
		}
	}

	return retReactant.get();
}

IReactant * PSIClusterReactionNetwork::getCompound(const std::string& type,
		const std::vector<int>& sizes) const {
	// Local Declarations
	static std::map<std::string, int> composition = { { heType, 0 },
			{ vType, 0 }, { iType, 0 }, { xeType, 0 } };
	std::shared_ptr<IReactant> retReactant;

	// Initialize the values because it's static
	composition[heType] = 0;
	composition[vType] = 0;
	composition[iType] = 0;

	// Only pull the reactant if the name is valid and there are enough sizes
	// to fill the composition.
	if ((type == heVType || type == heIType) && sizes.size() == 3) {
		composition[heType] = sizes[0];
		composition[vType] = sizes[1];
		composition[iType] = sizes[2];

		// Make sure the reactant is in the map
		std::string compStr = Reactant::toCanonicalString(type, composition);
		if (mixedSpeciesMap.count(compStr)) {
			retReactant = mixedSpeciesMap.at(compStr);
		}
	}

	return retReactant.get();
}

IReactant * PSIClusterReactionNetwork::getSuper(const std::string& type,
		const std::vector<int>& sizes) const {
	// Local Declarations
	static std::map<std::string, int> composition = { { heType, 0 },
			{ vType, 0 }, { iType, 0 }, { xeType, 0 } };
	std::shared_ptr<IReactant> retReactant;

	// Setup the composition map to default values
	composition[heType] = 0;
	composition[vType] = 0;
	composition[iType] = 0;

	// Only pull the reactant if the name is valid and there are enough sizes
	// to fill the composition.
	if (type == PSISuperType && sizes.size() == 3) {
		composition[heType] = sizes[0];
		composition[vType] = sizes[1];
		composition[iType] = sizes[2];
		// Make sure the reactant is in the map
		std::string compStr = Reactant::toCanonicalString(type, composition);
		if (superSpeciesMap.count(compStr)) {
			retReactant = superSpeciesMap.at(compStr);
		}
	}

	return retReactant.get();
}


void PSIClusterReactionNetwork::add(std::shared_ptr<IReactant> reactant) {
	// Local Declarations
	int numHe = 0, numV = 0, numI = 0;
	bool isMixed = false;
	int* numClusters = nullptr;
	int* maxClusterSize = nullptr;

	// Only add a complete reactant
	if (reactant != NULL) {
		// Get the composition
		auto composition = reactant->getComposition();
		std::string compStr = reactant->getCompositionString();
		// Get the species sizes
		numHe = composition.at(heType);
		numV = composition.at(vType);
		numI = composition.at(iType);

		// Determine if the cluster is a compound. If there is more than one
		// type, then the check below will sum to greater than one and we know
		// that we have a mixed cluster.
		isMixed = ((numHe > 0) + (numV > 0) + (numI > 0)) > 1;
		// Only add the element if we don't already have it
		// Add the compound or regular reactant.
		if (isMixed && mixedSpeciesMap.count(compStr) == 0) {
			// Put the compound in its map
			mixedSpeciesMap[compStr] = reactant;
			// Figure out whether we have HeV or HeI and set the keys
			if (numV > 0) {
				numClusters = &numHeVClusters;
				maxClusterSize = &maxHeVClusterSize;
			} else {
				numClusters = &numHeIClusters;
				maxClusterSize = &maxHeIClusterSize;
			}
		} else if (!isMixed && singleSpeciesMap.count(compStr) == 0) {
			/// Put the reactant in its map
			singleSpeciesMap[compStr] = reactant;

			// Figure out whether we have He, V or I and set the keys
			if (numHe > 0) {
				numClusters = &numHeClusters;
				maxClusterSize = &maxHeClusterSize;
			} else if (numV > 0) {
				numClusters = &numVClusters;
				maxClusterSize = &maxVClusterSize;
			} else {
				numClusters = &numIClusters;
				maxClusterSize = &maxIClusterSize;
			}
		} else {
			std::stringstream errStream;
			errStream << "PSIClusterReactionNetwork Message: "
					<< "Duplicate Reactant (He=" << numHe << ",V=" << numV
					<< ",I=" << numI << ") not added!" << std::endl;
			throw errStream.str();
		}

		// Increment the number of total clusters of this type
		(*numClusters)++;
		// Increment the max cluster size key
		int clusterSize = numHe + numV + numI;
		(*maxClusterSize) = std::max(clusterSize, (*maxClusterSize));
		// Update the size
		++networkSize;
		// Set the id for this cluster
		reactant->setId(networkSize);

        // Add reactant to our per-type map.
        clusterTypeMap.at(reactant->getType()).push_back(reactant);

		// Add the pointer to the list of all clusters
        allReactants.push_back(reactant.get());
	}

	return;
}

void PSIClusterReactionNetwork::addSuper(std::shared_ptr<IReactant> reactant) {
	// Local Declarations
	int numHe = 0, numV = 0, numI = 0;
	bool isMixed = false;
	int* numClusters = nullptr;

	// Only add a complete reactant
	if (reactant != NULL) {
		// Get the composition
		auto composition = reactant->getComposition();
		std::string compStr = reactant->getCompositionString();
		// Get the species sizes
		numHe = composition.at(heType);
		numV = composition.at(vType);
		numI = composition.at(iType);
		// Determine if the cluster is a compound. If there is more than one
		// type, then the check below will sum to greater than one and we know
		// that we have a mixed cluster.
		isMixed = ((numHe > 0) + (numV > 0) + (numI > 0)) > 1;
		// Only add the element if we don't already have it
		// Add the compound or regular reactant.
		if (isMixed && superSpeciesMap.count(compStr) == 0) {
			// Put the compound in its map
			superSpeciesMap[compStr] = reactant;
			// Set the key
			numClusters = &numSuperClusters;
		} else {
			std::stringstream errStream;
			errStream << "PSIClusterReactionNetwork Message: "
					<< "Duplicate Super Reactant (He=" << numHe << ",V=" << numV
					<< ",I=" << numI << ") not added!" << std::endl;
			throw errStream.str();
		}

		// Increment the number of total clusters of this type
		(*numClusters)++;
		// Update the size
		++networkSize;
		// Set the id for this cluster
		reactant->setId(networkSize);
		// Add cluster to our per-type map.
        clusterTypeMap.at(reactant->getType()).push_back(reactant);

		// Add the pointer to the list of all clusters
		allReactants.push_back(reactant.get());
	}

	return;
}

void PSIClusterReactionNetwork::removeReactants(
		const std::vector<IReactant*>& doomedReactants) {

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
	std::set<std::string> typesUsed;
	for (auto reactant : doomedReactants) {
		typesUsed.insert(reactant->getType());
	}

	// ...Next, examine each type's collection of clusters and remove the
	// doomed reactants.
	for (auto currType : typesUsed) {
		auto& clusters = clusterTypeMap[currType];
		auto citer = std::remove_if(clusters.begin(), clusters.end(),
				doomedReactantMatcher);
		clusters.erase(citer, clusters.end());
	}

	// Remove the doomed reactants from the SpeciesMap.
	// We cannot use std::remove_if and our ReactantMatcher here
	// because std::remove_if reorders the elements in the underlying
	// container to move the doomed elements to the end of the container,
	// but the std::map doesn't support reordering.
	for (auto reactant : doomedReactants) {
		if (reactant->isMixed())
			mixedSpeciesMap.erase(reactant->getCompositionString());
		else
			singleSpeciesMap.erase(reactant->getCompositionString());
	}

	return;
}

void PSIClusterReactionNetwork::reinitializeNetwork() {
	// Recount HeV clusters
	numHeVClusters = 0;
	// Reset the Ids
	int id = 0;
	for (auto it = allReactants.begin(); it != allReactants.end(); ++it) {
		id++;
		(*it)->setId(id);
		(*it)->setHeMomentumId(id);
		(*it)->setVMomentumId(id);

		if ((*it)->getType() == heVType)
			numHeVClusters++;
	}

	// Reset the network size
	networkSize = id;

	// Get all the super clusters and loop on them
    for (auto& currCluster : clusterTypeMap[PSISuperType]) {

		id++;
		currCluster->setHeMomentumId(id);
		id++;
		currCluster->setVMomentumId(id);

		// Update the HeV size
		auto cluster = (PSISuperCluster *) currCluster.get();
		auto bounds = cluster->getBoundaries();
		int clusterSize = bounds[1] + bounds[3];
		maxHeVClusterSize = max(maxHeVClusterSize, clusterSize);
	}

	return;
}

void PSIClusterReactionNetwork::reinitializeConnectivities() {
	// Loop on all the reactants to reset their connectivities
	for (auto it = allReactants.begin(); it != allReactants.end(); ++it) {
		(*it)->resetConnectivities();
	}

	return;
}

void PSIClusterReactionNetwork::updateConcentrationsFromArray(
		double * concentrations) {
	// Local Declarations
	int size = allReactants.size();
	int id = 0;

	// Set the concentrations
	concUpdateCounter->increment();	// increment the update concentration counter
	for (int i = 0; i < size; i++) {
		id = allReactants.at(i)->getId() - 1;
		allReactants.at(i)->setConcentration(concentrations[id]);
	}

	// Set the moments
	for (int i = size - numSuperClusters; i < size; i++) {
		auto cluster = (PSISuperCluster *) allReactants.at(i);
		id = cluster->getId() - 1;
		cluster->setZerothMomentum(concentrations[id]);
		id = cluster->getHeMomentumId() - 1;
		cluster->setHeMomentum(concentrations[id]);
		id = cluster->getVMomentumId() - 1;
		cluster->setVMomentum(concentrations[id]);
	}

	return;
}

void PSIClusterReactionNetwork::getDiagonalFill(int *diagFill) {
	// Get all the super clusters
	auto superClusters = getAll(PSISuperType);

	// Degrees of freedom is the total number of clusters in the network
	const int dof = getDOF();

	// Declarations for the loop
	std::vector<int> connectivity;
	int connectivityLength, id, index;

	// Get the connectivity for each reactant
	for (int i = 0; i < networkSize; i++) {
		// Get the reactant and its connectivity
		auto& reactant = allReactants.at(i);
		connectivity = reactant->getConnectivity();
		connectivityLength = connectivity.size();
		// Get the reactant id so that the connectivity can be lined up in
		// the proper column
		id = reactant->getId() - 1;
		// Create the vector that will be inserted into the dFill map
		std::vector<int> columnIds;
		// Add it to the diagonal fill block
		for (int j = 0; j < connectivityLength; j++) {
			// The id starts at j*connectivity length and is always offset
			// by the id, which denotes the exact column.
			index = id * dof + j;
			diagFill[index] = connectivity[j];
			// Add a column id if the connectivity is equal to 1.
			if (connectivity[j] == 1) {
				columnIds.push_back(j);
			}
		}
		// Update the map
		dFillMap[id] = columnIds;
	}
	// Get the connectivity for each moment
	for (int i = 0; i < superClusters.size(); i++) {
		// Get the reactant and its connectivity
		auto reactant = superClusters[i];
		connectivity = reactant->getConnectivity();
		connectivityLength = connectivity.size();
		// Get the helium momentum id so that the connectivity can be lined up in
		// the proper column
		id = reactant->getHeMomentumId() - 1;

		// Create the vector that will be inserted into the dFill map
		std::vector<int> columnIds;
		// Add it to the diagonal fill block
		for (int j = 0; j < connectivityLength; j++) {
			// The id starts at j*connectivity length and is always offset
			// by the id, which denotes the exact column.
			index = (id) * dof + j;
			diagFill[index] = connectivity[j];
			// Add a column id if the connectivity is equal to 1.
			if (connectivity[j] == 1) {
				columnIds.push_back(j);
			}
		}
		// Update the map
		dFillMap[id] = columnIds;

		// Get the vacancy momentum id so that the connectivity can be lined up in
		// the proper column
		id = reactant->getVMomentumId() - 1;

		// Add it to the diagonal fill block
		for (int j = 0; j < connectivityLength; j++) {
			// The id starts at j*connectivity length and is always offset
			// by the id, which denotes the exact column.
			index = (id) * dof + j;
			diagFill[index] = connectivity[j];
		}
		// Update the map
		dFillMap[id] = columnIds;
	}

	return;
}

double PSIClusterReactionNetwork::getTotalAtomConcentration() {
	// Initial declarations
	double heliumConc = 0.0;

	// Get all the He clusters
	auto heClusters = getAll(heType);
	// Loop on them
	for (int i = 0; i < heClusters.size(); i++) {
		// Get the cluster and its composition
		auto cluster = heClusters[i];
		double size = cluster->getSize();

		// Add the concentration times the He content to the total helium concentration
		heliumConc += cluster->getConcentration() * size;
	}

	// Get all the HeV clusters
	auto heVClusters = getAll(heVType);
	// Loop on them
	for (int i = 0; i < heVClusters.size(); i++) {
		// Get the cluster and its composition
		auto cluster = heVClusters[i];
		auto comp = cluster->getComposition();

		// Add the concentration times the He content to the total helium concentration
		heliumConc += cluster->getConcentration() * comp[heType];
	}

	// Get all the super clusters
	auto superClusters = getAll(PSISuperType);
	// Loop on them
	for (int i = 0; i < superClusters.size(); i++) {
		// Get the cluster
		auto cluster = (PSISuperCluster *) superClusters[i];

		// Add its total helium concentration helium concentration
		heliumConc += cluster->getTotalHeliumConcentration();
	}

	return heliumConc;
}

double PSIClusterReactionNetwork::getTotalTrappedAtomConcentration() {
	// Initial declarations
	double heliumConc = 0.0;

	// Get all the HeV clusters
	auto heVClusters = getAll(heVType);
	// Loop on them
	for (int i = 0; i < heVClusters.size(); i++) {
		// Get the cluster and its composition
		auto cluster = heVClusters[i];
		auto comp = cluster->getComposition();

		// Add the concentration times the He content to the total helium concentration
		heliumConc += cluster->getConcentration() * comp[heType];
	}

	// Get all the super clusters
	auto superClusters = getAll(PSISuperType);
	// Loop on them
	for (int i = 0; i < superClusters.size(); i++) {
		// Get the cluster
		auto cluster = (PSISuperCluster *) superClusters[i];

		// Add its total helium concentration
		heliumConc += cluster->getTotalHeliumConcentration();
	}

	return heliumConc;
}

double PSIClusterReactionNetwork::getTotalVConcentration() {
	// Initial declarations
	double vConc = 0.0;

	// Get all the V clusters
	auto vClusters = getAll(vType);
	// Loop on them
	for (int i = 0; i < vClusters.size(); i++) {
		// Get the cluster and its composition
		auto cluster = vClusters[i];
		double size = cluster->getSize();

		// Add the concentration times the V content to the total vacancy concentration
		vConc += cluster->getConcentration() * size;
	}

	// Get all the HeV clusters
	auto heVClusters = getAll(heVType);
	// Loop on them
	for (int i = 0; i < heVClusters.size(); i++) {
		// Get the cluster and its composition
		auto cluster = heVClusters[i];
		auto comp = cluster->getComposition();

		// Add the concentration times the V content to the total vacancy concentration
		vConc += cluster->getConcentration() * comp[vType];
	}

	// Get all the super clusters
	auto superClusters = getAll(PSISuperType);
	// Loop on them
	for (int i = 0; i < superClusters.size(); i++) {
		// Get the cluster
		auto cluster = (PSISuperCluster *) superClusters[i];

		// Add its total vacancy concentration
		vConc += cluster->getTotalVacancyConcentration();
	}

	return vConc;
}

double PSIClusterReactionNetwork::getTotalIConcentration() {
	// Initial declarations
	double iConc = 0.0;

	// Get all the V clusters
	auto iClusters = getAll(iType);
	// Loop on them
	for (int i = 0; i < iClusters.size(); i++) {
		// Get the cluster and its composition
		auto cluster = iClusters[i];
		double size = cluster->getSize();

		// Add the concentration times the I content to the total interstitial concentration
		iConc += cluster->getConcentration() * size;
	}

	return iConc;
}

void PSIClusterReactionNetwork::computeRateConstants() {
	// Local declarations
	double rate = 0.0;
	// Initialize the value for the biggest production rate
	double biggestProductionRate = 0.0;

	// Loop on all the production reactions
	for (auto iter = allProductionReactions.begin();
			iter != allProductionReactions.end(); iter++) {
		// Compute the rate
		rate = calculateReactionRateConstant(iter->get());
		// Set it in the reaction
		(*iter)->kConstant = rate;

		// Check if the rate is the biggest one up to now
		if (rate > biggestProductionRate)
			biggestProductionRate = rate;
	}

	// Loop on all the dissociation reactions
	for (auto iter = allDissociationReactions.begin();
			iter != allDissociationReactions.end(); iter++) {
		// Compute the rate
		rate = calculateDissociationConstant(iter->get());
		// Set it in the reaction
		(*iter)->kConstant = rate;
	}

	// Set the biggest rate
	biggestRate = biggestProductionRate;

	return;
}

void PSIClusterReactionNetwork::computeAllFluxes(double *updatedConcOffset) {
	// Initial declarations
	PSISuperCluster * superCluster;
	double flux = 0.0;
	int reactantIndex = 0;
	auto superClusters = getAll(PSISuperType);

	// ----- Compute all of the new fluxes -----
	for (int i = 0; i < networkSize; i++) {
		auto& cluster = allReactants.at(i);
		// Compute the flux
		flux = cluster->getTotalFlux();
		// Update the concentration of the cluster
		reactantIndex = cluster->getId() - 1;
		updatedConcOffset[reactantIndex] += flux;
	}

	// ---- Moments ----
	for (int i = 0; i < superClusters.size(); i++) {
		superCluster = (xolotlCore::PSISuperCluster *) superClusters[i];

		// Compute the helium momentum flux
		flux = superCluster->getHeMomentumFlux();
		// Update the concentration of the cluster
		reactantIndex = superCluster->getHeMomentumId() - 1;
		updatedConcOffset[reactantIndex] += flux;

		// Compute the vacancy momentum flux
		flux = superCluster->getVMomentumFlux();
		// Update the concentration of the cluster
		reactantIndex = superCluster->getVMomentumId() - 1;
		updatedConcOffset[reactantIndex] += flux;
	}

	return;
}

void PSIClusterReactionNetwork::computeAllPartials(double *vals, int *indices,
		int *size) {
	// Initial declarations
	int reactantIndex = 0, pdColIdsVectorSize = 0;
	const int dof = getDOF();
	std::vector<double> clusterPartials;
	clusterPartials.resize(dof, 0.0);
	// Get the super clusters
	auto superClusters = getAll(PSISuperType);

	// Update the column in the Jacobian that represents each normal reactant
	for (int i = 0; i < networkSize - superClusters.size(); i++) {
		auto& reactant = allReactants.at(i);
		// Get the reactant index
		reactantIndex = reactant->getId() - 1;

		// Get the partial derivatives
		reactant->getPartialDerivatives(clusterPartials);
		// Get the list of column ids from the map
		auto pdColIdsVector = dFillMap.at(reactantIndex);
		// Number of partial derivatives
		pdColIdsVectorSize = pdColIdsVector.size();
		size[reactantIndex] = pdColIdsVectorSize;

		// Loop over the list of column ids
		for (int j = 0; j < pdColIdsVectorSize; j++) {
			// Set the index
			indices[reactantIndex * dof + j] = pdColIdsVector[j];

			// Get the partial derivative from the array of all of the partials
			vals[reactantIndex * dof + j] = clusterPartials[pdColIdsVector[j]];

			// Reset the cluster partial value to zero. This is much faster
			// than using memset.
			clusterPartials[pdColIdsVector[j]] = 0.0;
		}
	}

	// Update the column in the Jacobian that represents the moment for the super clusters
	for (int i = 0; i < superClusters.size(); i++) {
		auto reactant = (PSISuperCluster *) superClusters[i];

		// Get the super cluster index
		reactantIndex = reactant->getId() - 1;

		// Get the partial derivatives
		reactant->getPartialDerivatives(clusterPartials);
		// Get the list of column ids from the map
		auto pdColIdsVector = dFillMap.at(reactantIndex);
		// Number of partial derivatives
		pdColIdsVectorSize = pdColIdsVector.size();
		size[reactantIndex] = pdColIdsVectorSize;

		// Loop over the list of column ids
		for (int j = 0; j < pdColIdsVectorSize; j++) {
			// Set the index
			indices[reactantIndex * dof + j] = pdColIdsVector[j];
			// Get the partial derivative from the array of all of the partials
			vals[reactantIndex * dof + j] = clusterPartials[pdColIdsVector[j]];

			// Reset the cluster partial value to zero. This is much faster
			// than using memset.
			clusterPartials[pdColIdsVector[j]] = 0.0;
		}

		// Get the helium momentum index
		reactantIndex = reactant->getHeMomentumId() - 1;

		// Get the partial derivatives
		reactant->getHeMomentPartialDerivatives(clusterPartials);
		// Get the list of column ids from the map
		pdColIdsVector = dFillMap.at(reactantIndex);
		// Number of partial derivatives
		pdColIdsVectorSize = pdColIdsVector.size();
		size[reactantIndex] = pdColIdsVectorSize;

		// Loop over the list of column ids
		for (int j = 0; j < pdColIdsVectorSize; j++) {
			// Set the index
			indices[reactantIndex * dof + j] = pdColIdsVector[j];
			// Get the partial derivative from the array of all of the partials
			vals[reactantIndex * dof + j] = clusterPartials[pdColIdsVector[j]];

			// Reset the cluster partial value to zero. This is much faster
			// than using memset.
			clusterPartials[pdColIdsVector[j]] = 0.0;
		}

		// Get the vacancy momentum index
		reactantIndex = reactant->getVMomentumId() - 1;

		// Get the partial derivatives
		reactant->getVMomentPartialDerivatives(clusterPartials);
		// Get the list of column ids from the map
		pdColIdsVector = dFillMap.at(reactantIndex);
		// Number of partial derivatives
		pdColIdsVectorSize = pdColIdsVector.size();
		size[reactantIndex] = pdColIdsVectorSize;

		// Loop over the list of column ids
		for (int j = 0; j < pdColIdsVectorSize; j++) {
			// Set the index
			indices[reactantIndex * dof + j] = pdColIdsVector[j];
			// Get the partial derivative from the array of all of the partials
			vals[reactantIndex * dof + j] = clusterPartials[pdColIdsVector[j]];

			// Reset the cluster partial value to zero. This is much faster
			// than using memset.
			clusterPartials[pdColIdsVector[j]] = 0.0;
		}
	}

	return;
}

double PSIClusterReactionNetwork::computeBindingEnergy(
		DissociationReaction * reaction) const {

	double bindingEnergy = 5.0;
	if (reaction->dissociating->getType() == heType
			&& reaction->first->getType() == heType) {
		if (reaction->dissociating->getSize() == 2)
			bindingEnergy = 0.5;
		else
			bindingEnergy = 1.0;
	}
	if (reaction->dissociating->getType() == vType
			&& reaction->first->getType() == vType) {
		int size = reaction->dissociating->getSize();
		bindingEnergy = 1.73
				- 2.59
						* (pow((double) size, 2.0 / 3.0)
								- pow((double) size - 1.0, 2.0 / 3.0));
	}
	if ((reaction->dissociating->getType() == heVType)
			&& (reaction->first->getType() == vType
					|| reaction->second->getType() == vType)) {
		auto comp = reaction->dissociating->getComposition();
		bindingEnergy = 1.73
				- 2.59
						* (pow((double) comp[vType], 2.0 / 3.0)
								- pow((double) comp[vType] - 1.0, 2.0 / 3.0))
				+ 2.5
						* log(
								1.0
										+ ((double) comp[heType]
												/ (double) comp[vType]));
	}
	if (reaction->dissociating->getType() == PSISuperType
			&& (reaction->first->getType() == vType
					|| reaction->second->getType() == vType)) {
		auto superCluster = (PSISuperCluster *) reaction->dissociating;
		auto comp = reaction->dissociating->getComposition();
		double numV = (double) comp[vType];
		double numHe = (double) comp[heType];
		bindingEnergy = 1.73
				- 2.59 * (pow(numV, 2.0 / 3.0) - pow(numV - 1.0, 2.0 / 3.0))
				+ 2.5 * log(1.0 + (numHe / numV));
	}
	if (reaction->first->getType() == iType
			|| reaction->second->getType() == iType) {
		if (reaction->dissociating->getType() == heVType) {
			auto comp = reaction->dissociating->getComposition();
			bindingEnergy =
					4.88
							+ 2.59
									* (pow((double) comp[vType], 2.0 / 3.0)
											- pow((double) comp[vType] - 1.0,
													2.0 / 3.0))
							- 2.5
									* log(
											1.0
													+ ((double) comp[heType]
															/ (double) comp[vType]));
		} else if (reaction->dissociating->getType() == PSISuperType) {
			auto superCluster = (PSISuperCluster *) reaction->dissociating;
			auto comp = reaction->dissociating->getComposition();
			double numV = (double) comp[vType];
			double numHe = (double) comp[heType];
			bindingEnergy = 4.88
					+ 2.59 * (pow(numV, 2.0 / 3.0) - pow(numV - 1.0, 2.0 / 3.0))
					- 2.5 * log(1.0 + (numHe / numV));
		} else if (reaction->dissociating->getType() == heType) {
			int size = reaction->dissociating->getSize();
			switch (size) {
			case 1:
				bindingEnergy = 4.31;
				break;
			case 2:
				bindingEnergy = 2.90;
				break;
			case 3:
				bindingEnergy = 2.02;
				break;
			case 4:
				bindingEnergy = 1.09;
				break;
			case 5:
				bindingEnergy = 0.58;
				break;
			case 6:
				bindingEnergy = 0.13;
				break;
			case 7:
				bindingEnergy = -0.25;
				break;
			case 8:
				bindingEnergy = -0.59;
				break;
			default:
				break;
			}
		}

	}

//	if (bindingEnergy < -5.0)
//	std::cout << "dissociation: " << reaction->dissociating->getName() << " -> "
//			<< reaction->first->getName() << " + "
//			<< reaction->second->getName() << " : " << bindingEnergy
//			<< std::endl;

	return max(bindingEnergy, -5.0);
}


IReactant * PSIClusterReactionNetwork::getSuperFromComp(int nHe, int nV) const {
	// Initial declarations
	IReactant * toReturn = nullptr;
	static std::vector<IReactant *> superClusters;
	int superSize = superClusters.size();
	if (superSize == 0) superClusters = getAll(PSISuperType);
	// Find the right indices for He and V
	int i = -1, j = -1;
	for (auto it = boundVector.begin(); it != boundVector.end(); it++) {
		if (nHe >= *it) i++;
		if (nV >= *it) j++;
	}

	// Compute the super index
	int index = i + (j * (boundVector.size() - 1)) - std::min(j+1, 3) * 3;
	// Get the super cluster
	if (index < superSize)
		toReturn = superClusters[index];

	return toReturn;
}
