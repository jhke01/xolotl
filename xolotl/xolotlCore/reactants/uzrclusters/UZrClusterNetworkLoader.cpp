#include <fstream>
#include <functional>
#include <cassert>
#include "UZrClusterNetworkLoader.h"
#include <UZrClusterReactionNetwork.h>
#include <UZrXeCluster.h>
#include <UZrVCluster.h>
#include <UZrXeVCluster.h>
#include <xolotlPerf.h>
#include <MathUtils.h>
#include <UZrSuperCluster.h>
#include "xolotlCore/io/XFile.h"


namespace xolotlCore {

std::unique_ptr<UZrCluster> UZrClusterNetworkLoader::createUZrCluster(int numXe,
		int numV, IReactionNetwork &network) const {

	// Local Declarations
	UZrCluster *cluster = nullptr;

	// Determine the type of the cluster given the number of each species.
	// Create a new cluster by that type and specify the names of the
	// property keys.
	/*
	if (numXe > 0 && numV > 0) {
		// Create a new XeVCluster
		cluster = new UZrXeVCluster(numXe, numV, network, handlerRegistry);
	} else if (numXe > 0) {
		// Create a new XeCluster
		cluster = new UZrXeCluster(numXe, network, handlerRegistry);
	} else if (numV > 0) {
		// Create a new VCluster
		cluster = new UZrVCluster(numV, network, handlerRegistry);
	}
	*/

  if (numXe > 0) {
		// Create a new XeCluster
		cluster = new UZrXeCluster(numXe, network, handlerRegistry);
	} else if (numV > 0) {
		// Create a new VCluster
		cluster = new UZrVCluster(numV, network, handlerRegistry);
	}


	assert(cluster != nullptr);

	return std::unique_ptr<UZrCluster>(cluster);
}

std::unique_ptr<UZrCluster> UZrClusterNetworkLoader::createUZrSuperCluster(
		int nTot, int maxXe, IReactionNetwork& network) const {
	// Create the cluster
	auto superCluster = new UZrSuperCluster(maxXe, nTot, network,
			handlerRegistry);

	// TODO when we have widespread C++14 support, use std::make_unique
	// and construct unique ptr and object pointed to in one memory operation.
	return std::unique_ptr<UZrCluster>(superCluster);
}

void UZrClusterNetworkLoader::pushUZrCluster(
		std::unique_ptr<UZrClusterReactionNetwork> & network,
		std::vector<std::reference_wrapper<Reactant> > & reactants,
		std::unique_ptr<UZrCluster> & cluster) {
	// Check if we want dummy reactions
	if (dummyReactions) {
		// Create a dummy cluster (Reactant) from the existing cluster
		auto dummyCluster = std::unique_ptr<Reactant>(new Reactant(*cluster));
		// Save access to it so we can trigger updates after
		// we add all to the network.
		reactants.emplace_back(*dummyCluster);

		// Give the cluster to the network
		network->add(std::move(dummyCluster));
	} else {
		// Save access to it so we can trigger updates after
		// we add all to the network.
		reactants.emplace_back(*cluster);

		// Give the cluster to the network
		network->add(std::move(cluster));
	}

	return;
}


UZrClusterNetworkLoader::UZrClusterNetworkLoader(
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) {
	networkStream = nullptr;
	handlerRegistry = registry;
	fileName = "";
	dummyReactions = false;
	//maxXe = -1;
	//maxV = -1;
	//xeMin = 1000000;
	//xeMax = -1;
	sectionWidth = 1;

	return;
}

UZrClusterNetworkLoader::UZrClusterNetworkLoader(
		const std::shared_ptr<std::istream> stream,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) {
	networkStream = stream;
	handlerRegistry = registry;
	fileName = "";
	dummyReactions = false;
	//maxXe = -1;
	//maxV = -1;
	//xeMin = 1000000;
	//xeMax = -1;
	sectionWidth = 1;

	return;
}

std::unique_ptr<IReactionNetwork> UZrClusterNetworkLoader::load(
		const IOptions &options) {
	// Get the dataset from the HDF5 files
	int normalSize = 0, superSize = 0;
	XFile networkFile(fileName);
	auto networkGroup = networkFile.getGroup<XFile::NetworkGroup>();
	assert(networkGroup);
	networkGroup->readNetworkSize(normalSize, superSize);

	// Initialization
	int numXe = 0, numV = 0;
	double formationEnergy = 0.0, migrationEnergy = 0.0;
	double diffusionFactor = 0.0;
	std::vector<std::reference_wrapper<Reactant> > reactants;

	// Prepare the network
	std::unique_ptr<UZrClusterReactionNetwork> network(
			new UZrClusterReactionNetwork(handlerRegistry));

	// Set the lattice parameter in the network
	double latticeParam = options.getLatticeParameter();
	if (!(latticeParam > 0.0))
		latticeParam = uraniumZirconiumLatticeConstant;
	network->setLatticeParameter(latticeParam);

	// Set the xenon radius in the network
	double radius = options.getImpurityRadius();
	if (!(radius > 0.0))
		radius = xenonRadius;
	network->setImpurityRadius(radius);

	// Loop on the clusters
	for (int i = 0; i < normalSize + superSize; i++) {
		// Open the cluster group
		XFile::ClusterGroup clusterGroup(*networkGroup, i);

		if (i < normalSize) {
			// Normal cluster
			// Read the composition
			auto comp = clusterGroup.readCluster(formationEnergy,
					migrationEnergy, diffusionFactor);
			numXe = comp[toCompIdx(Species::Xe)];
			numV = comp[toCompIdx(Species::V)];

			// Create the cluster
			auto nextCluster = createUZrCluster(numXe, numV, *network);

			// Set the formation energy
			nextCluster->setFormationEnergy(formationEnergy);
			// Set the diffusion factor and migration energy
			nextCluster->setMigrationEnergy(migrationEnergy);
			nextCluster->setDiffusionFactor(diffusionFactor);

			if (numXe == 1) {
				// If the diffusivity is given
				if (options.getXenonDiffusivity() > 0.0) {
					nextCluster->setDiffusionFactor(
							options.getXenonDiffusivity());
					nextCluster->setMigrationEnergy(-1.0);
				}
			}

			// Save access to it so we can trigger updates once
			// added to the network.
			//reactants.emplace_back(*nextCluster);

			// Give the cluster to the network
			//network->add(std::move(nextCluster));

			// Save it in the network
			pushUZrCluster(network, reactants, nextCluster);
		} else {
			// This is not happening yet
			// Super cluster
			int nTot = 0, maxXe = 0;
			clusterGroup.readUZrSuperCluster(nTot, maxXe);

			// Create the cluster
			auto nextCluster = createUZrSuperCluster(nTot, maxXe, *network);

			// Save it in the network
			pushUZrCluster(network, reactants, nextCluster);
		}
	}

	// Ask reactants to update now that they are in network.
	for (IReactant &currReactant : reactants) {
		currReactant.updateFromNetwork();
	}

	// Set the reactions
	networkGroup->readReactions(*network);

	// Recompute Ids and network size
	network->reinitializeNetwork();

	return std::move(network);
}

std::unique_ptr<IReactionNetwork> UZrClusterNetworkLoader::generate(
		const IOptions &options) {
	// Initial declarations
	maxXe = options.getMaxImpurity(), maxV = options.getMaxV();
	xeMax = options.getMaxImpurity();
	int numXe = 0, numV = 0;
	double formationEnergy = 0.0, migrationEnergy = 0.0;
	double diffusionFactor = 0.0;

	// Once we have C++14, use std::make_unique.
	std::unique_ptr<UZrClusterReactionNetwork> network(
			new UZrClusterReactionNetwork(handlerRegistry));
	std::vector<std::reference_wrapper<Reactant> > reactants;

	// Set the lattice parameter in the network
	double latticeParam = options.getLatticeParameter();
	if (!(latticeParam > 0.0))
		latticeParam = uraniumZirconiumLatticeConstant;
	network->setLatticeParameter(latticeParam);

	// Set the xenon radius in the network
	double radius = options.getImpurityRadius();
	if (!(radius > 0.0))
		radius = xenonRadius;
	network->setImpurityRadius(radius);

	// Set the density in a bubble
	network->setDensity(options.getDensity());

	// TODO: replace with the correct kinetics here and phase-space

	// Xe formation energies in eV
	std::vector<double> heFormationEnergies = { 0.0 };
	// Xe diffusion factors in nm^2/s
	std::vector<double> heDiffusion = { 1 };
	// Xe migration energies in eV
	std::vector<double> heMigration = { 0 };

	// V formation energies in eV
	std::vector<double> vFormationEnergies = { 0.0 };
	// V diffusion factors in nm^2/s
	std::vector<double> vDiffusion = { 1e14 };
	// V migration energies in eV
	std::vector<double> vMigration = { options.getVaMigrationE() };

	// Generate the Xe clusters
	for (int i = 1; i <= min(xeMin - 1,xeMax); ++i) {
		// Set the composition
		numXe = i;
		// Create the cluster
		auto nextCluster = createUZrCluster(numXe, numV, *network);

		// Set the other attributes
		nextCluster->setFormationEnergy(0.0);

		if (i <= 1){
			nextCluster->setFormationEnergy((-1*log(options.getXeSolubility()) * (xolotlCore::kBoltzmann * options.getConstTemperature()))); //1e-8 is the Xe solubility
		} else {
			//nextCluster->setFormationEnergy(pow(i,2.0/3.0)*0.6434*3*1.0); // 0.6434 is for 0.1 J/m^2 interface energy
			nextCluster->setFormationEnergy(pow(i,2.0/3.0)*6.4349535575*options.getXeInterfaceE());
			// pow(36*xolotlCore::pi/pow(options.getDensity(),2),1.0/3.0) = 6.4349535575
		}

		if (i <= heDiffusion.size()) {
			nextCluster->setDiffusionFactor(heDiffusion[i - 1]);
			nextCluster->setMigrationEnergy(heMigration[i - 1]);
			// If the diffusivity is given
			if (options.getXenonDiffusivity() > 0.0) {
				nextCluster->setDiffusionFactor(options.getXenonDiffusivity());
				nextCluster->setMigrationEnergy(-1.0);
			}
		} else {
			nextCluster->setDiffusionFactor(0.0);
			nextCluster->setMigrationEnergy(
					std::numeric_limits<double>::infinity());
		}

		// Save access to it so we can trigger updates once
		// added to the network.
		//reactants.emplace_back(*nextCluster);

		// Give the cluster to the network
		//network->add(std::move(nextCluster));

		// Save it in the network
		pushUZrCluster(network, reactants, nextCluster);
	}

	// Reset the Xe composition
	numXe = 0;
	double atomicVolume = 0.5 * pow(latticeParam, 3);

	//std::cout << "lattice para =  " << latticeParam << std::endl;

	// Generate the V clusters
	for (int i = 1; i <= maxV; ++i) {
		// Set the composition
		numV = i;
		// Create the cluster
		auto nextCluster = createUZrCluster(numXe, numV, *network);

		// Set the other attributes
		nextCluster->setFormationEnergy(0.0);

		if (i <= 1){
			nextCluster->setFormationEnergy(options.getVaFormationE()); //1.2 is vacancy formation energy
		} else {
			nextCluster->setFormationEnergy(pow(i,2.0/3.0)*pow(36*xolotlCore::pi
				*pow(atomicVolume*1e-27,2),1.0/3.0)/1.602e-19 * options.getVaInterfaceE());
				// pow(36*xolotlCore::pi*pow(atomicVolume*1e-27,2),1.0/3.0) = 6.4349535575
		}

		if (i <= vDiffusion.size()) {
			nextCluster->setDiffusionFactor(vDiffusion[i - 1]);
			nextCluster->setMigrationEnergy(vMigration[i - 1]);
			// If the diffusivity is given
			if (options.getXenonDiffusivity() > 0.0) {
				nextCluster->setDiffusionFactor(options.getXenonDiffusivity());
				nextCluster->setMigrationEnergy(-1.0);
			}
		} else {
			nextCluster->setDiffusionFactor(0.0);
			nextCluster->setMigrationEnergy(
					std::numeric_limits<double>::infinity());
		}

		// Save access to it so we can trigger updates once
		// added to the network.
		reactants.emplace_back(*nextCluster);

		// Give the cluster to the network
		network->add(std::move(nextCluster));
	}

	// Reset the V composition
	numV = 0;

/* Turn off mixed cluster for testing
	// Loop over vacancies in the outer loop.
	for (int i = 1; i <= maxV; ++i) {
		numV = i;

		// Loop on the xenon number
		for (int j = 1; j <= maxXe; j++) {
			numXe = j;

			// Create the cluster
			auto nextCluster = createUZrCluster(numXe, numV, *network);
			// Set its attributes
			nextCluster->setFormationEnergy(0.0);
			nextCluster->setDiffusionFactor(0.0);
			nextCluster->setMigrationEnergy(
					std::numeric_limits<double>::infinity());

			// Save access to it so we can trigger updates once
			// all are added to the network.
			reactants.emplace_back(*nextCluster);

			// Give the cluster to the network
			network->add(std::move(nextCluster));
		}
	}
	*/

	// Update reactants now that they are in network.
	for (Reactant &currCluster : reactants) {
		currCluster.updateFromNetwork();
	}

	// Check if we want dummy reactions
	if (!dummyReactions) {
		// Apply sectional grouping
		applyGrouping(*network);
	}

	// Create the reactions
	network->createReactionConnectivity();

	// Recompute Ids and network size and redefine the connectivities
	//network->reinitializeNetwork();

	return std::move(network);
}

void UZrClusterNetworkLoader::applyGrouping(IReactionNetwork& network) const {

	// Initialize variables for the loop
	int count = 0, superCount = 0, width = sectionWidth;
	int size = 0;

	int xeMax = maxXe;
	// Decide here which types will undergo grouping
	//std::vector<ReactantType> typeVec { ReactantType::Xe };

	std::vector<ReactantType> typeVec { ReactantType::Xe };
	// Loop on the xenon groups
	for (int k = xeMin; k < xeMax; k++) {

		// Increment the counter
		count++;

		// Track the size
		size = k;

		// Continue if we are not at the wanted width yet
		if (count < width && k < xeMax - 1)
			continue;

		// Create the cluster
		auto rawSuperCluster = new UZrSuperCluster(size, count, network,
				handlerRegistry);

//		std::cout << superCount << " " << count << " "
//				<< rawSuperCluster->getName() << std::endl;

		auto superCluster = std::unique_ptr<UZrSuperCluster>(rawSuperCluster);
		// Give the cluster to the network.
		network.add(std::move(superCluster));

		// Reinitialize everything
		size = 0;
		count = 0;
		superCount++;
		width += 1;
	}

	std::cout << "FFF = " << xeMin << " " << xeMax  << std::endl;

	if (xeMin < xeMax) {
		// Group the last one alone
		auto rawSuperCluster = new UZrSuperCluster(xeMax, 1, network,
				handlerRegistry);

//		std::cout << superCount << " last " << rawSuperCluster->getName()
//				<< std::endl;

		auto superCluster = std::unique_ptr<UZrSuperCluster>(rawSuperCluster);
		// Give the cluster to the network.
		network.add(std::move(superCluster));
	}

	// Recompute Ids and network size
	network.reinitializeNetwork();

	return;
}

} // namespace xolotlCore
