// Includes
#include "Reactant.h"
#include "ReactionNetwork.h"
#include <xolotlPerf.h>
#include <iostream>
#include <math.h>
#include <MathUtils.h>
#include <Constants.h>

namespace xolotlCore {

// TODO modify signature to take type as argument.
Reactant::Reactant(IReactionNetwork& _network,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry,
		const std::string& _name) :
		concentration(0.0), id(0), type(ReactantType::Invalid), network(
				_network), handlerRegistry(registry), size(0), formationEnergy(
				0.0), diffusionFactor(0.0), migrationEnergy(0.0), name(_name), reactionRadius(
				0.0) {

}

void Reactant::recomputeDiffusionCoefficient(double temp, int i) {
	// Return zero if the diffusion factor is zero.
	/*
	if (!xolotlCore::equal(diffusionFactor, 0.0)) {
		// Otherwise use the Arrhenius equation to compute the diffusion
		// coefficient
		double k_b = xolotlCore::kBoltzmann;
		double kernel = -1.0 * migrationEnergy / (k_b * temp);
		diffusionCoefficient[i] = diffusionFactor * exp(kernel);
	}
	*/

	/*
	for (auto const &currMapItem : network.getAll(ReactantType::V)) {
		// Get the cluster
    auto const &cluster = *(currMapItem.second);
		//double V_conc = cluster.getConcentration();
		if (cluster.getSize() == 1) {
			monomerConc1 = cluster.getConcentration();
			}
		}
	//std::cout << "monomerXv =  " << monomerConc << std::endl;
	*/
	double atomicVolume = 0.5 * pow(network.getLatticeParameter(), 3);
	//std::cout << "XxXXX =  " << atomicVolume << std::endl;

	auto singleXeCluster = network.get(Species::V, 1);
	auto monomerConc = singleXeCluster->getConcentration();
	//std::cout << "conc " << monomerConc  << std::endl;

	if (getType() == ReactantType::Xe) {
		// Intrinsic diffusion
		//double kernel_0 = -1.30 / (xolotlCore::kBoltzmann * temp);
		//double D0 = 7.6e8 * exp(kernel_0); // nm2/s

		double kernel = -0.5 / (xolotlCore::kBoltzmann * temp); // just migration barrier
		double D3_Xe = 1e-7 * exp(kernel) * monomerConc * atomicVolume * 1.0e18; // nm2/s

		// Athermal diffusion
		// We need the fission rate now
		double fissionRate = network.getFissionRate() * 1.0e27; // #/m3/s
		double D1_Xe = (8e-40 * fissionRate) * 1.0e18; // nm2/s

		diffusionCoefficient[i] = D3_Xe ;
		//std::cout << "DXe_athermal =  " << D1_Xe << std::endl;
		std::cout << "DXe_RED =  " << diffusionCoefficient[i] << std::endl;
	}

	if (getType() == ReactantType::UZrSuper) {
		diffusionCoefficient[i] = 0 ;
	}

	if (getType() == ReactantType::V) {
		// Vacancy diffusion
		// parameters from Smirnova, D.E., Kuksin, A.Y., Starikov, S.V. et al.
		// Atomistic modeling of the self-diffusion in γ-U and γ-U-Mo.
		// Phys. Metals Metallogr. 116, 445–455 (2015).
		// https://doi.org/10.1134/S0031918X1503014X
		//double kernel = -0.5 / (xolotlCore::kBoltzmann * temp); // just migration barrier
		//double D3 = 1.74e-7 * exp(kernel) * 1.0e18; // nm2/s
		double kernel = -1*migrationEnergy / (xolotlCore::kBoltzmann * temp); // just migration barrier
		double D3_Va = 1e-8 * exp(kernel) * 1.0e18; // nm2/s
		diffusionCoefficient[i] = D3_Va;
		//std::cout << "migrationEnergy =  " << migrationEnergy << std::endl;
		//std::cout << "DV =  " << diffusionCoefficient[i] << std::endl;
	}


	return;
}

void Reactant::addGridPoints(int i) {
	// Add grid points
	if (i > 0) {
		while (i > 0) {
			diffusionCoefficient.emplace(diffusionCoefficient.begin(), 0.0);
			temperature.emplace(temperature.begin(), 0.0);

			// Decrease i
			i--;
		}
	} else {
		diffusionCoefficient.erase(diffusionCoefficient.begin(),
				diffusionCoefficient.begin() - i);
		temperature.erase(temperature.begin(), temperature.begin() - i);
	}
	return;
}

std::vector<int> Reactant::getConnectivity() const {
	// The connectivity array by default is filled with
	// zeros.
	int connectivityLength = network.getDOF();
	std::vector<int> connectivity = std::vector<int>(connectivityLength, 0);

	// This reactant should be connected to itself
	connectivity[id - 1] = 1;

	return connectivity;
}

void Reactant::setTemperature(double temp, int i) {
	temperature[i] = temp;

	// Recompute the diffusion coefficient
	recomputeDiffusionCoefficient(temp, i);
}

void Reactant::setDiffusionFactor(const double factor) {
	// Set the diffusion factor
	diffusionFactor = factor;

	return;
}

void Reactant::setMigrationEnergy(const double energy) {
	// Set the migration energy
	migrationEnergy = energy;

	return;
}

std::ostream&
operator<<(std::ostream& os, const IReactant::Composition& comp) {
	std::vector<Species> compSpecies { Species::He, Species::D, Species::T,
			Species::I, Species::V, Species::Xe };
	for (auto const& currSpecies : compSpecies) {
		os << toString(currSpecies) << comp[toCompIdx(currSpecies)];
	}
	return os;
}

std::ostream&
operator<<(std::ostream& os, const IReactant& reactant) {
	os << "id: " << reactant.getId() << "; " << "type: "
			<< toString(reactant.getType()) << "; " << "comp: "
			<< reactant.getComposition();
	return os;
}

} // namespace xolotlCore
