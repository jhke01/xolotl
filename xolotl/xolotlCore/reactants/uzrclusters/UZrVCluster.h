#ifndef UZRVCLUSTER_H
#define UZRVCLUSTER_H

// Includes
#include "UZrCluster.h"
#include <xolotlPerf.h>
#include <Constants.h>

namespace xolotlCore {

/**
 * This class represents a cluster composed entirely of atomic vacancies.
 */
class UZrVCluster: public UZrCluster {

private:
	static std::string buildName(IReactant::SizeType nV) {
		std::stringstream nameStream;
		nameStream << "V_" << nV;
		return nameStream.str();
	}

	/**
	 * The sink strength
	 */
	double sinkStrength;

	/**
	 * The recombination term
	 */
	double recombStrength;
	double trapStrength;

public:

	/**
	 * Default constructor, deleted because we require info to construct.
	 */
	UZrVCluster() = delete;

	/**
	 * The constructor. All UZrVClusters must be initialized with a size.
	 *
	 * @param nV the number of atomic vacancies in the cluster
	 * @param registry The performance handler registry
	 */
	UZrVCluster(int nV, IReactionNetwork &_network,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
			UZrCluster(_network, registry, buildName(nV)) {

		// Set the size
		size = nV;
		// Set the typename appropriately
		type = ReactantType::V;

		// Update the composition map
		composition[toCompIdx(Species::V)] = size;

		// TODO: Compute the reaction radius
		// You can look at FeVCluster or PSIVCluster
		//reactionRadius = network.getLatticeParameter();

		// Compute the reaction radius
		double latticeParam = network.getLatticeParameter();
		// It is the same formula for HeV clusters
		reactionRadius = latticeParam
				* pow((3.0 * size) / xolotlCore::pi, (1.0 / 3.0)) * 0.5;

		//! Parameters for biased sink in the iron case
		double r0 = latticeParam * 0.75 * sqrt(2.0);
		double reactionRadius = latticeParam * cbrt(3.0 / xolotlCore::pi) * 0.5;
		constexpr double rho = 1e14 * 1e-18; //sink density? in a unit of nm^-2?

		sinkStrength = rho * 1.00;
		/*
		sinkStrength = -4.0 * xolotlCore::pi * rho
				/ log(
						xolotlCore::pi * rho * (reactionRadius + r0)
								* (reactionRadius + r0));
		*/

		recombStrength = 4.0 * xolotlCore::pi * reactionRadius * 1.00;



		trapStrength = 4.0 * xolotlCore::pi * reactionRadius * 0.00;


		return;
	}

	/**
	 * Copy constructor, deleted to prevent use.
	 */
	UZrVCluster(const UZrVCluster &other) = delete;

	//! Destructor
	~UZrVCluster() {
	}

	/**
	 * This operation returns the total change in this cluster due its
	 * own dissociation.
	 *
	 * @param i The location on the grid in the depth direction
	 * @return The flux due to its dissociation
	 */
	double getEmissionFlux(int i) override {
		// Initial declarations
		double flux = UZrCluster::getEmissionFlux(i);
		double k_b = xolotlCore::kBoltzmann;
		double temp = network.getTemperature();
		//double D_v = 1e-4*exp(-2.0/k_b/temp)*1e18; //D in nm^2/s
		double Di = 1e-8*exp(-0.6/k_b/temp)*1e18; //D in nm^2/s

		auto singleXeCluster = network.get(Species::Xe, 1);
		auto conc = singleXeCluster->getConcentration();
		//std::cout << "conc " << conc << std::endl;

		// Compute the loss to dislocation sinks
		if (size < 2) {
			// k^2 * D * C
			flux += sinkStrength * diffusionCoefficient[i] * concentration
					+ recombStrength * (Di + diffusionCoefficient[i])
						* diffusionCoefficient[i] / Di * concentration * concentration
					+ trapStrength * diffusionCoefficient[i] * concentration * conc;
		}

		return flux;
	}

	/**
	 * This operation computes the partial derivatives due to emission
	 * reactions.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 * @param i The location on the grid in the depth direction
	 */
	void getEmissionPartialDerivatives(std::vector<double> & partials,
			int i) const override {
		// Initial declarations
		UZrCluster::getEmissionPartialDerivatives(partials, i);
		double k_b = xolotlCore::kBoltzmann;
		double temp = network.getTemperature();
		//double D_v = 1e-4*exp(-2.0/k_b/temp)*1e18; //D in nm^2/s
		double Di = 1e-8*exp(-0.6/k_b/temp)*1e18; //D in nm^2/s

		auto singleXeCluster = network.get(Species::Xe, 1);
		auto conc = singleXeCluster->getConcentration();

		// Compute the loss to dislocation sinks
		if (size < 2) {
			// k^2 * D * C
			partials[id - 1] -= sinkStrength * diffusionCoefficient[i]
			+	recombStrength * (Di + diffusionCoefficient[i])
				* diffusionCoefficient[i] / Di * concentration * 2
			+ trapStrength * diffusionCoefficient[i] * conc;
		}

		return;
	}

	/**
	 * Add grid points to the vector of diffusion coefficients or remove
	 * them if the value is negative.
	 *
	 * @param i The number of grid point to add or remove
	 */
	void addGridPoints(int i) override {
		if (diffusionFactor > 0.0) {
			Reactant::addGridPoints(i);
		}

		// Don't do anything
		return;
	}

	/**
	 * This operation sets the temperature at which the reactant currently
	 * exists. Temperature-dependent quantities are recomputed when this
	 * operation is called, so the temperature should always be set first.
	 *
	 * @param temp The new cluster temperature
	 * @param i The location on the grid
	 */
	void setTemperature(double temp, int i) override {
		if (diffusionFactor > 0.0) {
			Reactant::setTemperature(temp, i);
		}

		// Don't do anything
		return;
	}

	/**
	 * This operation returns the diffusion coefficient for this reactant and is
	 * calculated from the diffusion factor.
	 *
	 * @param i The position on the grid
	 * @return The diffusion coefficient
	 */
	double getDiffusionCoefficient(int i) const override {
		if (diffusionFactor > 0.0) {
			return Reactant::getDiffusionCoefficient(i);
		}

		return 0.0;
	}

};
//end class UZrVCluster

} /* end namespace xolotlCore */

#endif
