#ifndef PULSEDFITFLUXHANDLER_H
#define PULSEDFITFLUXHANDLER_H

#include "FluxHandler.h"
#include <cmath>
#include <MathUtils.h>
#include <experimental/PSIReactionNetwork.h>

namespace xolotlCore {

/**
 * This class realizes the IFluxHandler interface to calculate the incident V and I flux
 * for tungsten material with a flux amplitude being intermittent.
 */
class PulsedFitFluxHandler: public FluxHandler {
private:

	/**
	 * Total time length of the pulse
	 */
	double deltaTime;

	/**
	 * Proportion of the total time where the flux amplitude is not 0.0
	 */
	double alpha;

	/**
	 * Parameters for the gaussian profile
	 */
	double mu = 2000.0;
	double sigma = 100.0;

	/**
	 * Function that calculate the flux at a given position x (in nm).
	 * This function is not normalized.
	 *
	 * @param x The position where to evaluate the fit
	 * @return The evaluated value
	 */
	double FitFunction(double x) {
		// Compute the polynomial fit
		double value = exp(-pow((x - mu) / (sqrt(2.0) * sigma), 2.0));

		return std::max(value, 0.0);
	}

public:

	/**
	 * The constructor
	 */
	PulsedFitFluxHandler() :
			deltaTime(0.0), alpha(0.0) {
	}

	/**
	 * The Destructor
	 */
	~PulsedFitFluxHandler() {
	}

	/**
	 * Compute and store the incident flux values at each grid point.
	 * \see IFluxHandler.h
	 */
	void initializeFluxHandler(experimental::IReactionNetwork& network,
			int surfacePos, std::vector<double> grid) {
		// Call the general method
		FluxHandler::initializeFluxHandler(network, surfacePos, grid);

		// Skip if the flux amplitude is 0.0 and we are not using a time profile
		if (equal(fluxAmplitude, 0.0) && !useTimeProfile)
			return;

		// Set the flux index corresponding the the single vacancy cluster here
		using NetworkType =
		experimental::PSIReactionNetwork<experimental::PSIFullSpeciesList>;
		auto psiNetwork = dynamic_cast<NetworkType*>(&network);

		// Set the flux index corresponding the the single helium cluster here
		NetworkType::Composition comp;
		// Initialize the composition
		for (auto i : psiNetwork->getSpeciesRange()) {
			comp[i] = 0;
		}
		comp[NetworkType::Species::V] = 1;
		auto cluster = psiNetwork->findCluster(comp, plsm::onHost);
		// Check that the vacancy cluster is present in the network
//		if (!cluster) {
//			throw std::string(
//					"\nThe single vacancy cluster is not present in the network, "
//							"cannot use the flux option!");
//		}
		fluxIndices.push_back(cluster.getId());

		// Set the flux index corresponding the the single interstitial cluster here
		// Initialize the composition
		for (auto i : psiNetwork->getSpeciesRange()) {
			comp[i] = 0;
		}
		comp[NetworkType::Species::I] = 1;
		cluster = psiNetwork->findCluster(comp, plsm::onHost);
		// Check that the interstitial cluster is present in the network
//		if (!cluster) {
//			throw std::string(
//					"\nThe single interstitial cluster is not present in the network, "
//							"cannot use the flux option!");
//		}
		fluxIndices.push_back(cluster.getId());

		return;
	}

	/**
	 * This operation computes the flux due to incoming particles at a given grid point.
	 * \see IFluxHandler.h
	 */
	void computeIncidentFlux(double currentTime, double *updatedConcOffset,
			int xi, int surfacePos) {
		// Check in which phase of the pulse we are
		int cycle = currentTime / deltaTime;
		// The flux is 0.0 after alpha * deltaTime
		if (currentTime - ((double) cycle * deltaTime) > alpha * deltaTime
				|| equal(deltaTime, 0.0) || equal(alpha, 0.0))
			return;

		// Update the concentration array
		updatedConcOffset[fluxIndices[0]] += incidentFluxVec[0][xi - surfacePos]; // V
		updatedConcOffset[fluxIndices[1]] += incidentFluxVec[0][xi - surfacePos]; // I

		return;
	}

	/**
	 * This operation sets the time of the pulse.
	 * \see IFluxHandler.h
	 */
	void setPulseTime(double time) {
		deltaTime = time;
		return;
	}

	/**
	 * This operation sets proportion of the pulse that is on.
	 * \see IFluxHandler.h
	 */
	void setProportion(double a) {
		alpha = a;
		return;
	}

};
//end class PulsedFitFluxHandler

}

#endif
