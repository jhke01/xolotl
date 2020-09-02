#ifndef UZRCLUSTERNETWORKLOADER_H_
#define UZRCLUSTERNETWORKLOADER_H_

//Includes
#include <UZrCluster.h>
#include <NetworkLoader.h>
#include <UZrClusterReactionNetwork.h>

namespace xolotlCore {

/**
 * This class will load a reaction network composed of UZrClusters from an
 * HDF5 file or generate it.
 *
 * The network will be returned as a ReactionNetwork of UZrClusters ordered with
 * single-species Xe, V clusters first and all mixed clusters coming
 * last. Each species is ordered from the smallest cluster size, (1), to the
 * maximum size for that cluster. Instances of the appropriate cluster type are
 * instantiated during the loading process, but returned as UZrClusters.
 */
class UZrClusterNetworkLoader: public NetworkLoader {

protected:

	/**
	 * The xenon size at which the grouping scheme starts
	 */
	int xeMin;

	/**
	 * The xenon size at which the grouping scheme ends
	 */
	int xeMax;

	/**
	 * The width of the group in the xenon direction.
	 */
	int sectionWidth;


	/**
	 * The maximum size for xenon clusters
	 */
	int maxXe;

	/**
	 * The maximum size for vacancy clusters
	 */
	int maxV;

	/**
	 * Private nullary constructor.
	 */
	UZrClusterNetworkLoader() :
			NetworkLoader(), maxXe(0), maxV(0), xeMin(-1), xeMax(-1), sectionWidth(-1) {
	}

	/**
	 * This operation creates a normal cluster.
	 *
	 * @param numXe The number of xenon atoms
	 * @param numV The number of atomic vacancies
	 * @return The new cluster
	 */
	std::unique_ptr<UZrCluster> createUZrCluster(int numXe, int numV,
			IReactionNetwork &network) const;

	/**
	 * This operation creates a super cluster
	 *
	 * @param nTot The total number of clusters
	 * @param maxXe The maximum number of xenon
	 * @return The new cluster
	 */
	std::unique_ptr<UZrCluster> createUZrSuperCluster(int nTot, int maxXe,
			IReactionNetwork& network) const;

	/**
	 * This operation will add the given cluster to the network and reactants vector
	 * as a standard cluster or a dummy one if we do not want the reactions to happen.
	 *
	 * @param network The network
	 * @param reactants The vector of reactants kept by the loader
	 * @param cluster The cluster to add to them
	 */
	virtual void pushUZrCluster(
			std::unique_ptr<UZrClusterReactionNetwork> & network,
			std::vector<std::reference_wrapper<Reactant> > & reactants,
			std::unique_ptr<UZrCluster> & cluster);

public:

	/**
	 * The default constructor. The setInputstream() operation must be called
	 * if this constructor is used.
	 *
	 * @param registry The performance handler registry
	 */
	UZrClusterNetworkLoader(
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * An alternative constructor provided for convenience.
	 *
	 * @param stream The inputstream from which the cluster data should be
	 * loaded
	 * @param registry The performance handler registry
	 */
	UZrClusterNetworkLoader(const std::shared_ptr<std::istream> stream,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * Destructor
	 */
	virtual ~UZrClusterNetworkLoader() {
	}

	/**
	 * This operation will load the reaction network from the inputstream in
	 * the format specified previously. The network will be empty if it can not
	 * be loaded.
	 *
	 * @return network The reaction network
	 */
	virtual std::unique_ptr<IReactionNetwork> load(const IOptions &options)
			override;

	/**
	 * This operation will generate the reaction network from options.
	 * The network will be empty if it can not be loaded.
	 *
	 * @param options The command line options
	 * @return network The reaction network
	 */
	virtual std::unique_ptr<IReactionNetwork> generate(const IOptions &options)
			override;

	/**
	 * This operation will apply a grouping method to the network.
	 *
	 * @param The network to be modified.
	 */
	void applyGrouping(IReactionNetwork& network) const;

	/**
	 * This operation will set the xenon size at which the grouping scheme starts.
	 *
	 * @param min The value for the size
	 */
	void setXeMin(int min) {
		xeMin = min;
	}

	/**
	 * This operation will set the xenon width for the grouping scheme.
	 *
	 * @param w The value of the width
	 */
	void setWidth(int w) {
		sectionWidth = w;
	}

};

} /* namespace xolotlCore */

#endif /* UZRCLUSTERNETWORKLOADER_H_ */
