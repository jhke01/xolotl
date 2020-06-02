// Includes
#include <petscts.h>
#include <petscsys.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <memory>
#include <xolotl/solver/PetscSolver.h>
#include <xolotl/perf/xolotlPerf.h>
#include <xolotl/factory/viz/VizHandlerRegistryFactory.h>
#include <xolotl/viz/PlotType.h>
#include <xolotl/viz/dataprovider/CvsXDataProvider.h>
#include <xolotl/viz/dataprovider/CvsXYDataProvider.h>
#include <xolotl/viz/LabelProvider.h>
#include <xolotl/core/Constants.h>
#include <xolotl/core/network/NEReactionNetwork.h>
#include <xolotl/core/network/PSIReactionNetwork.h>
#include <xolotl/core/MathUtils.h>
#include <xolotl/solver/handler/RandomNumberGenerator.h>
#include <xolotl/io/XFile.h>
#include <xolotl/solver/monitor/Monitor.h>

namespace xolotl {
namespace solver {
namespace monitor {

// Declaration of the functions defined in Monitor.cpp
extern PetscErrorCode checkTimeStep(TS ts);
extern PetscErrorCode monitorTime(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx);
extern PetscErrorCode computeFluence(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx);
extern PetscErrorCode monitorPerf(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx);

// Declaration of the variables defined in Monitor.cpp
extern std::shared_ptr<viz::IPlot> perfPlot;
extern double previousTime;
extern double timeStepThreshold;

//! The pointer to the plot used in monitorScatter1D.
std::shared_ptr<viz::IPlot> scatterPlot1D;
//! The pointer to the series plot used in monitorSeries1D.
std::shared_ptr<viz::IPlot> seriesPlot1D;
//! The variable to store the interstitial flux at the previous time step.
double previousIFlux1D = 0.0;
//! The variable to store the total number of interstitials going through the surface.
double nInterstitial1D = 0.0;
//! The variable to store the helium flux at the previous time step.
double previousHeFlux1D = 0.0;
//! The variable to store the total number of helium going through the bottom.
double nHelium1D = 0.0;
//! The variable to store the deuterium flux at the previous time step.
double previousDFlux1D = 0.0;
//! The variable to store the total number of deuterium going through the bottom.
double nDeuterium1D = 0.0;
//! The variable to store the tritium flux at the previous time step.
double previousTFlux1D = 0.0;
//! The variable to store the total number of tritium going through the bottom.
double nTritium1D = 0.0;
//! The variable to store the vacancy flux at the previous time step.
double previousVFlux1D = 0.0;
//! The variable to store the total number of vacancy going through the bottom.
double nVacancy1D = 0.0;
//! The variable to store the int flux at the previous time step.
double previousIBulkFlux1D = 0.0;
//! The variable to store the total number of int going through the bottom.
double nIBulk1D = 0.0;
//! The variable to store the sputtering yield at the surface.
double sputteringYield1D = 0.0;
//! The threshold for the negative concentration
double negThreshold1D = 0.0;
//! How often HDF5 file is written
PetscReal hdf5Stride1D = 0.0;
//! Previous time for HDF5
PetscInt hdf5Previous1D = 0;
//! HDF5 output file name
std::string hdf5OutputName1D = "xolotlStop.h5";
// The vector of depths at which bursting happens
std::vector<int> depthPositions1D;
// The vector of ids for diffusing interstitial clusters
std::vector<int> iClusterIds1D;

// Timers
std::shared_ptr<perf::ITimer> initTimer;
std::shared_ptr<perf::ITimer> checkNegativeTimer;
std::shared_ptr<perf::ITimer> tridynTimer;
std::shared_ptr<perf::ITimer> startStopTimer;
std::shared_ptr<perf::ITimer> heRetentionTimer;
std::shared_ptr<perf::ITimer> xeRetentionTimer;
std::shared_ptr<perf::ITimer> scatterTimer;
std::shared_ptr<perf::ITimer> seriesTimer;
std::shared_ptr<perf::ITimer> eventFuncTimer;
std::shared_ptr<perf::ITimer> postEventFuncTimer;

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "checkNegative1D")
/**
 * This is a monitoring method that looks at if there are negative concentrations at each time step.
 */
PetscErrorCode checkNegative1D(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void*) {

	perf::ScopedTimer myTimer(checkNegativeTimer);

	// Initial declaration
	PetscErrorCode ierr;
	double **solutionArray, *gridPointSolution;
	PetscInt xs, xm;

	PetscFunctionBeginUser;

	// Get the number of processes
	int worldSize;
	MPI_Comm_size(PETSC_COMM_WORLD, &worldSize);

	// Gets the process ID (important when it is running in parallel)
	int procId;
	MPI_Comm_rank(PETSC_COMM_WORLD, &procId);

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the solutionArray
	ierr = DMDAVecGetArrayDOF(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	CHKERRQ(ierr);

	// Get the solver handler
	auto &solverHandler = PetscSolver::getSolverHandler();

	// Get the network and dof
	auto &network = solverHandler.getNetwork();
	const int nClusters = network.getNumClusters();

	// Loop on the local grid
	for (PetscInt i = xs; i < xs + xm; i++) {
		// Get the pointer to the beginning of the solution data for this grid point
		gridPointSolution = solutionArray[i]; // Loop on the concentrations
		for (int l = 0; l < nClusters; l++) {
			if (gridPointSolution[l] < negThreshold1D
					&& gridPointSolution[l] > 0.0) {
				gridPointSolution[l] = negThreshold1D;
			} else if (gridPointSolution[l] > -negThreshold1D
					&& gridPointSolution[l] < 0.0) {
				gridPointSolution[l] = -negThreshold1D;
			}
		}
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOF(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "computeTRIDYN1D")
/**
 * This is a monitoring method that will compute the data to send to TRIDYN
 */
PetscErrorCode computeTRIDYN1D(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx) {

	perf::ScopedTimer myTimer(tridynTimer);

	// Initial declarations
	PetscErrorCode ierr;
	PetscInt xs, xm;

	PetscFunctionBeginUser;

	// Get the number of processes
	int worldSize;
	MPI_Comm_size(PETSC_COMM_WORLD, &worldSize);

	// Gets the process ID
	int procId;
	MPI_Comm_rank(PETSC_COMM_WORLD, &procId);

	// Get the solver handler
	auto &solverHandler = PetscSolver::getSolverHandler();

	// Get the network
	using NetworkType =
	core::network::PSIReactionNetwork<core::network::PSIFullSpeciesList>;
	using Spec = typename NetworkType::Species;
	auto &network = dynamic_cast<NetworkType&>(solverHandler.getNetwork());
	int dof = network.getDOF();

	// Get the position of the surface
	int surfacePos = solverHandler.getSurfacePosition();

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	CHKERRQ(ierr);

	// Get the total size of the grid
	PetscInt Mx;
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE);
	CHKERRQ(ierr);

	// Get the physical grid
	auto grid = solverHandler.getXGrid();

	// Get the array of concentration
	PetscReal **solutionArray;
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Save current concentrations as an HDF5 file.
	//
	// First create the file for parallel file access.
	std::ostringstream tdFileStr;
	tdFileStr << "TRIDYN_" << timestep << ".h5";
	io::HDF5File tdFile(tdFileStr.str(),
			io::HDF5File::AccessMode::CreateOrTruncateIfExists,
			PETSC_COMM_WORLD, true);

	// Define a dataset for concentrations.
	// Everyone must create the dataset with the same shape.
	constexpr auto numConcSpecies = 5;
	constexpr auto numValsPerGridpoint = numConcSpecies + 2;
	const auto firstIdxToWrite = (surfacePos + solverHandler.getLeftOffset());
	const auto numGridpointsWithConcs = (Mx - firstIdxToWrite);
	io::HDF5File::SimpleDataSpace<2>::Dimensions concsDsetDims = {
			(hsize_t) numGridpointsWithConcs, numValsPerGridpoint };
	io::HDF5File::SimpleDataSpace<2> concsDsetSpace(concsDsetDims);

	const std::string concsDsetName = "concs";
	io::HDF5File::DataSet<double> concsDset(tdFile, concsDsetName,
			concsDsetSpace);

	// Specify the concentrations we will write.
	// We only consider our own grid points.
	const auto myFirstIdxToWrite = std::max(xs, firstIdxToWrite);
	auto myEndIdx = (xs + xm);  // "end" in the C++ sense; i.e., one-past-last
	auto myNumPointsToWrite =
			(myEndIdx > myFirstIdxToWrite) ? (myEndIdx - myFirstIdxToWrite) : 0;
	io::HDF5File::DataSet<double>::DataType2D<numValsPerGridpoint> myConcs(
			myNumPointsToWrite);

	for (auto xi = myFirstIdxToWrite; xi < myEndIdx; ++xi) {

		if (xi >= firstIdxToWrite) {

			// Determine current gridpoint value.
			double x = (grid[xi] + grid[xi + 1]) / 2.0 - grid[1];

			// Access the solution data for this grid point.
			auto gridPointSolution = solutionArray[xi];
			using HostUnmanaged =
			Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
			auto hConcs = HostUnmanaged(gridPointSolution, dof);
			auto dConcs = Kokkos::View<double*>("Concentrations", dof);
			deep_copy(dConcs, hConcs);

			// Get the total concentrations at this grid point
			auto currIdx = xi - myFirstIdxToWrite;
			myConcs[currIdx][0] = (x - (grid[surfacePos + 1] - grid[1]));
			myConcs[currIdx][1] = network.getTotalAtomConcentration(dConcs,
					Spec::He, 0);
			myConcs[currIdx][2] = network.getTotalAtomConcentration(dConcs,
					Spec::D, 0);
			myConcs[currIdx][3] = network.getTotalAtomConcentration(dConcs,
					Spec::T, 0);
			myConcs[currIdx][4] = network.getTotalAtomConcentration(dConcs,
					Spec::V, 0);
			myConcs[currIdx][5] = network.getTotalAtomConcentration(dConcs,
					Spec::I, 0);
			myConcs[currIdx][6] = gridPointSolution[dof];
		}
	}

	// Write the concs dataset in parallel.
	// (We write only our part.)
	concsDset.parWrite2D<numValsPerGridpoint>(PETSC_COMM_WORLD,
			myFirstIdxToWrite - firstIdxToWrite, myConcs);

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "startStop1D")
/**
 * This is a monitoring method that update an hdf5 file at each time step.
 */
PetscErrorCode startStop1D(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void*) {

	perf::ScopedTimer myTimer(startStopTimer);

	// Initial declaration
	PetscErrorCode ierr;
	const double **solutionArray, *gridPointSolution;
	PetscInt xs, xm, Mx;

	PetscFunctionBeginUser;

	// Compute the dt
	double dt = time - previousTime;

	// Don't do anything if it is not on the stride
	if (((int) ((time + dt / 10.0) / hdf5Stride1D) <= hdf5Previous1D)
			&& timestep > 0) {
		PetscFunctionReturn(0);
	}

	// Update the previous time
	if ((int) ((time + dt / 10.0) / hdf5Stride1D) > hdf5Previous1D)
		hdf5Previous1D++;

	// Get the number of processes
	int worldSize;
	MPI_Comm_size(PETSC_COMM_WORLD, &worldSize);

	// Gets the process ID (important when it is running in parallel)
	int procId;
	MPI_Comm_rank(PETSC_COMM_WORLD, &procId);

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the solutionArray
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	CHKERRQ(ierr);
	// Get the size of the total grid
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE);
	CHKERRQ(ierr);

	// Get the solver handler
	auto &solverHandler = PetscSolver::getSolverHandler();

	// Get the network and dof
	auto &network = solverHandler.getNetwork();
	const int dof = network.getDOF();

	// Create an array for the concentration
	double concArray[dof][2];

	// Get the position of the surface
	int surfacePos = solverHandler.getSurfacePosition();

	// Open the existing HDF5 file
	io::XFile checkpointFile(hdf5OutputName1D, PETSC_COMM_WORLD,
			io::XFile::AccessMode::OpenReadWrite);

	// Get the current time step
	double currentTimeStep;
	ierr = TSGetTimeStep(ts, &currentTimeStep);
	CHKERRQ(ierr);

	// Add a concentration time step group for the current time step.
	auto concGroup = checkpointFile.getGroup<
			io::XFile::ConcentrationGroup>();
	assert(concGroup);
	auto tsGroup = concGroup->addTimestepGroup(timestep, time, previousTime,
			currentTimeStep);

	if (solverHandler.moveSurface()) {
		// Write the surface positions and the associated interstitial quantities
		// in the concentration sub group
		tsGroup->writeSurface1D(surfacePos, nInterstitial1D, previousIFlux1D);
	}

	// Write the bottom impurity information if the bottom is a free surface
	if (solverHandler.getRightOffset() == 1)
		tsGroup->writeBottom1D(nHelium1D, previousHeFlux1D, nDeuterium1D,
				previousDFlux1D, nTritium1D, previousTFlux1D, nVacancy1D,
				previousVFlux1D, nIBulk1D, previousIBulkFlux1D);

	// Determine the concentration values we will write.
	// We only examine and collect the grid points we own.
	// TODO measure impact of us building the flattened representation
	// rather than a ragged 2D representation.
    io::XFile::TimestepGroup::Concs1DType concs(xm);
	for (auto i = 0; i < xm; ++i) {

		// Access the solution data for the current grid point.
		auto gridPointSolution = solutionArray[xs + i];

		for (auto l = 0; l < dof + 1; ++l) {
			if (std::fabs(gridPointSolution[l]) > 1.0e-16) {
				concs[i].emplace_back(l, gridPointSolution[l]);
			}
		}
	}

	// Write our concentration data to the current timestep group
	// in the HDF5 file.
	// We only write the data for the grid points we own.
	tsGroup->writeConcentrations(checkpointFile, xs, concs);

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	ierr = computeTRIDYN1D(ts, timestep, time, solution, NULL);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "computeHeliumDesorption1D")
/**
 * This is a monitoring method that will compute the helium desorption at the surface
 */
PetscErrorCode computeHeliumDesorption1D(TS ts, PetscInt, PetscReal time,
		Vec solution, void*) {

	// Initial declarations
	PetscErrorCode ierr;
	PetscInt xs, xm;

	PetscFunctionBeginUser;

	// Get the solver handler
	auto &solverHandler = PetscSolver::getSolverHandler();

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	CHKERRQ(ierr);

	// Get the total size of the grid
	PetscInt Mx;
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE);
	CHKERRQ(ierr);

	// Get the physical grid
	auto grid = solverHandler.getXGrid();
	// Get the position of the surface
	int surfacePos = solverHandler.getSurfacePosition();

	// Get the network
	using NetworkType =
	core::network::PSIReactionNetwork<core::network::PSIFullSpeciesList>;
	using Spec = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;
	auto &network = dynamic_cast<NetworkType&>(solverHandler.getNetwork());
	const int dof = network.getDOF();

	// Get the array of concentration
	PetscReal **solutionArray;
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Store the He concentration at the surface
	double heConc = 0.0, diffCoeff = 0.0;
	Composition comp = Composition::zero();
	comp[Spec::He] = 1;
	auto heCluster = network.findCluster(comp, plsm::onHost);
	int heIndex = heCluster.getId();

	// Declare the pointer for the concentrations at a specific grid point
	PetscReal *gridPointSolution;

	// Temperature
	double localTemp = 0.0;

	// Loop on the grid
	for (PetscInt xi = xs; xi < xs + xm; xi++) {
		// Get the pointer to the beginning of the solution data for this grid point
		gridPointSolution = solutionArray[xi];

		// Check if we are next to the surface
		if (xi == surfacePos + 1) {
			heConc = gridPointSolution[heIndex];
			diffCoeff = heCluster.getDiffusionCoefficient(xi - xs);
			localTemp = gridPointSolution[dof];
		}
	}

	// Get the current process ID
	int procId;
	MPI_Comm_rank(PETSC_COMM_WORLD, &procId);

	// Send the concentration to proc Id 0
	double localFactor = heConc * diffCoeff, factor = 0.0;
	MPI_Reduce(&localFactor, &factor, 1, MPI_DOUBLE, MPI_SUM, 0,
			PETSC_COMM_WORLD);
	double temperature = 0.0;
	MPI_Reduce(&localTemp, &temperature, 1, MPI_DOUBLE, MPI_SUM, 0,
			PETSC_COMM_WORLD);

	// Master process
	if (procId == 0) {
		double hxLeft = 0.0;
		if (surfacePos < 0) {
			hxLeft = grid[surfacePos + 2] - grid[surfacePos + 1];
		} else {
			hxLeft = (grid[surfacePos + 2] - grid[surfacePos]) / 2.0;
		}
		double surfaceFlux = factor * hxLeft;
		// Write the flux at the boundary and temperature in a file
		std::ofstream outputFile;
		outputFile.open("thds.txt", std::ios::app);
		outputFile << temperature << " " << surfaceFlux << std::endl;
		outputFile.close();
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "computeHeliumRetention1D")
/**
 * This is a monitoring method that will compute the helium retention
 */
PetscErrorCode computeHeliumRetention1D(TS ts, PetscInt, PetscReal time,
		Vec solution, void*) {
	perf::ScopedTimer myTimer(heRetentionTimer);

	// Initial declarations
	PetscErrorCode ierr;
	PetscInt xs, xm;

	PetscFunctionBeginUser;

	// Get the solver handler
	auto &solverHandler = PetscSolver::getSolverHandler();

	// Get the flux handler that will be used to know the fluence
	auto fluxHandler = solverHandler.getFluxHandler();

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	CHKERRQ(ierr);

	// Get the total size of the grid
	PetscInt Mx;
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE);
	CHKERRQ(ierr);

	// Get the physical grid
	auto grid = solverHandler.getXGrid();
	// Get the position of the surface
	int surfacePos = solverHandler.getSurfacePosition();

	// Get the network
	using NetworkType =
	core::network::PSIReactionNetwork<core::network::PSIFullSpeciesList>;
	using Spec = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;
	auto &network = dynamic_cast<NetworkType&>(solverHandler.getNetwork());
	const int dof = network.getDOF();

	// Get the array of concentration
	PetscReal **solutionArray;
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Store the concentration over the grid
	double heConcentration = 0.0, dConcentration = 0.0, tConcentration = 0.0,
			vConcentration = 0.0, iConcentration = 0.0;

	// Declare the pointer for the concentrations at a specific grid point
	PetscReal *gridPointSolution;

	// Loop on the grid
	for (PetscInt xi = xs; xi < xs + xm; xi++) {

		// Boundary conditions
		if (xi < surfacePos + solverHandler.getLeftOffset()
				|| xi >= Mx - solverHandler.getRightOffset())
			continue;

		// Get the pointer to the beginning of the solution data for this grid point
		gridPointSolution = solutionArray[xi];

		double hx = grid[xi + 1] - grid[xi];

		using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
		auto hConcs = HostUnmanaged(gridPointSolution, dof);
		auto dConcs = Kokkos::View<double*>("Concentrations", dof);
		deep_copy(dConcs, hConcs);

		// Get the total concentrations at this grid point
		heConcentration += network.getTotalAtomConcentration(dConcs, Spec::He,
				0) * hx;
		dConcentration += network.getTotalAtomConcentration(dConcs, Spec::D, 0)
				* hx;
		tConcentration += network.getTotalAtomConcentration(dConcs, Spec::T, 0)
				* hx;
		vConcentration += network.getTotalAtomConcentration(dConcs, Spec::V, 0)
				* hx;
		iConcentration += network.getTotalAtomConcentration(dConcs, Spec::I, 0)
				* hx;
	}

	// Get the current process ID
	int procId;
	MPI_Comm_rank(PETSC_COMM_WORLD, &procId);

	// Determine total concentrations for He, D, T.
	std::array<double, 5> myConcData { heConcentration, dConcentration,
			tConcentration, vConcentration, iConcentration };
	std::array<double, 5> totalConcData;

	MPI_Reduce(myConcData.data(), totalConcData.data(), myConcData.size(),
	MPI_DOUBLE,
	MPI_SUM, 0, PETSC_COMM_WORLD);

	// Extract total He, D, T concentrations.  Values are valid only on rank 0.
	double totalHeConcentration = totalConcData[0];
	double totalDConcentration = totalConcData[1];
	double totalTConcentration = totalConcData[2];
	double totalVConcentration = totalConcData[3];
	double totalIConcentration = totalConcData[4];

	// TODO: will need to be refactored: getDiff ids from diffusion handler?
//	// Look at the fluxes going in the bulk if the bottom is a free surface
//	if (solverHandler.getRightOffset() == 1) {
//		// Set the bottom surface position
//		int xi = Mx - 2;
//
//		// Value to know on which processor is the bottom
//		int bottomProc = 0;
//
//		// Check we are on the right proc
//		if (xi >= xs && xi < xs + xm) {
//			// Get the delta time from the previous timestep to this timestep
//			double dt = time - previousTime;
//			// Compute the total number of impurities that went in the bulk
//			nHelium1D += previousHeFlux1D * dt;
//			nDeuterium1D += previousDFlux1D * dt;
//			nTritium1D += previousTFlux1D * dt;
//			nVacancy1D += previousVFlux1D * dt;
//			nIBulk1D += previousIBulkFlux1D * dt;
//
//			// Get the pointer to the beginning of the solution data for this grid point
//			gridPointSolution = solutionArray[xi];
//
//			// Factor for finite difference
//			double hxLeft = 0.0, hxRight = 0.0;
//			if (xi - 1 >= 0 && xi < Mx) {
//				hxLeft = (grid[xi + 1] - grid[xi - 1]) / 2.0;
//				hxRight = (grid[xi + 2] - grid[xi]) / 2.0;
//			} else if (xi - 1 < 0) {
//				hxLeft = grid[xi + 1] - grid[xi];
//				hxRight = (grid[xi + 2] - grid[xi]) / 2.0;
//			} else {
//				hxLeft = (grid[xi + 1] - grid[xi - 1]) / 2.0;
//				hxRight = grid[xi + 1] - grid[xi];
//			}
//			double factor = 2.0 / (hxRight * (hxLeft + hxRight));
//
//			// Initialize the value for the flux
//			double newFlux = 0.0;
//			// Consider each helium cluster.
//			for (auto const& heMapItem : network.getAll(ReactantType::He)) {
//				// Get the cluster
//				auto const& cluster = *(heMapItem.second);
//				// Get its id and concentration
//				int id = cluster.getId() - 1;
//				double conc = gridPointSolution[id];
//				// Get its size and diffusion coefficient
//				int size = cluster.getSize();
//				double coef = cluster.getDiffusionCoefficient(xi - xs);
//				// Compute the flux going to the right
//				newFlux += (double) size * factor * coef * conc * hxRight;
//			}
//			// Update the helium flux
//			previousHeFlux1D = newFlux;
//
//			// Initialize the value for the flux
//			newFlux = 0.0;
//			// Consider each deuterium cluster.
//			for (auto const& dMapItem : network.getAll(ReactantType::D)) {
//				// Get the cluster
//				auto const& cluster = *(dMapItem.second);
//				// Get its id and concentration
//				int id = cluster.getId() - 1;
//				double conc = gridPointSolution[id];
//				// Get its size and diffusion coefficient
//				int size = cluster.getSize();
//				double coef = cluster.getDiffusionCoefficient(xi - xs);
//				// Compute the flux going to the right
//				newFlux += (double) size * factor * coef * conc * hxRight;
//			}
//			// Update the deuterium flux
//			previousDFlux1D = newFlux;
//
//			// Initialize the value for the flux
//			newFlux = 0.0;
//			// Consider each tritium cluster.
//			for (auto const& tMapItem : network.getAll(ReactantType::T)) {
//				// Get the cluster
//				auto const& cluster = *(tMapItem.second);
//				// Get its id and concentration
//				int id = cluster.getId() - 1;
//				double conc = gridPointSolution[id];
//				// Get its size and diffusion coefficient
//				int size = cluster.getSize();
//				double coef = cluster.getDiffusionCoefficient(xi - xs);
//				// Compute the flux going to the right
//				newFlux += (double) size * factor * coef * conc * hxRight;
//			}
//			// Update the tritium flux
//			previousTFlux1D = newFlux;
//
//			// Initialize the value for the flux
//			newFlux = 0.0;
//			// Consider each vacancy cluster.
//			for (auto const& vMapItem : network.getAll(ReactantType::V)) {
//				// Get the cluster
//				auto const& cluster = *(vMapItem.second);
//				// Get it diffusion coefficient
//				double coef = cluster.getDiffusionCoefficient(xi - xs);
//				if (coef <= 0.0)
//					continue;
//				// Get its id and concentration
//				int id = cluster.getId() - 1;
//				double conc = gridPointSolution[id];
//				// Get its size
//				int size = cluster.getSize();
//				// Compute the flux going to the right
//				newFlux += (double) size * factor * coef * conc * hxRight;
//			}
//			// Update the tritium flux
//			previousVFlux1D = newFlux;
//
//			// Initialize the value for the flux
//			newFlux = 0.0;
//			// Consider each int cluster.
//			for (auto const& iMapItem : network.getAll(ReactantType::I)) {
//				// Get the cluster
//				auto const& cluster = *(iMapItem.second);
//				// Get its id and concentration
//				int id = cluster.getId() - 1;
//				double conc = gridPointSolution[id];
//				// Get its size and diffusion coefficient
//				int size = cluster.getSize();
//				double coef = cluster.getDiffusionCoefficient(xi - xs);
//				// Compute the flux going to the right
//				newFlux += (double) size * factor * coef * conc * hxRight;
//			}
//			// Update the tritium flux
//			previousIBulkFlux1D = newFlux;
//
//			// Set the bottom processor
//			bottomProc = procId;
//		}
//
//		// Get which processor will send the information
//		// TODO do we need to do this allreduce just to figure out
//		// who owns the data?
//		// And is it supposed to be a sum?   Why not a min?
//		int bottomId = 0;
//		MPI_Allreduce(&bottomProc, &bottomId, 1, MPI_INT, MPI_SUM,
//				PETSC_COMM_WORLD);
//
//		// Send the information about impurities
//		// to the other processes
//		std::array<double, 10> countFluxData { nHelium1D, previousHeFlux1D,
//				nDeuterium1D, previousDFlux1D, nTritium1D, previousTFlux1D,
//				nVacancy1D, previousVFlux1D, nIBulk1D, previousIBulkFlux1D };
//		MPI_Bcast(countFluxData.data(), countFluxData.size(), MPI_DOUBLE,
//				bottomId, PETSC_COMM_WORLD);
//
//		// Extract inpurity data from broadcast buffer.
//		nHelium1D = countFluxData[0];
//		previousHeFlux1D = countFluxData[1];
//		nDeuterium1D = countFluxData[2];
//		previousDFlux1D = countFluxData[3];
//		nTritium1D = countFluxData[4];
//		previousTFlux1D = countFluxData[5];
//		nVacancy1D = countFluxData[6];
//		previousVFlux1D = countFluxData[7];
//		nIBulk1D = countFluxData[8];
//		previousIBulkFlux1D = countFluxData[9];
//	}

	// Master process
	if (procId == 0) {
		// Get the fluence
		double fluence = fluxHandler->getFluence();

		// Print the result
		std::cout << "\nTime: " << time << std::endl;
		std::cout << "Helium content = " << totalHeConcentration << std::endl;
		std::cout << "Deuterium content = " << totalDConcentration << std::endl;
		std::cout << "Tritium content = " << totalTConcentration << std::endl;
		std::cout << "Vacancy content = " << totalVConcentration << std::endl;
		std::cout << "Interstitial content = " << totalIConcentration
				<< std::endl;
		std::cout << "Fluence = " << fluence << "\n" << std::endl;

		// Uncomment to write the retention and the fluence in a file
		std::ofstream outputFile;
		outputFile.open("retentionOut.txt", std::ios::app);
		outputFile << fluence << " " << totalHeConcentration << " "
				<< totalDConcentration << " " << totalTConcentration << " "
				<< totalVConcentration << " " << totalIConcentration << " "
				<< nHelium1D << " " << nDeuterium1D << " " << nTritium1D << " "
				<< nVacancy1D << " " << nIBulk1D << std::endl;
		outputFile.close();
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "computeXenonRetention1D")
/**
 * This is a monitoring method that will compute the xenon retention
 */
PetscErrorCode computeXenonRetention1D(TS ts, PetscInt, PetscReal time,
		Vec solution, void*) {

	perf::ScopedTimer myTimer(xeRetentionTimer);

	// Initial declarations
	PetscErrorCode ierr;
	PetscInt xs, xm;

	PetscFunctionBeginUser;

	// Get the solver handler
	auto &solverHandler = PetscSolver::getSolverHandler();

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	CHKERRQ(ierr);

	// Get the total size of the grid
	PetscInt Mx;
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE);
	CHKERRQ(ierr);

	// Get the physical grid
	auto grid = solverHandler.getXGrid();

	using NetworkType =
	core::network::NEReactionNetwork;
	using Spec = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;

	// Degrees of freedom is the total number of clusters in the network
	auto &network = dynamic_cast<NetworkType&>(solverHandler.getNetwork());
	const int dof = network.getDOF();

	// Get the array of concentration
	PetscReal **solutionArray;
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Store the concentration and other values over the grid
	double xeConcentration = 0.0, bubbleConcentration = 0.0, radii = 0.0,
			partialBubbleConcentration = 0.0, partialRadii = 0.0, partialSize =
					0.0;

	// Declare the pointer for the concentrations at a specific grid point
	PetscReal *gridPointSolution;

	// Get the minimum size for the radius
	auto minSizes = solverHandler.getMinSizes();

	// Loop on the grid
	for (PetscInt xi = xs; xi < xs + xm; xi++) {

		// Get the pointer to the beginning of the solution data for this grid point
		gridPointSolution = solutionArray[xi];

		using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
		auto hConcs = HostUnmanaged(gridPointSolution, dof);
		auto dConcs = Kokkos::View<double*>("Concentrations", dof);
		deep_copy(dConcs, hConcs);

		double hx = grid[xi + 1] - grid[xi];

		// Get the concentrations
		xeConcentration += network.getTotalAtomConcentration(dConcs, Spec::Xe,
				0) * hx;
		bubbleConcentration += network.getTotalConcentration(dConcs, Spec::Xe,
				0) * hx;
		radii += network.getTotalRadiusConcentration(dConcs, Spec::Xe, 0) * hx;
		partialBubbleConcentration = network.getTotalConcentration(dConcs,
				Spec::Xe, minSizes[0]) * hx;
		partialRadii += network.getTotalRadiusConcentration(dConcs, Spec::Xe,
				minSizes[0]) * hx;
		partialSize += network.getTotalAtomConcentration(dConcs, Spec::Xe,
				minSizes[0]) * hx;

//		if (xi == 1) {
//				std::ofstream outFile;
//				outFile.open("size.txt");
//				for (unsigned int i = 0; i < network.getNumClusters(); i++) {
//					// Add the current concentration times the number of xenon in the cluster
//					// (from the weight vector)
//					double conc = gridPointSolution[i];
//					const auto& clReg = network.getCluster(i, plsm::onHost).getRegion();
//					Composition lo = clReg.getOrigin();
//					outFile << lo[Spec::Xe] << " " << conc << std::endl;
//				}
//				outFile.close();
//		}
	}

	// Get the current process ID
	int procId;
	MPI_Comm_rank(PETSC_COMM_WORLD, &procId);

	// Sum all the concentrations through MPI reduce
	std::array<double, 6> myConcData { xeConcentration, bubbleConcentration,
			radii, partialBubbleConcentration, partialRadii, partialSize };
	std::array<double, 6> totalConcData { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	MPI_Reduce(myConcData.data(), totalConcData.data(), myConcData.size(),
	MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);

	// Master process
	if (procId == 0) {
		// Print the result
		std::cout << "\nTime: " << time << std::endl;
		std::cout << "Xenon concentration = " << totalConcData[0] << std::endl
				<< std::endl;

		// Make sure the average partial radius makes sense
		double averagePartialRadius = 0.0, averagePartialSize = 0.0;
		if (totalConcData[3] > 1.e-16) {
			averagePartialRadius = totalConcData[4] / totalConcData[3];
			averagePartialSize = totalConcData[5] / totalConcData[3];
		}

		// Uncomment to write the content in a file
		std::ofstream outputFile;
		outputFile.open("retentionOut.txt", std::ios::app);
		outputFile << time << " " << totalConcData[0] << " "
				<< totalConcData[2] / totalConcData[1] << " "
				<< averagePartialRadius << " " << totalConcData[3] << " "
				<< averagePartialSize << std::endl;
		outputFile.close();
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "profileTemperature1D")
/**
 * This is a monitoring method that will store the temperature profile
 */
PetscErrorCode profileTemperature1D(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx) {
	// Initial declarations
	PetscErrorCode ierr;
	PetscInt xs, xm;

	PetscFunctionBeginUser;

	// Get the number of processes
	int worldSize;
	MPI_Comm_size(PETSC_COMM_WORLD, &worldSize);

	// Gets the process ID
	int procId;
	MPI_Comm_rank(PETSC_COMM_WORLD, &procId);

	// Get the solver handler
	auto &solverHandler = PetscSolver::getSolverHandler();

	// Get the network and dof
	auto &network = solverHandler.getNetwork();
	const int dof = network.getDOF();

	// Get the position of the surface
	int surfacePos = solverHandler.getSurfacePosition();

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	CHKERRQ(ierr);

	// Get the total size of the grid
	PetscInt Mx;
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE);
	CHKERRQ(ierr);

	// Get the physical grid
	auto grid = solverHandler.getXGrid();

	// Get the array of concentration
	PetscReal **solutionArray;
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Declare the pointer for the concentrations at a specific grid point
	PetscReal *gridPointSolution;

	// Create the output file
	std::ofstream outputFile;
	if (procId == 0) {
		outputFile.open("tempProf.txt", std::ios::app);
		outputFile << time;
	}

	// Loop on the entire grid
	for (int xi = surfacePos + solverHandler.getLeftOffset();
			xi < Mx - solverHandler.getRightOffset(); xi++) {
		// Set x
		double x = (grid[xi] + grid[xi + 1]) / 2.0 - grid[1];

		double localTemp = 0.0;
		// Check if this process is in charge of xi
		if (xi >= xs && xi < xs + xm) {
			// Get the pointer to the beginning of the solution data for this grid point
			gridPointSolution = solutionArray[xi];

			// Get the local temperature
			localTemp = gridPointSolution[dof];
		}

		// Get the value on procId = 0
		double temperature = 0.0;
		MPI_Reduce(&localTemp, &temperature, 1, MPI_DOUBLE,
		MPI_SUM, 0, PETSC_COMM_WORLD);

		// The master process writes in the file
		if (procId == 0) {
			outputFile << " " << temperature;
		}
	}

	// Close the file
	if (procId == 0) {
		outputFile << std::endl;
		outputFile.close();
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "computeAlloy1D")
/**
 * This is a monitoring method that will compute average density and diameter
 */
PetscErrorCode computeAlloy1D(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx) {

	// Initial declarations
	PetscErrorCode ierr;
	PetscInt xs, xm;

	PetscFunctionBeginUser;

//	// Get the number of processes
//	int worldSize;
//	MPI_Comm_size(PETSC_COMM_WORLD, &worldSize);
//
//	// Gets the process ID
//	int procId;
//	MPI_Comm_rank(PETSC_COMM_WORLD, &procId);
//
//	// Get the solver handler
//	auto &solverHandler = PetscSolver::getSolverHandler();
//
//	// Get the physical grid and its length
//	auto grid = solverHandler.getXGrid();
//	int xSize = grid.size();
//
//	// Get the position of the surface
//	int surfacePos = solverHandler.getSurfacePosition();
//
//	// Get the da from ts
//	DM da;
//	ierr = TSGetDM(ts, &da);
//	CHKERRQ(ierr);
//
//	// Get the corners of the grid
//	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
//	CHKERRQ(ierr);
//
//	// Get the total size of the grid
//	PetscInt Mx;
//	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE,
//	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
//	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
//	PETSC_IGNORE);
//	CHKERRQ(ierr);
//
//	// Get the array of concentration
//	PetscReal **solutionArray;
//	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
//	CHKERRQ(ierr);
//
//	// Get the network
//	auto &network = solverHandler.getNetwork();
//
//	// Initial declarations for the density and diameter
//	double iDensity = 0.0, vDensity = 0.0, voidDensity = 0.0,
//			frankDensity = 0.0, faultedDensity = 0.0, perfectDensity = 0.0,
//			voidPartialDensity = 0.0, frankPartialDensity = 0.0,
//			faultedPartialDensity = 0.0, perfectPartialDensity = 0.0,
//			iDiameter = 0.0, vDiameter = 0.0, voidDiameter = 0.0,
//			frankDiameter = 0.0, faultedDiameter = 0.0, perfectDiameter = 0.0,
//			voidPartialDiameter = 0.0, frankPartialDiameter = 0.0,
//			faultedPartialDiameter = 0.0, perfectPartialDiameter = 0.0;
//
//	// Get the minimum size for the loop densities and diameters
//	auto minSizes = solverHandler.getMinSizes();
//
//	// Declare the pointer for the concentrations at a specific grid point
//	PetscReal *gridPointSolution;
//
//	// Loop on the grid
//	for (PetscInt xi = xs; xi < xs + xm; xi++) {
//
//		// Boundary conditions
//		if (xi < surfacePos + solverHandler.getLeftOffset()
//				|| xi == Mx - solverHandler.getRightOffset())
//			continue;
//
//		// Get the pointer to the beginning of the solution data for this grid point
//		gridPointSolution = solutionArray[xi];
//
//		// Update the concentration in the network
//		network.updateConcentrationsFromArray(gridPointSolution);
//
//		// Loop on I
//		for (auto const &iMapItem : network.getAll(ReactantType::I)) {
//			// Get the cluster
//			auto const &cluster = *(iMapItem.second);
//			iDensity += gridPointSolution[cluster.getId() - 1];
//			iDiameter += gridPointSolution[cluster.getId() - 1]
//					* cluster.getReactionRadius() * 2.0;
//		}
//
//		// Loop on V
//		for (auto const &vMapItem : network.getAll(ReactantType::V)) {
//			// Get the cluster
//			auto const &cluster = *(vMapItem.second);
//			vDensity += gridPointSolution[cluster.getId() - 1];
//			vDiameter += gridPointSolution[cluster.getId() - 1]
//					* cluster.getReactionRadius() * 2.0;
//		}
//
//		// Loop on Void
//		for (auto const &voidMapItem : network.getAll(ReactantType::Void)) {
//			// Get the cluster
//			auto const &cluster = *(voidMapItem.second);
//			voidDensity += gridPointSolution[cluster.getId() - 1];
//			voidDiameter += gridPointSolution[cluster.getId() - 1]
//					* cluster.getReactionRadius() * 2.0;
//			if (cluster.getSize() >= minSizes[0]) {
//				voidPartialDensity += gridPointSolution[cluster.getId() - 1];
//				voidPartialDiameter += gridPointSolution[cluster.getId() - 1]
//						* cluster.getReactionRadius() * 2.0;
//			}
//		}
//		for (auto const &voidMapItem : network.getAll(ReactantType::VoidSuper)) {
//			// Get the cluster
//			auto const &cluster =
//					static_cast<AlloySuperCluster&>(*(voidMapItem.second));
//			voidDensity += cluster.getTotalConcentration();
//			voidDiameter += cluster.getTotalConcentration()
//					* cluster.getReactionRadius() * 2.0;
//			if (cluster.getSize() >= minSizes[0]) {
//				voidPartialDensity += cluster.getTotalConcentration();
//				voidPartialDiameter += cluster.getTotalConcentration()
//						* cluster.getReactionRadius() * 2.0;
//			}
//		}
//
//		// Loop on Faulted
//		for (auto const &faultedMapItem : network.getAll(ReactantType::Faulted)) {
//			// Get the cluster
//			auto const &cluster = *(faultedMapItem.second);
//			faultedDensity += gridPointSolution[cluster.getId() - 1];
//			faultedDiameter += gridPointSolution[cluster.getId() - 1]
//					* cluster.getReactionRadius() * 2.0;
//			if (cluster.getSize() >= minSizes[1]) {
//				faultedPartialDensity += gridPointSolution[cluster.getId() - 1];
//				faultedPartialDiameter += gridPointSolution[cluster.getId() - 1]
//						* cluster.getReactionRadius() * 2.0;
//			}
//		}
//		for (auto const &faultedMapItem : network.getAll(
//				ReactantType::FaultedSuper)) {
//			// Get the cluster
//			auto const &cluster =
//					static_cast<AlloySuperCluster&>(*(faultedMapItem.second));
//			faultedDensity += cluster.getTotalConcentration();
//			faultedDiameter += cluster.getTotalConcentration()
//					* cluster.getReactionRadius() * 2.0;
//			if (cluster.getSize() >= minSizes[1]) {
//				faultedPartialDensity += cluster.getTotalConcentration();
//				faultedPartialDiameter += cluster.getTotalConcentration()
//						* cluster.getReactionRadius() * 2.0;
//			}
//		}
//
//		// Loop on Perfect
//		for (auto const &perfectMapItem : network.getAll(ReactantType::Perfect)) {
//			// Get the cluster
//			auto const &cluster = *(perfectMapItem.second);
//			perfectDensity += gridPointSolution[cluster.getId() - 1];
//			perfectDiameter += gridPointSolution[cluster.getId() - 1]
//					* cluster.getReactionRadius() * 2.0;
//			if (cluster.getSize() >= minSizes[2]) {
//				perfectPartialDensity += gridPointSolution[cluster.getId() - 1];
//				perfectPartialDiameter += gridPointSolution[cluster.getId() - 1]
//						* cluster.getReactionRadius() * 2.0;
//			}
//		}
//		for (auto const &perfectMapItem : network.getAll(
//				ReactantType::PerfectSuper)) {
//			// Get the cluster
//			auto const &cluster =
//					static_cast<AlloySuperCluster&>(*(perfectMapItem.second));
//			perfectDensity += cluster.getTotalConcentration();
//			perfectDiameter += cluster.getTotalConcentration()
//					* cluster.getReactionRadius() * 2.0;
//			if (cluster.getSize() >= minSizes[2]) {
//				perfectPartialDensity += cluster.getTotalConcentration();
//				perfectPartialDiameter += cluster.getTotalConcentration()
//						* cluster.getReactionRadius() * 2.0;
//			}
//		}
//
//		// Loop on Frank
//		for (auto const &frankMapItem : network.getAll(ReactantType::Frank)) {
//			// Get the cluster
//			auto const &cluster = *(frankMapItem.second);
//			frankDensity += gridPointSolution[cluster.getId() - 1];
//			frankDiameter += gridPointSolution[cluster.getId() - 1]
//					* cluster.getReactionRadius() * 2.0;
//			if (cluster.getSize() >= minSizes[3]) {
//				frankPartialDensity += gridPointSolution[cluster.getId() - 1];
//				frankPartialDiameter += gridPointSolution[cluster.getId() - 1]
//						* cluster.getReactionRadius() * 2.0;
//			}
//		}
//		for (auto const &frankMapItem : network.getAll(ReactantType::FrankSuper)) {
//			// Get the cluster
//			auto const &cluster =
//					static_cast<AlloySuperCluster&>(*(frankMapItem.second));
//			frankDensity += cluster.getTotalConcentration();
//			frankDiameter += cluster.getTotalConcentration()
//					* cluster.getReactionRadius() * 2.0;
//			if (cluster.getSize() >= minSizes[3]) {
//				frankPartialDensity += cluster.getTotalConcentration();
//				frankPartialDiameter += cluster.getTotalConcentration()
//						* cluster.getReactionRadius() * 2.0;
//			}
//		}
//	}
//
//	// Sum all the concentrations through MPI reduce
//	double iTotalDensity = 0.0, vTotalDensity = 0.0, voidTotalDensity = 0.0,
//			frankTotalDensity = 0.0, faultedTotalDensity = 0.0,
//			perfectTotalDensity = 0.0, voidPartialTotalDensity = 0.0,
//			frankPartialTotalDensity = 0.0, faultedPartialTotalDensity = 0.0,
//			perfectPartialTotalDensity = 0.0, iTotalDiameter = 0.0,
//			vTotalDiameter = 0.0, voidTotalDiameter = 0.0, frankTotalDiameter =
//					0.0, faultedTotalDiameter = 0.0, perfectTotalDiameter = 0.0,
//			voidPartialTotalDiameter = 0.0, frankPartialTotalDiameter = 0.0,
//			faultedPartialTotalDiameter = 0.0,
//			perfectPartialTotalDiameter = 0.0;
//	MPI_Reduce(&iDensity, &iTotalDensity, 1, MPI_DOUBLE, MPI_SUM, 0,
//			PETSC_COMM_WORLD);
//	MPI_Reduce(&vDensity, &vTotalDensity, 1, MPI_DOUBLE, MPI_SUM, 0,
//			PETSC_COMM_WORLD);
//	MPI_Reduce(&voidDensity, &voidTotalDensity, 1, MPI_DOUBLE, MPI_SUM, 0,
//			PETSC_COMM_WORLD);
//	MPI_Reduce(&perfectDensity, &perfectTotalDensity, 1, MPI_DOUBLE, MPI_SUM, 0,
//			PETSC_COMM_WORLD);
//	MPI_Reduce(&frankDensity, &frankTotalDensity, 1, MPI_DOUBLE, MPI_SUM, 0,
//			PETSC_COMM_WORLD);
//	MPI_Reduce(&faultedDensity, &faultedTotalDensity, 1, MPI_DOUBLE, MPI_SUM, 0,
//			PETSC_COMM_WORLD);
//	MPI_Reduce(&voidPartialDensity, &voidPartialTotalDensity, 1, MPI_DOUBLE,
//	MPI_SUM, 0, PETSC_COMM_WORLD);
//	MPI_Reduce(&perfectPartialDensity, &perfectPartialTotalDensity, 1,
//	MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
//	MPI_Reduce(&frankPartialDensity, &frankPartialTotalDensity, 1, MPI_DOUBLE,
//	MPI_SUM, 0, PETSC_COMM_WORLD);
//	MPI_Reduce(&faultedPartialDensity, &faultedPartialTotalDensity, 1,
//	MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
//	MPI_Reduce(&iDiameter, &iTotalDiameter, 1, MPI_DOUBLE, MPI_SUM, 0,
//			PETSC_COMM_WORLD);
//	MPI_Reduce(&vDiameter, &vTotalDiameter, 1, MPI_DOUBLE, MPI_SUM, 0,
//			PETSC_COMM_WORLD);
//	MPI_Reduce(&voidDiameter, &voidTotalDiameter, 1, MPI_DOUBLE, MPI_SUM, 0,
//			PETSC_COMM_WORLD);
//	MPI_Reduce(&perfectDiameter, &perfectTotalDiameter, 1, MPI_DOUBLE, MPI_SUM,
//			0, PETSC_COMM_WORLD);
//	MPI_Reduce(&frankDiameter, &frankTotalDiameter, 1, MPI_DOUBLE, MPI_SUM, 0,
//			PETSC_COMM_WORLD);
//	MPI_Reduce(&faultedDiameter, &faultedTotalDiameter, 1, MPI_DOUBLE, MPI_SUM,
//			0, PETSC_COMM_WORLD);
//	MPI_Reduce(&voidPartialDiameter, &voidPartialTotalDiameter, 1, MPI_DOUBLE,
//	MPI_SUM, 0, PETSC_COMM_WORLD);
//	MPI_Reduce(&perfectPartialDiameter, &perfectPartialTotalDiameter, 1,
//	MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
//	MPI_Reduce(&frankPartialDiameter, &frankPartialTotalDiameter, 1, MPI_DOUBLE,
//	MPI_SUM, 0, PETSC_COMM_WORLD);
//	MPI_Reduce(&faultedPartialDiameter, &faultedPartialTotalDiameter, 1,
//	MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
//
//	// Average the data
//	if (procId == 0) {
//		iTotalDensity = iTotalDensity / (grid[Mx] - grid[surfacePos + 1]);
//		vTotalDensity = vTotalDensity / (grid[Mx] - grid[surfacePos + 1]);
//		voidTotalDensity = voidTotalDensity / (grid[Mx] - grid[surfacePos + 1]);
//		perfectTotalDensity = perfectTotalDensity
//				/ (grid[Mx] - grid[surfacePos + 1]);
//		faultedTotalDensity = faultedTotalDensity
//				/ (grid[Mx] - grid[surfacePos + 1]);
//		frankTotalDensity = frankTotalDensity
//				/ (grid[Mx] - grid[surfacePos + 1]);
//		voidPartialTotalDensity = voidPartialTotalDensity
//				/ (grid[Mx] - grid[surfacePos + 1]);
//		perfectPartialTotalDensity = perfectPartialTotalDensity
//				/ (grid[Mx] - grid[surfacePos + 1]);
//		faultedPartialTotalDensity = faultedPartialTotalDensity
//				/ (grid[Mx] - grid[surfacePos + 1]);
//		frankPartialTotalDensity = frankPartialTotalDensity
//				/ (grid[Mx] - grid[surfacePos + 1]);
//		iTotalDiameter = iTotalDiameter
//				/ (iTotalDensity * (grid[Mx] - grid[surfacePos + 1]));
//		vTotalDiameter = vTotalDiameter
//				/ (vTotalDensity * (grid[Mx] - grid[surfacePos + 1]));
//		voidTotalDiameter = voidTotalDiameter
//				/ (voidTotalDensity * (grid[Mx] - grid[surfacePos + 1]));
//		perfectTotalDiameter = perfectTotalDiameter
//				/ (perfectTotalDensity * (grid[Mx] - grid[surfacePos + 1]));
//		faultedTotalDiameter = faultedTotalDiameter
//				/ (faultedTotalDensity * (grid[Mx] - grid[surfacePos + 1]));
//		frankTotalDiameter = frankTotalDiameter
//				/ (frankTotalDensity * (grid[Mx] - grid[surfacePos + 1]));
//		voidPartialTotalDiameter = voidPartialTotalDiameter
//				/ (voidPartialTotalDensity * (grid[Mx] - grid[surfacePos + 1]));
//		perfectPartialTotalDiameter = perfectPartialTotalDiameter
//				/ (perfectPartialTotalDensity
//						* (grid[Mx] - grid[surfacePos + 1]));
//		faultedPartialTotalDiameter = faultedPartialTotalDiameter
//				/ (faultedPartialTotalDensity
//						* (grid[Mx] - grid[surfacePos + 1]));
//		frankPartialTotalDiameter =
//				frankPartialTotalDiameter
//						/ (frankPartialTotalDensity
//								* (grid[Mx] - grid[surfacePos + 1]));
//
//		// Set the output precision
//		const int outputPrecision = 5;
//
//		// Open the output file
//		std::fstream outputFile;
//		outputFile.open("Alloy.dat", std::fstream::out | std::fstream::app);
//		outputFile << std::setprecision(outputPrecision);
//
//		// Output the data
//		outputFile << timestep << " " << time << " " << iTotalDensity << " "
//				<< iTotalDiameter << " " << vTotalDensity << " "
//				<< vTotalDiameter << " " << voidTotalDensity << " "
//				<< voidTotalDiameter << " " << faultedTotalDensity << " "
//				<< faultedTotalDiameter << " " << perfectTotalDensity << " "
//				<< perfectTotalDiameter << " " << frankTotalDensity << " "
//				<< frankTotalDiameter << " " << voidPartialTotalDensity << " "
//				<< voidPartialTotalDiameter << " " << faultedPartialTotalDensity
//				<< " " << faultedPartialTotalDiameter << " "
//				<< perfectPartialTotalDensity << " "
//				<< perfectPartialTotalDiameter << " "
//				<< frankPartialTotalDensity << " " << frankPartialTotalDiameter
//				<< std::endl;
//
//		// Close the output file
//		outputFile.close();
//	}
//
//	// Restore the PETSC solution array
//	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
//	CHKERRQ(ierr);

	PetscFunctionReturn(0);

}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "monitorScatter1D")
/**
 * This is a monitoring method that will save 1D plots of the xenon concentration
 * distribution at the middle of the grid.
 */
PetscErrorCode monitorScatter1D(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void*) {
	perf::ScopedTimer myTimer(scatterTimer);

	// Initial declarations
	PetscErrorCode ierr;
	double **solutionArray, *gridPointSolution;
	PetscInt xs, xm, xi, Mx;

	PetscFunctionBeginUser;

	// Don't do anything if it is not on the stride
	if (timestep % 200 != 0)
		PetscFunctionReturn(0);

	// Get the number of processes
	int worldSize;
	MPI_Comm_size(PETSC_COMM_WORLD, &worldSize);

	// Gets the process ID (important when it is running in parallel)
	int procId;
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the solutionArray
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	CHKERRQ(ierr);

	// Get the size of the total grid
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE);
	CHKERRQ(ierr);

	// Get the solver handler
	auto &solverHandler = PetscSolver::getSolverHandler();

	// Get the network and its size
	using NetworkType =
	core::network::NEReactionNetwork;
	using Spec = typename NetworkType::Species;
	using Region = typename NetworkType::Region;
	auto &network = dynamic_cast<NetworkType&>(solverHandler.getNetwork());
	int networkSize = network.getNumClusters();

	// Get the index of the middle of the grid
	PetscInt ix = Mx / 2;

	// If the middle is on this process
	if (ix >= xs && ix < xs + xm) {
		// Create a Point vector to store the data to give to the data provider
		// for the visualization
		auto myPoints = std::make_shared<std::vector<viz::dataprovider::Point> >();

		// Get the pointer to the beginning of the solution data for this grid point
		gridPointSolution = solutionArray[ix];

		for (int i = 0; i < networkSize; i++) {
			// Create a Point with the concentration[i] as the value
			// and add it to myPoints
			auto cluster = network.getCluster(i);
			const Region &clReg = cluster.getRegion();
			for (std::size_t j : makeIntervalRange(clReg[Spec::Xe])) {
                viz::dataprovider::Point aPoint;
				aPoint.value = gridPointSolution[i];
				aPoint.t = time;
				aPoint.x = (double) j;
				myPoints->push_back(aPoint);
			}
		}

		// Get the data provider and give it the points
		scatterPlot1D->getDataProvider()->setPoints(myPoints);

		// Change the title of the plot and the name of the data
		std::stringstream title;
		title << "Size Distribution";
		scatterPlot1D->getDataProvider()->setDataName(title.str());
		scatterPlot1D->plotLabelProvider->titleLabel = title.str();
		// Give the time to the label provider
		std::stringstream timeLabel;
		timeLabel << "time: " << std::setprecision(4) << time << "s";
		scatterPlot1D->plotLabelProvider->timeLabel = timeLabel.str();
		// Get the current time step
		PetscReal currentTimeStep;
		ierr = TSGetTimeStep(ts, &currentTimeStep);
		CHKERRQ(ierr);
		// Give the timestep to the label provider
		std::stringstream timeStepLabel;
		timeStepLabel << "dt: " << std::setprecision(4) << currentTimeStep
				<< "s";
		scatterPlot1D->plotLabelProvider->timeStepLabel = timeStepLabel.str();

		// Render and save in file
		std::stringstream fileName;
		fileName << "Scatter_TS" << timestep << ".png";
		scatterPlot1D->write(fileName.str());
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "monitorSeries1D")
/**
 * This is a monitoring method that will save 1D plots of many concentrations
 */
PetscErrorCode monitorSeries1D(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void*) {

	perf::ScopedTimer myTimer(seriesTimer);

	// Initial declarations
	PetscErrorCode ierr;
	const double **solutionArray, *gridPointSolution;
	PetscInt xs, xm, xi;
	double x = 0.0;

	PetscFunctionBeginUser;

	// Don't do anything if it is not on the stride
	if (timestep % 10 != 0)
		PetscFunctionReturn(0);

	// Get the number of processes
	int worldSize;
	MPI_Comm_size(PETSC_COMM_WORLD, &worldSize);

	// Gets the process ID (important when it is running in parallel)
	int procId;
	MPI_Comm_rank(PETSC_COMM_WORLD, &procId);

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the solutionArray
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	CHKERRQ(ierr);

	// Get the solver handler
	auto &solverHandler = PetscSolver::getSolverHandler();

	// Get the network and its size
	auto &network = solverHandler.getNetwork();
	const int networkSize = network.getNumClusters();

	// Get the physical grid
	auto grid = solverHandler.getXGrid();

	// To plot a maximum of 18 clusters of the whole benchmark
	const int loopSize = std::min(18, networkSize);

	if (procId == 0) {
		// Create a Point vector to store the data to give to the data provider
		// for the visualization
		std::vector<std::vector<viz::dataprovider::Point> > myPoints(loopSize);

		// Loop on the grid
		for (xi = xs; xi < xs + xm; xi++) {
			// Get the pointer to the beginning of the solution data for this grid point
			gridPointSolution = solutionArray[xi];

			for (int i = 0; i < loopSize; i++) {
				// Create a Point with the concentration[i] as the value
				// and add it to myPoints
                viz::dataprovider::Point aPoint;
				aPoint.value = gridPointSolution[i];
				aPoint.t = time;
				aPoint.x = (grid[xi] + grid[xi + 1]) / 2.0 - grid[1];
				myPoints[i].push_back(aPoint);
			}
		}

		// Loop on the other processes
		for (int i = 1; i < worldSize; i++) {
			// Get the size of the local grid of that process
			int localSize = 0;
			MPI_Recv(&localSize, 1, MPI_INT, i, 20, PETSC_COMM_WORLD,
					MPI_STATUS_IGNORE);

			// Loop on their grid
			for (int k = 0; k < localSize; k++) {
				// Get the position
				MPI_Recv(&x, 1, MPI_DOUBLE, i, 21, PETSC_COMM_WORLD,
						MPI_STATUS_IGNORE);

				for (int j = 0; j < loopSize; j++) {
					// and the concentrations
					double conc = 0.0;
					MPI_Recv(&conc, 1, MPI_DOUBLE, i, 22, PETSC_COMM_WORLD,
							MPI_STATUS_IGNORE);

					// Create a Point with the concentration[i] as the value
					// and add it to myPoints
                    viz::dataprovider::Point aPoint;
					aPoint.value = conc;
					aPoint.t = time;
					aPoint.x = x;
					myPoints[j].push_back(aPoint);
				}
			}
		}

		for (int i = 0; i < loopSize; i++) {
			// Get the data provider and give it the points
			auto thePoints = std::make_shared<std::vector<viz::dataprovider::Point> >(
					myPoints[i]);
			seriesPlot1D->getDataProvider(i)->setPoints(thePoints);
			// TODO: get the name or comp of the cluster
			seriesPlot1D->getDataProvider(i)->setDataName(std::to_string(i));
		}

		// Change the title of the plot
		std::stringstream title;
		title << "Concentrations";
		seriesPlot1D->plotLabelProvider->titleLabel = title.str();
		// Give the time to the label provider
		std::stringstream timeLabel;
		timeLabel << "time: " << std::setprecision(4) << time << "s";
		seriesPlot1D->plotLabelProvider->timeLabel = timeLabel.str();
		// Get the current time step
		PetscReal currentTimeStep;
		ierr = TSGetTimeStep(ts, &currentTimeStep);
		CHKERRQ(ierr);
		// Give the timestep to the label provider
		std::stringstream timeStepLabel;
		timeStepLabel << "dt: " << std::setprecision(4) << currentTimeStep
				<< "s";
		seriesPlot1D->plotLabelProvider->timeStepLabel = timeStepLabel.str();

		// Render and save in file
		std::stringstream fileName;
		fileName << "log_series_TS" << timestep << ".ppm";
		seriesPlot1D->write(fileName.str());
	}

	else {
		// Send the value of the local grid size to the master process
		MPI_Send(&xm, 1, MPI_DOUBLE, 0, 20, PETSC_COMM_WORLD);

		// Loop on the grid
		for (xi = xs; xi < xs + xm; xi++) {
			// Dump x
			x = (grid[xi] + grid[xi + 1]) / 2.0 - grid[1];

			// Get the pointer to the beginning of the solution data for this grid point
			gridPointSolution = solutionArray[xi];

			// Send the value of the local position to the master process
			MPI_Send(&x, 1, MPI_DOUBLE, 0, 21, PETSC_COMM_WORLD);

			for (int i = 0; i < loopSize; i++) {
				// Send the value of the concentrations to the master process
				MPI_Send(&gridPointSolution[i], 1, MPI_DOUBLE, 0, 22,
						PETSC_COMM_WORLD);
			}
		}
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "eventFunction1D")
/**
 * This is a method that checks if the surface should move or bursting happen
 */
PetscErrorCode eventFunction1D(TS ts, PetscReal time, Vec solution,
		PetscScalar *fvalue, void*) {
	perf::ScopedTimer myTimer(eventFuncTimer);

	// Initial declaration
	PetscErrorCode ierr;
	double **solutionArray, *gridPointSolution;
	PetscInt xs, xm, xi, Mx;
	depthPositions1D.clear();
	fvalue[0] = 1.0, fvalue[1] = 1.0, fvalue[2] = 1.0;

	PetscFunctionBeginUser;

	// Gets the process ID
	int procId;
	MPI_Comm_rank(PETSC_COMM_WORLD, &procId);

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the solutionArray
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	CHKERRQ(ierr);

	// Get the size of the total grid
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE);
	CHKERRQ(ierr);

	// Get the solver handler
	auto &solverHandler = PetscSolver::getSolverHandler();

	// Get the position of the surface
	int surfacePos = solverHandler.getSurfacePosition();
	xi = surfacePos + solverHandler.getLeftOffset();

	// Get the network
	auto &network = solverHandler.getNetwork();

	// Get the physical grid
	auto grid = solverHandler.getXGrid();

	// Get the flux handler to know the flux amplitude.
	auto fluxHandler = solverHandler.getFluxHandler();
	double heliumFluxAmplitude = fluxHandler->getFluxAmplitude();

	// Get the delta time from the previous timestep to this timestep
	double dt = time - previousTime;

	// Work of the moving surface first
	if (solverHandler.moveSurface()) {
		// Write the initial surface position
		if (procId == 0 && core::equal(time, 0.0)) {
			std::ofstream outputFile;
			outputFile.open("surface.txt", std::ios::app);
			outputFile << time << " " << grid[surfacePos + 1] - grid[1]
					<< std::endl;
			outputFile.close();
		}

		// Value to know on which processor is the location of the surface,
		// for MPI usage
		int surfaceProc = 0;

		// if xi is on this process
		if (xi >= xs && xi < xs + xm) {
			// Get the concentrations at xi = surfacePos + 1
			gridPointSolution = solutionArray[xi];

			// Compute the total density of intersitials that escaped from the
			// surface since last timestep using the stored flux
			nInterstitial1D += previousIFlux1D * dt;

			// Remove the sputtering yield since last timestep
			nInterstitial1D -= sputteringYield1D * heliumFluxAmplitude * dt;

			// Initialize the value for the flux
			double newFlux = 0.0;

			// Consider each interstitial cluster.
			for (int i = 0; i < iClusterIds1D.size(); i++) {
				auto currId = iClusterIds1D[i];
				// Get the cluster
				auto cluster = network.getClusterCommon(currId);
				// Get its concentration
				double conc = gridPointSolution[currId];
				// Get its size and diffusion coefficient
				int size = i + 1;
				double coef = cluster.getDiffusionCoefficient(xi - xs);

				// Factor for finite difference
				double hxLeft = 0.0, hxRight = 0.0;
				if (xi - 1 >= 0 && xi < Mx) {
					hxLeft = (grid[xi + 1] - grid[xi - 1]) / 2.0;
					hxRight = (grid[xi + 2] - grid[xi]) / 2.0;
				} else if (xi - 1 < 0) {
					hxLeft = (grid[xi + 1] + grid[xi]) / 2.0;
					hxRight = (grid[xi + 2] - grid[xi]) / 2.0;
				} else {
					hxLeft = (grid[xi + 1] - grid[xi - 1]) / 2.0;
					hxRight = (grid[xi + 1] - grid[xi]) / 2;
				}
				double factor = 2.0 / (hxLeft * (hxLeft + hxRight));
				// Compute the flux going to the left
				newFlux += (double) size * factor * coef * conc * hxLeft;
			}

			// Update the previous flux
			previousIFlux1D = newFlux;

			// Set the surface processor
			surfaceProc = procId;
		}

		// Get which processor will send the information
		int surfaceId = 0;
		MPI_Allreduce(&surfaceProc, &surfaceId, 1, MPI_INT, MPI_SUM,
				PETSC_COMM_WORLD);

		// Send the information about nInterstitial1D and previousFlux1D
		// to the other processes
		MPI_Bcast(&nInterstitial1D, 1, MPI_DOUBLE, surfaceId, PETSC_COMM_WORLD);
		MPI_Bcast(&previousIFlux1D, 1, MPI_DOUBLE, surfaceId, PETSC_COMM_WORLD);

		// Now that all the processes have the same value of nInterstitials, compare
		// it to the threshold to now if we should move the surface

		// Get the initial vacancy concentration
		double initialVConc = solverHandler.getInitialVConc();

		// The density of tungsten is 62.8 atoms/nm3, thus the threshold is
		double threshold = (62.8 - initialVConc) * (grid[xi] - grid[xi - 1]);
		if (nInterstitial1D > threshold) {
			// The surface is moving
			fvalue[0] = 0.0;
		}

		// Moving the surface back
		else if (nInterstitial1D < -threshold / 10.0) {
			// The surface is moving
			fvalue[1] = 0.0;
		}
	}

	// Now work on the bubble bursting
	if (solverHandler.burstBubbles()) {
		using NetworkType =
		core::network::PSIReactionNetwork<core::network::PSIFullSpeciesList>;
		using Spec = typename NetworkType::Species;
		auto psiNetwork = dynamic_cast<NetworkType*>(&network);
		auto dof = network.getDOF();

		// Compute the prefactor for the probability (arbitrary)
		double prefactor = heliumFluxAmplitude * dt * 0.1;

		// The depth parameter to know where the bursting should happen
		double depthParam = solverHandler.getTauBursting();			// nm

		// For now we are not bursting
		bool burst = false;

		// Loop on the full grid of interest
		for (xi = surfacePos + solverHandler.getLeftOffset();
				xi < Mx - solverHandler.getRightOffset(); xi++) {

			// If this is the locally owned part of the grid
			if (xi >= xs && xi < xs + xm) {

				// Get the pointer to the beginning of the solution data for this grid point
				gridPointSolution = solutionArray[xi];

				using HostUnmanaged =
				Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
				auto hConcs = HostUnmanaged(gridPointSolution, dof);
				auto dConcs = Kokkos::View<double*>("Concentrations", dof);
				deep_copy(dConcs, hConcs);

				// Get the distance from the surface
				double distance = (grid[xi] + grid[xi + 1]) / 2.0
						- grid[surfacePos + 1];

				// Compute the helium density at this grid point
				double heDensity = psiNetwork->getTotalAtomConcentration(dConcs,
						Spec::He, 0);

				// Compute the radius of the bubble from the number of helium
				double nV = heDensity * (grid[xi + 1] - grid[xi]) / 4.0;
//			double nV = pow(heDensity / 5.0, 1.163) * (grid[xi + 1] - grid[xi]);
				double latticeParam = network.getLatticeParameter();
				double tlcCubed = latticeParam * latticeParam * latticeParam;
				double radius = (sqrt(3.0) / 4) * latticeParam
						+ cbrt((3.0 * tlcCubed * nV) / (8.0 * core::pi))
						- cbrt((3.0 * tlcCubed) / (8.0 * core::pi));

				// If the radius is larger than the distance to the surface, burst
				if (radius > distance) {
					burst = true;
					depthPositions1D.push_back(xi);
					// Exit the loop
					continue;
				}
				// Add randomness
				double prob = prefactor * (1.0 - (distance - radius) / distance)
						* std::min(1.0,
								exp(
										-(distance - depthParam)
												/ (depthParam * 2.0)));
				double test = solverHandler.getRNG().GetRandomDouble();

				if (prob > test) {
					burst = true;
					depthPositions1D.push_back(xi);
				}
			}
		}

		// If at least one grid point is bursting
		if (burst) {
			// The event is happening
			fvalue[2] = 0.0;
		}
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "postEventFunction1D")
/**
 * This is a method that moves the surface or burst bubbles
 */
PetscErrorCode postEventFunction1D(TS ts, PetscInt nevents,
		PetscInt eventList[], PetscReal time, Vec solution, PetscBool, void*) {

	perf::ScopedTimer myTimer(postEventFuncTimer);

	// Initial declaration
	PetscErrorCode ierr;
	double **solutionArray, *gridPointSolution;
	PetscInt xs, xm, xi;

	PetscFunctionBeginUser;

	// Call monitor time hear because it is skipped when post event is used
	ierr = computeFluence(ts, 0, time, solution, NULL);
	CHKERRQ(ierr);
	ierr = monitorTime(ts, 0, time, solution, NULL);
	CHKERRQ(ierr);

	// Check if the surface has moved
	if (nevents == 0) {
		PetscFunctionReturn(0);
	}

	// Check if both events happened
	if (nevents == 3)
		throw std::string(
				"\nxolotlSolver::Monitor1D: This is not supposed to happen, the surface cannot "
						"move in both directions at the same time!!");

	// Gets the process ID
	int procId;
	MPI_Comm_rank(PETSC_COMM_WORLD, &procId);

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the solutionArray
	ierr = DMDAVecGetArrayDOF(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	CHKERRQ(ierr);

	// Get the solver handler
	auto &solverHandler = PetscSolver::getSolverHandler();

	// Get the position of the surface
	int surfacePos = solverHandler.getSurfacePosition();

	// Get the network
	auto &network = solverHandler.getNetwork();
	int dof = network.getDOF();

	// Get the physical grid
	auto grid = solverHandler.getXGrid();

	// Take care of bursting
	using NetworkType =
	core::network::PSIReactionNetwork<core::network::PSIFullSpeciesList>;
	using Spec = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;
	auto psiNetwork = dynamic_cast<NetworkType*>(&network);

	// Loop on each bursting depth
	for (int i = 0; i < depthPositions1D.size(); i++) {
		// Get the pointer to the beginning of the solution data for this grid point
		gridPointSolution = solutionArray[depthPositions1D[i]];

		// Get the distance from the surface
		double distance = (grid[depthPositions1D[i]]
				+ grid[depthPositions1D[i] + 1]) / 2.0 - grid[surfacePos + 1];

		// Write the bursting information
		std::ofstream outputFile;
		outputFile.open("bursting.txt", std::ios::app);
		outputFile << time << " " << distance << std::endl;
		outputFile.close();

		// Pinhole case
		// Loop on every cluster
		for (unsigned int i = 0; i < network.getNumClusters(); i++) {
			const auto &clReg =
					psiNetwork->getCluster(i, plsm::onHost).getRegion();
			// Non-grouped clusters
			if (clReg.isSimplex()) {
				// Get the composition
				Composition comp = clReg.getOrigin();
				// Pure He, D, or T case
				if (comp.isOnAxis(Spec::He) || comp.isOnAxis(Spec::D)
						|| comp.isOnAxis(Spec::T)) {
					// Reset concentration
					gridPointSolution[i] = 0.0;
				}
				// Mixed cluster case
				else if (!comp.isOnAxis(Spec::V) && !comp.isOnAxis(Spec::I)) {
					// Transfer concentration to V of the same size
					Composition vComp = Composition::zero();
					vComp[Spec::V] = comp[Spec::V];
					auto vCluster = psiNetwork->findCluster(vComp,
							plsm::onHost);
					gridPointSolution[vCluster.getId()] += gridPointSolution[i];
					gridPointSolution[i] = 0.0;
				}
			}
			// Grouped clusters
			else {
				// Get the factor
				double concFactor = clReg.volume() / clReg[Spec::V].length();
				// Loop on the Vs
				for (auto j : makeIntervalRange(clReg[Spec::V])) {
					// Transfer concentration to V of the same size
					Composition vComp = Composition::zero();
					vComp[Spec::V] = j;
					auto vCluster = psiNetwork->findCluster(vComp,
							plsm::onHost);
					// TODO: refine formula with V moment
					gridPointSolution[vCluster.getId()] += gridPointSolution[i]
							* concFactor;
				}

				// Reset the concentration and moments
				gridPointSolution[i] = 0.0;
				auto momentIds =
						psiNetwork->getCluster(i, plsm::onHost).getMomentIds();
				for (std::size_t j = 0; j < momentIds.extent(0); j++) {
					gridPointSolution[momentIds(j)] = 0.0;
				}
			}
		}
	}

	// Now takes care of moving surface
	bool moving = false;
	bool movingUp = false;
	for (int i = 0; i < nevents; i++) {
		if (eventList[i] < 2)
			moving = true;
		if (eventList[i] == 0)
			movingUp = true;
	}

	// Skip if nothing is moving
	if (!moving) {
		// Restore the solutionArray
		ierr = DMDAVecRestoreArrayDOF(da, solution, &solutionArray);
		CHKERRQ(ierr);

		PetscFunctionReturn(0);
	}

	// Set the surface position
	xi = surfacePos + solverHandler.getLeftOffset();

	// Get the initial vacancy concentration
	double initialVConc = solverHandler.getInitialVConc();

	// The density of tungsten is 62.8 atoms/nm3, thus the threshold is
	double threshold = (62.8 - initialVConc) * (grid[xi] - grid[xi - 1]);

	if (movingUp) {
		int nGridPoints = 0;
		// Move the surface up until it is smaller than the next threshold
		while (nInterstitial1D > threshold) {
			// Move the surface higher
			surfacePos--;
			xi = surfacePos + solverHandler.getLeftOffset();
			nGridPoints++;
			// Update the number of interstitials
			nInterstitial1D -= threshold;
			// Update the thresold
			threshold = (62.8 - initialVConc) * (grid[xi] - grid[xi - 1]);
		}

		// Throw an exception if the position is negative
		if (surfacePos < 0) {
			PetscBool flagCheck;
			ierr = PetscOptionsHasName(NULL, NULL, "-check_collapse",
					&flagCheck);
			CHKERRQ(ierr);
			if (flagCheck) {
				// Write the convergence reason
				std::ofstream outputFile;
				outputFile.open("solverStatus.txt");
				outputFile << "overgrid" << std::endl;
				outputFile.close();
			}
			throw std::string(
					"\nxolotlSolver::Monitor1D: The surface is trying to go outside of the grid!!");
		}

		// Printing information about the extension of the material
		if (procId == 0) {
			std::cout << "Adding " << nGridPoints
					<< " points to the grid at time: " << time << " s."
					<< std::endl;
		}

		// Set it in the solver
		solverHandler.setSurfacePosition(surfacePos);

		// Initialize the vacancy concentration and the temperature on the new grid points
		// Get the single vacancy ID
		auto singleVacancyCluster = network.getSingleVacancy();
		auto vacancyIndex = core::network::IReactionNetwork::invalidIndex();
		if (singleVacancyCluster.getId()
				!= core::network::IReactionNetwork::invalidIndex())
			vacancyIndex = singleVacancyCluster.getId();
		// Get the surface temperature
		double temp = 0.0;
		if (xi >= xs && xi < xs + xm) {
			temp = solutionArray[xi][dof];
		}
		double surfTemp = 0.0;
		MPI_Allreduce(&temp, &surfTemp, 1, MPI_DOUBLE, MPI_SUM,
				PETSC_COMM_WORLD);

		// Loop on the new grid points
		while (nGridPoints >= 0) {
			// Position of the newly created grid point
			xi = surfacePos + nGridPoints;

			// If xi is on this process
			if (xi >= xs && xi < xs + xm) {
				// Get the concentrations
				gridPointSolution = solutionArray[xi];

				// Set the new surface temperature
				gridPointSolution[dof] = surfTemp;

				if (vacancyIndex
						!= core::network::IReactionNetwork::invalidIndex()
						&& nGridPoints > 0) {
					// Initialize the vacancy concentration
					gridPointSolution[vacancyIndex] = initialVConc;
				}
			}

			// Decrease the number of grid points
			--nGridPoints;
		}
	}

	// Moving the surface back
	else {
		// Move it back as long as the number of interstitials in negative
		while (nInterstitial1D < 0.0) {
			// Compute the threshold to a deeper grid point
			threshold = (62.8 - initialVConc) * (grid[xi + 1] - grid[xi]);
			// Set all the concentrations to 0.0 at xi = surfacePos + 1
			// if xi is on this process
			if (xi >= xs && xi < xs + xm) {
				// Get the concentrations at xi = surfacePos + 1
				gridPointSolution = solutionArray[xi];
				// Loop on DOF
				for (int i = 0; i < dof; i++) {
					gridPointSolution[i] = 0.0;
				}
			}

			// Move the surface deeper
			surfacePos++;
			xi = surfacePos + solverHandler.getLeftOffset();
			// Update the number of interstitials
			nInterstitial1D += threshold;
		}

		// Printing information about the extension of the material
		if (procId == 0) {
			std::cout << "Removing grid points to the grid at time: " << time
					<< " s." << std::endl;
		}

		// Set it in the solver
		solverHandler.setSurfacePosition(surfacePos);
	}

	// Set the new surface location in the surface advection handler
	auto advecHandler = solverHandler.getAdvectionHandler();
	advecHandler->setLocation(grid[surfacePos + 1] - grid[1]);

	// Set the new surface in the temperature handler
	auto tempHandler = solverHandler.getTemperatureHandler();
	tempHandler->updateSurfacePosition(surfacePos);

	// Get the flux handler to reinitialize it
	auto fluxHandler = solverHandler.getFluxHandler();
	fluxHandler->initializeFluxHandler(solverHandler.getNetwork(), surfacePos,
			grid);

	// Get the modified trap-mutation handler to reinitialize it
	auto mutationHandler = solverHandler.getMutationHandler();
	auto advecHandlers = solverHandler.getAdvectionHandlers();
	mutationHandler->initializeIndex1D(surfacePos, solverHandler.getNetwork(),
			advecHandlers, grid, xm, xs);

	// Write the updated surface position
	if (procId == 0) {
		std::ofstream outputFile;
		outputFile.open("surface.txt", std::ios::app);
		outputFile << time << " " << grid[surfacePos + 1] - grid[1]
				<< std::endl;
		outputFile.close();
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOF(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

/**
 * This operation sets up different monitors
 *  depending on the options.
 * @param ts The time stepper
 * @return A standard PETSc error code
 */
PetscErrorCode setupPetsc1DMonitor(TS ts,
		std::shared_ptr<perf::IHandlerRegistry> handlerRegistry) {

	PetscErrorCode ierr;

	// Initialize the timers, including the one for this function.
	initTimer = handlerRegistry->getTimer("monitor1D:init");
	perf::ScopedTimer myTimer(initTimer);
	checkNegativeTimer = handlerRegistry->getTimer("monitor1D:checkNeg");
	tridynTimer = handlerRegistry->getTimer("monitor1D:tridyn");
	startStopTimer = handlerRegistry->getTimer("monitor1D:startStop");
	heRetentionTimer = handlerRegistry->getTimer("monitor1D:heRet");
	xeRetentionTimer = handlerRegistry->getTimer("monitor1D:xeRet");
	scatterTimer = handlerRegistry->getTimer("monitor1D:scatter");
	seriesTimer = handlerRegistry->getTimer("monitor1D:series");
	eventFuncTimer = handlerRegistry->getTimer("monitor1D:event");
	postEventFuncTimer = handlerRegistry->getTimer("monitor1D:postEvent");

	// Get the process ID
	int procId;
	MPI_Comm_rank(PETSC_COMM_WORLD, &procId);

	// Get xolotlViz handler registry
	auto vizHandlerRegistry = factory::viz::getVizHandlerRegistry();

	// Flags to launch the monitors or not
	PetscBool flagNeg, flagCollapse, flag2DPlot, flag1DPlot, flagSeries,
			flagPerf, flagHeDesorption, flagHeRetention, flagStatus,
			flagXeRetention, flagTRIDYN, flagAlloy, flagTemp;

	// Check the option -check_negative
	ierr = PetscOptionsHasName(NULL, NULL, "-check_negative", &flagNeg);
	checkPetscError(ierr,
			"setupPetsc1DMonitor: PetscOptionsHasName (-check_negative) failed.");

	// Check the option -check_collapse
	ierr = PetscOptionsHasName(NULL, NULL, "-check_collapse", &flagCollapse);
	checkPetscError(ierr,
			"setupPetsc1DMonitor: PetscOptionsHasName (-check_collapse) failed.");

	// Check the option -plot_perf
	ierr = PetscOptionsHasName(NULL, NULL, "-plot_perf", &flagPerf);
	checkPetscError(ierr,
			"setupPetsc1DMonitor: PetscOptionsHasName (-plot_perf) failed.");

	// Check the option -plot_series
	ierr = PetscOptionsHasName(NULL, NULL, "-plot_series", &flagSeries);
	checkPetscError(ierr,
			"setupPetsc1DMonitor: PetscOptionsHasName (-plot_series) failed.");

	// Check the option -plot_1d
	ierr = PetscOptionsHasName(NULL, NULL, "-plot_1d", &flag1DPlot);
	checkPetscError(ierr,
			"setupPetsc1DMonitor: PetscOptionsHasName (-plot_1d) failed.");

	// Check the option -helium_desorption
	ierr = PetscOptionsHasName(NULL, NULL, "-helium_desorption",
			&flagHeDesorption);
	checkPetscError(ierr,
			"setupPetsc1DMonitor: PetscOptionsHasName (-helium_desorption) failed.");

	// Check the option -helium_retention
	ierr = PetscOptionsHasName(NULL, NULL, "-helium_retention",
			&flagHeRetention);
	checkPetscError(ierr,
			"setupPetsc1DMonitor: PetscOptionsHasName (-helium_retention) failed.");

	// Check the option -xenon_retention
	ierr = PetscOptionsHasName(NULL, NULL, "-xenon_retention",
			&flagXeRetention);
	checkPetscError(ierr,
			"setupPetsc1DMonitor: PetscOptionsHasName (-xenon_retention) failed.");

	// Check the option -start_stop
	ierr = PetscOptionsHasName(NULL, NULL, "-start_stop", &flagStatus);
	checkPetscError(ierr,
			"setupPetsc1DMonitor: PetscOptionsHasName (-start_stop) failed.");

	// Check the option -tridyn
	ierr = PetscOptionsHasName(NULL, NULL, "-tridyn", &flagTRIDYN);
	checkPetscError(ierr,
			"setupPetsc1DMonitor: PetscOptionsHasName (-tridyn) failed.");

	// Check the option -alloy
	ierr = PetscOptionsHasName(NULL, NULL, "-alloy", &flagAlloy);
	checkPetscError(ierr,
			"setupPetsc1DMonitor: PetscOptionsHasName (-alloy) failed.");

	// Check the option -temp_profile
	ierr = PetscOptionsHasName(NULL, NULL, "-temp_profile", &flagTemp);
	checkPetscError(ierr,
			"setupPetsc1DMonitor: PetscOptionsHasName (-temp_profile) failed.");

	// Get the solver handler
	auto &solverHandler = PetscSolver::getSolverHandler();

	// Get the network and its size
	auto &network = solverHandler.getNetwork();
	const int networkSize = network.getNumClusters();

	// Determine if we have an existing restart file,
	// and if so, it it has had timesteps written to it.
	std::unique_ptr<io::XFile> networkFile;
	std::unique_ptr<io::XFile::TimestepGroup> lastTsGroup;
	std::string networkName = solverHandler.getNetworkName();
	bool hasConcentrations = false;
	if (not networkName.empty()) {
		networkFile.reset(new io::XFile(networkName));
		auto concGroup = networkFile->getGroup<
				io::XFile::ConcentrationGroup>();
		hasConcentrations = (concGroup and concGroup->hasTimesteps());
		if (hasConcentrations) {
			lastTsGroup = concGroup->getLastTimestepGroup();
		}
	}

	// Set the post step processing to stop the solver if the time step collapses
	if (flagCollapse) {
		// Find the threshold
		PetscBool flag;
		ierr = PetscOptionsGetReal(NULL, NULL, "-check_collapse",
				&timeStepThreshold, &flag);
		checkPetscError(ierr,
				"setupPetsc1DMonitor: PetscOptionsGetInt (-check_collapse) failed.");
		if (!flag)
			timeStepThreshold = 1.0e-16;

		// Set the post step process that tells the solver when to stop if the time step collapse
		ierr = TSSetPostStep(ts, checkTimeStep);
		checkPetscError(ierr,
				"setupPetsc1DMonitor: TSSetPostStep (checkTimeStep) failed.");
	}

	// Set the monitor to check the negative concentrations
	if (flagNeg) {
		// Find the stride to know how often we want to check
		PetscBool flag;
		ierr = PetscOptionsGetReal(NULL, NULL, "-check_negative",
				&negThreshold1D, &flag);
		checkPetscError(ierr,
				"setupPetsc1DMonitor: PetscOptionsGetInt (-check_negative) failed.");
		if (!flag)
			negThreshold1D = 1.0e-30;

		// checkNegative1D will be called at each timestep
		ierr = TSMonitorSet(ts, checkNegative1D, NULL, NULL);
		checkPetscError(ierr,
				"setupPetsc1DMonitor: TSMonitorSet (checkNegative1D) failed.");
	}

	// Set the monitor to save the status of the simulation in hdf5 file
	if (flagStatus) {
		// Find the stride to know how often the HDF5 file has to be written
		PetscBool flag;
		ierr = PetscOptionsGetReal(NULL, NULL, "-start_stop", &hdf5Stride1D,
				&flag);
		checkPetscError(ierr,
				"setupPetsc1DMonitor: PetscOptionsGetInt (-start_stop) failed.");
		if (!flag)
			hdf5Stride1D = 1.0;

		// Compute the correct hdf5Previous1D for a restart
		// Get the last time step written in the HDF5 file
		if (hasConcentrations) {

			assert(lastTsGroup);

			// Get the previous time from the HDF5 file
			previousTime = lastTsGroup->readPreviousTime();
			hdf5Previous1D = (int) (previousTime / hdf5Stride1D);
		}

		// Don't do anything if both files have the same name
		if (hdf5OutputName1D != solverHandler.getNetworkName()) {

			PetscInt Mx;
			PetscErrorCode ierr;

			// Get the da from ts
			DM da;
			ierr = TSGetDM(ts, &da);
			checkPetscError(ierr, "setupPetsc1DMonitor: TSGetDM failed.");

			// Get the size of the total grid
			ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE);
			checkPetscError(ierr, "setupPetsc1DMonitor: DMDAGetInfo failed.");

			// Get the solver handler
			auto &solverHandler = PetscSolver::getSolverHandler();

			// Get the physical grid
			auto grid = solverHandler.getXGrid();

			// Create and initialize a checkpoint file.
			// We do this in its own scope so that the file
			// is closed when the file object goes out of scope.
			// We want it to close before we (potentially) copy
			// the network from another file using a single-process
			// MPI communicator.
			{
				io::XFile checkpointFile(hdf5OutputName1D, grid,
						PETSC_COMM_WORLD);
			}

			// Copy the network group from the given file (if it has one).
			// We open the files using a single-process MPI communicator
			// because it is faster for a single process to do the
			// copy with HDF5's H5Ocopy implementation than it is
			// when all processes call the copy function.
			// The checkpoint file must be closed before doing this.
			writeNetwork(PETSC_COMM_WORLD, solverHandler.getNetworkName(),
					hdf5OutputName1D, network);
		}

		// startStop1D will be called at each timestep
		ierr = TSMonitorSet(ts, startStop1D, NULL, NULL);
		checkPetscError(ierr,
				"setupPetsc1DMonitor: TSMonitorSet (startStop1D) failed.");
	}

// If the user wants the surface to be able to move or bursting
	if (solverHandler.moveSurface() || solverHandler.burstBubbles()) {
		// Surface
		if (solverHandler.moveSurface()) {

			using NetworkType =
			core::network::PSIReactionNetwork<core::network::PSIFullSpeciesList>;
			auto psiNetwork = dynamic_cast<NetworkType*>(&network);

			// Initialize the composition
			NetworkType::Composition comp = NetworkType::Composition::zero();

			// Loop on interstital clusters
			bool iClusterExists = true;
			int iSize = 1;
			while (iClusterExists) {
				comp[NetworkType::Species::I] = iSize;
				auto cluster = psiNetwork->findCluster(comp, plsm::onHost);
				// Check that the helium cluster is present in the network
				if (cluster.getId() != NetworkType::invalidIndex()) {
					iClusterIds1D.push_back(cluster.getId());
					iSize++;
				} else
					iClusterExists = false;
			}

			// Get the interstitial information at the surface if concentrations were stored
			if (hasConcentrations) {

				assert(lastTsGroup);

				// Get the interstitial quantity from the HDF5 file
				nInterstitial1D = lastTsGroup->readData1D("nInterstitial");
				// Get the previous I flux from the HDF5 file
				previousIFlux1D = lastTsGroup->readData1D("previousIFlux");
				// Get the previous time from the HDF5 file
				previousTime = lastTsGroup->readPreviousTime();
			}

			// Get the sputtering yield
			sputteringYield1D = solverHandler.getSputteringYield();

			// Master process
			if (procId == 0) {
				// Clear the file where the surface will be written
				std::ofstream outputFile;
				outputFile.open("surface.txt");
				outputFile.close();
			}
		}

		// Bursting
		if (solverHandler.burstBubbles()) {
			// No need to seed the random number generator here.
			// The solver handler has already done it.
		}

		// Set directions and terminate flags for the surface event
		PetscInt direction[3];
		PetscBool terminate[3];
		direction[0] = 0, direction[1] = 0, direction[2] = 0;
		terminate[0] = PETSC_FALSE, terminate[1] = PETSC_FALSE, terminate[2] =
				PETSC_FALSE;
		// Set the TSEvent
		ierr = TSSetEventHandler(ts, 3, direction, terminate, eventFunction1D,
				postEventFunction1D, NULL);
		checkPetscError(ierr,
				"setupPetsc1DMonitor: TSSetEventHandler (eventFunction1D) failed.");

		// Master process
		if (procId == 0) {
			// Uncomment to clear the file where the bursting info will be written
			std::ofstream outputFile;
			outputFile.open("bursting.txt");
			outputFile.close();
		}
	}

// Set the monitor to save 1D plot of xenon distribution
	if (flag1DPlot) {
		// Only the master process will create the plot
		if (procId == 0) {
			// Create a ScatterPlot
			scatterPlot1D = vizHandlerRegistry->getPlot("scatterPlot1D",
					viz::PlotType::SCATTER);

			scatterPlot1D->setLogScale();

			// Create and set the label provider
			auto labelProvider = std::make_shared<viz::LabelProvider>(
					"labelProvider");
			labelProvider->axis1Label = "Xenon Size";
			labelProvider->axis2Label = "Concentration";

			// Give it to the plot
			scatterPlot1D->setLabelProvider(labelProvider);

			// Create the data provider
			auto dataProvider = std::make_shared<viz::dataprovider::CvsXDataProvider>(
					"dataProvider");

			// Give it to the plot
			scatterPlot1D->setDataProvider(dataProvider);
		}

		// monitorScatter1D will be called at each timestep
		ierr = TSMonitorSet(ts, monitorScatter1D, NULL, NULL);
		checkPetscError(ierr,
				"setupPetsc1DMonitor: TSMonitorSet (monitorScatter1D) failed.");
	}

// Set the monitor to save 1D plot of many concentrations
	if (flagSeries) {
		// Only the master process will create the plot
		if (procId == 0) {
			// Create a ScatterPlot
			seriesPlot1D = vizHandlerRegistry->getPlot("seriesPlot1D",
					viz::PlotType::SERIES);

			// set the log scale
//			seriesPlot1D->setLogScale();

			// Create and set the label provider
			auto labelProvider = std::make_shared<viz::LabelProvider>(
					"labelProvider");
			labelProvider->axis1Label = "x Position on the Grid";
			labelProvider->axis2Label = "Concentration";

			// Give it to the plot
			seriesPlot1D->setLabelProvider(labelProvider);

			// To plot a maximum of 18 clusters of the whole benchmark
			const int loopSize = std::min(18, networkSize);

			// Create a data provider for each cluster in the network
			for (int i = 0; i < loopSize; i++) {
				// Set the name for Identifiable
				std::stringstream dataProviderName;
				dataProviderName << "dataprovider" << i;
				// Create the data provider
				auto dataProvider =
						std::make_shared<viz::dataprovider::CvsXDataProvider>(
								dataProviderName.str());

				// Give it to the plot
				seriesPlot1D->addDataProvider(dataProvider);
			}
		}

		// monitorSeries1D will be called at each timestep
		ierr = TSMonitorSet(ts, monitorSeries1D, NULL, NULL);
		checkPetscError(ierr,
				"setupPetsc1DMonitor: TSMonitorSet (monitorSeries1D) failed.");
	}

// Set the monitor to save performance plots (has to be in parallel)
	if (flagPerf) {
		// Only the master process will create the plot
		if (procId == 0) {
			// Create a ScatterPlot
			perfPlot = vizHandlerRegistry->getPlot("perfPlot",
					viz::PlotType::SCATTER);

			// Create and set the label provider
			auto labelProvider = std::make_shared<viz::LabelProvider>(
					"labelProvider");
			labelProvider->axis1Label = "Process ID";
			labelProvider->axis2Label = "Solver Time";

			// Give it to the plot
			perfPlot->setLabelProvider(labelProvider);

			// Create the data provider
			auto dataProvider = std::make_shared<viz::dataprovider::CvsXDataProvider>(
					"dataProvider");

			// Give it to the plot
			perfPlot->setDataProvider(dataProvider);
		}

		// monitorPerf will be called at each timestep
		ierr = TSMonitorSet(ts, monitorPerf, NULL, NULL);
		checkPetscError(ierr,
				"setupPetsc1DMonitor: TSMonitorSet (monitorPerf) failed.");
	}

	// Set the monitor to compute the helium desorption
	if (flagHeDesorption) {
		// computeHeliumDesorption1D will be called at each timestep
		ierr = TSMonitorSet(ts, computeHeliumDesorption1D, NULL, NULL);
		checkPetscError(ierr,
				"setupPetsc1DMonitor: TSMonitorSet (computeHeliumDesorption1D) failed.");

		// Master process
		if (procId == 0) {
			// Uncomment to clear the file where the desorption
			std::ofstream outputFile;
			outputFile.open("thds.txt");
			outputFile.close();
		}
	}

// Set the monitor to compute the helium fluence and the retention
// for the retention calculation
	if (flagHeRetention) {

		// Get the previous time if concentrations were stored and initialize the fluence
		if (hasConcentrations) {

			assert(lastTsGroup);

			// Get the previous time from the HDF5 file
			double time = lastTsGroup->readPreviousTime();
			// Initialize the fluence
			auto fluxHandler = solverHandler.getFluxHandler();
			// The length of the time step
			double dt = time;
			// Increment the fluence with the value at this current timestep
			fluxHandler->incrementFluence(dt);
			// Get the previous time from the HDF5 file
			// TODO isn't this the same as 'time' above?
			previousTime = lastTsGroup->readPreviousTime();

			// If the bottom is a free surface
			if (solverHandler.getRightOffset() == 1) {
				// Read about the impurity fluxes in the bulk
				nHelium1D = lastTsGroup->readData1D("nHelium");
				previousHeFlux1D = lastTsGroup->readData1D("previousHeFlux");
				nDeuterium1D = lastTsGroup->readData1D("nDeuterium");
				previousDFlux1D = lastTsGroup->readData1D("previousDFlux");
				nTritium1D = lastTsGroup->readData1D("nTritium");
				previousTFlux1D = lastTsGroup->readData1D("previousTFlux");
				nVacancy1D = lastTsGroup->readData1D("nVacancy");
				previousVFlux1D = lastTsGroup->readData1D("previousVFlux");
				nIBulk1D = lastTsGroup->readData1D("nIBulk");
				previousIBulkFlux1D = lastTsGroup->readData1D(
						"previousIBulkFlux");
			}
		}

		// computeFluence will be called at each timestep
		ierr = TSMonitorSet(ts, computeFluence, NULL, NULL);
		checkPetscError(ierr,
				"setupPetsc1DMonitor: TSMonitorSet (computeFluence) failed.");

		// computeHeliumRetention1D will be called at each timestep
		ierr = TSMonitorSet(ts, computeHeliumRetention1D, NULL, NULL);
		checkPetscError(ierr,
				"setupPetsc1DMonitor: TSMonitorSet (computeHeliumRetention1D) failed.");

		// Master process
		if (procId == 0) {
			// Uncomment to clear the file where the retention will be written
			std::ofstream outputFile;
			outputFile.open("retentionOut.txt");
			outputFile.close();
		}
	}

// Set the monitor to compute the xenon fluence and the retention
// for the retention calculation
	if (flagXeRetention) {
		// Get the previous time if concentrations were stored and initialize the fluence
		if (hasConcentrations) {

			assert(lastTsGroup);

			// Get the previous time from the HDF5 file
			double time = lastTsGroup->readPreviousTime();
			// Initialize the fluence
			auto fluxHandler = solverHandler.getFluxHandler();
			// The length of the time step
			double dt = time;
			// Increment the fluence with the value at this current timestep
			fluxHandler->incrementFluence(dt);
			// Get the previous time from the HDF5 file
			// TODO isn't this the same as 'time' above?
			previousTime = lastTsGroup->readPreviousTime();
		}

		// computeFluence will be called at each timestep
		ierr = TSMonitorSet(ts, computeFluence, NULL, NULL);
		checkPetscError(ierr,
				"setupPetsc1DMonitor: TSMonitorSet (computeFluence) failed.");

		// computeXenonRetention1D will be called at each timestep
		ierr = TSMonitorSet(ts, computeXenonRetention1D, NULL, NULL);
		checkPetscError(ierr,
				"setupPetsc1DMonitor: TSMonitorSet (computeXenonRetention1D) failed.");

		// Master process
		if (procId == 0) {
			// Uncomment to clear the file where the retention will be written
			std::ofstream outputFile;
			outputFile.open("retentionOut.txt");
			outputFile.close();
		}
	}

	// Set the monitor to output data for TRIDYN
	if (flagTRIDYN) {
		// computeTRIDYN1D will be called at each timestep
		ierr = TSMonitorSet(ts, computeTRIDYN1D, NULL, NULL);
		checkPetscError(ierr,
				"setupPetsc1DMonitor: TSMonitorSet (computeTRIDYN1D) failed.");
	}

	// Set the monitor to output data for Alloy
	if (flagAlloy) {
		if (procId == 0) {
			// Create/open the output files
			std::fstream outputFile;
			outputFile.open("Alloy.dat", std::fstream::out);
			outputFile.close();
		}

		// computeAlloy1D will be called at each timestep
		ierr = TSMonitorSet(ts, computeAlloy1D, NULL, NULL);
		checkPetscError(ierr,
				"setupPetsc1DMonitor: TSMonitorSet (computeAlloy1D) failed.");
	}

	// Set the monitor to compute the temperature profile
	if (flagTemp) {

		if (procId == 0) {
			// Uncomment to clear the file where the retention will be written
			std::ofstream outputFile;
			outputFile.open("tempProf.txt");

			// Get the da from ts
			DM da;
			ierr = TSGetDM(ts, &da);
			checkPetscError(ierr, "setupPetsc1DMonitor: TSGetDM failed.");

			// Get the total size of the grid
			PetscInt Mx;
			ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE);
			checkPetscError(ierr, "setupPetsc1DMonitor: DMDAGetInfo failed.");

			// Get the physical grid
			auto grid = solverHandler.getXGrid();
			// Get the position of the surface
			int surfacePos = solverHandler.getSurfacePosition();

			// Loop on the entire grid
			for (int xi = surfacePos + solverHandler.getLeftOffset();
					xi < Mx - solverHandler.getRightOffset(); xi++) {
				// Set x
				double x = (grid[xi] + grid[xi + 1]) / 2.0 - grid[1];
				outputFile << x << " ";

			}
			outputFile << std::endl;
			outputFile.close();
		}

		// computeCumulativeHelium1D will be called at each timestep
		ierr = TSMonitorSet(ts, profileTemperature1D, NULL, NULL);
		checkPetscError(ierr,
				"setupPetsc1DMonitor: TSMonitorSet (profileTemperature1D) failed.");
	}

	// Set the monitor to simply change the previous time to the new time
	// monitorTime will be called at each timestep
	ierr = TSMonitorSet(ts, monitorTime, NULL, NULL);
	checkPetscError(ierr,
			"setupPetsc1DMonitor: TSMonitorSet (monitorTime) failed.");

	PetscFunctionReturn(0);
}

} /* end namespace monitor */
} /* end namespace solver */
} /* end namespace xolotl */