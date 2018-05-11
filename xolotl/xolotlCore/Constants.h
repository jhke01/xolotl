/*
 * Constants.h
 *
 *  Created on: May 6, 2013
 *      Author: bkj
 */
#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include <memory>
#include <cmath>
#include "MathUtils.h"

namespace xolotlCore {

//! Definitions of fundamental constants used in Xolotl

//! The Boltzmann constant in units of eV K^-1.
constexpr double kBoltzmann = 8.61733240000000000E-5;

//! Pi, taken from "100000 digits of Pi,"
//! at http://www.geom.uiuc.edu/~huberty/math5337/groupe/digits.html
constexpr double pi = 3.1415926535897932;

//! Lattice Parameter. Equal to 3.17 Angstroms, taken from Becquart et. al.
//! Journal of Nuclear Materials 403 (2010) 75–88. Given in units here of nm.
constexpr double tungstenLatticeConstant = 0.31700000000000000;

//! Lattice Parameter for UO2
constexpr double uraniumDioxydeLatticeConstant = 0.57400000000000000;

//! Lattice Parameter for Iron
constexpr double ironLatticeConstant = 0.28700000000000000;

// Tungsten heat coefficient = lambda / (rho * C) in nm2 s-1
constexpr double tungstenHeatCoefficient = 6.835e13;

// UO2 heat coefficient = lambda / (rho * C) in nm2 s-1
constexpr double uo2HeatCoefficient = 0.0;

// Iron heat coefficient = lambda / (rho * C) in nm2 s-1
constexpr double feHeatCoefficient = 0.0;

//! Parameters for biased sink in the iron case
static const double reactionRadius = ironLatticeConstant
		* std::cbrt(3.0 / pi) * 0.5;
static const double r0 = ironLatticeConstant * 0.75 * std::sqrt(3.0);
constexpr double rho = 0.0003;
static const double sinkStrength = -4.0 * pi * rho
		/ log(pi * rho * ipow<2>(reactionRadius + r0));
constexpr double sinkBias = 1.05;

} /* end namespace xolotlCore */
#endif /* CONSTANTS_H_ */
