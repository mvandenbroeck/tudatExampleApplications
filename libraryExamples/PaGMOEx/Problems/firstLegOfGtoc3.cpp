/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <limits>

#include <Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/convertMeanToEccentricAnomalies.h>
#include <Tudat/Astrodynamics/MissionSegments/multiRevolutionLambertTargeterIzzo.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>

#include <Thesis/constants.h>
#include <Thesis/taylorSeriesIntegrator.h>

#include "PaGMOEx/Problems/firstLegOfGtoc3.h"

namespace pagmo { namespace problem {

FirstLegOfGtoc3::FirstLegOfGtoc3(
    const std::vector< std::vector< double > > problemBounds ) :
    base( problemBounds[ 0 ], problemBounds[ 1 ], 0, 1 ),
    problemBounds_( problemBounds )
{ }

//! Clone method.
base_ptr FirstLegOfGtoc3::clone( ) const {
        return base_ptr( new FirstLegOfGtoc3( *this ) );
}

//! Descriptive name of the problem
std::string FirstLegOfGtoc3::get_name() const {
    return "First leg (Earth - Asteroid 49) of GTOC3 problem";
}

//! Implementation of the objective function.
void FirstLegOfGtoc3::objfun_impl( fitness_vector &f, const decision_vector &xv ) const{
    using taylorSeriesIntegration ;

    // Gravitational parameter of the Sun
    double mu = 1.32712440018e+20;

    // Set initial and final position as those of Earth and Mars at
    // departure and arrival respectively.
    StateType initialState = getPlanetPosition( xv[0], "Earth" );
    StateType finalState   = getPlanetPosition( xv[0] + xv[1], "Asteroid 49" );

    MultiRevolutionLambertTargeterIzzo lambertTargeter( initialState.segment(0,3),
	finalState.segment(0,3), xv[1]*86400, mu );
    double deltaV = std::numeric_limits<double>::infinity();
    unsigned int maxrev = lambertTargeter.getMaximumNumberOfRevolutions( );

    // Go through all multi-revolution solutions and select the one
    // with the lowest delta-V
    for( unsigned int i = 0; i <= maxrev; ++i)
    {
        lambertTargeter.computeForRevolutionsAndBranch( i, false );
        deltaV = std::min( deltaV, ( initialState.segment(3,3)
            - lambertTargeter.getInertialVelocityAtDeparture( )).norm() +
            + ( finalState.segment(3,3)
            - lambertTargeter.getInertialVelocityAtArrival( )).norm());
    }
    f[0] = deltaV;
}

//! Function to obtain position of Earth and Asteroid
StateType FirstLegOfGtoc3::getPlanetPosition( const double date,
                                                const std::string planetName ) const {
    using tudat::orbital_element_conversions::convertKeplerianToCartesianElements;
    using tudat::orbital_element_conversions::convertMeanAnomalyToEccentricAnomaly;
    using tudat::orbital_element_conversions::convertEccentricAnomalyToTrueAnomaly;

    using namespace tudat;
    using namespace unit_conversions;

    // Create constants object
    ConstantsPointer constantsPointer;

    // Gravitational parameter of the Sun
    double mu = 1.32712440018e20;

    StateType stateKepl, stateCart;
    double n, jd0;
    if( planetName == "Earth" )
    {
        jd0 = 54000;
        stateKepl << 0.999988049532578 * constantsPointer->astronomicalUnit_,
                1.671681163160e-2,
                tudat::unit_conversions::convertDegreesToRadians( 0.8854353079654e-3 ),
                tudat::unit_conversions::convertDegreesToRadians( 287.61577546182 ),
                tudat::unit_conversions::convertDegreesToRadians( 175.40647696473 ),
                tudat::unit_conversions::convertDegreesToRadians( 257.606837075163 );
        n   = std::sqrt( mu / std::pow( stateKepl( 1 ), 3 ) ); // mean motion
    }
    else
    {
        jd0 = 54200;
        stateKepl << 0.9774002 * constantsPointer->astronomicalUnit_,
                0.06697124,
                tudat::unit_conversions::convertDegreesToRadians( 0.11024 ),
                tudat::unit_conversions::convertDegreesToRadians( 274.92230 ),
                tudat::unit_conversions::convertDegreesToRadians( 192.31139 ),
                tudat::unit_conversions::convertDegreesToRadians( 180.331421177838 );
        n   = std::sqrt( mu / std::pow( stateKepl( 1 ), 3 ) ); // mean motion
    }
    stateKepl( 5 ) = convertMeanAnomalyToEccentricAnomaly( stateKepl( 1 ),
        fmod( stateKepl( 5 ) + ( date - jd0 ) * 86400. * n, 2.*M_PI ) );
    stateKepl( 5 ) = convertEccentricAnomalyToTrueAnomaly( stateKepl( 5 ), stateKepl( 1 ) );
    stateCart = convertKeplerianToCartesianElements( stateKepl , mu );
    return stateCart;
}

}} // namespace problem; namespace pagmo

BOOST_CLASS_EXPORT_IMPLEMENT( pagmo::problem::FirstLegOfGtoc3 );
