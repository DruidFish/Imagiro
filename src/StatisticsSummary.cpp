/**
  @class StatisticsSummary

  A class that maintains a statistical summary of events passed to it

  @author Benjamin M Wynne bwynne@cern.ch
  @date 09-03-2011
 */

#include "StatisticsSummary.h"
#include <cmath>
#include <iostream>

//Default constructor
StatisticsSummary::StatisticsSummary() : freshStart( true )
{
}

//Destructor
StatisticsSummary::~StatisticsSummary()
{
}

//Add a new event to the summary - update values accordingly
void StatisticsSummary::StoreEvent( double Value, double Weight )
{
	//Update the mean
	meanNumerator += Value * Weight;
	meanDenominator += Weight;

	//Update the mean of the value squared
	meanSquaredNumerator += Value * Value * Weight * Weight;
	meanSquaredDenominator += Weight * Weight;

	if ( freshStart )
	{
		//Populate the maximum and minimum
		currentMaximum = Value;
		currentMinimum = Value;
		freshStart = false;
	}
	else
	{
		//Update the maximum and minimum
		if ( Value > currentMaximum )
		{
			currentMaximum = Value;
		}

		if ( Value < currentMinimum )
		{
			currentMinimum = Value;
		}
	}
}

//Return the mean of a vector of doubles
double StatisticsSummary::Mean()
{
	return meanNumerator / meanDenominator;
}

//Return the variance of a vector of doubles
double StatisticsSummary::Variance()
{
	return ( meanSquaredNumerator / meanSquaredDenominator ) - ( Mean() * Mean() );
}

//Return the ideal number of bins for a histogram of a vector of doubles
//Uses D. Scott's method, published 1979
double StatisticsSummary::OptimumBinWidth()
{
	return 3.49 * sqrt( Variance() ) * pow( meanDenominator, -( 1.0 / 3.0 ) );
}
int StatisticsSummary::OptimumBinNumber( double Minimum, double Maximum )
{
	double width = OptimumBinWidth();
	double range = Maximum - Minimum;
	return (int)ceil( range / width );
}

//Returns the maximum and minimum of a vector of doubles 
double StatisticsSummary::Maximum()
{
	return currentMaximum;
}
double StatisticsSummary::Minimum()
{
	return currentMinimum;	
}
