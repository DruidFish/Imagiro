/**
  @class UniformPrior

  Unused class that gives a uniform 1D distribution, rather than using the MC truth as a prior input

  @author Benjamin M Wynne bwynne@cern.ch
  @date 17-06-2010
 */

#include "UniformPrior.h"

//Default constructor
UniformPrior::UniformPrior()
{
}

//Uniform distribution
UniformPrior::UniformPrior( double InputIntegral, Indices * InputIndices )
{
	int binNumber = InputIndices->GetBinNumber();
	binValues = vector<double>( binNumber, InputIntegral / (double)binNumber );
	integral = InputIntegral;
	indexCalculator = InputIndices;
}

//Destructor
UniformPrior::~UniformPrior()
{
}
