/**
  @class UniformPrior

  Unused class that gives a uniform 1D distribution, rather than using the MC truth as a prior input

  @author Benjamin M Wynne bwynne@cern.ch
  @date 17-06-2010
 */


#ifndef UNIFORM_PRIOR_H
#define UNIFORM_PRIOR_H

#include "Distribution.h"

class UniformPrior : public Distribution
{
	public:
		UniformPrior();
		UniformPrior( double InputIntegral, IIndexCalculator * InputIndices );
		~UniformPrior();
};

#endif
