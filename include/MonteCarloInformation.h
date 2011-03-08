/**
  @class MonteCarloInformation

  A simple class that holds lists of file locations and formatting instructions for input Monte Carlo samples

  @author Benjamin M Wynne bwynne@cern.ch
  @date 06-01-2011
 */

#ifndef MONTE_CARLO_INFORMATION_H
#define MONTE_CARLO_INFORMATION_H

#include "Rtypes.h"
#include <string>
#include <vector>

using namespace std;

class MonteCarloInformation
{
	public:
		MonteCarloInformation();
		~MonteCarloInformation();

		int NumberOfSources();
		string Description( int Index );
		string TruthFilePath( int Index );
		string ReconstructedFilePath( int Index );
		Color_t LineColour( int Index );
		Style_t LineStyle( int Index );

	private:
		vector<string> descriptions, truthPaths, recoPaths;
		vector< Color_t > colours;
		vector< Style_t > styles;

};

#endif
