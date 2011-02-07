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
