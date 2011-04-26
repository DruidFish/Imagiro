/**
  @class MonteCarloInformation

  A class that holds lists of file locations and formatting instructions for input Monte Carlo samples
  It can instantiate objects to read these files

  @author Benjamin M Wynne bwynne@cern.ch
  @date 06-01-2011
 */

#ifndef MONTE_CARLO_INFORMATION_H
#define MONTE_CARLO_INFORMATION_H

#include "Rtypes.h"
#include <string>
#include <vector>
#include "IFileInput.h"

using namespace std;

class MonteCarloInformation
{
	public:
		MonteCarloInformation();
		~MonteCarloInformation();

		//Return the number of different MC samples defined here
		int NumberOfSources();

		//Make an object that will read the truth or reconstructed files(s) for this MC sample
		IFileInput * MakeTruthInput( unsigned int Index );
		IFileInput * MakeReconstructedInput( unsigned int Index );

		//Return the colour and style of the line to plot for this MC sample
		Color_t LineColour( unsigned int Index );
		Style_t LineStyle( unsigned int Index );

		//Return the human-readable description of this MC sample
		string Description( unsigned int Index );

	private:
		//Make a file reader of the correct type
		IFileInput * InstantiateSingleInput( string FilePath, string InternalPath, int Index );

		//Return the index of the additional files to combine for this MC sample
		int FindExtraIndex( unsigned int RegularIndex );

		vector< vector< string > > extraTruthPaths, extraRecoPaths;
		vector< bool > combineFiles;
		vector< vector< double > > inputWeights;
		vector<string> descriptions, truthPaths, recoPaths, inputTypes, internalTruth, internalReco;
		vector< Color_t > colours;
		vector< Style_t > styles;

};

#endif
