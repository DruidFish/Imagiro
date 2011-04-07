/**
  @class MonteCarloInformation

  A class that holds lists of file locations and formatting instructions for input Monte Carlo samples
  It can instantiate objects to read these files

  @author Benjamin M Wynne bwynne@cern.ch
  @date 06-01-2011
 */


#include "MonteCarloInformation.h"
#include "CombinedFileInput.h"
#include "InputNtuple.h"
#include "InputUETree.h"
#include <cstdlib>
#include <iostream>

const string NTUPLE_TYPE_STRING = "InputNtuple";
const string UE_TREE_TYPE_STRING = "InputUETree";

MonteCarloInformation::MonteCarloInformation()
{
	////////////////////////////////////////////////////////////
	//                                                        //
	// Define the file locations of the MC events             //
	// Also set how they appear on plots                      //
	//                                                        //
	////////////////////////////////////////////////////////////

	//Pythia J3
	combineFiles.push_back( false );
	truthPaths.push_back( "data/L1_EM5/D3PDpythiaJ3/combinedTruth.root" );
	recoPaths.push_back( "data/L1_EM5/D3PDpythiaJ3/combinedReconstructed.root" );
	descriptions.push_back( "PYTHIA J3" );
	colours.push_back( kMagenta );
	styles.push_back( 2 );
	inputTypes.push_back( "InputUETree" );
	internalPaths.push_back( "benTuple" );

	//Alpgen
	combineFiles.push_back( true );

	truthPaths.push_back( "data/L1_EM5/D3PDalpgenJimmyJ2/combinedTruth.root" );
	vector< string > alpgenExtraTruth;
	alpgenExtraTruth.push_back( "data/L1_EM5/D3PDalpgenJimmyJ3/combinedTruth.root" );
	alpgenExtraTruth.push_back( "data/L1_EM5/D3PDalpgenJimmyJ4/combinedTruth.root" );
	extraTruthPaths.push_back( alpgenExtraTruth );

	recoPaths.push_back( "data/L1_EM5/D3PDalpgenJimmyJ2/combinedReconstructed.root" );
	vector< string > alpgenExtraReco;
	alpgenExtraReco.push_back( "data/L1_EM5/D3PDalpgenJimmyJ3/combinedReconstructed.root" );
	alpgenExtraReco.push_back( "data/L1_EM5/D3PDalpgenJimmyJ4/combinedReconstructed.root" );
	extraRecoPaths.push_back( alpgenExtraReco );

	vector< double > alpgenWeights;
	alpgenWeights.push_back( 1.0 );
	alpgenWeights.push_back( 1.0 );
	alpgenWeights.push_back( 1.0 );
	inputWeights.push_back( alpgenWeights );

	descriptions.push_back( "ALPGEN" );
	colours.push_back( kRed );
	styles.push_back( 1 );
	inputTypes.push_back( "InputUETree" );
	internalPaths.push_back( "benTuple" );
}

MonteCarloInformation::~MonteCarloInformation()
{
}

int MonteCarloInformation::NumberOfSources()
{
	return descriptions.size();
}

string MonteCarloInformation::Description( unsigned int Index )
{
	if ( Index < 0 || Index >= descriptions.size() )
	{
		cerr << "Index out of range" << endl;
		exit(1);
	}
	else
	{
		return descriptions[Index];
	}
}

Color_t MonteCarloInformation::LineColour( unsigned int Index )
{
	if ( Index < 0 || Index >= colours.size() )
	{
		cerr << "Index out of range" << endl;
		exit(1);
	}
	else
	{
		return colours[Index];
	}
}

Style_t MonteCarloInformation::LineStyle( unsigned int Index )
{
	if ( Index < 0 || Index >= styles.size() )
	{
		cerr << "Index out of range" << endl;
		exit(1);
	}
	else
	{
		return styles[Index];
	}
}

//For a given regular index, find out the corresponding index in the "extra" vectors
int MonteCarloInformation::FindExtraIndex( unsigned int RegularIndex )
{
	//Check that there should be an extra index
	if ( combineFiles[ RegularIndex ] )
	{
		int extraIndex = 0;

		for ( unsigned int searchIndex = 0; searchIndex < RegularIndex; searchIndex++ )
		{
			if ( combineFiles[ searchIndex ] )
			{
				extraIndex++;
			}
		}

		return extraIndex;
	}
	else
	{
		return -1;
	}
}

IFileInput * MonteCarloInformation::MakeTruthInput( unsigned int Index )
{
	if (Index < 0 || Index >= styles.size() )
	{
		cerr << "Index out of range" << endl;
		exit(1);
	}
	else
	{
		vector< IFileInput* > inputs;
		inputs.push_back( InstantiateSingleInput( truthPaths[ Index ], Index ) );

		//Find out whether to combine files or not
		int extraIndex = FindExtraIndex( Index );
		if ( extraIndex >= 0 )
		{
			//Combined object
			for ( unsigned int pathIndex = 0; pathIndex < extraTruthPaths[ extraIndex ].size(); pathIndex++ )
			{
				inputs.push_back( InstantiateSingleInput( extraTruthPaths[ extraIndex ][ pathIndex ], Index ) );
			}

			//Check there's the same number of weights as files
			if ( inputs.size() != inputWeights[ extraIndex ].size() )
			{
				cerr << "Not the same number of truth inputs as weights: " << inputs.size() << " vs " << inputWeights[ extraIndex ].size() << endl;
				exit(1);
			}

			return new CombinedFileInput( inputs, inputWeights[ extraIndex ], descriptions[ Index ], Index );
		}
		else
		{
			//Just a single file input
			return inputs[0];
		}
	}
}

IFileInput * MonteCarloInformation::MakeReconstructedInput( unsigned int Index )
{
	if (Index < 0 || Index >= styles.size() )
	{
		cerr << "Index out of range" << endl;
		exit(1);
	}
	else
	{
		vector< IFileInput* > inputs;
		inputs.push_back( InstantiateSingleInput( recoPaths[ Index ], Index ) );

		//Find out whether to combine files or not
		int extraIndex = FindExtraIndex( Index );
		if ( extraIndex >= 0 )
		{
			//Combined object
			for ( unsigned int pathIndex = 0; pathIndex < extraRecoPaths[ extraIndex ].size(); pathIndex++ )
			{
				inputs.push_back( InstantiateSingleInput( extraRecoPaths[ extraIndex ][ pathIndex ], Index ) );
			}

			//Check there's the same number of weights as files
			if ( inputs.size() != inputWeights[ extraIndex ].size() )
			{
				cerr << "Not the same number of reco inputs as weights: " << inputs.size() << " vs " << inputWeights[ extraIndex ].size() << endl;
				exit(1);
			}

			return new CombinedFileInput( inputs, inputWeights[ extraIndex ], descriptions[ Index ], Index );
		}
		else
		{
			//Just a single file input
			return inputs[0];
		}
	}
}

IFileInput * MonteCarloInformation::InstantiateSingleInput( string Path, int Index )
{
	if ( inputTypes[ Index ] == NTUPLE_TYPE_STRING )
	{
		return new InputNtuple( Path, internalPaths[ Index ], descriptions[ Index ], Index );
	}
	else if ( inputTypes[ Index ] == UE_TREE_TYPE_STRING )
	{
		return new InputUETree( Path, internalPaths[ Index ], descriptions[ Index ], Index );
	}
	else
	{
		cerr << "Unrecognised input type: " << inputTypes[ Index ] << endl;
		exit(1);
	}
}
