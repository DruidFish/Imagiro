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

	//Pythia
	combineFiles.push_back( true );

	truthPaths.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/L1_J5.v2/J0_pythia/combinedTruth.root" );
	vector< string > pythiaExtraTruth;
	pythiaExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/L1_J5.v2/J1_pythia/combinedTruth.root" );
	pythiaExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/L1_J5.v2/J2_pythia/combinedTruth.root" );
	pythiaExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/L1_J5.v2/J3_pythia/combinedTruth.root" );
	//pythiaExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/L1_J5.v2/J4_pythia/combinedTruth.root" );
	//pythiaExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/L1_J5.v2/J5_pythia/combinedTruth.root" );
	//pythiaExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/L1_J5.v2/J6_pythia/combinedTruth.root" );
	//pythiaExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/L1_J5.v2/J7_pythia/combinedTruth.root" );
	//pythiaExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/L1_J5.v2/J8_pythia/combinedTruth.root" );
	extraTruthPaths.push_back( pythiaExtraTruth );

	recoPaths.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/L1_J5.v2/J0_pythia/combinedReco.root" );
	vector< string > pythiaExtraReco;
	pythiaExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/L1_J5.v2/J1_pythia/combinedReco.root" );
	pythiaExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/L1_J5.v2/J2_pythia/combinedReco.root" );
	pythiaExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/L1_J5.v2/J3_pythia/combinedReco.root" );
	//pythiaExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/L1_J5.v2/J4_pythia/combinedReco.root" );
	//pythiaExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/L1_J5.v2/J5_pythia/combinedReco.root" );
	//pythiaExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/L1_J5.v2/J6_pythia/combinedReco.root" );
	//pythiaExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/L1_J5.v2/J7_pythia/combinedReco.root" );
	//pythiaExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/L1_J5.v2/J8_pythia/combinedReco.root" );
	extraRecoPaths.push_back( pythiaExtraReco );

	vector< double > pythiaWeights;
	pythiaWeights.push_back( 9.8608E+06 / 199973.0 );
	pythiaWeights.push_back( 6.7818E+05 / 199968.0 );
	pythiaWeights.push_back( 4.0982E+04 / 199940.0 );
	pythiaWeights.push_back( 2.1929E+03 / 199929.0 );
	//pythiaWeights.push_back( 87.701 / 199837.0 );
	//pythiaWeights.push_back( 2.3501 / 199640.0 );
	//pythiaWeights.push_back( 0.033614 / 199206.0 );
	//pythiaWeights.push_back( 1.3744E-04 / 198767.0 );
	//pythiaWeights.push_back( 6.2099E-09 / 197077.0 );
	inputWeights.push_back( pythiaWeights );

	descriptions.push_back( "PYTHIA" );
	colours.push_back( kMagenta );
	styles.push_back( 2 );
	inputTypes.push_back( "InputUETree" );
	internalTruth.push_back( "benTuple" );
	internalReco.push_back( "benTuple" );

	//Herwig
/*	combineFiles.push_back( true );

	truthPaths.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/L1_J5.v2/J0_herwig/combinedTruth.root" );
	vector< string > herwigExtraTruth;
	herwigExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/L1_J5.v2/J1_herwig/combinedTruth.root" );
	herwigExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/L1_J5.v2/J2_herwig/combinedTruth.root" );
	herwigExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/L1_J5.v2/J3_herwig/combinedTruth.root" );
	//herwigExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/L1_J5.v2/J4_herwig/combinedTruth.root" );
	//herwigExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/L1_J5.v2/J5_herwig/combinedTruth.root" );
	//herwigExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/L1_J5.v2/J6_herwig/combinedTruth.root" );
	extraTruthPaths.push_back( herwigExtraTruth );

	recoPaths.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/L1_J5.v2/J0_herwig/combinedReco.root" );
	vector< string > herwigExtraReco;
	herwigExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/L1_J5.v2/J1_herwig/combinedReco.root" );
	herwigExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/L1_J5.v2/J2_herwig/combinedReco.root" );
	herwigExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/L1_J5.v2/J3_herwig/combinedReco.root" );
	//herwigExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/L1_J5.v2/J4_herwig/combinedReco.root" );
	//herwigExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/L1_J5.v2/J5_herwig/combinedReco.root" );
	//herwigExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/L1_J5.v2/J6_herwig/combinedReco.root" );
	extraRecoPaths.push_back( herwigExtraReco );

	vector< double > herwigWeights;
	herwigWeights.push_back( 9.6139E+06 / 398799.0 );
	herwigWeights.push_back( 7.4366E+05 / 397897.0 );
	herwigWeights.push_back( 4.4307E+04 / 398498.0 );
	herwigWeights.push_back( 2357.6 / 399598.0 );
	//herwigWeights.push_back( 94.236 / 397443.0 );
	//herwigWeights.push_back( 2.5813 / 398094.0 );
	//herwigWeights.push_back( 0.039439 / 397597.0 );
	inputWeights.push_back( herwigWeights );

	descriptions.push_back( "HERWIG" );
	colours.push_back( kRed + 3 );
	styles.push_back( 8 );
	inputTypes.push_back( "InputUETree" );
	internalTruth.push_back( "benTuple" );
	internalReco.push_back( "benTuple" );*/

	//Pythia 6 (MC09?)
	/*combineFiles.push_back( false );
	truthPaths.push_back( "data/user.bwynne.LeadingJetModifiedv3.Pythia6.MC.TruthJet/mergedFile.root" );
	recoPaths.push_back( "data/user.bwynne.LeadingJetModifiedv3.Pythia6.MC.CaloJet/mergedFile.root" );
	descriptions.push_back( "PYTHIA6 ATLAS MC09" );
	colours.push_back( kMagenta );
	styles.push_back( 2 );
	inputTypes.push_back( "InputNtuple" );
	internalTruth.push_back( "benTuple" );
	internalReco.push_back( "benTuple" );

	//AMBT1
	combineFiles.push_back( false );
	truthPaths.push_back( "data/user.bwynne.LeadingJetModifiedv3.AMBT1.MC.TruthJet/mergedFile.root" );
	recoPaths.push_back( "data/user.bwynne.LeadingJetModifiedv3.AMBT1.MC.CaloJet/mergedFile.root" );
	descriptions.push_back( "PYTHIA6 AMBT1" );
	colours.push_back( kRed );
	styles.push_back( 1 );
	inputTypes.push_back( "InputNtuple" );
	internalTruth.push_back( "benTuple" );
	internalReco.push_back( "benTuple" );

	//DW
	combineFiles.push_back( false );
	truthPaths.push_back( "data/user.bwynne.LeadingJetModifiedv3.DW.MC.TruthJet/mergedFile.root" );
	recoPaths.push_back( "data/user.bwynne.LeadingJetModifiedv3.DW.MC.CaloJet/mergedFile.root" );
	descriptions.push_back( "PYTHIA6 DW" );
	colours.push_back( kBlue );
	styles.push_back( 7 );
	inputTypes.push_back( "InputNtuple" );
	internalTruth.push_back( "benTuple" );
	internalReco.push_back( "benTuple" );

	//Herwig++
	combineFiles.push_back( false );
	truthPaths.push_back( "data/user.bwynne.LeadingJetModifiedv3.HerwigPP.MC.TruthJet/mergedFile.root" );
	recoPaths.push_back( "data/user.bwynne.LeadingJetModifiedv3.HerwigPP.MC.CaloJet/mergedFile.root" );
	descriptions.push_back( "Herwig++ ATLAS MC09" );
	colours.push_back( kRed + 3 );
	styles.push_back( 8 );
	inputTypes.push_back( "InputNtuple" );
	internalTruth.push_back( "benTuple" );
	internalReco.push_back( "benTuple" );

	//Pythia 8 non diffractive
	combineFiles.push_back( false );
	truthPaths.push_back( "data/user.bwynne.LeadingJetModifiedv3.Pythia8NonDiff.MC.TruthJet/mergedFile.root" );
	recoPaths.push_back( "data/user.bwynne.LeadingJetModifiedv3.Pythia8NonDiff.MC.CaloJet/mergedFile.root" );
	descriptions.push_back( "PYTHIA8 Non-Diffractive" );
	colours.push_back( kGreen + 2 );
	styles.push_back( 5 );
	inputTypes.push_back( "InputNtuple" );
	internalTruth.push_back( "benTuple" );
	internalReco.push_back( "benTuple" );*/
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
		inputs.push_back( InstantiateSingleInput( truthPaths[ Index ], internalTruth[ Index ], Index ) );

		//Find out whether to combine files or not
		int extraIndex = FindExtraIndex( Index );
		if ( extraIndex >= 0 )
		{
			//Combined object
			for ( unsigned int pathIndex = 0; pathIndex < extraTruthPaths[ extraIndex ].size(); pathIndex++ )
			{
				inputs.push_back( InstantiateSingleInput( extraTruthPaths[ extraIndex ][ pathIndex ], internalTruth[ Index ], Index ) );
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
		inputs.push_back( InstantiateSingleInput( recoPaths[ Index ], internalReco[ Index ], Index ) );

		//Find out whether to combine files or not
		int extraIndex = FindExtraIndex( Index );
		if ( extraIndex >= 0 )
		{
			//Combined object
			for ( unsigned int pathIndex = 0; pathIndex < extraRecoPaths[ extraIndex ].size(); pathIndex++ )
			{
				inputs.push_back( InstantiateSingleInput( extraRecoPaths[ extraIndex ][ pathIndex ], internalReco[ Index ], Index ) );
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

IFileInput * MonteCarloInformation::InstantiateSingleInput( string FilePath, string InternalPath, int Index )
{
	if ( inputTypes[ Index ] == NTUPLE_TYPE_STRING )
	{
		return new InputNtuple( FilePath, InternalPath, descriptions[ Index ], Index );
	}
	else if ( inputTypes[ Index ] == UE_TREE_TYPE_STRING )
	{
		return new InputUETree( FilePath, InternalPath, descriptions[ Index ], Index );
	}
	else
	{
		cerr << "Unrecognised input type: " << inputTypes[ Index ] << endl;
		exit(1);
	}
}
