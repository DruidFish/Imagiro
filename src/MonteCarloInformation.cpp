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
#include "TriggerChoosingInput.h"
#include <cstdlib>
#include <iostream>

const string NTUPLE_TYPE_STRING = "InputNtuple";
const string UE_TREE_TYPE_STRING = "InputUETree";
const string TRIGGER_CHOOSING_TYPE_STRING = "TriggerChoosingInput";

MonteCarloInformation::MonteCarloInformation()
{
	////////////////////////////////////////////////////////////
	//                                                        //
	// Define the file locations of the MC events             //
	// Also set how they appear on plots                      //
	//                                                        //
	////////////////////////////////////////////////////////////

	//Pythia 6 AMBT1
	combineFiles.push_back( true );

	truthPaths.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J0_AMBT1/combined.Truth.AntiKt4TopoEM.root" );
	vector< string > pythiaExtraTruth;
	pythiaExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J1_AMBT1/combined.Truth.AntiKt4TopoEM.root" );
	pythiaExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J2_AMBT1/combined.Truth.AntiKt4TopoEM.root" );
	pythiaExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J3_AMBT1/combined.Truth.AntiKt4TopoEM.root" );
	pythiaExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J4_AMBT1/combined.Truth.AntiKt4TopoEM.root" );
	pythiaExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J5_AMBT1/combined.Truth.AntiKt4TopoEM.root" );
	pythiaExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J6_AMBT1/combined.Truth.AntiKt4TopoEM.root" );
	//pythiaExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J7_AMBT1/combined.Truth.AntiKt4TopoEM.root" );
	//pythiaExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J8_AMBT1/combined.Truth.AntiKt4TopoEM.root" );
	extraTruthPaths.push_back( pythiaExtraTruth );

	recoPaths.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J0_AMBT1/combined.TriggerName.AntiKt4TopoEM.root" );
	vector< string > pythiaExtraReco;
	pythiaExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J1_AMBT1/combined.TriggerName.AntiKt4TopoEM.root" );
	pythiaExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J2_AMBT1/combined.TriggerName.AntiKt4TopoEM.root" );
	pythiaExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J3_AMBT1/combined.TriggerName.AntiKt4TopoEM.root" );
	pythiaExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J4_AMBT1/combined.TriggerName.AntiKt4TopoEM.root" );
	pythiaExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J5_AMBT1/combined.TriggerName.AntiKt4TopoEM.root" );
	pythiaExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J6_AMBT1/combined.TriggerName.AntiKt4TopoEM.root" );
	//pythiaExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J7_AMBT1/combined.TriggerName.AntiKt4TopoEM.root" );
	//pythiaExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J8_AMBT1/combined.TriggerName.AntiKt4TopoEM.root" );
	extraRecoPaths.push_back( pythiaExtraReco );

	vector< double > pythiaWeights;
	pythiaWeights.push_back( 9.8608E+06 / 199973.0 );
	pythiaWeights.push_back( 6.7818E+05 / 199968.0 );
	pythiaWeights.push_back( 4.0982E+04 / 199940.0 );
	pythiaWeights.push_back( 2.1929E+03 / 199929.0 );
	pythiaWeights.push_back( 87.701 / 199837.0 );
	pythiaWeights.push_back( 2.3501 / 199640.0 );
	pythiaWeights.push_back( 0.033614 / 199206.0 );
	//pythiaWeights.push_back( 1.3744E-04 / 198767.0 );
	//pythiaWeights.push_back( 6.2099E-09 / 197077.0 );
	inputWeights.push_back( pythiaWeights );

	descriptions.push_back( "PYTHIA6 AMBT1" );
	colours.push_back( kMagenta );
	styles.push_back( 2 );
	inputTypes.push_back( "TriggerChoosingInput" );
	internalTruth.push_back( "benTuple" );
	internalReco.push_back( "benTuple" );

	//Herwig++
/*	combineFiles.push_back( true );

	truthPaths.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J0_Herwig/combined.Truth.AntiKt4TopoEM.root" );
	vector< string > herwigExtraTruth;
	herwigExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J1_Herwig/combined.Truth.AntiKt4TopoEM.root" );
	herwigExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J2_Herwig/combined.Truth.AntiKt4TopoEM.root" );
	herwigExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J3_Herwig/combined.Truth.AntiKt4TopoEM.root" );
	herwigExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J4_Herwig/combined.Truth.AntiKt4TopoEM.root" );
	herwigExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J5_Herwig/combined.Truth.AntiKt4TopoEM.root" );
	herwigExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J6_Herwig/combined.Truth.AntiKt4TopoEM.root" );
	extraTruthPaths.push_back( herwigExtraTruth );

	recoPaths.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J0_Herwig/combined.TriggerName.AntiKt4TopoEM.root" );
	vector< string > herwigExtraReco;
	herwigExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J1_Herwig/combined.TriggerName.AntiKt4TopoEM.root" );
	herwigExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J2_Herwig/combined.TriggerName.AntiKt4TopoEM.root" );
	herwigExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J3_Herwig/combined.TriggerName.AntiKt4TopoEM.root" );
	herwigExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J4_Herwig/combined.TriggerName.AntiKt4TopoEM.root" );
	herwigExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J5_Herwig/combined.TriggerName.AntiKt4TopoEM.root" );
	herwigExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J6_Herwig/combined.TriggerName.AntiKt4TopoEM.root" );
	extraRecoPaths.push_back( herwigExtraReco );

	vector< double > herwigWeights;
	herwigWeights.push_back( 9.6139E+06 / 398799.0 );
	herwigWeights.push_back( 7.4366E+05 / 397897.0 );
	herwigWeights.push_back( 4.4307E+04 / 398498.0 );
	herwigWeights.push_back( 2357.6 / 399598.0 );
	herwigWeights.push_back( 94.236 / 397443.0 );
	herwigWeights.push_back( 2.5813 / 398094.0 );
	herwigWeights.push_back( 0.039439 / 397597.0 );
	inputWeights.push_back( herwigWeights );

	descriptions.push_back( "HERWIG++" );
	colours.push_back( kRed + 3 );
	styles.push_back( 8 );
	inputTypes.push_back( "TriggerChoosingInput" );
	internalTruth.push_back( "benTuple" );
	internalReco.push_back( "benTuple" );*/

	//Pythia 6 Perugia 2010
	combineFiles.push_back( true );

	truthPaths.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J0_Perugia/combined.Truth.AntiKt4TopoEM.root" );
	vector< string > perugiaExtraTruth;
	perugiaExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J1_Perugia/combined.Truth.AntiKt4TopoEM.root" );
	perugiaExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J2_Perugia/combined.Truth.AntiKt4TopoEM.root" );
	perugiaExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J3_Perugia/combined.Truth.AntiKt4TopoEM.root" );
	perugiaExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J4_Perugia/combined.Truth.AntiKt4TopoEM.root" );
	perugiaExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J5_Perugia/combined.Truth.AntiKt4TopoEM.root" );
	perugiaExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J6_Perugia/combined.Truth.AntiKt4TopoEM.root" );
	//perugiaExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J7_Perugia/combined.Truth.AntiKt4TopoEM.root" );
	//perugiaExtraTruth.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J8_Perugia/combined.Truth.AntiKt4TopoEM.root" );
	extraTruthPaths.push_back( perugiaExtraTruth );

	recoPaths.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J0_Perugia/combined.TriggerName.AntiKt4TopoEM.root" );
	vector< string > perugiaExtraReco;
	perugiaExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J1_Perugia/combined.TriggerName.AntiKt4TopoEM.root" );
	perugiaExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J2_Perugia/combined.TriggerName.AntiKt4TopoEM.root" );
	perugiaExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J3_Perugia/combined.TriggerName.AntiKt4TopoEM.root" );
	perugiaExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J4_Perugia/combined.TriggerName.AntiKt4TopoEM.root" );
	perugiaExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J5_Perugia/combined.TriggerName.AntiKt4TopoEM.root" );
	perugiaExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J6_Perugia/combined.TriggerName.AntiKt4TopoEM.root" );
	//perugiaExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J7_Perugia/combined.TriggerName.AntiKt4TopoEM.root" );
	//perugiaExtraReco.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version7/J8_Perugia/combined.TriggerName.AntiKt4TopoEM.root" );
	extraRecoPaths.push_back( perugiaExtraReco );

	vector< double > perugiaWeights;
	perugiaWeights.push_back( 7.7714E+06 / 398849.0 );
	perugiaWeights.push_back( 5.0385E+05 / 1498393.0 );
	perugiaWeights.push_back( 2.9358E+04 / 988144.0 );
	perugiaWeights.push_back( 1.5600E+03 / 399497.0 );
	perugiaWeights.push_back( 64.393 / 398199.0 );
	perugiaWeights.push_back( 1.8764 / 398046.0 );
	perugiaWeights.push_back( 3.0412E-02 / 397900.0 );
	//perugiaWeights.push_back( 1.3212E-04 / 399145.0 );
	//perugiaWeights.push_back( 4.9978E-09 / 396886.0 );
	inputWeights.push_back( perugiaWeights );

	descriptions.push_back( "PYTHIA6 Perugia2010" );
	colours.push_back( kGreen );
	styles.push_back( 3 );
	inputTypes.push_back( "TriggerChoosingInput" );
	internalTruth.push_back( "benTuple" );
	internalReco.push_back( "benTuple" );

	//Sherpa - this one doesn't need jet slice combining, so don't merge the files
/*	combineFiles.push_back( true );

	char filePath[1000];
	//truthPaths.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version4/user.bwynne.D3PDtoUE.mc10_7TeV.114680.SherpaQCD6jets20GeV.merge.NTUP_JETMET.e719_s933_s946_r2302_r2300_p515.v4.110713075526/user.bwynne.002024._00001.truthOutput.AntiKt4TopoEM.root" );
	truthPaths.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version6/Sherpa/user.bwynne.002076._00001.truthOutput.AntiKt4TopoEM.root" );
	//recoPaths.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version4/user.bwynne.D3PDtoUE.mc10_7TeV.114680.SherpaQCD6jets20GeV.merge.NTUP_JETMET.e719_s933_s946_r2302_r2300_p515.v4.110713075526/user.bwynne.002024._00001.reconstructedOutput.TriggerName.AntiKt4TopoEM.root" );
	recoPaths.push_back( "/Disk/speyside7/Grid/grid-files/bwynne/Version6/Sherpa/user.bwynne.002076._00001.reconstructedOutput.TriggerName.AntiKt4TopoEM.root" );
	vector< string > extraSherpaTruth, extraSherpaReco;
	//for ( unsigned int fileIndex = 2; fileIndex <= 72; fileIndex++ )
	for ( unsigned int fileIndex = 2; fileIndex <= 20; fileIndex++ )
	{
		//Extra truth
		//sprintf( filePath, "/Disk/speyside7/Grid/grid-files/bwynne/Version4/user.bwynne.D3PDtoUE.mc10_7TeV.114680.SherpaQCD6jets20GeV.merge.NTUP_JETMET.e719_s933_s946_r2302_r2300_p515.v4.110713075526/user.bwynne.002024._000%02d.truthOutput.AntiKt4TopoEM.root", fileIndex );
		sprintf( filePath, "/Disk/speyside7/Grid/grid-files/bwynne/Version6/Sherpa/user.bwynne.002076._000%02d.truthOutput.AntiKt4TopoEM.root", fileIndex );
		extraSherpaTruth.push_back( filePath );

		//Extra reco
		//sprintf( filePath, "/Disk/speyside7/Grid/grid-files/bwynne/Version4/user.bwynne.D3PDtoUE.mc10_7TeV.114680.SherpaQCD6jets20GeV.merge.NTUP_JETMET.e719_s933_s946_r2302_r2300_p515.v4.110713075526/user.bwynne.002024._000%02d.reconstructedOutput.TriggerName.AntiKt4TopoEM.root", fileIndex );
		sprintf( filePath, "/Disk/speyside7/Grid/grid-files/bwynne/Version6/Sherpa/user.bwynne.002076._000%02d.reconstructedOutput.TriggerName.AntiKt4TopoEM.root", fileIndex );
		extraSherpaReco.push_back( filePath );
	}
	extraTruthPaths.push_back( extraSherpaTruth );
	extraRecoPaths.push_back( extraSherpaReco );
	inputWeights.push_back( vector< double >( 20, 1.0 ) );

	descriptions.push_back( "SHERPA" );
	colours.push_back( kOrange );
	styles.push_back( 4 );
	inputTypes.push_back( "TriggerChoosingInput" );
	internalTruth.push_back( "benTuple" );
	internalReco.push_back( "benTuple" );*/
}

MonteCarloInformation::~MonteCarloInformation()
{
}

unsigned int MonteCarloInformation::NumberOfSources()
{
	return descriptions.size();
}

string MonteCarloInformation::Description( unsigned int Index )
{
	if ( Index >= descriptions.size() )
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
	if ( Index >= colours.size() )
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
	if ( Index >= styles.size() )
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

IFileInput * MonteCarloInformation::MakeTruthInput( unsigned int Index, ObservableList * RelevanceChecker )
{
	if ( Index >= styles.size() )
	{
		cerr << "Index out of range" << endl;
		exit(1);
	}
	else
	{
		//Find out whether to combine files or not
		int extraIndex = FindExtraIndex( Index );
		if ( extraIndex >= 0 )
		{
			//Make a vector of all the file paths
			vector<string> allTruthPaths( 1, truthPaths[ Index ] );
			allTruthPaths.insert( allTruthPaths.end(), extraTruthPaths[ extraIndex ].begin(), extraTruthPaths[ extraIndex ].end() );

			//Check there's the same number of weights as files
			if ( allTruthPaths.size() != inputWeights[ extraIndex ].size() )
			{
				cerr << "Not the same number of truth inputs as weights: " << allTruthPaths.size() << " vs " << inputWeights[ extraIndex ].size() << endl;
				exit(1);
			}

			return new CombinedFileInput( allTruthPaths, inputWeights[ extraIndex ], internalTruth[ Index ], inputTypes[ Index ], descriptions[ Index ], Index, RelevanceChecker );
		}
		else
		{
			//Just a single file input
			return InstantiateSingleInput( truthPaths[ Index ], internalTruth[ Index ], Index, RelevanceChecker );
		}
	}
}

IFileInput * MonteCarloInformation::MakeReconstructedInput( unsigned int Index, ObservableList * RelevanceChecker )
{
	if ( Index >= styles.size() )
	{
		cerr << "Index out of range" << endl;
		exit(1);
	}
	else
	{
		//Find out whether to combine files or not
		int extraIndex = FindExtraIndex( Index );
		if ( extraIndex >= 0 )
		{
			//Make a vector of all the file paths
			vector<string> allRecoPaths( 1, recoPaths[ Index ] );
			allRecoPaths.insert( allRecoPaths.end(), extraRecoPaths[ extraIndex ].begin(), extraRecoPaths[ extraIndex ].end() );

			//Check there's the same number of weights as files
			if ( allRecoPaths.size() != inputWeights[ extraIndex ].size() )
			{
				cerr << "Not the same number of reco inputs as weights: " << allRecoPaths.size() << " vs " << inputWeights[ extraIndex ].size() << endl;
				exit(1);
			}

			return new CombinedFileInput( allRecoPaths, inputWeights[ extraIndex ], internalReco[ Index ], inputTypes[ Index ], descriptions[ Index ], Index, RelevanceChecker );
		}
		else
		{
			//Just a single file input
			return InstantiateSingleInput( recoPaths[ Index ], internalReco[ Index ], Index, RelevanceChecker );
		}
	}
}

IFileInput * MonteCarloInformation::InstantiateSingleInput( string FilePath, string InternalPath, unsigned int Index, ObservableList * RelevanceChecker )
{
	if ( inputTypes[ Index ] == NTUPLE_TYPE_STRING )
	{
		return new InputNtuple( FilePath, InternalPath, descriptions[ Index ], Index, RelevanceChecker );
	}
	else if ( inputTypes[ Index ] == UE_TREE_TYPE_STRING )
	{
		return new InputUETree( FilePath, InternalPath, descriptions[ Index ], Index, RelevanceChecker );
	}
	else if ( inputTypes[ Index ] == TRIGGER_CHOOSING_TYPE_STRING )
	{
		return new TriggerChoosingInput( FilePath, InternalPath, descriptions[ Index ], Index, RelevanceChecker );
	}
	else
	{
		cerr << "Unrecognised input type: " << inputTypes[ Index ] << endl;
		exit(1);
	}
}
