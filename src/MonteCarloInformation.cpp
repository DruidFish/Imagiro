#include "MonteCarloInformation.h"
#include <cstdlib>
#include <iostream>

MonteCarloInformation::MonteCarloInformation()
{
	////////////////////////////////////////////////////////////
        //                                                        //
	// Define the file locations of the MC events             //
	// Also set how they appear on plots                      //
	//                                                        //
	////////////////////////////////////////////////////////////

	//Pythia 6 (MC09?)
	truthPaths.push_back( "data/user.bwynne.LeadingJetModifiedv3.Pythia6.MC.TruthJet/mergedFile.root" );
	recoPaths.push_back( "data/user.bwynne.LeadingJetModifiedv3.Pythia6.MC.CaloJet/mergedFile.root" );
	descriptions.push_back( "Pythia6" );
	colours.push_back( kMagenta );
	styles.push_back( 2 );

	//AMBT1
	truthPaths.push_back( "data/user.bwynne.LeadingJetModifiedv3.AMBT1.MC.TruthJet/mergedFile.root" );
	recoPaths.push_back( "data/user.bwynne.LeadingJetModifiedv3.AMBT1.MC.CaloJet/mergedFile.root" );
	descriptions.push_back( "AMBT1" );
	colours.push_back( kRed );
	styles.push_back( 1 );

	//DW
	truthPaths.push_back( "data/user.bwynne.LeadingJetModifiedv3.DW.MC.TruthJet/mergedFile.root" );
	recoPaths.push_back( "data/user.bwynne.LeadingJetModifiedv3.DW.MC.CaloJet/mergedFile.root" );
	descriptions.push_back( "DW" );
	colours.push_back( kBlue );
	styles.push_back( 7 );

	//Herwig++
	truthPaths.push_back( "data/user.bwynne.LeadingJetModifiedv3.HerwigPP.MC.TruthJet/mergedFile.root" );
	recoPaths.push_back( "data/user.bwynne.LeadingJetModifiedv3.HerwigPP.MC.CaloJet/mergedFile.root" );
	descriptions.push_back( "HerwigPP" );
	colours.push_back( kRed + 3 );
	styles.push_back( 8 );

	//Pythia 8 non diffractive
	truthPaths.push_back( "data/user.bwynne.LeadingJetModifiedv3.Pythia8NonDiff.MC.TruthJet/mergedFile.root" );
	recoPaths.push_back( "data/user.bwynne.LeadingJetModifiedv3.Pythia8NonDiff.MC.CaloJet/mergedFile.root" );
	descriptions.push_back( "Pythia8NonDiff" );
	colours.push_back( kGreen + 2 );
	styles.push_back( 5 );
}

MonteCarloInformation::~MonteCarloInformation()
{
}

int MonteCarloInformation::NumberOfSources()
{
	return descriptions.size();
}

string MonteCarloInformation::Description( int Index )
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

string MonteCarloInformation::TruthFilePath( int Index )
{
	if ( Index < 0 || Index >= truthPaths.size() )
	{
		cerr << "Index out of range" << endl;
		exit(1);
	}
	else
	{
		return truthPaths[Index];
	}
}

string MonteCarloInformation::ReconstructedFilePath( int Index )
{
	if ( Index < 0 || Index >= recoPaths.size() )
	{
		cerr << "Index out of range" << endl;
		exit(1);
	}
	else
	{
		return recoPaths[Index];
	}
}

Color_t MonteCarloInformation::LineColour( int Index )
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

Style_t MonteCarloInformation::LineStyle( int Index )
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
