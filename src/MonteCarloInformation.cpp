/**
  @class MonteCarloInformation

  A simple class that holds lists of file locations and formatting instructions for input Monte Carlo samples

  @author Benjamin M Wynne bwynne@cern.ch
  @date 06-01-2011
 */


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
	descriptions.push_back( "PYTHIA6 ATLAS MC09" );
	colours.push_back( kMagenta );
	styles.push_back( 2 );

	//AMBT1
	truthPaths.push_back( "data/user.bwynne.LeadingJetModifiedv3.AMBT1.MC.TruthJet/mergedFile.root" );
	recoPaths.push_back( "data/user.bwynne.LeadingJetModifiedv3.AMBT1.MC.CaloJet/mergedFile.root" );
	descriptions.push_back( "PYTHIA6 AMBT1" );
	colours.push_back( kRed );
	styles.push_back( 1 );

	//DW
	truthPaths.push_back( "data/user.bwynne.LeadingJetModifiedv3.DW.MC.TruthJet/mergedFile.root" );
	recoPaths.push_back( "data/user.bwynne.LeadingJetModifiedv3.DW.MC.CaloJet/mergedFile.root" );
	descriptions.push_back( "PYTHIA6 DW" );
	colours.push_back( kBlue );
	styles.push_back( 7 );

	//Herwig++
	truthPaths.push_back( "data/user.bwynne.LeadingJetModifiedv3.HerwigPP.MC.TruthJet/mergedFile.root" );
	recoPaths.push_back( "data/user.bwynne.LeadingJetModifiedv3.HerwigPP.MC.CaloJet/mergedFile.root" );
	descriptions.push_back( "Herwig++ ATLAS MC09" );
	colours.push_back( kRed + 3 );
	styles.push_back( 8 );

	//Pythia 8 non diffractive
	truthPaths.push_back( "data/user.bwynne.LeadingJetModifiedv3.Pythia8NonDiff.MC.TruthJet/mergedFile.root" );
	recoPaths.push_back( "data/user.bwynne.LeadingJetModifiedv3.Pythia8NonDiff.MC.CaloJet/mergedFile.root" );
	descriptions.push_back( "PYTHIA8 Non-Diffractive" );
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

string MonteCarloInformation::TruthFilePath( unsigned int Index )
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

string MonteCarloInformation::ReconstructedFilePath( unsigned int Index )
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
