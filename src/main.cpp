/**
  @file main.cpp

  The entry point for Imagiro

  @author Benjamin M Wynne
  @date 08-03-2011
 */

#include "IFileInput.h"
#include "CombinedFileInput.h"
#include "InputUETree.h"
#include "XvsYNormalisedPlotMaker.h"
#include "XvsYNormalisedFolding.h"
#include "MonteCarloSummaryPlotMaker.h"
#include "MonteCarloInformation.h"
#include "XPlotMaker.h"
#include "XFolding.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include <ctime>
#include <iostream>
#include <string>
#include <cmath>

using namespace std;

//Method declarations
void MakeSmearingMatrices( IFileInput * TruthInput, IFileInput * ReconstructedInput );
void DoTheUnfolding( IFileInput * DataInput );

//The plotmakers
vector< MonteCarloSummaryPlotMaker* > allPlotMakers;

//Time
time_t timeNow;

////////////////////////////////////////////////////////////
//                                                        //
// Combine all MC into one smearing matrix,               //
// or use separate matrices? (recommend true)             //
//                                                        //
////////////////////////////////////////////////////////////
const bool COMBINE_MC = true;

////////////////////////////////////////////////////////////
//                                                        //
// Set the error calculation mode:                        //
// 0) Bin-by-bin scaling of uncorrected errors (fast)     //
// 1) Variances calculated using D'Agostini's method      //
// 2) Full D'Agositini moethd covariance matrix (slow)    //
//                                                        //
////////////////////////////////////////////////////////////
const int ERROR_MODE = 2;

////////////////////////////////////////////////////////////
//                                                        //
// Set whether to smooth the prior distribution           //
// before applying it each iteration                      //
// (recommend false)                                      //
//                                                        //
////////////////////////////////////////////////////////////
const bool WITH_SMOOTHING = false;

////////////////////////////////////////////////////////////
//                                                        //
// Set the output file name                               //
//                                                        //
////////////////////////////////////////////////////////////
const string OUTPUT_FILE_NAME = "UnfoldedFinal.Data.root";

int main ( int argc, char * argv[] )
{
	//Status message
	time(&timeNow);
	cout << endl << "Imagiro started: " << ctime( &timeNow ) << endl;

	////////////////////////////////////////////////////////////
	//                                                        //
	// Load the MC - Set this up in MonteCarloInformation.cpp //
	//                                                        //
	////////////////////////////////////////////////////////////
	MonteCarloInformation * mcInfo = new MonteCarloInformation();
	vector< IFileInput* > truthInputs, reconstructedInputs;
	for ( int sourceIndex = 0; sourceIndex < mcInfo->NumberOfSources(); sourceIndex++ )
	{
		truthInputs.push_back( mcInfo->MakeTruthInput( sourceIndex ) );
		reconstructedInputs.push_back( mcInfo->MakeReconstructedInput( sourceIndex ) );
	}

	////////////////////////////////////////////////////////////
	//                                                        //
	// Load the data - Again, set this up yourself            //
	//                                                        //
	////////////////////////////////////////////////////////////
	//IFileInput * dataInput = new InputUETree( "/media/My Book/UEdata/L1_J5.v1/combinedJetTauEtmiss/combinedJetTauEtmiss.root", "benTuple", "Data D3PD", mcInfo->NumberOfSources() );
	//IFileInput * dataInput = new InputUETree( "/media/My Book/UEdata/L1_EM5.v1/D3PDjetEtMissData/combinedData.root", "benTuple", "Data D3PD", mcInfo->NumberOfSources() );
	int sourceIndex = 0;
	IFileInput * dataInput = mcInfo->MakeReconstructedInput( sourceIndex );

	////////////////////////////////////////////////////////////
	//                                                        //
	// Initialise the plot makers                             //
	//                                                        //
	// In each case, choose the type of plot you want to make //
	// Instantiate one of these                               //
	// Then put it into a MonteCarloSummaryPlotMaker          //
	//                                                        //
	////////////////////////////////////////////////////////////
/*	double scaleFactor = 3.0 / ( 10.0 * M_PI );
	int jetPtBins = 100;
	double jetPtMin = 0.0;
	double jetPtMax = 200000.0;
	int nChargeBins = 500;
	double nChargeMin = 0.5;
	double nChargeMax = 500.5;

	//1D Unfolding

	//Lead jet pT
	XPlotMaker * leadJetPtPlot = new XPlotMaker( "LeadJetPt", "PYTHIA", jetPtBins, jetPtMin, jetPtMax, 1.0, true );
	MonteCarloSummaryPlotMaker * leadJetPtSummary = new MonteCarloSummaryPlotMaker( leadJetPtPlot, mcInfo, COMBINE_MC );
	leadJetPtSummary->SetYRange( 0.0, 0.4 );
	allPlotMakers.push_back( leadJetPtSummary );

	//N charge towards
	XPlotMaker * nChargedTowardsPlot = new XPlotMaker( "NChargedTowards", "PYTHIA", nChargeBins, nChargeMin, nChargeMax, 1.0 );
	MonteCarloSummaryPlotMaker * nChargedTowardsSummary = new MonteCarloSummaryPlotMaker( nChargedTowardsPlot, mcInfo, COMBINE_MC );
	allPlotMakers.push_back( nChargedTowardsSummary );

	//2D Unfolding

	//Make a plot of number of charged particles in the toward region vs lead jet pT
	XvsYNormalisedPlotMaker * pTvsNChargedTowardPlot = new XvsYNormalisedPlotMaker( "LeadJetPt", "NChargedTowards", "PYTHIA",
			jetPtBins, jetPtMin, jetPtMax, nChargeBins, nChargeMin, nChargeMax, scaleFactor );
	MonteCarloSummaryPlotMaker * pTvsNChargedTowardSummary = new MonteCarloSummaryPlotMaker( pTvsNChargedTowardPlot, mcInfo, COMBINE_MC );
	pTvsNChargedTowardSummary->SetYRange( 0.1, 2.9 );
	pTvsNChargedTowardSummary->SetAxisLabels( "p_{T}^{lead} [GeV]", "<d^{2}N_{ch}/d#etad#phi>" );
	allPlotMakers.push_back( pTvsNChargedTowardSummary );

	//Make a plot of number of charged particles in the away region vs lead jet pT
	XvsYNormalisedPlotMaker * pTvsNChargedAwayPlot = new XvsYNormalisedPlotMaker( "LeadJetPt", "NChargedAway", "PYTHIA",
			jetPtBins, jetPtMin, jetPtMax, nChargeBins, nChargeMin, nChargeMax, scaleFactor );
	MonteCarloSummaryPlotMaker * pTvsNChargedAwaySummary = new MonteCarloSummaryPlotMaker( pTvsNChargedAwayPlot, mcInfo, COMBINE_MC );
	pTvsNChargedAwaySummary->SetYRange( 0.1, 2.9 );
	pTvsNChargedAwaySummary->SetAxisLabels( "p_{T}^{lead} [GeV]", "<d^{2}N_{ch}/d#etad#phi>" );
	allPlotMakers.push_back( pTvsNChargedAwaySummary );

	//Make a plot of number of charged particles in the transverse region vs lead jet pT
	XvsYNormalisedPlotMaker * pTvsNChargedTransPlot = new XvsYNormalisedPlotMaker( "LeadJetPt", "NChargedTransverse", "PYTHIA",
			jetPtBins, jetPtMin, jetPtMax, nChargeBins, nChargeMin, nChargeMax, scaleFactor );
	MonteCarloSummaryPlotMaker * pTvsNChargedTransSummary = new MonteCarloSummaryPlotMaker( pTvsNChargedTransPlot, mcInfo, COMBINE_MC );
	pTvsNChargedTransSummary->SetYRange( 0.1, 2.9 );
	pTvsNChargedTransSummary->SetAxisLabels( "p_{T}^{lead} [GeV]", "<d^{2}N_{ch}/d#etad#phi>" );
	allPlotMakers.push_back( pTvsNChargedTransSummary );*/

	double scaleFactor = 3.0 / ( 10.0 * M_PI );
	int jetPtBins = 30;
	double jetPtMin = 20.0;
	double jetPtMax = 50.0;
	int nChargeBins = 50;
	double nChargeMin = 0.5;
	double nChargeMax = 50.5;

/*	//2D Unfolding

	//Make a plot of number of charged particles in the toward region vs lead jet pT
	XvsYNormalisedPlotMaker * pTvsNChargedTowardPlot = new XvsYNormalisedPlotMaker( "MaxJetPt", "NChargeToward", "Pythia6",
			jetPtBins, jetPtMin, jetPtMax, nChargeBins, nChargeMin, nChargeMax, scaleFactor );
	MonteCarloSummaryPlotMaker * pTvsNChargedTowardSummary = new MonteCarloSummaryPlotMaker( pTvsNChargedTowardPlot, mcInfo, COMBINE_MC );
	pTvsNChargedTowardSummary->SetYRange( 0.1, 2.9 );
	pTvsNChargedTowardSummary->SetAxisLabels( "p_{T}^{lead} [GeV]", "<d^{2}N_{ch}/d#etad#phi>" );
	allPlotMakers.push_back( pTvsNChargedTowardSummary );

	//Make a plot of number of charged particles in the away region vs lead jet pT
	XvsYNormalisedPlotMaker * pTvsNChargedAwayPlot = new XvsYNormalisedPlotMaker( "MaxJetPt", "NChargeAway", "Pythia6",
			jetPtBins, jetPtMin, jetPtMax, nChargeBins, nChargeMin, nChargeMax, scaleFactor );
	MonteCarloSummaryPlotMaker * pTvsNChargedAwaySummary = new MonteCarloSummaryPlotMaker( pTvsNChargedAwayPlot, mcInfo, COMBINE_MC );
	pTvsNChargedAwaySummary->SetYRange( 0.1, 2.9 );
	pTvsNChargedAwaySummary->SetAxisLabels( "p_{T}^{lead} [GeV]", "<d^{2}N_{ch}/d#etad#phi>" );
	allPlotMakers.push_back( pTvsNChargedAwaySummary );

	//Make a plot of number of charged particles in the transverse region vs lead jet pT
	XvsYNormalisedPlotMaker * pTvsNChargedTransPlot = new XvsYNormalisedPlotMaker( "MaxJetPt", "NChargeTrans", "Pythia6",
			jetPtBins, jetPtMin, jetPtMax, nChargeBins, nChargeMin, nChargeMax, scaleFactor );
	MonteCarloSummaryPlotMaker * pTvsNChargedTransSummary = new MonteCarloSummaryPlotMaker( pTvsNChargedTransPlot, mcInfo, COMBINE_MC );
	pTvsNChargedTransSummary->SetYRange( 0.1, 2.9 );
	pTvsNChargedTransSummary->SetAxisLabels( "p_{T}^{lead} [GeV]", "<d^{2}N_{ch}/d#etad#phi>" );
	allPlotMakers.push_back( pTvsNChargedTransSummary );*/

	//1D Unfolding

	//Make a plot of lead jet pT in each event
	XPlotMaker * pTPlot = new XPlotMaker( "MaxJetPt", "Pythia6", jetPtBins, jetPtMin, jetPtMax, 1.0, true );
	MonteCarloSummaryPlotMaker * pTSummary = new MonteCarloSummaryPlotMaker( pTPlot, mcInfo, COMBINE_MC );
	pTSummary->SetYRange( 0.01, 0.19 );
	pTSummary->SetAxisLabels( "p_{T}^{lead} [GeV]", "Events" );
	pTSummary->UseLogScale();
	allPlotMakers.push_back( pTSummary );

	//Make a plot of number of charged particles in the toward region in each event
	XPlotMaker * nChargedTowardPlot = new XPlotMaker( "NChargeToward", "Pythia6", nChargeBins, nChargeMin, nChargeMax, 1.0, true );
	MonteCarloSummaryPlotMaker * nChargeTowardSummary = new MonteCarloSummaryPlotMaker( nChargedTowardPlot, mcInfo, COMBINE_MC );
	nChargeTowardSummary->SetYRange( 0.01, 0.09 );
	nChargeTowardSummary->SetAxisLabels( "N_{ch}", "Events" );
	allPlotMakers.push_back( nChargeTowardSummary );

	//Make a plot of number of charged particles in the away region in each event
	XPlotMaker * nChargedAwayPlot = new XPlotMaker( "NChargeAway", "Pythia6", nChargeBins, nChargeMin, nChargeMax, 1.0, true );
	MonteCarloSummaryPlotMaker * nChargeAwaySummary = new MonteCarloSummaryPlotMaker( nChargedAwayPlot, mcInfo, COMBINE_MC );
	nChargeAwaySummary->SetYRange( 0.01, 0.09 );
	nChargeAwaySummary->SetAxisLabels( "N_{ch}", "Events" );
	allPlotMakers.push_back( nChargeAwaySummary );

	//Make a plot of number of charged particles in the transverse region in each event
	XPlotMaker * nChargedTransPlot = new XPlotMaker( "NChargeTrans", "Pythia6", nChargeBins, nChargeMin, nChargeMax, 1.0, true );
	MonteCarloSummaryPlotMaker * nChargeTransSummary = new MonteCarloSummaryPlotMaker( nChargedTransPlot, mcInfo, COMBINE_MC );
	nChargeTransSummary->SetYRange( 0.01, 0.09 );
	nChargeTransSummary->SetAxisLabels( "N_{ch}", "Events" );
	allPlotMakers.push_back( nChargeTransSummary );

	/////////////////////////////////////////////////////////////
	//                                                         //
	// No further user edits required                          //
	//                                                         //
	/////////////////////////////////////////////////////////////

	//Populate the smearing matrices
	for ( unsigned int mcIndex = 0; mcIndex < truthInputs.size(); mcIndex++ )
	{
		MakeSmearingMatrices( truthInputs[ mcIndex ], reconstructedInputs[ mcIndex ] );
	}

	//Unfold!
	DoTheUnfolding( dataInput );

	//Status message
	time(&timeNow);
	cout << endl << "Imagiro finished: " << ctime( &timeNow ) << endl;
}

//Match up event numbers between truth and reco inputs
void MakeSmearingMatrices( IFileInput * TruthInput, IFileInput * ReconstructedInput )
{
	//Initialise
	long matchedEvents = 0;
	long fakeEvents = 0;
	long missedEvents = 0;
	vector<bool> recoWasMatched( ReconstructedInput->NumberOfRows(), false );
	cout << endl << "Loading " << *( TruthInput->Description() ) << " events" << endl;

	//Loop over all truth events to try and match to reconstructed events
	for ( long truthIndex = 0; truthIndex < TruthInput->NumberOfRows(); truthIndex++ )
	{
		//Get the event number for each truth event
		TruthInput->ReadRow( truthIndex );
		float truthEventNumber = TruthInput->EventNumber();
		long truthEventFile = TruthInput->CurrentFile();

		//Try to find a reconstructed event with that number
		if ( ReconstructedInput->ReadEvent( truthEventNumber, truthEventFile ) )
		{
			//Add the match to all plot makers
			for ( unsigned int plotIndex = 0; plotIndex < allPlotMakers.size(); plotIndex++ )
			{
				allPlotMakers[ plotIndex ]->StoreMatch( TruthInput, ReconstructedInput );
			}

			//Record the match
			recoWasMatched[ ReconstructedInput->CurrentRow() ] = true;
			matchedEvents++;
		}
		else
		{
			//Add the miss to all plot makers
			for ( unsigned int plotIndex = 0; plotIndex < allPlotMakers.size(); plotIndex++ )
			{
				allPlotMakers[ plotIndex ]->StoreMiss( TruthInput );
			}

			//Record the miss
			missedEvents++;
		}
	}

	//Loop over all reconstructed events looking for fakes
	for ( long reconstructedIndex = 0; reconstructedIndex < ReconstructedInput->NumberOfRows(); reconstructedIndex++ )
	{
		//Use the cached search success/fail rather than search again (half the disk io: big time saver)
		if ( !recoWasMatched[ reconstructedIndex ] )
		{
			//Update the line being accessed
			ReconstructedInput->ReadRow( reconstructedIndex );

			//Add the fake to all plotmakers
			for ( unsigned int plotIndex = 0; plotIndex < allPlotMakers.size(); plotIndex++ )
			{
				allPlotMakers[plotIndex]->StoreFake( ReconstructedInput );
			}

			//Record the fake
			fakeEvents++;
		}
	}

	//Debug
	cout << "Matched: " << matchedEvents << endl;
	cout << "Missed: " << missedEvents << endl;
	cout << "Fake: " << fakeEvents << endl;
}

void DoTheUnfolding( IFileInput * DataInput )
{
	//Populate the data distribution
	cout << endl << "Loading " << *( DataInput->Description() ) << " events" << endl;
	for ( long rowIndex = 0; rowIndex < DataInput->NumberOfRows(); rowIndex++ )
	{
		//Load the row into memory
		DataInput->ReadRow( rowIndex );

		//Store the row in all plot makers
		for ( unsigned int plotIndex = 0; plotIndex < allPlotMakers.size(); plotIndex++ )
		{
			allPlotMakers[ plotIndex ]->StoreData( DataInput );
		}
	}
	cout << "Total: " << DataInput->NumberOfRows() << endl;

	//Status message
	time(&timeNow);
	cout << endl << "Loading finished: " << ctime( &timeNow ) << endl;

	//Produce the plots
	TFile * OutputFile = new TFile( OUTPUT_FILE_NAME.c_str(), "RECREATE" );
	for ( unsigned int plotIndex = 0; plotIndex < allPlotMakers.size(); plotIndex++ )
	{
		//Unfold the data
		allPlotMakers[ plotIndex ]->Process( ERROR_MODE, WITH_SMOOTHING );

		//Write the result to file
		allPlotMakers[ plotIndex ]->SaveResult( OutputFile );

		//Free up some memory
		delete allPlotMakers[ plotIndex ];
	}
	OutputFile->Close();
}
