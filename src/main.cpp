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
const bool COMBINE_MC = false;

////////////////////////////////////////////////////////////
//                                                        //
// Set the error calculation mode:                        //
// 0) Bin-by-bin scaling of uncorrected errors (fast)     //
// 1) Variances calculated using D'Agostini's method      //
// 2) Full D'Agositini moethd covariance matrix (slow)    //
//                                                        //
////////////////////////////////////////////////////////////
const int ERROR_MODE = 0;

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
	IFileInput * dataInput = new InputUETree( "/Disk/speyside7/Grid/grid-files/bwynne/L1_J5.v2/JetTauEtmiss/combined.root", "benTuple", "JetTauEtmiss Data (2010)", mcInfo->NumberOfSources() );
	//IFileInput * dataInput = mcInfo->MakeReconstructedInput( 0 );

	////////////////////////////////////////////////////////////
	//                                                        //
	// Initialise the plot makers                             //
	//                                                        //
	// In each case, choose the type of plot you want to make //
	// Instantiate one of these                               //
	// Then put it into a MonteCarloSummaryPlotMaker          //
	//                                                        //
	////////////////////////////////////////////////////////////
	double scaleFactor = 3.0 / ( 10.0 * M_PI );
	int jetPtBins = 100;
	double jetPtMin = 0.0;
	double jetPtMax = 200000.0;
	int nChargeBins = 60;
	double nChargeMin = 0.5;
	double nChargeMax = 60.5;
	int XnChargeBins = 30;
	double XnChargeMin = 0.5;
	double XnChargeMax = 30.5;
	int meanPtBins = 1500;
	double meanPtMin = 0.0;
	double meanPtMax = 150000.0;
	int sumPtBins = 5000;
	double sumPtMin = 0.0;
	double sumPtMax = 5000000.0;

	//1D Unfolding

	//Lead jet pT
	XPlotMaker * leadJetPtPlot = new XPlotMaker( "LeadJetPt", "PYTHIA", jetPtBins, jetPtMin, jetPtMax, 1.0, true );
	MonteCarloSummaryPlotMaker * leadJetPtSummary = new MonteCarloSummaryPlotMaker( leadJetPtPlot, mcInfo, COMBINE_MC );
	leadJetPtSummary->SetYRange( 1E-7, 1.0 );
	leadJetPtSummary->UseLogScale();
	allPlotMakers.push_back( leadJetPtSummary );

	//N charge towards
	XPlotMaker * nChargedTowardsPlot = new XPlotMaker( "NChargedTowards", "PYTHIA", nChargeBins, nChargeMin, nChargeMax, 1.0, true );
	MonteCarloSummaryPlotMaker * nChargedTowardsSummary = new MonteCarloSummaryPlotMaker( nChargedTowardsPlot, mcInfo, COMBINE_MC );
	allPlotMakers.push_back( nChargedTowardsSummary );

	//2D Unfolding

	//Make a plot of number of charged particles in the toward region vs lead jet pT
	XvsYNormalisedPlotMaker * pTvsNChargedTowardPlot = new XvsYNormalisedPlotMaker( "LeadJetPt", "NChargedTowards", "PYTHIA",
			jetPtBins, jetPtMin, jetPtMax, nChargeBins, nChargeMin, nChargeMax, scaleFactor );
	MonteCarloSummaryPlotMaker * pTvsNChargedTowardSummary = new MonteCarloSummaryPlotMaker( pTvsNChargedTowardPlot, mcInfo, COMBINE_MC );
	pTvsNChargedTowardSummary->SetYRange( 0.1, 5.9 );
	pTvsNChargedTowardSummary->SetAxisLabels( "p_{T}^{lead} [GeV]", "<d^{2}N_{ch}/d#etad#phi>" );
	allPlotMakers.push_back( pTvsNChargedTowardSummary );

	//Make a plot of number of charged particles in the away region vs lead jet pT
	XvsYNormalisedPlotMaker * pTvsNChargedAwayPlot = new XvsYNormalisedPlotMaker( "LeadJetPt", "NChargedAway", "PYTHIA",
			jetPtBins, jetPtMin, jetPtMax, nChargeBins, nChargeMin, nChargeMax, scaleFactor );
	MonteCarloSummaryPlotMaker * pTvsNChargedAwaySummary = new MonteCarloSummaryPlotMaker( pTvsNChargedAwayPlot, mcInfo, COMBINE_MC );
	pTvsNChargedAwaySummary->SetYRange( 0.1, 5.9 );
	pTvsNChargedAwaySummary->SetAxisLabels( "p_{T}^{lead} [GeV]", "<d^{2}N_{ch}/d#etad#phi>" );
	allPlotMakers.push_back( pTvsNChargedAwaySummary );

	//Make a plot of number of charged particles in the transverse region vs lead jet pT
	XvsYNormalisedPlotMaker * pTvsNChargedTransPlot = new XvsYNormalisedPlotMaker( "LeadJetPt", "NChargedTransverse", "PYTHIA",
			jetPtBins, jetPtMin, jetPtMax, nChargeBins, nChargeMin, nChargeMax, scaleFactor );
	MonteCarloSummaryPlotMaker * pTvsNChargedTransSummary = new MonteCarloSummaryPlotMaker( pTvsNChargedTransPlot, mcInfo, COMBINE_MC );
	pTvsNChargedTransSummary->SetYRange( 0.1, 2.9 );
	pTvsNChargedTransSummary->SetAxisLabels( "p_{T}^{lead} [GeV]", "<d^{2}N_{ch}/d#etad#phi>" );
	allPlotMakers.push_back( pTvsNChargedTransSummary );

	//Make a plot of mean pT of charged particles in the toward region vs lead jet pT
	XvsYNormalisedPlotMaker * pTvsMeanPtTowardPlot = new XvsYNormalisedPlotMaker( "LeadJetPt", "MeanPtTowards", "PYTHIA",
			jetPtBins, jetPtMin, jetPtMax, meanPtBins, meanPtMin, meanPtMax, scaleFactor );
	MonteCarloSummaryPlotMaker * pTvsMeanPtTowardSummary = new MonteCarloSummaryPlotMaker( pTvsMeanPtTowardPlot, mcInfo, COMBINE_MC );
//	pTvsMeanPtTowardSummary->SetYRange( 0.1, 5.9 );
	pTvsMeanPtTowardSummary->SetAxisLabels( "p_{T}^{lead} [GeV]", "<d^{2}<p_{T}>/d#etad#phi>" );
	allPlotMakers.push_back( pTvsMeanPtTowardSummary );

	//Make a plot of mean pT of charged particles in the away region vs lead jet pT
	XvsYNormalisedPlotMaker * pTvsMeanPtAwayPlot = new XvsYNormalisedPlotMaker( "LeadJetPt", "MeanPtAway", "PYTHIA",
			jetPtBins, jetPtMin, jetPtMax, meanPtBins, meanPtMin, meanPtMax, scaleFactor );
	MonteCarloSummaryPlotMaker * pTvsMeanPtAwaySummary = new MonteCarloSummaryPlotMaker( pTvsMeanPtAwayPlot, mcInfo, COMBINE_MC );
//	pTvsMeanPtAwaySummary->SetYRange( 0.1, 5.9 );
	pTvsMeanPtAwaySummary->SetAxisLabels( "p_{T}^{lead} [GeV]", "<d^{2}<p_{T}>/d#etad#phi>" );
	allPlotMakers.push_back( pTvsMeanPtAwaySummary );

	//Make a plot of mean pT of charged particles in the transverse region vs lead jet pT
	XvsYNormalisedPlotMaker * pTvsMeanPtTransPlot = new XvsYNormalisedPlotMaker( "LeadJetPt", "MeanPtTransverse", "PYTHIA",
			jetPtBins, jetPtMin, jetPtMax, meanPtBins, meanPtMin, meanPtMax, scaleFactor );
	MonteCarloSummaryPlotMaker * pTvsMeanPtTransSummary = new MonteCarloSummaryPlotMaker( pTvsMeanPtTransPlot, mcInfo, COMBINE_MC );
//	pTvsMeanPtTransSummary->SetYRange( 0.1, 2.9 );
	pTvsMeanPtTransSummary->SetAxisLabels( "p_{T}^{lead} [GeV]", "<d^{2}<p_{T}>/d#etad#phi>" );
	allPlotMakers.push_back( pTvsMeanPtTransSummary );

	//Make a plot of sum pT of charged particles in the toward region vs lead jet pT
	XvsYNormalisedPlotMaker * pTvsPtSumTowardPlot = new XvsYNormalisedPlotMaker( "LeadJetPt", "PtSumTowards", "PYTHIA",
			jetPtBins, jetPtMin, jetPtMax, sumPtBins, sumPtMin, sumPtMax, scaleFactor );
	MonteCarloSummaryPlotMaker * pTvsPtSumTowardSummary = new MonteCarloSummaryPlotMaker( pTvsPtSumTowardPlot, mcInfo, COMBINE_MC );
//	pTvsPtSumTowardSummary->SetYRange( 0.1, 5.9 );
	pTvsPtSumTowardSummary->SetAxisLabels( "p_{T}^{lead} [GeV]", "<d^{2}#Sigmap_{T}/d#etad#phi>" );
	allPlotMakers.push_back( pTvsPtSumTowardSummary );

	//Make a plot of sum pT of charged particles in the away region vs lead jet pT
	XvsYNormalisedPlotMaker * pTvsPtSumAwayPlot = new XvsYNormalisedPlotMaker( "LeadJetPt", "PtSumAway", "PYTHIA",
			jetPtBins, jetPtMin, jetPtMax, sumPtBins, sumPtMin, sumPtMax, scaleFactor );
	MonteCarloSummaryPlotMaker * pTvsPtSumAwaySummary = new MonteCarloSummaryPlotMaker( pTvsPtSumAwayPlot, mcInfo, COMBINE_MC );
//	pTvsPtSumAwaySummary->SetYRange( 0.1, 5.9 );
	pTvsPtSumAwaySummary->SetAxisLabels( "p_{T}^{lead} [GeV]", "<d^{2}#Sigmap_{T}/d#etad#phi>" );
	allPlotMakers.push_back( pTvsPtSumAwaySummary );

	//Make a plot of sum pT of charged particles in the transverse region vs lead jet pT
	XvsYNormalisedPlotMaker * pTvsPtSumTransPlot = new XvsYNormalisedPlotMaker( "LeadJetPt", "PtSumTransverse", "PYTHIA",
			jetPtBins, jetPtMin, jetPtMax, sumPtBins, sumPtMin, sumPtMax, scaleFactor );
	MonteCarloSummaryPlotMaker * pTvsPtSumTransSummary = new MonteCarloSummaryPlotMaker( pTvsPtSumTransPlot, mcInfo, COMBINE_MC );
//	pTvsPtSumTransSummary->SetYRange( 0.1, 2.9 );
	pTvsPtSumTransSummary->SetAxisLabels( "p_{T}^{lead} [GeV]", "<d^{2}#Sigmap_{T}/d#etad#phi>" );
	allPlotMakers.push_back( pTvsPtSumTransSummary );

	//Make a plot of number of charged particles in the toward region vs lead jet pT
	XvsYNormalisedPlotMaker * pTmeanvsNChargedTowardPlot = new XvsYNormalisedPlotMaker( "NChargedTowards", "MeanPtTowards", "PYTHIA",
			XnChargeBins, XnChargeMin, XnChargeMax, meanPtBins, meanPtMin, meanPtMax, 1.0 );
	MonteCarloSummaryPlotMaker * pTmeanvsNChargedTowardSummary = new MonteCarloSummaryPlotMaker( pTmeanvsNChargedTowardPlot, mcInfo, COMBINE_MC );
	pTmeanvsNChargedTowardSummary->SetYRange( 1000.0, 5000.0 );
	pTmeanvsNChargedTowardSummary->SetAxisLabels( "N_{ch}", "<p_{T}>" );
	allPlotMakers.push_back( pTmeanvsNChargedTowardSummary );

	//Make a plot of number of charged particles in the away region vs lead jet pT
	XvsYNormalisedPlotMaker * pTmeanvsNChargedAwayPlot = new XvsYNormalisedPlotMaker( "NChargedAway", "MeanPtAway", "PYTHIA",
			XnChargeBins, XnChargeMin, XnChargeMax, meanPtBins, meanPtMin, meanPtMax, 1.0 );
	MonteCarloSummaryPlotMaker * pTmeanvsNChargedAwaySummary = new MonteCarloSummaryPlotMaker( pTmeanvsNChargedAwayPlot, mcInfo, COMBINE_MC );
	pTmeanvsNChargedAwaySummary->SetYRange( 1000.0, 2000.0 );
	pTmeanvsNChargedAwaySummary->SetAxisLabels( "N_{ch}", "<p_{T}>" );
	allPlotMakers.push_back( pTmeanvsNChargedAwaySummary );

	//Make a plot of number of charged particles in the transverse region vs lead jet pT
	XvsYNormalisedPlotMaker * pTmeanvsNChargedTransPlot = new XvsYNormalisedPlotMaker( "NChargedTransverse", "MeanPtTransverse", "PYTHIA",
			XnChargeBins, XnChargeMin, XnChargeMax, meanPtBins, meanPtMin, meanPtMax, 1.0 );
	MonteCarloSummaryPlotMaker * pTmeanvsNChargedTransSummary = new MonteCarloSummaryPlotMaker( pTmeanvsNChargedTransPlot, mcInfo, COMBINE_MC );
	pTmeanvsNChargedTransSummary->SetYRange( 800.0, 1600.0 );
	pTmeanvsNChargedTransSummary->SetAxisLabels( "N_{ch}", "<p_{T}>" );
	allPlotMakers.push_back( pTmeanvsNChargedTransSummary );

	//Make a plot of mean pT of charged particles in the toward region vs lead jet pT
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
	for ( unsigned long truthIndex = 0; truthIndex < TruthInput->NumberOfRows(); truthIndex++ )
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
	for ( unsigned long reconstructedIndex = 0; reconstructedIndex < ReconstructedInput->NumberOfRows(); reconstructedIndex++ )
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
	for ( unsigned long rowIndex = 0; rowIndex < DataInput->NumberOfRows(); rowIndex++ )
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
