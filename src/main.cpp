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
#include <cstdlib>

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

	////////////////////////////////////////////////////////////
	//                                                        //
	// Load the data - Again, set this up yourself            //
	//                                                        //
	////////////////////////////////////////////////////////////
	//IFileInput * dataInput = new InputUETree( "/Disk/speyside7/Grid/grid-files/bwynne/L1_J5.v2/JetTauEtmiss/combined.root", "benTuple", "JetTauEtmiss Data (2010)", mcInfo->NumberOfSources() );
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
	XPlotMaker * nChargedTowardsPlot = new XPlotMaker( "NChargedTowards500", "PYTHIA", nChargeBins, nChargeMin, nChargeMax, 1.0, true );
	MonteCarloSummaryPlotMaker * nChargedTowardsSummary = new MonteCarloSummaryPlotMaker( nChargedTowardsPlot, mcInfo, COMBINE_MC );
	allPlotMakers.push_back( nChargedTowardsSummary );

	//2D Unfolding

	//Make a plot of number of charged particles in the toward region vs lead jet pT
	XvsYNormalisedPlotMaker * pTvsNChargedTowardPlot = new XvsYNormalisedPlotMaker( "LeadJetPt", "NChargedTowards500", "PYTHIA6 AMBT1",
			jetPtBins, jetPtMin, jetPtMax, nChargeBins, nChargeMin, nChargeMax, scaleFactor );
	MonteCarloSummaryPlotMaker * pTvsNChargedTowardSummary = new MonteCarloSummaryPlotMaker( pTvsNChargedTowardPlot, mcInfo, COMBINE_MC );
	pTvsNChargedTowardSummary->SetYRange( 0.1, 5.9 );
	pTvsNChargedTowardSummary->SetAxisLabels( "p_{T}^{lead} [GeV]", "<d^{2}N_{ch}/d#etad#phi>" );
	allPlotMakers.push_back( pTvsNChargedTowardSummary );

	//Make a plot of number of charged particles in the away region vs lead jet pT
	XvsYNormalisedPlotMaker * pTvsNChargedAwayPlot = new XvsYNormalisedPlotMaker( "LeadJetPt", "NChargedAway500", "PYTHIA",
			jetPtBins, jetPtMin, jetPtMax, nChargeBins, nChargeMin, nChargeMax, scaleFactor );
	MonteCarloSummaryPlotMaker * pTvsNChargedAwaySummary = new MonteCarloSummaryPlotMaker( pTvsNChargedAwayPlot, mcInfo, COMBINE_MC );
	pTvsNChargedAwaySummary->SetYRange( 0.1, 5.9 );
	pTvsNChargedAwaySummary->SetAxisLabels( "p_{T}^{lead} [GeV]", "<d^{2}N_{ch}/d#etad#phi>" );
	allPlotMakers.push_back( pTvsNChargedAwaySummary );

	//Make a plot of number of charged particles in the transverse region vs lead jet pT
	XvsYNormalisedPlotMaker * pTvsNChargedTransPlot = new XvsYNormalisedPlotMaker( "LeadJetPt", "NChargedTransverse500", "PYTHIA",
			jetPtBins, jetPtMin, jetPtMax, nChargeBins, nChargeMin, nChargeMax, scaleFactor );
	MonteCarloSummaryPlotMaker * pTvsNChargedTransSummary = new MonteCarloSummaryPlotMaker( pTvsNChargedTransPlot, mcInfo, COMBINE_MC );
	pTvsNChargedTransSummary->SetYRange( 0.1, 2.9 );
	pTvsNChargedTransSummary->SetAxisLabels( "p_{T}^{lead} [GeV]", "<d^{2}N_{ch}/d#etad#phi>" );
	allPlotMakers.push_back( pTvsNChargedTransSummary );

	//Make a plot of mean pT of charged particles in the toward region vs lead jet pT
	XvsYNormalisedPlotMaker * pTvsMeanPtTowardPlot = new XvsYNormalisedPlotMaker( "LeadJetPt", "MeanPtTowards500", "PYTHIA",
			jetPtBins, jetPtMin, jetPtMax, meanPtBins, meanPtMin, meanPtMax, scaleFactor );
	MonteCarloSummaryPlotMaker * pTvsMeanPtTowardSummary = new MonteCarloSummaryPlotMaker( pTvsMeanPtTowardPlot, mcInfo, COMBINE_MC );
	//	pTvsMeanPtTowardSummary->SetYRange( 0.1, 5.9 );
	pTvsMeanPtTowardSummary->SetAxisLabels( "p_{T}^{lead} [GeV]", "<d^{2}<p_{T}>/d#etad#phi>" );
	allPlotMakers.push_back( pTvsMeanPtTowardSummary );

	//Make a plot of mean pT of charged particles in the away region vs lead jet pT
	XvsYNormalisedPlotMaker * pTvsMeanPtAwayPlot = new XvsYNormalisedPlotMaker( "LeadJetPt", "MeanPtAway500", "PYTHIA",
			jetPtBins, jetPtMin, jetPtMax, meanPtBins, meanPtMin, meanPtMax, scaleFactor );
	MonteCarloSummaryPlotMaker * pTvsMeanPtAwaySummary = new MonteCarloSummaryPlotMaker( pTvsMeanPtAwayPlot, mcInfo, COMBINE_MC );
	//	pTvsMeanPtAwaySummary->SetYRange( 0.1, 5.9 );
	pTvsMeanPtAwaySummary->SetAxisLabels( "p_{T}^{lead} [GeV]", "<d^{2}<p_{T}>/d#etad#phi>" );
	allPlotMakers.push_back( pTvsMeanPtAwaySummary );

	//Make a plot of mean pT of charged particles in the transverse region vs lead jet pT
	XvsYNormalisedPlotMaker * pTvsMeanPtTransPlot = new XvsYNormalisedPlotMaker( "LeadJetPt", "MeanPtTransverse500", "PYTHIA",
			jetPtBins, jetPtMin, jetPtMax, meanPtBins, meanPtMin, meanPtMax, scaleFactor );
	MonteCarloSummaryPlotMaker * pTvsMeanPtTransSummary = new MonteCarloSummaryPlotMaker( pTvsMeanPtTransPlot, mcInfo, COMBINE_MC );
	//	pTvsMeanPtTransSummary->SetYRange( 0.1, 2.9 );
	pTvsMeanPtTransSummary->SetAxisLabels( "p_{T}^{lead} [GeV]", "<d^{2}<p_{T}>/d#etad#phi>" );
	allPlotMakers.push_back( pTvsMeanPtTransSummary );

	//Make a plot of sum pT of charged particles in the toward region vs lead jet pT
	XvsYNormalisedPlotMaker * pTvsPtSumTowardPlot = new XvsYNormalisedPlotMaker( "LeadJetPt", "PtSumTowards500", "PYTHIA",
			jetPtBins, jetPtMin, jetPtMax, sumPtBins, sumPtMin, sumPtMax, scaleFactor );
	MonteCarloSummaryPlotMaker * pTvsPtSumTowardSummary = new MonteCarloSummaryPlotMaker( pTvsPtSumTowardPlot, mcInfo, COMBINE_MC );
	//	pTvsPtSumTowardSummary->SetYRange( 0.1, 5.9 );
	pTvsPtSumTowardSummary->SetAxisLabels( "p_{T}^{lead} [GeV]", "<d^{2}#Sigmap_{T}/d#etad#phi>" );
	allPlotMakers.push_back( pTvsPtSumTowardSummary );

	//Make a plot of sum pT of charged particles in the away region vs lead jet pT
	XvsYNormalisedPlotMaker * pTvsPtSumAwayPlot = new XvsYNormalisedPlotMaker( "LeadJetPt", "PtSumAway500", "PYTHIA",
			jetPtBins, jetPtMin, jetPtMax, sumPtBins, sumPtMin, sumPtMax, scaleFactor );
	MonteCarloSummaryPlotMaker * pTvsPtSumAwaySummary = new MonteCarloSummaryPlotMaker( pTvsPtSumAwayPlot, mcInfo, COMBINE_MC );
	//	pTvsPtSumAwaySummary->SetYRange( 0.1, 5.9 );
	pTvsPtSumAwaySummary->SetAxisLabels( "p_{T}^{lead} [GeV]", "<d^{2}#Sigmap_{T}/d#etad#phi>" );
	allPlotMakers.push_back( pTvsPtSumAwaySummary );

	//Make a plot of sum pT of charged particles in the transverse region vs lead jet pT
	XvsYNormalisedPlotMaker * pTvsPtSumTransPlot = new XvsYNormalisedPlotMaker( "LeadJetPt", "PtSumTransverse500", "PYTHIA",
			jetPtBins, jetPtMin, jetPtMax, sumPtBins, sumPtMin, sumPtMax, scaleFactor );
	MonteCarloSummaryPlotMaker * pTvsPtSumTransSummary = new MonteCarloSummaryPlotMaker( pTvsPtSumTransPlot, mcInfo, COMBINE_MC );
	//	pTvsPtSumTransSummary->SetYRange( 0.1, 2.9 );
	pTvsPtSumTransSummary->SetAxisLabels( "p_{T}^{lead} [GeV]", "<d^{2}#Sigmap_{T}/d#etad#phi>" );
	allPlotMakers.push_back( pTvsPtSumTransSummary );

	//Make a plot of number of charged particles in the toward region vs lead jet pT
	XvsYNormalisedPlotMaker * pTmeanvsNChargedTowardPlot = new XvsYNormalisedPlotMaker( "NChargedTowards500", "MeanPtTowards500", "PYTHIA",
			XnChargeBins, XnChargeMin, XnChargeMax, meanPtBins, meanPtMin, meanPtMax, 1.0 );
	MonteCarloSummaryPlotMaker * pTmeanvsNChargedTowardSummary = new MonteCarloSummaryPlotMaker( pTmeanvsNChargedTowardPlot, mcInfo, COMBINE_MC );
	pTmeanvsNChargedTowardSummary->SetYRange( 1000.0, 5000.0 );
	pTmeanvsNChargedTowardSummary->SetAxisLabels( "N_{ch}", "<p_{T}>" );
	allPlotMakers.push_back( pTmeanvsNChargedTowardSummary );

	//Make a plot of number of charged particles in the away region vs lead jet pT
	XvsYNormalisedPlotMaker * pTmeanvsNChargedAwayPlot = new XvsYNormalisedPlotMaker( "NChargedAway500", "MeanPtAway500", "PYTHIA",
			XnChargeBins, XnChargeMin, XnChargeMax, meanPtBins, meanPtMin, meanPtMax, 1.0 );
	MonteCarloSummaryPlotMaker * pTmeanvsNChargedAwaySummary = new MonteCarloSummaryPlotMaker( pTmeanvsNChargedAwayPlot, mcInfo, COMBINE_MC );
	pTmeanvsNChargedAwaySummary->SetYRange( 1000.0, 2000.0 );
	pTmeanvsNChargedAwaySummary->SetAxisLabels( "N_{ch}", "<p_{T}>" );
	allPlotMakers.push_back( pTmeanvsNChargedAwaySummary );

	//Make a plot of number of charged particles in the transverse region vs lead jet pT
	XvsYNormalisedPlotMaker * pTmeanvsNChargedTransPlot = new XvsYNormalisedPlotMaker( "NChargedTransverse500", "MeanPtTransverse500", "PYTHIA",
			XnChargeBins, XnChargeMin, XnChargeMax, meanPtBins, meanPtMin, meanPtMax, 1.0 );
	MonteCarloSummaryPlotMaker * pTmeanvsNChargedTransSummary = new MonteCarloSummaryPlotMaker( pTmeanvsNChargedTransPlot, mcInfo, COMBINE_MC );
	pTmeanvsNChargedTransSummary->SetYRange( 800.0, 1600.0 );
	pTmeanvsNChargedTransSummary->SetAxisLabels( "N_{ch}", "<p_{T}>" );
	allPlotMakers.push_back( pTmeanvsNChargedTransSummary );

	/////////////////////////////////////////////////////////////
	//                                                         //
	// No further user edits required                          //
	//                                                         //
	/////////////////////////////////////////////////////////////

	//Populate the smearing matrices
	for ( unsigned int mcIndex = 0; mcIndex < mcInfo->NumberOfSources(); mcIndex++ )
	{
		IFileInput * truthInput = mcInfo->MakeTruthInput( mcIndex );
		IFileInput * reconstructedInput = mcInfo->MakeReconstructedInput( mcIndex );
		MakeSmearingMatrices( truthInput, reconstructedInput );
		delete truthInput;
		delete reconstructedInput;
	}

	//Unfold!
	IFileInput * dataInput = mcInfo->MakeReconstructedInput( 0 );
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
	cout << endl << "Loading " << *( TruthInput->Description() ) << " events" << endl;

	//Check that truth and reco comprise the same number of files
	if ( TruthInput->NumberOfFiles() != ReconstructedInput->NumberOfFiles() )
	{
		cerr << "Smearing matrix construction given " << TruthInput->NumberOfFiles() << " truth files and " << ReconstructedInput->NumberOfFiles() << " reconstructed" << endl;
		cerr << "These numbers must be the same" << endl;
		exit(1);
	}

	//Loop over each truth-reco file pair
	for ( unsigned int fileIndex = 0; fileIndex < TruthInput->NumberOfFiles(); fileIndex++ )
	{
		//Force loading the file, so that NumberOfRows is accurate
		TruthInput->ReadRow( 0, fileIndex );
		ReconstructedInput->ReadRow( 0, fileIndex );

		//Cache to store results of event number searches
		vector< bool > recoMatched( ReconstructedInput->NumberOfRows(), false );

		//Attempt to pair each truth event with a reco event
		for ( unsigned long truthIndex = 0; truthIndex < TruthInput->NumberOfRows(); truthIndex++ )
		{
			//Read the row, find the event number;
			if ( TruthInput->ReadRow( truthIndex, fileIndex ) )
			{
				UInt_t truthEventNumber = TruthInput->EventNumber();

				//Look for that event number in the reco file
				if ( ReconstructedInput->ReadEvent( truthEventNumber, fileIndex ) )
				{
					//Matched event

					//Add the match to all plot makers
					for ( unsigned int plotIndex = 0; plotIndex < allPlotMakers.size(); plotIndex++ )
					{
						allPlotMakers[ plotIndex ]->StoreMatch( TruthInput, ReconstructedInput );
					}

					//Count the match
					recoMatched[ ReconstructedInput->CurrentRow() ] = true;
					matchedEvents++;
				}
				else
				{
					//Missed event

					//Add the miss to all plot makers
					for ( unsigned int plotIndex = 0; plotIndex < allPlotMakers.size(); plotIndex++ )
					{
						allPlotMakers[ plotIndex ]->StoreMiss( TruthInput );
					}

					//Count the miss
					missedEvents++;
				}
			}
			else
			{
				cerr << "Stupidity fail" << endl;
				exit(1);
			}
		}

		//Any reconstructed events not paired are fake
		for ( unsigned long recoIndex = 0; recoIndex < ReconstructedInput->NumberOfRows(); recoIndex++ )
		{
			//Find the rows that weren't matched
			if ( !recoMatched[ recoIndex ] )
			{
				//Fake event

				//Read the event from disk
				if ( ReconstructedInput->ReadRow( recoIndex, fileIndex ) )
				{
					//Add the fake to all plotmakers
					for ( unsigned int plotIndex = 0; plotIndex < allPlotMakers.size(); plotIndex++ )
					{
						allPlotMakers[plotIndex]->StoreFake( ReconstructedInput );
					}

					//Record the fake
					fakeEvents++;
				}
				else
				{
					cerr << "Stupidity fail" << endl;
					exit(1);
				}
			}
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
	long dataTotal = 0;
	for ( unsigned int fileIndex = 0; fileIndex < DataInput->NumberOfFiles(); fileIndex++ )
	{
		//Force loading the file, so that NumberOfRows is accurate
		DataInput->ReadRow( 0, fileIndex );

		//Loop over each row in the file
		for ( unsigned long dataIndex = 0; dataIndex < DataInput->NumberOfRows(); dataIndex++ )
		{
			//Read the row from disk
			if ( DataInput->ReadRow( dataIndex, fileIndex ) )
			{
				//Store the row in all plot makers
				for ( unsigned int plotIndex = 0; plotIndex < allPlotMakers.size(); plotIndex++ )
				{
					allPlotMakers[ plotIndex ]->StoreData( DataInput );
				}

				//Count the data
				dataTotal++;
			}
			else
			{
				cerr << "Stupidity fail" << endl;
				exit(1);
			}
		}
	}
	cout << "Total: " << dataTotal << endl;
	delete DataInput;

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
