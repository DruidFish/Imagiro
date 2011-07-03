/**
  @file main.cpp

  The entry point for Imagiro

  @author Benjamin M Wynne
  @date 08-03-2011
 */

#include "IFileInput.h"
#include "CombinedFileInput.h"
#include "TriggerChoosingInput.h"
#include "InputUETree.h"
#include "XvsYNormalisedPlotMaker.h"
#include "XvsYNormalisedFolding.h"
#include "XsubEventPlotMaker.h"
#include "MonteCarloSummaryPlotMaker.h"
#include "MonteCarloInformation.h"
#include "ObservableList.h"
#include "XPlotMaker.h"
#include "XFolding.h"
#include "TFile.h"
#include "TROOT.h"
#include "TStyle.h"
#include <ctime>
#include <fstream>
#include <unistd.h>
#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>

using namespace std;

//Method declarations
void MakeSmearingMatrices( IFileInput * TruthInput, IFileInput * ReconstructedInput );
void DoTheUnfolding( IFileInput * DataInput );
void TimeAndMemory();
TStyle * PlotStyle( string StyleName );

//The plotmakers
vector< MonteCarloSummaryPlotMaker* > allPlotMakers;

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
const unsigned int ERROR_MODE = 0;

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
	cout << endl << "Imagiro started" << endl;
	TimeAndMemory();

	////////////////////////////////////////////////////////////
	//                                                        //
	// Load the MC - Set this up in MonteCarloInformation.cpp //
	//                                                        //
	////////////////////////////////////////////////////////////
	MonteCarloInformation * mcInfo = new MonteCarloInformation();

	////////////////////////////////////////////////////////////
	//                                                        //
	// Set the plot style - customise this below              //
	//                                                        //
	////////////////////////////////////////////////////////////
	TStyle * rootPlotStyle = PlotStyle( "ATLAS" );
        gROOT->SetStyle( "ATLAS" );
	gROOT->ForceStyle();

	////////////////////////////////////////////////////////////
	//                                                        //
	// Initialise the plot makers                             //
	//                                                        //
	// In each case, choose the type of plot you want to make //
	// Instantiate one of these                               //
	// Then put it into a MonteCarloSummaryPlotMaker          //
	//                                                        //
	////////////////////////////////////////////////////////////

	unsigned int phiBins = 100;
	double phiMin = -M_PI;
	double phiMax = M_PI;
	//Phi distribution
	XsubEventPlotMaker * phiPlot = new XsubEventPlotMaker( "Track500Phis", "PYTHIA6 AMBT1", phiBins, phiMin, phiMax, vector< string >( 1, "Track500Pts" ), 1.0, true );
	MonteCarloSummaryPlotMaker * phiSummary = new MonteCarloSummaryPlotMaker( phiPlot, mcInfo, COMBINE_MC );
	allPlotMakers.push_back( phiSummary );

	/*double scaleFactor = 3.0 / ( 10.0 * M_PI );
	unsigned int jetPtBins = 30;
	double jetPtMin = 0.0;
	double jetPtMax = 200000.0;
	unsigned int nChargeBins = 60;
	double nChargeMin = 0.5;
	double nChargeMax = 60.5;
	unsigned int XnChargeBins = 30;
	double XnChargeMin = 0.5;
	double XnChargeMax = 30.5;
	unsigned int meanPtBins = 200;
	double meanPtMin = 0.0;
	double meanPtMax = 10000.0;
	unsigned int sumPtBins = 100;
	double sumPtMin = 0.0;
	double sumPtMax = 100000.0;

	vector< double > jetPtBinEdges;
	for ( unsigned int binIndex = 0; binIndex <= jetPtBins; binIndex++ )
	{
		jetPtBinEdges.push_back( jetPtMin + ( ( jetPtMax - jetPtMin ) * (double)binIndex / (double) jetPtBins ) );
	}
	vector< double > nChargeBinEdges;
	for ( unsigned int binIndex = 0; binIndex <= nChargeBins; binIndex++ )
	{
		nChargeBinEdges.push_back( nChargeMin + ( ( nChargeMax - nChargeMin ) * (double)binIndex / (double) nChargeBins ) );
	}
	vector< double > XnChargeBinEdges;
	for ( unsigned int binIndex = 0; binIndex <= XnChargeBins; binIndex++ )
	{
		XnChargeBinEdges.push_back( XnChargeMin + ( ( XnChargeMax - XnChargeMin ) * (double)binIndex / (double) XnChargeBins ) );
	}
	vector< double > meanPtBinEdges;
	for ( unsigned int binIndex = 0; binIndex <= meanPtBins; binIndex++ )
	{
		meanPtBinEdges.push_back( meanPtMin + ( ( meanPtMax - meanPtMin ) * (double)binIndex / (double) meanPtBins ) );
	}
	vector< double > sumPtBinEdges;
	for ( unsigned int binIndex = 0; binIndex <= sumPtBins; binIndex++ )
	{
		sumPtBinEdges.push_back( sumPtMin + ( ( sumPtMax - sumPtMin ) * (double)binIndex / (double) sumPtBins ) );
	}

	//1D Unfolding

	//Lead jet pT
	XPlotMaker * leadJetPtPlot = new XPlotMaker( "LeadJetPt", "PYTHIA6-AMBT1", jetPtBinEdges, 1.0, true );
	MonteCarloSummaryPlotMaker * leadJetPtSummary = new MonteCarloSummaryPlotMaker( leadJetPtPlot, mcInfo, COMBINE_MC );
	leadJetPtSummary->SetYRange( 1E-13, 1.0 );
	leadJetPtSummary->UseLogScale();
	allPlotMakers.push_back( leadJetPtSummary );

	//N charge towards
	XPlotMaker * nChargedTowardsPlot = new XPlotMaker( "NChargedTowards500", "PYTHIA6-AMBT1", nChargeBinEdges, 1.0, true );
	MonteCarloSummaryPlotMaker * nChargedTowardsSummary = new MonteCarloSummaryPlotMaker( nChargedTowardsPlot, mcInfo, COMBINE_MC );
	allPlotMakers.push_back( nChargedTowardsSummary );

	//mean pT towards
	XPlotMaker * meanPtTowardsPlot = new XPlotMaker( "MeanPtTowards500", "PYTHIA6-AMBT1", meanPtBinEdges, 1.0, true );
	MonteCarloSummaryPlotMaker * meanPtTowardsSummary = new MonteCarloSummaryPlotMaker( meanPtTowardsPlot, mcInfo, COMBINE_MC );
	allPlotMakers.push_back( meanPtTowardsSummary );

	//pT sum towards
	XPlotMaker * pTSumTowardsPlot = new XPlotMaker( "PtSumTowards500", "PYTHIA6-AMBT1", sumPtBinEdges, 1.0, true );
	MonteCarloSummaryPlotMaker * pTSumTowardsSummary = new MonteCarloSummaryPlotMaker( pTSumTowardsPlot, mcInfo, COMBINE_MC );
	allPlotMakers.push_back( pTSumTowardsSummary );

	//2D Unfolding

	//Make a plot of number of charged particles in the toward region vs lead jet pT
	XvsYNormalisedPlotMaker * pTvsNChargedTowardPlot = new XvsYNormalisedPlotMaker( "LeadJetPt", "NChargedTowards500", "PYTHIA6-AMBT1",
			jetPtBinEdges, nChargeBinEdges, scaleFactor );
	MonteCarloSummaryPlotMaker * pTvsNChargedTowardSummary = new MonteCarloSummaryPlotMaker( pTvsNChargedTowardPlot, mcInfo, COMBINE_MC );
	pTvsNChargedTowardSummary->SetYRange( 0.1, 5.9 );
	pTvsNChargedTowardSummary->SetAxisLabels( "p_{T}^{lead} [MeV]", "<d^{2}N_{ch}/d#etad#phi>" );
	allPlotMakers.push_back( pTvsNChargedTowardSummary );

	//Make a plot of number of charged particles in the away region vs lead jet pT
	XvsYNormalisedPlotMaker * pTvsNChargedAwayPlot = new XvsYNormalisedPlotMaker( "LeadJetPt", "NChargedAway500", "PYTHIA6-AMBT1",
			jetPtBinEdges, nChargeBinEdges, scaleFactor );
	MonteCarloSummaryPlotMaker * pTvsNChargedAwaySummary = new MonteCarloSummaryPlotMaker( pTvsNChargedAwayPlot, mcInfo, COMBINE_MC );
	pTvsNChargedAwaySummary->SetYRange( 0.1, 5.9 );
	pTvsNChargedAwaySummary->SetAxisLabels( "p_{T}^{lead} [MeV]", "<d^{2}N_{ch}/d#etad#phi>" );
	allPlotMakers.push_back( pTvsNChargedAwaySummary );

	//Make a plot of number of charged particles in the transverse region vs lead jet pT
	XvsYNormalisedPlotMaker * pTvsNChargedTransPlot = new XvsYNormalisedPlotMaker( "LeadJetPt", "NChargedTransverse500", "PYTHIA6-AMBT1",
			jetPtBinEdges, nChargeBinEdges, scaleFactor );
	MonteCarloSummaryPlotMaker * pTvsNChargedTransSummary = new MonteCarloSummaryPlotMaker( pTvsNChargedTransPlot, mcInfo, COMBINE_MC );
	pTvsNChargedTransSummary->SetYRange( 0.1, 2.9 );
	pTvsNChargedTransSummary->SetAxisLabels( "p_{T}^{lead} [MeV]", "<d^{2}N_{ch}/d#etad#phi>" );
	allPlotMakers.push_back( pTvsNChargedTransSummary );

	//Make a plot of mean pT of charged particles in the toward region vs lead jet pT
	XvsYNormalisedPlotMaker * pTvsMeanPtTowardPlot = new XvsYNormalisedPlotMaker( "LeadJetPt", "MeanPtTowards500", "PYTHIA6-AMBT1",
			jetPtBinEdges, meanPtBinEdges, 1.0 );
	MonteCarloSummaryPlotMaker * pTvsMeanPtTowardSummary = new MonteCarloSummaryPlotMaker( pTvsMeanPtTowardPlot, mcInfo, COMBINE_MC );
	//	pTvsMeanPtTowardSummary->SetYRange( 0.1, 5.9 );
	pTvsMeanPtTowardSummary->SetAxisLabels( "p_{T}^{lead} [MeV]", "<p_{T}> [MeV]" );
	allPlotMakers.push_back( pTvsMeanPtTowardSummary );

	//Make a plot of mean pT of charged particles in the away region vs lead jet pT
	XvsYNormalisedPlotMaker * pTvsMeanPtAwayPlot = new XvsYNormalisedPlotMaker( "LeadJetPt", "MeanPtAway500", "PYTHIA6-AMBT1",
			jetPtBinEdges, meanPtBinEdges, 1.0 );
	MonteCarloSummaryPlotMaker * pTvsMeanPtAwaySummary = new MonteCarloSummaryPlotMaker( pTvsMeanPtAwayPlot, mcInfo, COMBINE_MC );
	//	pTvsMeanPtAwaySummary->SetYRange( 0.1, 5.9 );
	pTvsMeanPtAwaySummary->SetAxisLabels( "p_{T}^{lead} [MeV]", "<p_{T}> [MeV]" );
	allPlotMakers.push_back( pTvsMeanPtAwaySummary );

	//Make a plot of mean pT of charged particles in the transverse region vs lead jet pT
	XvsYNormalisedPlotMaker * pTvsMeanPtTransPlot = new XvsYNormalisedPlotMaker( "LeadJetPt", "MeanPtTransverse500", "PYTHIA6-AMBT1",
			jetPtBinEdges, meanPtBinEdges, 1.0 );
	MonteCarloSummaryPlotMaker * pTvsMeanPtTransSummary = new MonteCarloSummaryPlotMaker( pTvsMeanPtTransPlot, mcInfo, COMBINE_MC );
	//	pTvsMeanPtTransSummary->SetYRange( 0.1, 2.9 );
	pTvsMeanPtTransSummary->SetAxisLabels( "p_{T}^{lead} [MeV]", "<p_{T}> [MeV]" );
	allPlotMakers.push_back( pTvsMeanPtTransSummary );

	//Make a plot of sum pT of charged particles in the toward region vs lead jet pT
	XvsYNormalisedPlotMaker * pTvsPtSumTowardPlot = new XvsYNormalisedPlotMaker( "LeadJetPt", "PtSumTowards500", "PYTHIA6-AMBT1",
			jetPtBinEdges, sumPtBinEdges, scaleFactor );
	MonteCarloSummaryPlotMaker * pTvsPtSumTowardSummary = new MonteCarloSummaryPlotMaker( pTvsPtSumTowardPlot, mcInfo, COMBINE_MC );
	//	pTvsPtSumTowardSummary->SetYRange( 0.1, 5.9 );
	pTvsPtSumTowardSummary->SetAxisLabels( "p_{T}^{lead} [MeV]", "<d^{2}#Sigmap_{T}/d#etad#phi> [MeV]" );
	allPlotMakers.push_back( pTvsPtSumTowardSummary );

	//Make a plot of sum pT of charged particles in the away region vs lead jet pT
	XvsYNormalisedPlotMaker * pTvsPtSumAwayPlot = new XvsYNormalisedPlotMaker( "LeadJetPt", "PtSumAway500", "PYTHIA6-AMBT1",
			jetPtBinEdges, sumPtBinEdges, scaleFactor );
	MonteCarloSummaryPlotMaker * pTvsPtSumAwaySummary = new MonteCarloSummaryPlotMaker( pTvsPtSumAwayPlot, mcInfo, COMBINE_MC );
	//	pTvsPtSumAwaySummary->SetYRange( 0.1, 5.9 );
	pTvsPtSumAwaySummary->SetAxisLabels( "p_{T}^{lead} [MeV]", "<d^{2}#Sigmap_{T}/d#etad#phi> [MeV]" );
	allPlotMakers.push_back( pTvsPtSumAwaySummary );

	//Make a plot of sum pT of charged particles in the transverse region vs lead jet pT
	XvsYNormalisedPlotMaker * pTvsPtSumTransPlot = new XvsYNormalisedPlotMaker( "LeadJetPt", "PtSumTransverse500", "PYTHIA6-AMBT1",
			jetPtBinEdges, sumPtBinEdges, scaleFactor );
	MonteCarloSummaryPlotMaker * pTvsPtSumTransSummary = new MonteCarloSummaryPlotMaker( pTvsPtSumTransPlot, mcInfo, COMBINE_MC );
	//	pTvsPtSumTransSummary->SetYRange( 0.1, 2.9 );
	pTvsPtSumTransSummary->SetAxisLabels( "p_{T}^{lead} [MeV]", "<d^{2}#Sigmap_{T}/d#etad#phi> [MeV]" );
	allPlotMakers.push_back( pTvsPtSumTransSummary );

	//Make a plot of number of charged particles in the toward region vs lead jet pT
	XvsYNormalisedPlotMaker * pTmeanvsNChargedTowardPlot = new XvsYNormalisedPlotMaker( "NChargedTowards500", "MeanPtTowards500", "PYTHIA6-AMBT1",
			XnChargeBinEdges, meanPtBinEdges, 1.0 );
	MonteCarloSummaryPlotMaker * pTmeanvsNChargedTowardSummary = new MonteCarloSummaryPlotMaker( pTmeanvsNChargedTowardPlot, mcInfo, COMBINE_MC );
	pTmeanvsNChargedTowardSummary->SetYRange( 1000.0, 5000.0 );
	pTmeanvsNChargedTowardSummary->SetAxisLabels( "N_{ch}", "<p_{T}> [MeV]" );
	allPlotMakers.push_back( pTmeanvsNChargedTowardSummary );

	//Make a plot of number of charged particles in the away region vs lead jet pT
	XvsYNormalisedPlotMaker * pTmeanvsNChargedAwayPlot = new XvsYNormalisedPlotMaker( "NChargedAway500", "MeanPtAway500", "PYTHIA6-AMBT1",
			XnChargeBinEdges, meanPtBinEdges, 1.0 );
	MonteCarloSummaryPlotMaker * pTmeanvsNChargedAwaySummary = new MonteCarloSummaryPlotMaker( pTmeanvsNChargedAwayPlot, mcInfo, COMBINE_MC );
	pTmeanvsNChargedAwaySummary->SetYRange( 1000.0, 2000.0 );
	pTmeanvsNChargedAwaySummary->SetAxisLabels( "N_{ch}", "<p_{T}> [MeV]" );
	allPlotMakers.push_back( pTmeanvsNChargedAwaySummary );

	//Make a plot of number of charged particles in the transverse region vs lead jet pT
	XvsYNormalisedPlotMaker * pTmeanvsNChargedTransPlot = new XvsYNormalisedPlotMaker( "NChargedTransverse500", "MeanPtTransverse500", "PYTHIA6-AMBT1",
			XnChargeBinEdges, meanPtBinEdges, 1.0 );
	MonteCarloSummaryPlotMaker * pTmeanvsNChargedTransSummary = new MonteCarloSummaryPlotMaker( pTmeanvsNChargedTransPlot, mcInfo, COMBINE_MC );
	pTmeanvsNChargedTransSummary->SetYRange( 800.0, 1600.0 );
	pTmeanvsNChargedTransSummary->SetAxisLabels( "N_{ch}", "<p_{T}> [MeV]" );
	allPlotMakers.push_back( pTmeanvsNChargedTransSummary );*/

	///////////////////////////////////////////////////////////
	//                                                       //
	//  End of plot area                                     //
	//                                                       //
	///////////////////////////////////////////////////////////

	//Make an object to keep track of which observables we actually need
	ObservableList * relevanceChecker = new ObservableList( allPlotMakers );

	//Populate the smearing matrices
	for ( unsigned int mcIndex = 0; mcIndex < mcInfo->NumberOfSources(); mcIndex++ )
	{
		IFileInput * truthInput = mcInfo->MakeTruthInput( mcIndex, relevanceChecker );
		IFileInput * reconstructedInput = mcInfo->MakeReconstructedInput( mcIndex, relevanceChecker );
		MakeSmearingMatrices( truthInput, reconstructedInput );
		delete truthInput;
		delete reconstructedInput;
	}

	////////////////////////////////////////////////////////////
	//                                                        //
	// Load the data - Again, set this up yourself            //
	//                                                        //
	////////////////////////////////////////////////////////////
	IFileInput * dataInput = mcInfo->MakeReconstructedInput( 0, relevanceChecker );
	//IFileInput * dataInput = new TriggerChoosingInput( "/Disk/speyside7/Grid/grid-files/bwynne/Version4/JetTauEtmiss/PeriodGtoI/combined.TriggerName.AntiKt4TopoEM.root",
	//		"benTuple", "JetTauEtmiss Data (2010)", mcInfo->NumberOfSources(), relevanceChecker );

	//Unfold!
	DoTheUnfolding( dataInput );
	delete rootPlotStyle;

	//Status message
	cout << endl << "Imagiro finished" << endl;
	TimeAndMemory();
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
		//Status message
		cout << "File " << fileIndex << " - ";
		TimeAndMemory();

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

		recoMatched.clear();
	}

	//Debug
	cout << "Matched: " << matchedEvents << endl;
	cout << "Missed: " << missedEvents << endl;
	cout << "Fake: " << fakeEvents << endl;
}

void DoTheUnfolding( IFileInput * DataInput )
{
	//Status message
	cout << endl << "Loading " << *( DataInput->Description() ) << " events" << endl;

	//Populate the data distribution
	long dataTotal = 0;
	for ( unsigned int fileIndex = 0; fileIndex < DataInput->NumberOfFiles(); fileIndex++ )
	{
		//Status message
		cout << "File " << fileIndex << " ";
		TimeAndMemory();

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
	cout << endl << "Loading finished" << endl;
	TimeAndMemory();

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

//Ouput the current memory usage and time to stdout
void TimeAndMemory()
{
	//Connect to the file that has the memory usage for this process
	static std::ifstream meminfo( "/proc/self/statm" );
	meminfo.seekg(0);
	meminfo.sync();

	//Read the file
	unsigned long current_virt;
	meminfo >> current_virt;

	//Convert from pages to megabytes
	double memoryUsage = (double)current_virt * (double)getpagesize() / ( 1024.0 * 1024.0 );

	//Find out the time
	time_t timeNow;
	time( &timeNow );

	//Status message
	cout << "Memory usage: " << memoryUsage << " MB at time " << ctime( &timeNow ) << endl;
}

//Create an ATLAS style object
TStyle * PlotStyle( string StyleName )
{
	TStyle * atlasStyle = new TStyle( StyleName.c_str(), "Atlas style" );

	// use plain black on white colors
	Int_t icol=0; // WHITE
	atlasStyle->SetFrameBorderMode(icol);
	atlasStyle->SetFrameFillColor(icol);
	atlasStyle->SetCanvasBorderMode(icol);
	atlasStyle->SetCanvasColor(icol);
	atlasStyle->SetPadBorderMode(icol);
	atlasStyle->SetPadColor(icol);
	atlasStyle->SetStatColor(icol);
	//atlasStyle->SetFillColor(icol); // don't use: white fill color floa *all* objects

	// set the paper & margin sizes
	atlasStyle->SetPaperSize(20,26);

	// set margin sizes
	atlasStyle->SetPadTopMargin(0.05);
	atlasStyle->SetPadRightMargin(0.05);
	atlasStyle->SetPadBottomMargin(0.16);
	atlasStyle->SetPadLeftMargin(0.16);

	// set title offsets (for axis label)
	atlasStyle->SetTitleXOffset(1.4);
	atlasStyle->SetTitleYOffset(1.4);

	// use large fonts
	//Int_t font=72; // Helvetica italics
	Int_t font=42; // Helvetica
	Double_t tsize=0.05;
	atlasStyle->SetTextFont(font);

	atlasStyle->SetTextSize(tsize);
	atlasStyle->SetLabelFont(font,"x");
	atlasStyle->SetTitleFont(font,"x");
	atlasStyle->SetLabelFont(font,"y");
	atlasStyle->SetTitleFont(font,"y");
	atlasStyle->SetLabelFont(font,"z");
	atlasStyle->SetTitleFont(font,"z");

	atlasStyle->SetLabelSize(tsize,"x");
	atlasStyle->SetTitleSize(tsize,"x");
	atlasStyle->SetLabelSize(tsize,"y");
	atlasStyle->SetTitleSize(tsize,"y");
	atlasStyle->SetLabelSize(tsize,"z");
	atlasStyle->SetTitleSize(tsize,"z");

	// use bold lines and markers
	atlasStyle->SetMarkerStyle(20);
	atlasStyle->SetMarkerSize(1.2);
	atlasStyle->SetHistLineWidth(2);
	atlasStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

	// get rid of X error bars and y error bar caps
	//atlasStyle->SetErrorX(0.001);

	// do not display any of the standard histogram decorations
	atlasStyle->SetOptTitle(0);
	//atlasStyle->SetOptStat(1111);
	atlasStyle->SetOptStat(0);
	//atlasStyle->SetOptFit(1111);
	atlasStyle->SetOptFit(0);

	// put tick marks on top and RHS of plots
	atlasStyle->SetPadTickX(1);
	atlasStyle->SetPadTickY(1);

	return atlasStyle;
}
