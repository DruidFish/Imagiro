#include "InputNtuple.h"
#include "XvsYNormalisedPlotMaker.h"
#include "XvsYNormalisedFolding.h"
#include "MonteCarloSummaryPlotMaker.h"
#include "MonteCarloInformation.h"
#include "BasicPlotMaker.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TRandom3.h"
#include <ctime>
#include <iostream>
#include <string>
#include <cmath>

using namespace std;

//Method declarations
void MakeSmearingMatrices( InputNtuple * TruthNtuple, InputNtuple * ReconstructedNtuple );
void DoTheUnfolding( InputNtuple * DataNtuple );

//The plotmakers
vector< MonteCarloSummaryPlotMaker* > allPlotMakers;

//Time
time_t timeNow;

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
	vector< InputNtuple* > truthNtuples, reconstructedNtuples;
	for ( int sourceIndex = 0; sourceIndex < mcInfo->NumberOfSources(); sourceIndex++ )
	{
		truthNtuples.push_back( new InputNtuple( mcInfo->TruthFilePath(sourceIndex), "benTuple", mcInfo->Description( sourceIndex ), sourceIndex ) );
		reconstructedNtuples.push_back( new InputNtuple( mcInfo->ReconstructedFilePath(sourceIndex), "benTuple", mcInfo->Description( sourceIndex ), sourceIndex ) );
	}

	////////////////////////////////////////////////////////////
	//                                                        //
	// Load the data - Again, set this up yourself            //
	//                                                        //
	////////////////////////////////////////////////////////////
	InputNtuple * dataNtuple = new InputNtuple( "data/user.bwynne.LeadingJetModifiedv3.Data.CaloJet/mergedFile.root", "benTuple", "7TeVData", mcInfo->NumberOfSources() );
	//int sourceIndex = 4;
	//InputNtuple * dataNtuple = new InputNtuple( mcInfo->TruthFilePath(sourceIndex), "benTuple", mcInfo->Description( sourceIndex ), sourceIndex );
	//InputNtuple * dataNtuple = new InputNtuple( mcInfo->ReconstructedFilePath(sourceIndex), "benTuple", mcInfo->Description( sourceIndex ), sourceIndex );

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
	int jetPtBins = 30;
	double jetPtMin = 20.0;
	double jetPtMax = 50.0;
	int nChargeBins = 50;
	double nChargeMin = 0.5;
	double nChargeMax = 50.5;
	XvsYNormalisedPlotMaker * pTvsNChargedTowardPlot = new XvsYNormalisedPlotMaker( "MaxJetPt", "NChargeToward", "Pythia6",
			jetPtBins, jetPtMin, jetPtMax, nChargeBins, nChargeMin, nChargeMax, scaleFactor, allPlotMakers.size() );
	allPlotMakers.push_back( new MonteCarloSummaryPlotMaker( pTvsNChargedTowardPlot, mcInfo, 0.0, 3.0 ) );

	XvsYNormalisedPlotMaker * pTvsNChargedAwayPlot = new XvsYNormalisedPlotMaker( "MaxJetPt", "NChargeAway", "Pythia6",
			jetPtBins, jetPtMin, jetPtMax, nChargeBins, nChargeMin, nChargeMax, scaleFactor, allPlotMakers.size() );
	allPlotMakers.push_back( new MonteCarloSummaryPlotMaker( pTvsNChargedAwayPlot, mcInfo, 0.0, 3.0 ) );

	XvsYNormalisedPlotMaker * pTvsNChargedTransPlot = new XvsYNormalisedPlotMaker( "MaxJetPt", "NChargeTrans", "Pythia6",
			jetPtBins, jetPtMin, jetPtMax, nChargeBins, nChargeMin, nChargeMax, scaleFactor, allPlotMakers.size() );
	allPlotMakers.push_back( new MonteCarloSummaryPlotMaker( pTvsNChargedTransPlot, mcInfo, 0.0, 3.0 ) );

	/*BasicPlotMaker * pTvsNChargedTowardBasic = new BasicPlotMaker( "MaxJetPt", "NChargeToward", "Pythia6",
			jetPtBins, jetPtMin, jetPtMax, scaleFactor, allPlotMakers.size() );
	allPlotMakers.push_back( new MonteCarloSummaryPlotMaker( pTvsNChargedTowardBasic, mcInfo, 0.0, 3.0 ) );

	BasicPlotMaker * pTvsNChargedAwayBasic = new BasicPlotMaker( "MaxJetPt", "NChargeAway", "Pythia6",
			jetPtBins, jetPtMin, jetPtMax, scaleFactor, allPlotMakers.size() );
	allPlotMakers.push_back( new MonteCarloSummaryPlotMaker( pTvsNChargedAwayBasic, mcInfo, 0.0, 3.0 ) );

	BasicPlotMaker * pTvsNChargedTransBasic = new BasicPlotMaker( "MaxJetPt", "NChargeTrans", "Pythia6",
			jetPtBins, jetPtMin, jetPtMax, scaleFactor, allPlotMakers.size() );
	allPlotMakers.push_back( new MonteCarloSummaryPlotMaker( pTvsNChargedTransBasic, mcInfo, 0.0, 3.0 ) );*/

	/*XvsYNormalisedFolding * pTvsNChargedTowardFold = new XvsYNormalisedFolding( "MaxJetPt", "NChargeToward", "Pythia6",
			jetPtBins, jetPtMin, jetPtMax, nChargeBins, nChargeMin, nChargeMax, scaleFactor, allPlotMakers.size() );
	allPlotMakers.push_back( new MonteCarloSummaryPlotMaker( pTvsNChargedTowardFold, mcInfo, 0.0, 3.0 ) );

	XvsYNormalisedFolding * pTvsNChargedAwayFold = new XvsYNormalisedFolding( "MaxJetPt", "NChargeAway", "Pythia6",
			jetPtBins, jetPtMin, jetPtMax, nChargeBins, nChargeMin, nChargeMax, scaleFactor, allPlotMakers.size() );
	allPlotMakers.push_back( new MonteCarloSummaryPlotMaker( pTvsNChargedAwayFold, mcInfo, 0.0, 3.0 ) );

	XvsYNormalisedFolding * pTvsNChargedTransFold = new XvsYNormalisedFolding( "MaxJetPt", "NChargeTrans", "Pythia6",
			jetPtBins, jetPtMin, jetPtMax, nChargeBins, nChargeMin, nChargeMax, scaleFactor, allPlotMakers.size() );
	allPlotMakers.push_back( new MonteCarloSummaryPlotMaker( pTvsNChargedTransFold, mcInfo, 0.0, 3.0 ) );*/


	/////////////////////////////////////////////////////////////
	//                                                         //
	// No further user edits required                          //
	//                                                         //
	/////////////////////////////////////////////////////////////
	
	//Populate the smearing matrices
	for ( int mcIndex = 0; mcIndex < truthNtuples.size(); mcIndex++ )
	{
		MakeSmearingMatrices( truthNtuples[mcIndex], reconstructedNtuples[mcIndex] );
	}

	//Unfold!
	DoTheUnfolding( dataNtuple );

	//Status message
	time(&timeNow);
	cout << endl << "Imagiro finished: " << ctime( &timeNow ) << endl;
}

//Match up event numbers between truth and reco inputs
void MakeSmearingMatrices( InputNtuple * TruthNtuple, InputNtuple * ReconstructedNtuple )
{
	//Initialise
	long matchedEvents = 0;
	long fakeEvents = 0;
	long missedEvents = 0;
	vector<bool> recoWasMatched( ReconstructedNtuple->NumberOfRows(), false );
	cout << endl << "Loading " << *( TruthNtuple->Description() ) << " events" << endl;

	//Loop over all truth events to try and match to reconstructed events
	for ( long truthIndex = 0; truthIndex < TruthNtuple->NumberOfRows(); truthIndex++ )
	{
		//Get the event number for each truth event
		TruthNtuple->ReadRow(truthIndex);
		float truthEventNumber = TruthNtuple->EventNumber();

		//Try to find a reconstructed event with that number
		if ( ReconstructedNtuple->ReadEvent(truthEventNumber) )
		{
			//Add the match to all plot makers
			for ( int plotIndex = 0; plotIndex < allPlotMakers.size(); plotIndex++ )
			{
				allPlotMakers[plotIndex]->StoreMatch( TruthNtuple, ReconstructedNtuple );
			}

			//Record the match
			recoWasMatched[ ReconstructedNtuple->CurrentRow() ] = true;
			matchedEvents++;
		}
		else
		{
			//Add the miss to all plot makers
			for ( int plotIndex = 0; plotIndex < allPlotMakers.size(); plotIndex++ )
			{
				allPlotMakers[plotIndex]->StoreMiss( TruthNtuple );
			}

			//Record the miss
			missedEvents++;
		}
	}

	//Loop over all reconstructed events looking for fakes
	for ( long reconstructedIndex = 0; reconstructedIndex < ReconstructedNtuple->NumberOfRows(); reconstructedIndex++ )
	{
		//Use the cached search success/fail rather than search again (half the disk io: big time saver)
		if ( !recoWasMatched[reconstructedIndex] )
		{
			//Update the line being accessed
			ReconstructedNtuple->ReadRow(reconstructedIndex);

			//Add the fake to all plotmakers
			for ( int plotIndex = 0; plotIndex < allPlotMakers.size(); plotIndex++ )
			{
				allPlotMakers[plotIndex]->StoreFake( ReconstructedNtuple );
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

void DoTheUnfolding( InputNtuple * DataNtuple )
{
	//Populate the data distribution
	cout << endl << "Loading " << *( DataNtuple->Description() ) << " events" << endl;
	for ( long rowIndex = 0; rowIndex < DataNtuple->NumberOfRows(); rowIndex++ )
	{
		//Load the row into memory
		DataNtuple->ReadRow( rowIndex );

		//Store the row in all plot makers
		for ( int plotIndex = 0; plotIndex < allPlotMakers.size(); plotIndex++ )
		{
			allPlotMakers[ plotIndex ]->StoreData( DataNtuple );
		}
	}
	cout << "Total: " << DataNtuple->NumberOfRows() << endl;

	//Status message
	time(&timeNow);
	cout << endl << "Loading finished: " << ctime( &timeNow ) << endl;

	//Produce the plots
	TFile * OutputFile = new TFile( OUTPUT_FILE_NAME.c_str(), "RECREATE" );
	for ( int plotIndex = 0; plotIndex < allPlotMakers.size(); plotIndex++ )
	{
		//Unfold the data
		allPlotMakers[ plotIndex ]->Unfold();

		//Save the output
		allPlotMakers[ plotIndex ]->ResultPlot()->Write();
		allPlotMakers[ plotIndex ]->SmearingMatrix()->Write();
		OutputFile->Flush();

		//Free up some memory
		delete allPlotMakers[ plotIndex ];
	}
	OutputFile->Close();
}
