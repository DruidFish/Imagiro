/**
  @class UseSeparateUnfoldingInputs

  The class that performs the unfolding by using the ROOT objects made by MakeSeparateUnfoldingInputs

  @author Benjamin M Wynne bwynne@cern.ch
  @date 17-06-2010
 */

#ifndef USE_SEPARATE_UNFOLDING_INPUTS_H
#define USE_SEPARATE_UNFOLDING_INPUTS_H

#include "SmearingMatrix.h"
#include "Indices.h"
#include "Distribution.h"
#include <vector>
#include "TH1F.h"
#include "TH2F.h"

using namespace std;

class UseSeparateUnfoldingInputs
{
	public:
		//Default constructor - useless
		UseSeparateUnfoldingInputs();

		//Constructor taking the required bin number,
		//minimum and maximum of the output distribution as arguments
		//NB: the unfolding scales roughly with bin number ^ 2, the
		//error calculation scales roughly with bin number ^ 3.
		UseSeparateUnfoldingInputs( int BinNumber, double Minimum, double Maximum, bool DebugMode = false );

		//N-Dimensional version
		UseSeparateUnfoldingInputs( vector<int> BinNumbers, vector<double> Minima, vector<double> Maxima, bool DebugMode = false );


		//Destructor
		~UseSeparateUnfoldingInputs();

		//If you want to make a smearing matrix from a bunch of
		//other ones stored as (un-normalised) TH2Fs, use this
		void StoreSmearingMatrix( TH2F* );

		//Store distribution(s) to use as a prior for the unfolding
		void StorePriorDistribution( TH1F* );

		//Store data distribution(s) to combine and unfold
		void StoreDataDistribution( TH1F* );

		//Once all data is stored, run the unfolding
		//You can specify when the iterations should end,
		//with an upper limit on iteration number, or by
		//comparing the results from the last two iterations.
		//Iteration ends if the chi squared comparison value of
		//the two results is lower than the threshold, or if
		//the Kolmogorov-Smirnof comparison value is higher
		//Set WithSmoothing = true to smooth the prior distribution
		//before each iteration, as it might reduce statistical
		//fluctuations when convergence is slow
		void Unfold( int MostIterations = 20, double ChiSquaredThreshold = 10.0, double KolmogorovThreshold = 0.1, bool WithSmoothing = false );

		//Retrieve a TH1F* containing the unfolded data
		//distribution, with or without errors
		//NB: the error calculation is only performed
		//when you run the method with errors for the first time
		TH1F * UnfoldedDistribution( bool WithErrors = false );

		//Retrieve the unfolded distribution from each iteration
		//as a vector of TH1F*.
		//The 0th entry will be the prior (MC truth) distirbution,
		//and the last entry will be the same as that returned
		//by UnfoldedDistribution()
		vector< TH1F* > AllIterationResults();

	private:
		vector< TH1F* > allResults, inputPriorDistributions, inputDataDistributions;
		SmearingMatrix * inputSmearing;
		Indices * indexCalculator;
		Distribution * lastResult;
		bool debug;
};

#endif
