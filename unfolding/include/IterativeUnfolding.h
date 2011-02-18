/**
  @class IterativeUnfolding

  The class that performs unfolding on a data distribution in the simple case when all events are available

  @author Benjamin M Wynne bwynne@cern.ch
  @date 17-06-2010
 */

#ifndef ITERATIVE_UNFOLDING_H
#define ITERATIVE_UNFOLDING_H

#include "Comparison.h"
#include "SmearingMatrix.h"
#include "Indices.h"
#include "Distribution.h"
#include <vector>
#include <string>
#include "TH1F.h"

using namespace std;

class IterativeUnfolding
{
	public:
		//Default constructor - useless
		IterativeUnfolding();

		//Constructor taking the required bin number,
		//minimum and maximum of the output distribution as arguments
		//NB: the unfolding scales roughly with bin number ^ 2, the
		//error calculation scales roughly with bin number ^ 3.
		IterativeUnfolding( int BinNumber, double Minimum, double Maximum, string Name = "", int UniqueID = 0, bool DebugMode = false );

		//N-Dimensional version
		IterativeUnfolding( vector<int> BinNumbers, vector<double> Minima, vector<double> Maxima, string Name = "", int UniqueID = 0, bool DebugMode = false );

		//Destructor
		~IterativeUnfolding();


		//Use this method to supply a value from the truth
		//distribution, and the corresponding reconstructed
		//value
		//NB: These values must both come from the SAME
		//Monte Carlo event, or the whole process is meaningless
		void StoreTruthRecoPair( double Truth, double Reco, double TruthWeight = 1.0, double RecoWeight = 1.0, bool UseInPrior = true );

		//N-Dimensional version
		void StoreTruthRecoPair( vector<double> Truth, vector<double> Reco, double TruthWeight = 1.0, double RecoWeight = 1.0, bool UseInPrior = true );


		//If an MC event is not reconstructed at all, use this
		//method to store the truth value alone
		void StoreUnreconstructedTruth( double Truth, double Weight = 1.0, bool UseInPrior = true );

		//N-Dimensional version
		void StoreUnreconstructedTruth( vector<double> Truth, double Weight = 1.0, bool UseInPrior = true );


		//If there is a fake reconstructed event with no
		//corresponding truth, use this method
		void StoreReconstructedFake( double Reco, double Weight = 1.0, bool UseInPrior = true );

		//N-Dimensional version
		void StoreReconstructedFake( vector<double> Reco, double Weight = 1.0, bool UseInPrior = true );


		//Store a value from the uncorrected data distribution
		void StoreDataValue( double Data, double Weight = 1.0 );

		//N-Dimensional version
		void StoreDataValue( vector<double> Data, double Weight = 1.0 );


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

		//Perform a closure test
		//Unfold the MC reco distribution with the corresponding truth information as a prior
		//It should give the truth information back exactly...
		void ClosureTest();

		//Perform an unfolding cross-check
		//Use MC truth A as a prior to unfold MC reco B
		//Iterations cease when result is sufficiently close to MC truth B (passed as argument)
		//Returns the number of iterations required. Convergence criteria as output arguments
		int MonteCarloCrossCheck( Distribution * ReferenceDistribution, double & ChiSquaredThreshold, double & KolmogorovThreshold, bool WithSmoothing = false );

		//Retrieve a TH1F* containing the unfolded data
		//distribution, with or without errors
		//NB: the error calculation is only performed
		//when you run the method with errors for the first time
		TH1F * GetUnfoldedHistogram( string Name = "unfolded", string Title = "Unfolded distribution", bool WithErrors = false );

		//Retrieve the smearing matrix used
		TH2F * GetSmearingMatrix( string Name, string Title );

		//Retrieve the truth distribution
		TH1F * GetTruthHistogram( string Name, string Title );
		Distribution * GetTruthDistribution();

		//Retrieve the uncorrected data distribution
		TH1F * GetUncorrectedDataHistogram( string Name, string Title );

		//Handy for error calculation
		vector<double> SumOfDataWeightSquares();

	private:
		Comparison * distributionComparison;
		int uniqueID;
		string name;
		vector< double > sumOfDataWeightSquares;
		Indices * indexCalculator;
		Distribution *dataDistribution, *unfoldedDistribution, *truthDistribution, *reconstructedDistribution;
		SmearingMatrix * inputSmearing;
		bool debug;
};

#endif
