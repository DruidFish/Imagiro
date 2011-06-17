/**
  @class IterativeUnfolding

  The class that manages unfolding a distribution, and tests the quality of the unfolding

  @author Benjamin M Wynne bwynne@cern.ch
  @date 17-06-2010
 */

#ifndef ITERATIVE_UNFOLDING_H
#define ITERATIVE_UNFOLDING_H

#include "CovarianceMatrix.h"
#include "Comparison.h"
#include "SmearingMatrix.h"
#include "IIndexCalculator.h"
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

		//Constructor taking an IIndexCalculator to define the bins
                IterativeUnfolding( IIndexCalculator * DistributionIndices, string Name, unsigned int UniqueID );

		//Destructor
		~IterativeUnfolding();


		//Use this method to supply a value from the truth
		//distribution, and the corresponding reconstructed
		//value
		//NB: These values must both come from the SAME
		//Monte Carlo event, or the whole process is meaningless
		void StoreTruthRecoPair( double Truth, double Reco, double TruthWeight = 1.0, double RecoWeight = 1.0, bool UseInPrior = true );

		//N-Dimensional version
		void StoreTruthRecoPair( vector< double > Truth, vector<double> Reco, double TruthWeight = 1.0, double RecoWeight = 1.0, bool UseInPrior = true );


		//If an MC event is not reconstructed at all, use this
		//method to store the truth value alone
		void StoreUnreconstructedTruth( double Truth, double Weight = 1.0, bool UseInPrior = true );

		//N-Dimensional version
		void StoreUnreconstructedTruth( vector< double > Truth, double Weight = 1.0, bool UseInPrior = true );


		//If there is a fake reconstructed event with no
		//corresponding truth, use this method
		void StoreReconstructedFake( double Reco, double Weight = 1.0, bool UseInPrior = true );

		//N-Dimensional version
		void StoreReconstructedFake( vector< double > Reco, double Weight = 1.0, bool UseInPrior = true );


		//Store a value from the uncorrected data distribution
		void StoreDataValue( double Data, double Weight = 1.0 );

		//N-Dimensional version
		void StoreDataValue( vector< double > Data, double Weight = 1.0 );


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
		void Unfold( unsigned int MostIterations, double ChiSquaredThreshold, double KolmogorovThreshold, unsigned int ErrorMode = 0, bool WithSmoothing = false );

		//Perform a closure test
		//Unfold the MC reco distribution with the corresponding truth information as a prior
		//It should give the truth information back exactly...
		//Return true if test passed
		bool ClosureTest( unsigned int MostIterations, double ChiSquaredThreshold, double KolmogorovThreshold, bool WithSmoothing = false );

		//Perform an unfolding cross-check
		//Use MC truth A as a prior to unfold MC reco B
		//Iterations cease when result is sufficiently close to MC truth B (passed as argument)
		//Returns the number of iterations required. Convergence criteria as output arguments
		unsigned int MonteCarloCrossCheck( Distribution * ReferenceDistribution, double & ChiSquaredThreshold, double & KolmogorovThreshold, bool WithSmoothing = false );

		//Retrieve a TH1F* containing the unfolded data
		//distribution, with or without errors
		//NB: the error calculation is only performed
		//when you run the method with errors for the first time
		TH1F * GetUnfoldedHistogram( string Name, string Title, bool Normalise = false );

		//Retrieve the smearing matrix used
		TH2F * GetSmearingMatrix( string Name, string Title );

		//Retrieve the truth distribution
		TH1F * GetTruthHistogram( string Name, string Title, bool Normalise = false );
		Distribution * GetTruthDistribution();

		//Retrieve the uncorrected data distribution
		TH1F * GetUncorrectedDataHistogram( string Name, string Title, bool Normalise = false );

		//Handy for error calculation
		vector< double > SumOfDataWeightSquares();
		vector< double > DAgostiniVariance();
		TH2F * DAgostiniCovariance( string Name, string Title );

	private:
		Comparison * distributionComparison;
		unsigned int uniqueID;
		string name;
		vector< double > sumOfDataWeightSquares, dagostiniVariance;
		IIndexCalculator * indexCalculator;
		Distribution *dataDistribution, *unfoldedDistribution, *truthDistribution, *reconstructedDistribution;
		CovarianceMatrix * fullErrors;
		SmearingMatrix * inputSmearing;
		double totalPaired, totalFake, totalMissed;
		bool debug;
};

#endif
