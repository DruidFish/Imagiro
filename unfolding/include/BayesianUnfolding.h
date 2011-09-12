/**
  @class BayesianUnfolding

  The class that manages unfolding a distribution, and tests the quality of the unfolding

  @author Benjamin M Wynne bwynne@cern.ch
  @date 17-06-2010
 */

#ifndef BAYESIAN_UNFOLDING_H
#define BAYESIAN_UNFOLDING_H

#include "ICorrection.h"
#include "CovarianceMatrix.h"
#include "Comparison.h"
#include "SmearingMatrix.h"
#include "IIndexCalculator.h"
#include "Distribution.h"
#include <vector>
#include <string>
#include "TH1F.h"

using namespace std;

class BayesianUnfolding : public ICorrection
{
	public:
		//Default constructor - useless
		BayesianUnfolding();

		//Constructor taking an IIndexCalculator to define the bins
                BayesianUnfolding( IIndexCalculator * DistributionIndices, string Name, unsigned int UniqueID );

		//Destructor
		virtual ~BayesianUnfolding();


		//Use this method to supply a value from the truth
		//distribution, and the corresponding reconstructed
		//value
		//NB: These values must both come from the SAME
		//Monte Carlo event, or the whole process is meaningless
		virtual void StoreTruthRecoPair( vector< double > Truth, vector<double> Reco, double TruthWeight = 1.0, double RecoWeight = 1.0, bool UseInPrior = true );

		//If an MC event is not reconstructed at all, use this
		//method to store the truth value alone
		virtual void StoreUnreconstructedTruth( vector< double > Truth, double Weight = 1.0, bool UseInPrior = true );

		//If there is a fake reconstructed event with no
		virtual void StoreReconstructedFake( vector< double > Reco, double Weight = 1.0, bool UseInPrior = true );

		//Store a value from the uncorrected data distribution
		virtual void StoreDataValue( vector< double > Data, double Weight = 1.0 );

		//Once all data is stored, run the unfolding
		//You can specify when the iterations should end,
		//with an upper limit on iteration number
		//Set WithSmoothing = true to smooth the prior distribution
		//before each iteration, as it might reduce statistical
		//fluctuations when convergence is slow
		virtual void Correct( unsigned int MostIterations, unsigned int ErrorMode = 0, bool WithSmoothing = false );

		//Perform a closure test
		//Unfold the MC reco distribution with the corresponding truth information as a prior
		//It should give the truth information back exactly...
		//Return true if test passed
		virtual bool ClosureTest( unsigned int MostIterations, bool WithSmoothing = false );

		//Perform an unfolding cross-check
		//Use MC truth A as a prior to unfold MC reco B
		//Iterations cease when result is sufficiently close to MC truth B (passed as argument)
		//Returns the number of iterations required
		virtual unsigned int MonteCarloCrossCheck( Distribution * InputPriorDistribution, SmearingMatrix * InputSmearing, bool WithSmoothing = false );

		//Retrieve a TH1F* containing the corrected data distribution
		virtual TH1F * GetCorrectedHistogram( string Name, string Title, bool Normalise = false );

		//Retrieve the smearing matrix used
		virtual TH2F * GetSmearingHistogram( string Name, string Title );
		virtual SmearingMatrix * GetSmearingMatrix();

		//Retrieve the truth distribution
		virtual TH1F * GetTruthHistogram( string Name, string Title, bool Normalise = false );
		virtual Distribution * GetTruthDistribution();

		//Retrieve the uncorrected data distribution
		virtual TH1F * GetUncorrectedHistogram( string Name, string Title, bool Normalise = false );

		//Handy for error calculation
		virtual vector< double > Variances();
		virtual TH2F * DAgostiniCovariance( string Name, string Title );

		//Make another instance of the ICorrection which shares the smearing matrix
                virtual BayesianUnfolding * CloneShareSmearingMatrix();

	private:
		//For use with Clone
		BayesianUnfolding( IIndexCalculator * DistributionIndices, string Name, unsigned int UniqueID,
				Comparison * SharedComparison, Distribution * SharedTruthDistribution, SmearingMatrix * SharedSmearingMatrix );

		Comparison * distributionComparison;
		unsigned int uniqueID;
		string name;
		vector< double > sumOfDataWeightSquares, dagostiniVariance;
		IIndexCalculator * indexCalculator;
		Distribution *dataDistribution, *unfoldedDistribution, *truthDistribution, *reconstructedDistribution;
		CovarianceMatrix * fullErrors;
		SmearingMatrix * inputSmearing;
		bool isClone;
};

#endif
