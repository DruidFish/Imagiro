/**
  @interface ICorrection

  General interface for the different ways of correcting a plot

  @author Benjamin M Wynne bwynne@cern.ch
  @date 04-07-2011
 */

#ifndef I_CORRECTION_H
#define I_CORRECTION_H

#include "SmearingMatrix.h"
#include "Distribution.h"
#include <vector>
#include <string>
#include "TH1F.h"
#include "TH2F.h"

using namespace std;

class ICorrection
{
	public:
		//Destructor
		virtual ~ICorrection()
		{
		}

		//Use this method to supply a value from the truth
		//distribution, and the corresponding reconstructed
		//value
		//NB: These values must both come from the SAME
		//Monte Carlo event, or the whole process is meaningless
		virtual void StoreTruthRecoPair( vector< double > Truth, vector<double> Reco, double TruthWeight = 1.0, double RecoWeight = 1.0, bool UseInPrior = true ) = 0;

		//If an MC event is not reconstructed at all, use this
		//method to store the truth value alone
		virtual void StoreUnreconstructedTruth( vector< double > Truth, double Weight = 1.0, bool UseInPrior = true ) = 0;

		//If there is a fake reconstructed event with no
		//corresponding truth, use this method
		virtual void StoreReconstructedFake( vector< double > Reco, double Weight = 1.0, bool UseInPrior = true ) = 0;

		//Store a value from the uncorrected data distribution
		virtual void StoreDataValue( vector< double > Data, double Weight = 1.0 ) = 0;

		//Once all data is stored, run the correction
		//You can specify when the iterations should end,
		//with an upper limit on iteration number
		//Set WithSmoothing = true to smooth the prior distribution
		//before each iteration, as it might reduce statistical
		//fluctuations when convergence is slow
		virtual void Correct( unsigned int MostIterations, unsigned int ErrorMode = 0, bool WithSmoothing = false ) = 0;

		//Perform a closure test
		//Unfold the MC reco distribution with the corresponding truth information as a prior
		//It should give the truth information back exactly...
		//Return true if test passed
		virtual bool ClosureTest( unsigned int MostIterations, bool WithSmoothing = false ) = 0;

		//Perform an unfolding cross-check
		//Use MC truth A as a prior to unfold MC reco B
		//Iterations cease when result is sufficiently close to MC truth B (passed as argument)
		//Returns the number of iterations required
		virtual unsigned int MonteCarloCrossCheck( Distribution * InputPriorDistribution, SmearingMatrix * InputSmearing, bool WithSmoothing = false ) = 0;

		//Retrieve a TH1F* containing the corrected data distribution
		virtual TH1F * GetCorrectedHistogram( string Name, string Title, bool Normalise = false ) = 0;

		//Retrieve the smearing matrix used
		virtual TH2F * GetSmearingHistogram( string Name, string Title ) = 0;
		virtual SmearingMatrix * GetSmearingMatrix() = 0;

		//Retrieve the truth distribution
		virtual TH1F * GetTruthHistogram( string Name, string Title, bool Normalise = false ) = 0;
		virtual Distribution * GetTruthDistribution() = 0;

		//Retrieve the uncorrected data distribution
		virtual TH1F * GetUncorrectedHistogram( string Name, string Title, bool Normalise = false ) = 0;

		//Handy for error calculation
		virtual vector< double > Variances() = 0;
		virtual TH2F * DAgostiniCovariance( string Name, string Title ) = 0;

		//Make another instance of the ICorrection which shares the smearing matrix
		virtual ICorrection * CloneShareSmearingMatrix() = 0;
};

#endif
