/**
  @class NoCorrection

  Just creates distributions from the given data, with no attempt at correction

  @author Benjamin M Wynne bwynne@cern.ch
  @date 07-07-2011
 */

#ifndef NO_CORRECTION_H
#define NO_CORRECTION_H

#include "ICorrection.h"
#include "IIndexCalculator.h"

using namespace std;

class NoCorrection : public ICorrection
{
	public:
		//Default constructor - useless
		NoCorrection();

		//Constructor taking an IIndexCalculator to define the bins
		NoCorrection( IIndexCalculator * DistributionIndices, string Name, unsigned int UniqueID );

		//Destructor
		virtual ~NoCorrection();


		//Use this method to supply a value from the truth
		//distribution, and the corresponding reconstructed
		//value
		//NB: These values must both come from the SAME
		//Monte Carlo event, or the whole process is meaningless
		virtual void StoreTruthRecoPair( vector< double > Truth, vector< double > Reco, double TruthWeight = 1.0, double RecoWeight = 1.0, bool UseInPrior = true );


		//If an MC event is not reconstructed at all, use this
		//method to store the truth value alone
		virtual void StoreUnreconstructedTruth( vector< double > Truth, double Weight = 1.0, bool UseInPrior = true );


		//If there is a fake reconstructed event with no
		//corresponding truth, use this method
		virtual void StoreReconstructedFake( vector< double > Reco, double Weight = 1.0, bool UseInPrior = true );


		//Store a value from the distribution to be smeared
		virtual void StoreDataValue( vector< double > ToFold, double Weight = 1.0 );

		//Once all data is stored, run the correction
		//The arguments are dummies in this case - folding is a simple process
		virtual void Correct( unsigned int MostIterations, unsigned int ErrorMode = 0, bool WithSmoothing = false );

		//Perform a closure test
		//Fold the MC truth information
		//It should give the reco information back exactly...
		//Return true if test passed
		virtual bool ClosureTest( unsigned int MostIterations, bool WithSmoothing = false );

		//Perform an unfolding cross-check
		//Dummy - no iterations
		virtual unsigned int MonteCarloCrossCheck( Distribution * ReferenceDistribution, SmearingMatrix * InputSmearing, bool WithSmoothing = false );

		//Retrieve a TH1F* containing the unfolded data
		//distribution, with or without errors
		//NB: the error calculation is only performed
		//when you run the method with errors for the first time
		virtual TH1F * GetCorrectedHistogram( string Name, string Title, bool Normalise = false );

		//Retrieve the smearing matrix used
		virtual	TH2F * GetSmearingHistogram( string Name, string Title );
		virtual SmearingMatrix * GetSmearingMatrix();

		//Retrieve the reconstructed distribution (since the "true" folded value is the reconstructed one)
		virtual TH1F * GetTruthHistogram( string Name, string Title, bool Normalise = false );
		virtual Distribution * GetTruthDistribution();

		//Retrieve the reconstructed distribution
                virtual TH1F * GetReconstructedHistogram( string Name, string Title, bool Normalise = false );

		//Retrieve the input distribution to fold
		virtual TH1F * GetUncorrectedHistogram( string Name, string Title, bool Normalise = false );

		//Handy for error calculation
		virtual vector< double > Variances();
                virtual TH2F * DAgostiniCovariance( string Name, string Title );

		//Make another instance of the ICorrection which shares the smearing matrix
                virtual NoCorrection * CloneShareSmearingMatrix();

	private:
		//For use with Clone
		NoCorrection( IIndexCalculator * DistributionIndices, string Name, unsigned int UniqueID,
				Distribution * SharedReconstructed, SmearingMatrix * SharedSmearing );

		bool isClone;
		unsigned int uniqueID;
		string name;
		vector< double > sumOfInputWeightSquares;
		IIndexCalculator * indexCalculator;
		Distribution *inputDistribution, *reconstructedDistribution, *truthDistribution;
		SmearingMatrix * inputSmearing;
};

#endif
