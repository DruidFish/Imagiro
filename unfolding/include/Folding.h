/**
  @class Folding

  The class that applies a smearing matrix to a "truth" distribution, to simulate the detector effects

  @author Benjamin M Wynne bwynne@cern.ch
  @date 11-02-2011
 */

#ifndef FOLDING_H
#define FOLDING_H

#include "ICorrection.h"
#include "IIndexCalculator.h"
#include "Distribution.h"
#include "Comparison.h"

using namespace std;

class Folding : public ICorrection
{
	public:
		//Default constructor - useless
		Folding();

		//Constructor taking an IIndexCalculator to define the bins
		Folding( IIndexCalculator * DistributionIndices, string Name, unsigned int UniqueID );

		//Destructor
		virtual ~Folding();


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
		virtual unsigned int MonteCarloCrossCheck( Distribution * ReferenceDistribution, bool WithSmoothing = false );

		//Retrieve a TH1F* containing the unfolded data
		//distribution, with or without errors
		//NB: the error calculation is only performed
		//when you run the method with errors for the first time
		virtual TH1F * GetCorrectedHistogram( string Name, string Title, bool Normalise = false );

		//Retrieve the smearing matrix used
		virtual	TH2F * GetSmearingMatrix( string Name, string Title );

		//Retrieve the reconstructed distribution (since the "true" folded value is the reconstructed one)
		virtual TH1F * GetTruthHistogram( string Name, string Title, bool Normalise = false );
		virtual Distribution * GetTruthDistribution();

		//Retrieve the input distribution to fold
		virtual TH1F * GetUncorrectedHistogram( string Name, string Title, bool Normalise = false );

		//Handy for error calculation
		virtual vector< double > Variances();
                virtual TH2F * DAgostiniCovariance( string Name, string Title );

		//Make another instance of the ICorrection which shares the smearing matrix
                virtual Folding * CloneShareSmearingMatrix();

	private:
		//For use with Clone
		Folding( IIndexCalculator * DistributionIndices, string Name, unsigned int UniqueID,
				Comparison * SharedComparison, Distribution * SharedReconstructed, SmearingMatrix * SharedSmearing, double PairedMC, double MissedMC, double FakeMC );

		bool isClone;
		Comparison * distributionComparison;
		unsigned int uniqueID;
		string name;
		vector< double > sumOfInputWeightSquares;
		IIndexCalculator * indexCalculator;
		Distribution *inputDistribution, *smearedDistribution, *truthDistribution, *reconstructedDistribution;
		double totalPaired, totalFake, totalMissed;
		SmearingMatrix * inputSmearing;
};

#endif
