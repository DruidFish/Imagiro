/**
  @class Folding

  The class that applies a smearing matrix to a "truth" distribution, to simulate the detector effects

  @author Benjamin M Wynne bwynne@cern.ch
  @date 11-02-2011
 */

#ifndef FOLDING_H
#define FOLDING_H

#include "SmearingMatrix.h"
#include "Indices.h"
#include "Distribution.h"
#include "Comparison.h"
#include <vector>
#include <string>
#include "TH1F.h"

using namespace std;

class Folding
{
	public:
		//Default constructor - useless
		Folding();

		//Constructor taking the required bin number,
		//minimum and maximum of the output distribution as arguments
		//NB: the unfolding scales roughly with bin number ^ 2, the
		//error calculation scales roughly with bin number ^ 3.
		Folding( int BinNumber, double Minimum, double Maximum, string Name = "", int UniqueID = 0 );

		//N-Dimensional version
		Folding( vector<int> BinNumbers, vector<double> Minima, vector<double> Maxima, string Name = "", int UniqueID = 0 );

		//Destructor
		~Folding();


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


		//Store a value from the distribution to be smeared
		void StoreValueToFold( double ToFold, double Weight = 1.0 );

		//N-Dimensional version
		void StoreValueToFold( vector<double> ToFold, double Weight = 1.0 );

		//Smear the input distribution
		void Fold();

		//Perform a closure test
                //Fold the MC truth information - should give the MC reco exactly
		//Return true if test passed
                bool ClosureTest();

		//Retrieve a TH1F* containing the unfolded data
		//distribution, with or without errors
		//NB: the error calculation is only performed
		//when you run the method with errors for the first time
		TH1F * GetFoldedHistogram( string Name, string Title, bool Normalise = false );

		//Retrieve the smearing matrix used
		TH2F * GetSmearingMatrix( string Name, string Title );

		//Retrieve the reconstructed distribution
		TH1F * GetReconstructedHistogram( string Name, string Title, bool Normalise = false );

                //Retrieve the input distribution to fold
                TH1F * GetInputHistogram( string Name, string Title, bool Normalise = false );

                //Handy for error calculation
                vector<double> SumOfInputWeightSquares();

	private:
		Comparison * distributionComparison;
		int uniqueID;
		string name;
		vector< double > sumOfInputWeightSquares;
		Indices * indexCalculator;
		Distribution *inputDistribution, *smearedDistribution, *truthDistribution, *reconstructedDistribution;
		double totalPaired, totalFake, totalMissed;
		SmearingMatrix * inputSmearing;
};

#endif
