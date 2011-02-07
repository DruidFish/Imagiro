/**
  @class MakeSeparateUnfoldingInputs

  The class to make the inputs for an unfolding process in situations where only some fraction of the events are available

  @author Benjamin M Wynne bwynne@cern.ch
  @date 17-06-2010
 */

#ifndef MAKE_SEPARATE_UNFOLDING_INPUTS_H
#define MAKE_SEPARATE_UNFOLDING_INPUTS_H

#include "Distribution.h"
#include "SmearingMatrix.h"
#include "Indices.h"
#include <vector>
#include "TH1F.h"
#include "TH2F.h"

using namespace std;

class MakeSeparateUnfoldingInputs
{
	public:
		//Default constructor - useless
		MakeSeparateUnfoldingInputs();

		//Constructor taking the required bin number,
		//minimum and maximum of the output distribution as arguments
		//NB: the unfolding scales roughly with bin number ^ 2, the
		//error calculation scales roughly with bin number ^ 3.
		MakeSeparateUnfoldingInputs( int BinNumber, double Minimum, double Maximum );

		//N-Dimensional version
		MakeSeparateUnfoldingInputs( vector<int> BinNumbers, vector<double> Minima, vector<double> Maxima );


		//Destructor
		~MakeSeparateUnfoldingInputs();

		//Use this method to supply a value from the truth
		//distribution, and the corresponding reconstructed
		//value
		//NB: These values must both come from the SAME
		//Monte Carlo event, or the whole process is meaningless
		void StoreTruthRecoPair( double Truth, double Reco, double Weight = 1.0 );

		//N-Dimensional version
		void StoreTruthRecoPair( vector<double> Truth, vector<double> Reco, double Weight = 1.0 );


		//If an MC event is not reconstructed at all, use this
		//method to store the truth value alone
		void StoreUnreconstructedTruth( double Truth, double Weight = 1.0 );

		//N-Dimensional version
		void StoreUnreconstructedTruth( vector<double> Truth, double Weight = 1.0 );


		//If there is a fake reconstructed event with no
		//corresponding truth, use this method
		void StoreReconstructedFake( double Reco, double Weight = 1.0 );

		//N-Dimensional version
		void StoreReconstructedFake( vector<double> Reco, double Weight = 1.0 );


		//Store a value from the uncorrected data distribution
		void StoreDataValue( double Data, double Weight = 1.0 );

		//N-Dimensional version
		void StoreDataValue( vector<double> Data, double Weight = 1.0 );


		//Retrieve an (un-normalised) smearing matrix to use
		//in a later calculation
		TH2F * UnnormalisedSmearingMatrix();

		//Retrieve an (un-normalised) truth distribution
		//to use in a later calculation
		TH1F * UnnormalisedTruthDistribution();

		//Retrieve an (un-normalised) data distribution
		//to use in a later calculation
		TH1F * UnnormalisedDataDistribution();

	private:
		Distribution * truth;
	       	Distribution * data;
		Indices * indexCalculator;
		SmearingMatrix * inputSmearing;
};

#endif
