#ifndef SMEARING_COVARIANCE_H
#define SMEARING_COVARIANCE_H

#include <map>
#include <vector>
#include "SmearingMatrix.h"
#include "UnfoldingMatrix.h"

using namespace std;

class SmearingCovariance
{
	public:
		SmearingCovariance();
		SmearingCovariance( SmearingMatrix * InputSmearing, UnfoldingMatrix * InputUnfolding );
		~SmearingCovariance();

		//For full covariance
		double ThisContribution( int I, int J, int K, int L );

	private:
		SmearingMatrix * smearing;
		UnfoldingMatrix * unfolding;
		vector< double > smearingErrors;
                vector< int > rIndices, sIndices, uIndices;

		//For fast lookups
		vector< vector< int > > u_to_entries;
		vector< vector< vector< int > > > r_to_S_to_entries, r_to_U_to_entries, u_to_S_to_entries;

		//Simple calculation caching
		vector< double > oneOverEfficiency;
		map< pair< int, int >, double > oneOverSmearing, deltaUR;

		//Caching for the u==k==l case
		vector< double > sumOverAll_RNotI_SNotJ;
		vector< vector< double > > sumOverThis_SNotJ, sumOverThis_SIsJ, sumOverThis_RNotI;
		vector< map< pair< int, int >, double > > correctionU_RS;
};

#endif
