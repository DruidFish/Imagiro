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
		double ThisContribution( unsigned int I, unsigned int J, unsigned int K, unsigned int L );

	private:
		SmearingMatrix * smearing;
		UnfoldingMatrix * unfolding;
		vector< double > smearingErrors;
                vector< unsigned int > rIndices, sIndices, uIndices;

		//For fast lookups
		vector< vector< unsigned int > > u_to_entries;
		vector< vector< vector< unsigned int > > > r_to_S_to_entries, r_to_U_to_entries, u_to_S_to_entries;

		//Simple calculation caching
		vector< double > oneOverEfficiency;
		map< pair< unsigned int, unsigned int >, double > oneOverSmearing, deltaUR;

		//Caching for the u==k==l case
		vector< double > sumOverAll_RNotI_SNotJ;
		vector< vector< double > > sumOverThis_SNotJ, sumOverThis_SIsJ, sumOverThis_RNotI;
		vector< map< pair< unsigned int, unsigned int >, double > > correctionU_RS;
};

#endif
