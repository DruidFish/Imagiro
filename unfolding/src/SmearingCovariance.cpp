#include "SmearingCovariance.h"
#include <iostream>

SmearingCovariance::SmearingCovariance()
{
}

SmearingCovariance::SmearingCovariance( SmearingMatrix * InputSmearing, UnfoldingMatrix * InputUnfolding )
{
	smearing = InputSmearing;
	unfolding = InputUnfolding;

	//Initialise lookups - apologies for horrible syntax
	unsigned int binNumber = InputSmearing->GetBinNumber();
	r_to_S_to_entries = vector< vector< vector< unsigned int > > >( binNumber, vector< vector< unsigned int > >( binNumber, vector< unsigned int >() ) );
	r_to_U_to_entries = vector< vector< vector< unsigned int > > >( binNumber, vector< vector< unsigned int > >( binNumber, vector< unsigned int >() ) );
	u_to_S_to_entries = vector< vector< vector< unsigned int > > >( binNumber, vector< vector< unsigned int > >( binNumber, vector< unsigned int >() ) );
	u_to_entries = vector< vector< unsigned int > >( binNumber, vector< unsigned int >() );

	unsigned int firstEntryNumber = InputSmearing->GetEntryNumberAndResetIterator();
	for ( unsigned int firstEntryIndex = 0; firstEntryIndex < firstEntryNumber; firstEntryIndex++ )
	{
		//Get a non-zero smearing matrix entry
		unsigned int r, u;
		double firstEntryValue = InputSmearing->GetNextEntry( u, r );

		//Cache the inverse of the entry
		oneOverSmearing[ pair< unsigned int, unsigned int >( u, r ) ] = 1.0 / firstEntryValue;

		//Cache part of the delta calculation
		deltaUR[ pair< unsigned int, unsigned int >( u, r ) ] = -unfolding->GetElement( u, r ) * smearing->GetEfficiency( u ) * oneOverSmearing[ pair< unsigned int, unsigned int >( u, r ) ];

		//Get all entries with the same cause index
		vector< unsigned int > * theseSIndices = InputSmearing->GetIndicesWithFirstIndex( u );
		vector< double > * secondEntryValues = InputSmearing->GetEntriesWithFirstIndex( u );

		//Loop over these entries
		for ( unsigned int secondEntryIndex = 0; secondEntryIndex < theseSIndices->size(); secondEntryIndex++ )
		{
			unsigned int s = ( *theseSIndices )[ secondEntryIndex ];
			double smearingError;
			double truthNumber = InputSmearing->GetTruthTotal( u );

			if ( r == s )
			{
				smearingError = firstEntryValue * ( 1.0 - firstEntryValue ) / truthNumber;
			}
			else
			{
				smearingError = -1.0 * firstEntryValue * ( *secondEntryValues )[ secondEntryIndex ] / truthNumber;
			}

			//Store the sparse matrix entry
			smearingErrors.push_back( smearingError );
			rIndices.push_back( r );
			sIndices.push_back( s );
			uIndices.push_back( u );
			unsigned int entryIndex = smearingErrors.size() - 1;

			//Store the quick lookups
			r_to_S_to_entries[ r ][ s ].push_back( entryIndex );
			r_to_U_to_entries[ r ][ u ].push_back( entryIndex );
			u_to_S_to_entries[ u ][ s ].push_back( entryIndex );
			u_to_entries[ u ].push_back( entryIndex );
		}
	}

	//Cache the inverse of the efficiencies
	oneOverEfficiency = vector< double >( binNumber, 0.0 );
	for ( unsigned int binIndex = 0; binIndex < binNumber; binIndex++ )
	{
		oneOverEfficiency[ binIndex ] = 1.0 / smearing->GetEfficiency( binIndex );
	}

	//Cache some results for the u==k==l case
	sumOverAll_RNotI_SNotJ = vector< double >( binNumber, 0.0 );
	sumOverThis_SNotJ = vector< vector< double > >( binNumber, vector< double >( binNumber, 0.0 ) );
	sumOverThis_SIsJ = vector< vector< double > >( binNumber, vector< double >( binNumber, 0.0 ) );
	sumOverThis_RNotI = vector< vector< double > >( binNumber, vector< double >( binNumber, 0.0 ) );
	for ( unsigned int uIndex = 0; uIndex < binNumber; uIndex++ )
	{
		//Make the correction factor map
		map< pair< unsigned int, unsigned int >, double > correctionsForThisU;

		//Loop over all entries with this u value
		for ( unsigned int searchIndex = 0; searchIndex < u_to_entries[ uIndex ].size(); searchIndex++ )
		{
			unsigned int entryIndex = u_to_entries[ uIndex ][ searchIndex ];
			unsigned int s = sIndices[ entryIndex ];
			unsigned int r = rIndices[ entryIndex ];

			//Work out the KIRU component
			double deltaKIRU = -oneOverEfficiency[ uIndex ];

			//Work out the LJSU component for s != j
			double deltaLJSU = -oneOverEfficiency[ uIndex ];

			//Store result for s != j
			double simpleResult = deltaKIRU * deltaLJSU * smearingErrors[ entryIndex ];
			sumOverThis_SNotJ[ uIndex ][ s ] += simpleResult;
			sumOverThis_RNotI[ uIndex ][ r ] += simpleResult;

			//Work out the LJSU component for s == j
			deltaLJSU += oneOverSmearing[ pair< unsigned int, unsigned int >( uIndex, s ) ];
			deltaLJSU += deltaUR[ pair< unsigned int, unsigned int >( uIndex, s ) ];

			//Store result for s == j
			double complexResult = deltaKIRU * deltaLJSU * smearingErrors[ entryIndex ];
			sumOverThis_SIsJ[ uIndex ][ s ] += complexResult;

			//Store the corection factor for r == i and s == j
			correctionsForThisU[ pair< unsigned int, unsigned int >( r, s ) ] = simpleResult - complexResult;
		}

		//Store the correction map
		correctionU_RS.push_back( correctionsForThisU );

		//Now sum the sums
		for ( unsigned int sIndex = 0; sIndex < binNumber; sIndex++ )
		{
			sumOverAll_RNotI_SNotJ[ uIndex ] += sumOverThis_SNotJ[ uIndex ][ sIndex ];
		}
	}
}

SmearingCovariance::~SmearingCovariance()
{
	smearingErrors.clear();
	rIndices.clear();
	sIndices.clear();
	uIndices.clear();
	r_to_S_to_entries.clear();
	r_to_U_to_entries.clear();
	u_to_S_to_entries.clear();
	u_to_entries.clear();
}

double SmearingCovariance::ThisContribution( unsigned int I, unsigned int J, unsigned int K, unsigned int L )
{
	double returnValue = 0.0;

	//Keep track of entries added already
	vector< bool > usedAlready( smearingErrors.size(), false );

	//Look for entries with r == i and s == j
	for ( unsigned int searchIndex = 0; searchIndex < r_to_S_to_entries[ I ][ J ].size(); searchIndex++ )
	{
		unsigned int entryIndex = r_to_S_to_entries[ I ][ J ][ searchIndex ];
		usedAlready[ entryIndex ] = true;

		//U information
		unsigned int u = uIndices[ entryIndex ];

		//Work out the KIRU component
		double deltaKIRU = deltaUR[ pair< unsigned int, unsigned int >( u, I ) ];
		if ( K == u )
		{
			deltaKIRU -= oneOverEfficiency[ u ];
			deltaKIRU += oneOverSmearing[ pair< unsigned int, unsigned int >( u, I ) ];
		}

		//Work out the LJSU component
		double deltaLJSU = deltaUR[ pair< unsigned int, unsigned int >( u, J ) ];
		if ( L == u )
		{
			deltaLJSU -= oneOverEfficiency[ u ];
			deltaLJSU += oneOverSmearing[ pair< unsigned int, unsigned int >( u, J ) ];
		}

		//Store result
		returnValue += deltaKIRU * deltaLJSU * smearingErrors[ entryIndex ];
	}

	//Look for entries with r == i and u == l
	for ( unsigned int searchIndex = 0; searchIndex < r_to_U_to_entries[ I ][ L ].size(); searchIndex++ )
	{
		unsigned int entryIndex = r_to_U_to_entries[ I ][ L ][ searchIndex ];
		if ( !usedAlready[ entryIndex ] )
		{
			//Not used already, therefore s != j
			usedAlready[ entryIndex ] = true;

			//Work out the KIRU component
			double deltaKIRU = deltaUR[ pair< unsigned int, unsigned int >( L, I ) ];
			if ( K == L )
			{
				deltaKIRU -= oneOverEfficiency[ L ];
				deltaKIRU += oneOverSmearing[ pair< unsigned int, unsigned int >( L, I ) ];
			}

			//Work out the LJSU component
			double deltaLJSU = -oneOverEfficiency[ L ];

			//Store result
			returnValue += deltaKIRU * deltaLJSU * smearingErrors[ entryIndex ];
		}
	}

	//Look for entries with u == k == l
	if ( K == L )
	{
		//Use cached sums over all entries with u == k
		double centralCorrection;
		map< pair< unsigned int, unsigned int >, double >::iterator searchResult = correctionU_RS[ K ].find( pair< unsigned int, unsigned int >( I, J ) );
		if ( searchResult == correctionU_RS[ K ].end() )
		{
			centralCorrection = 0.0;
		}
		else
		{
			centralCorrection = searchResult->second;
		}
		returnValue += sumOverAll_RNotI_SNotJ[ K ] - sumOverThis_SNotJ[ K ][ J ] + sumOverThis_SIsJ[ K ][ J ] - sumOverThis_RNotI[ K ][ I ] + centralCorrection;
	}
	else
	{
		//Look for entries with u == k and s == j
		for ( unsigned int searchIndex = 0; searchIndex < u_to_S_to_entries[ K ][ J ].size(); searchIndex++ )
		{
			unsigned int entryIndex = u_to_S_to_entries[ K ][ J ][ searchIndex ];
			if ( !usedAlready[ entryIndex ] )
			{
				//Not used already, therefore r != i
				usedAlready[ entryIndex ] = true;

				//Work out the KIRU component
				double deltaKIRU = -oneOverEfficiency[ K ];

				//Work out the LJSU component
				double deltaLJSU = deltaUR[ pair< unsigned int, unsigned int >( K, J ) ];

				//Store result
				returnValue += deltaKIRU * deltaLJSU * smearingErrors[ entryIndex ];
			}
		}
	}

	return returnValue;
}
