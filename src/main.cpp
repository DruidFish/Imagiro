#include "CustomIndices.h"
#include "UniformIndices.h"
#include "TRandom3.h"
#include <iostream>
#include <vector>

using namespace std;

int main( int argc, char * argv[] )
{
	//Uniform bins
	unsigned int binNumber = 98;
	double binMinimum = 1.0;
	double binMaximum = 99.0;
	UniformIndices * uniform = new UniformIndices( vector< unsigned int >( 1, binNumber ), vector< double >( 1, binMinimum ), vector< double >( 1, binMaximum ) );

	//Custom bins
	vector< double > binLowEdges;
	for ( unsigned int binIndex = 1; binIndex < 100; binIndex++ )
	{
		binLowEdges.push_back( (double) binIndex );
	}
	CustomIndices * custom = new CustomIndices( vector< vector< double > >( 1, binLowEdges ) );

	//Random number generator
	TRandom3 * generator = new TRandom3(0);

	//Test loop
	unsigned long testFails = 0;
	for ( unsigned int testIndex = 0; testIndex < 10000000; testIndex++ )
	{
		vector< double > testNumber( 1, generator->Rndm() * 100.0 );
		unsigned int uniformIndex = uniform->GetIndex( testNumber );
	       	unsigned int customIndex = custom->GetIndex( testNumber );

		if ( uniformIndex != customIndex )
		{
			cout << uniformIndex << " vs " << customIndex << endl;
			testFails++;
		}
	}

	cout << testFails << endl;
}
