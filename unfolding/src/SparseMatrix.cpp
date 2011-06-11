/**
  @class SparseMatrix

  A matrix with many zero values, stored as a list of the non-zero entries
  Searchable or readable

  @author Benjamin M Wynne bwynne@cern.ch
  @date 09-04-2011
 */

#include "SparseMatrix.h"
#include <iostream>
#include <cstdlib>
#include <cassert>

SparseMatrix::SparseMatrix()
{
	nextEntry = matrix.begin();
	otherNextEntry = matrix.begin();
	entryNumber = 0;
	vectorsMade = false;
}

SparseMatrix::~SparseMatrix()
{
	matrix.clear();
	matrixValues.clear();
	secondIndices.clear();
}

//Add to the existing entry at these indices, or create a new entry if one does not exist
void SparseMatrix::AddToEntry( unsigned int FirstIndex, unsigned int SecondIndex, double Value )
{
	if ( vectorsMade )
	{
		cerr << "Trying to add elements to finalised sparse matrix" << endl;
		exit(1);
	}

	//Find if this matrix element already exists
	const pair< unsigned int, unsigned int > searchPair( FirstIndex, SecondIndex );
	const map< pair< unsigned int, unsigned int >, double >::iterator searchResult = matrix.find( searchPair );

	//Either create a new element or add to the existing one
	if ( searchResult == matrix.end() )
	{
		matrix[ searchPair ] = Value;
		entryNumber++;
	}
	else
	{
		matrix[ searchPair ] += Value;
	}
}

//Get the next non-zero entry in an iteration through them
double SparseMatrix::GetNextEntry( unsigned int & FirstIndex, unsigned int & SecondIndex, bool UseSecondIterator )
{
	if ( UseSecondIterator )
	{
		if ( otherNextEntry == matrix.end() )
		{
			cerr << "Acessing beyond end of sparse matrix: reset iterator" << endl;
			exit(1);
		}

		FirstIndex = otherNextEntry->first.first;
		SecondIndex = otherNextEntry->first.second;
		double value = otherNextEntry->second;

		otherNextEntry++;
		return value;
	}
	else
	{
		if ( nextEntry == matrix.end() )
		{
			cerr << "Acessing beyond end of sparse matrix: reset iterator" << endl;
			exit(1);
		}

		FirstIndex = nextEntry->first.first;
		SecondIndex = nextEntry->first.second;
		double value = nextEntry->second;

		nextEntry++;
		return value;
	}
}

//Get any element of the matrix
double SparseMatrix::GetElement( unsigned int FirstIndex, unsigned int SecondIndex )
{
	//Find if this matrix element already exists
	pair< unsigned int, unsigned int > searchPair( FirstIndex, SecondIndex );
	map< pair< unsigned int, unsigned int >, double >::iterator searchResult;
	searchResult = matrix.find( searchPair );

	if ( searchResult == matrix.end() )
	{
		return 0.0;
	}
	else
	{
		return searchResult->second;
	}
}

void SparseMatrix::VectorsFromMap( unsigned int BinNumber )
{
	//Set up the data structures
	matrixValues = vector< vector< double > >( BinNumber, vector< double >() );
	secondIndices = vector< vector< unsigned int > >( BinNumber, vector< unsigned int >() );

	//Read out the map
	map< pair< unsigned int, unsigned int >, double >::iterator matrixIterator;
	for ( matrixIterator = matrix.begin(); matrixIterator != matrix.end(); matrixIterator++ )
	{
		secondIndices[ matrixIterator->first.first ].push_back( matrixIterator->first.second );

		matrixValues[ matrixIterator->first.first ].push_back( matrixIterator->second );
	}

	nextEntry = matrix.begin();
	vectorsMade = true;
}

//Return a root 2D histogram containing the smearing matrix
TH2F * SparseMatrix::MakeRootHistogram( string Name, string Title )
{
	unsigned int binNumber = matrixValues.size();

	//Create the histogram object
	TH2F * outputHistogram = new TH2F( Name.c_str(), Title.c_str(), binNumber, 0.0, (double)binNumber, binNumber, 0.0, (double)binNumber );

	//Loop over all filled entries
	for ( unsigned int firstIndex = 0; firstIndex < matrixValues.size(); firstIndex++ )
	{
		for ( unsigned int entryIndex = 0; entryIndex < matrixValues[ firstIndex ].size(); entryIndex++ )
		{
			unsigned int outputBin = outputHistogram->GetBin( firstIndex + 1, secondIndices[ firstIndex ][ entryIndex ] + 1, 0 );

			outputHistogram->SetBinContent( outputBin, matrixValues[ firstIndex ][ entryIndex ] );
		}
	}

	return outputHistogram;
}

unsigned int SparseMatrix::GetBinNumber()
{
	return matrixValues.size();
}

unsigned int SparseMatrix::GetEntryNumberAndResetIterator( bool UseSecondIterator )
{
	if ( UseSecondIterator )
	{
		otherNextEntry = matrix.begin();
	}
	else
	{
		nextEntry = matrix.begin();
	}
	return entryNumber;
}

//Get all non-zero entries of the matrix with the given FirstIndex
vector< double > * SparseMatrix::GetEntriesWithFirstIndex( unsigned int FirstIndex )
{
	return &( matrixValues[ FirstIndex ] );
}
vector< unsigned int > * SparseMatrix::GetIndicesWithFirstIndex( unsigned int FirstIndex )
{
	return &( secondIndices[ FirstIndex ] );
}
