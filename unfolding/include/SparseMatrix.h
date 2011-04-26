/**
  @class SparseMatrix

  A matrix with many zero values, stored as a list of the non-zero entries
  Searchable or readable

  @author Benjamin M Wynne bwynne@cern.ch
  @date 09-04-2011
  */

#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <map>
#include <vector>
#include <string>
#include "TH2F.h"

using namespace std;

class SparseMatrix
{
	public:
		SparseMatrix();
		~SparseMatrix();

		//Get the number of non-zero entries
		int GetEntryNumberAndResetIterator( bool UseSecondIterator = false );

		//Get any element of the matrix
		double GetElement( int FirstIndex, int SecondIndex );

		//Get the next non-zero entry in an iteration through them
		double GetNextEntry( int & FirstIndex, int & SecondIndex, bool UseSecondIterator = false );

		//Get all non-zero entries of the matrix with the given FirstIndex
		vector< double > * GetEntriesWithFirstIndex( int FirstIndex );
		vector< int > * GetIndicesWithFirstIndex( int FirstIndex );

		//Return a root histogram containing the matrix
		TH2F * MakeRootHistogram( string Name, string Title );

		//Get the number of bins along one side of the matrix
		int GetBinNumber();

	protected:
		//Add to the existing entry at these indices, or create a new entry if one does not exist
		void AddToEntry( int FirstIndex, int SecondIndex, double Value );

		//Make the vectors from the map
		void VectorsFromMap( int BinNumber );

		map< pair< int, int >, double > matrix;
		map< pair< int, int >, double >::iterator nextEntry;
		map< pair< int, int >, double >::iterator otherNextEntry;
		vector< vector< double > > matrixValues;
		vector< vector< int > > secondIndices;

	private:
		bool vectorsMade;
		int entryNumber;
};

#endif
