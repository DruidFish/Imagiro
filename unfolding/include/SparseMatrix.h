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
		unsigned int GetEntryNumberAndResetIterator( bool UseSecondIterator = false );

		//Get any element of the matrix
		double GetElement( unsigned int FirstIndex, unsigned int SecondIndex );

		//Get the next non-zero entry in an iteration through them
		double GetNextEntry( unsigned int & FirstIndex, unsigned int & SecondIndex, bool UseSecondIterator = false );

		//Get all non-zero entries of the matrix with the given FirstIndex
		vector< double > * GetEntriesWithFirstIndex( unsigned int FirstIndex );
		vector< unsigned int > * GetIndicesWithFirstIndex( unsigned int FirstIndex );

		//Return a root histogram containing the matrix
		TH2F * MakeRootHistogram( string Name, string Title );

		//Get the number of bins along one side of the matrix
		unsigned int GetBinNumber();

	protected:
		//Add to the existing entry at these indices, or create a new entry if one does not exist
		void AddToEntry( unsigned int FirstIndex, unsigned int SecondIndex, double Value );

		//Make the vectors from the map
		void VectorsFromMap( unsigned int BinNumber );

		map< pair< unsigned int, unsigned int >, double > matrix;
		map< pair< unsigned int, unsigned int >, double >::iterator nextEntry;
		map< pair< unsigned int, unsigned int >, double >::iterator otherNextEntry;
		vector< vector< double > > matrixValues;
		vector< vector< unsigned int > > secondIndices;

	private:
		bool vectorsMade;
		unsigned int entryNumber;
};

#endif
