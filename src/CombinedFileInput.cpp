/**
  @class CombinedFileInput

  Combines input files on the fly without affecting their structure on disk

  @author Benjamin M Wynne bwynne@cern.ch
  @date 06-04-2011
 */

#include "InputUETree.h"
#include "CombinedFileInput.h"
#include <iostream>
#include <cstdlib>

CombinedFileInput::CombinedFileInput()
{
}

CombinedFileInput::CombinedFileInput( vector< IFileInput* > FileInputs, vector< double > FileWeights, string Description, unsigned int DescriptionIndex )
{
	sourceDescription = Description;
	sourceIndex = DescriptionIndex;
	fileWeights = FileWeights;
	totalRows = 0;
	currentFile = 0;

	//Check the input
	if ( FileInputs.size() < 1 )
	{
		cerr << "CombinedFileInput requires at least one input file" << endl;
		exit(1);
	}
	else if ( FileInputs.size() != FileWeights.size() )
	{
		cerr << "Mismatch in number of input files and number of file weights" << endl;
		exit(1);
	}
	else
	{
		//Set up each slice input
		fileInputs = FileInputs;
		for ( unsigned int fileIndex = 0; fileIndex < fileInputs.size(); fileIndex++ )
		{
			currentFile = fileIndex;
			totalRows += fileInputs[ fileIndex ]->NumberOfRows();
		}
	}
}

CombinedFileInput::~CombinedFileInput()
{
	for ( unsigned int fileIndex = 0; fileIndex < fileInputs.size(); fileIndex++ )
	{
		delete fileInputs[ fileIndex ];
	}
}

//Access a particular event, return false if the event is not found
bool CombinedFileInput::ReadRow( unsigned long RowNumber )
{
	//Check for out of range
	if ( RowNumber < 0 || RowNumber >= totalRows )
	{
		return false;
	}
	else
	{
		//Find the slice containing this row
		for ( unsigned int fileIndex = 0; fileIndex < fileInputs.size(); fileIndex++ )
		{
			if ( RowNumber >= fileInputs[ fileIndex ]->NumberOfRows() )
			{
				//The row is not in this slice
				RowNumber -= fileInputs[ fileIndex ]->NumberOfRows();
			}
			else
			{
				//Found the slice containing the row
				currentFile = fileIndex;
				return fileInputs[ fileIndex ]->ReadRow( RowNumber );
			}
		}

		//To stop -Wall complaining
		return false;
	}
}

//Access a particular event, return false if the event is not found
bool CombinedFileInput::ReadEvent( UInt_t EventNumber )
{
	bool found = false;

	//Look for the event in each slice
	for ( unsigned int fileIndex = 0; fileIndex < fileInputs.size(); fileIndex++ )
	{
		if ( fileInputs[ fileIndex ]->ReadEvent( EventNumber ) )
		{
			//Check there are no duplicate event numbers - should use the other ReadEvent method if so
			if ( found )
			{
				cerr << "Duplicate event numbers!" << endl;
				exit(1);
			}
			else
			{
				currentFile = fileIndex;
				found = true;
			}
		}
	}

	return found;
}
bool CombinedFileInput::ReadEvent( UInt_t EventNumber, unsigned int FileIndex )
{
	//Constrain the search to a particular file to avoid event number conflicts
	currentFile = FileIndex;
	return fileInputs[ FileIndex ]->ReadEvent( EventNumber );
}

//Get the standard event number and weight information
UInt_t CombinedFileInput::EventNumber()
{
	return fileInputs[ currentFile ]->EventNumber();
}

double CombinedFileInput::EventWeight()
{
	return fileInputs[ currentFile ]->EventWeight() * fileWeights[ currentFile ];
}

//Get any other column value by name
double CombinedFileInput::GetValue( string VariableName )
{
	return fileInputs[ currentFile ]->GetValue( VariableName );
}

//Get the number of rows
unsigned long CombinedFileInput::NumberOfRows()
{
	return totalRows;
}

unsigned long CombinedFileInput::CurrentRow()
{
	long currentRow = 0;

	//Add up the rows in earlier slices
	for ( int fileIndex = 0; fileIndex < currentFile; fileIndex++ )
	{
		currentRow += fileInputs[ fileIndex ]->NumberOfRows();
	}

	//Then add the row of the active slice
	currentRow += fileInputs[ currentFile ]->CurrentRow();

	return currentRow;
}

unsigned int CombinedFileInput::CurrentFile()
{
	return currentFile;
}

//Get the description of the source
string * CombinedFileInput::Description()
{
	return &sourceDescription;
}

unsigned int CombinedFileInput::DescriptionIndex()
{
	return sourceIndex;
}
