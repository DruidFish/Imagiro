/**
  @class InputUETree

  Accesses a slightly less general file format I use for the underlying event analysis

  @author Benjamin M Wynne bwynne@cern.ch
  @date 06-04-2011
 */

#include "InputUETree.h"
#include "TLeaf.h"
#include <iostream>
#include <cstdlib>

const string EVENT_NUMBER_COLUMN_NAME = "EventNumber";
const string EVENT_WEIGHT_COLUMN_NAME = "EventWeight";

//Default constructor - useless
InputUETree::InputUETree()
{
}

//Constructor taking arguments pointing to a particular Ntuple in a root file
InputUETree::InputUETree( string FilePath, string NtuplePath, string Description, int DescriptionIndex )
{
	sourceDescription = Description;
	sourceDescriptionIndex = DescriptionIndex;

	//Get the Ntuple from the file
	inputFile = new TFile( FilePath.c_str(), "READ" );
	wrappedNtuple = ( TTree* )inputFile->Get( NtuplePath.c_str() );

	//Find out the number of rows
	numberOfRows = wrappedNtuple->GetEntries();

	//Loop over all columns/branches to find their names
	TIter branchIterator( wrappedNtuple->GetListOfLeaves() );
	TLeaf * nextLeaf;
	bool eventNumberFound = false;
	bool eventWeightFound = false;
	vector<string> otherColumnNames;
	while ( ( nextLeaf = ( TLeaf* )branchIterator() ) )
	{
		string leafName = nextLeaf->GetName();

		//Find the special columns for EventNumber and EventWeight
		if ( leafName == EVENT_NUMBER_COLUMN_NAME )
		{
			wrappedNtuple->SetBranchAddress( EVENT_NUMBER_COLUMN_NAME.c_str(), &currentEventNumber, &eventNumberBranch );
			eventNumberFound = true;
		}
		else if ( leafName == EVENT_WEIGHT_COLUMN_NAME )
		{
			wrappedNtuple->SetBranchAddress( EVENT_WEIGHT_COLUMN_NAME.c_str(), &currentEventWeight, &eventWeightBranch );
			eventWeightFound = true;
		}
		else
		{
			if ( leafName != "DeltaPhi" )
			{
				otherColumnNames.push_back(leafName);
			}
		}
	}

	//Check we found the special columns
	if ( !eventNumberFound )
	{
		cerr << "Ntuple contains no column called \"" << EVENT_NUMBER_COLUMN_NAME << "\"" << endl;
		exit(1);
	}
	if ( !eventWeightFound )
	{
		cerr << "Ntuple contains no column called \"" << EVENT_WEIGHT_COLUMN_NAME << "\"" << endl;
		exit(1);
	}

	//Set up access to the other columns
	//You mustn't edit the lengths of these two vectors once you start using SetBranchAddress with them
	currentValues = vector< double >( otherColumnNames.size() );
	branches = vector< TBranch* >( otherColumnNames.size() );
	for ( unsigned int columnIndex = 0; columnIndex < otherColumnNames.size(); columnIndex++ )
	{
		string leafName = otherColumnNames[columnIndex];

		//Map the column name
		columnNameToIndex[ leafName ] = columnIndex;

		//Do the nasty root Ntuple access thing
		wrappedNtuple->SetBranchAddress( leafName.c_str(), &( currentValues[ columnIndex ] ), &( branches[ columnIndex ] ) );
	}

	//Loop over all rows to map event number to row index
	for ( long rowIndex = 0; rowIndex < numberOfRows; rowIndex++ )
	{
		//Load the row
		currentRowNumber = rowIndex;
		wrappedNtuple->GetEvent( rowIndex );

		//Map the event number
		eventNumberToRow[ currentEventNumber ] = rowIndex;
	}
}

//Destructor
InputUETree::~InputUETree()
{
	inputFile->Close();
	delete inputFile;
	delete wrappedNtuple;
}

//Change the Ntuple row being examined
bool InputUETree::ReadRow( long RowIndex )
{
	//Check if we're already there
	if ( currentRowNumber == RowIndex )
	{
		return true;
	}
	else
	{
		//Check if the row index is in range
		if ( RowIndex < 0 || RowIndex >= numberOfRows )
		{
			return false;
		}
		else
		{
			//Load the corresponding row from the file
			currentRowNumber = RowIndex;
			wrappedNtuple->GetEvent( RowIndex );
			return true;
		}
	}
}
bool InputUETree::ReadEvent( UInt_t EventNumber )
{
	//Check if we're already there
	if ( currentEventNumber == EventNumber )
	{
		return true;
	}
	else
	{
		//Look for an event with this number
		eventIterator = eventNumberToRow.find( EventNumber );

		//Check if the event exists
		if ( eventIterator == eventNumberToRow.end() )
		{
			return false;
		}
		else
		{
			//Load the corresponding row from the file
			currentRowNumber = eventIterator->second;
			wrappedNtuple->GetEvent( eventIterator->second );
			return true;
		}
	}
}
bool InputUETree::ReadEvent( UInt_t EventNumber, int FileIndex )
{
        return ReadEvent( EventNumber );
}

//Get the standard event number and weight information
UInt_t InputUETree::EventNumber()
{
	return currentEventNumber;
}
double InputUETree::EventWeight()
{
	return currentEventWeight;
}

//Get any other column value by name
double InputUETree::GetValue( string VariableName )
{
	//Look for the column
	columnIterator = columnNameToIndex.find( VariableName );

	//Check if the column exists
	if ( columnIterator == columnNameToIndex.end() )
	{
		//Not found
		cerr << "Column named \"" << VariableName << "\" not found in Ntuple" << endl;
		exit(1);
	}
	else
	{
		//Found - return corresponding value
		return currentValues[ columnIterator->second ];
	}
}

//Get the number of rows
long InputUETree::NumberOfRows()
{
	return numberOfRows;
}
long InputUETree::CurrentRow()
{
	return currentRowNumber;
}
long InputUETree::CurrentFile()
{
        return 0;
}

//Get the description of the source
string * InputUETree::Description()
{
	return &sourceDescription;
}

//Get the index of the source
int InputUETree::DescriptionIndex()
{
	return sourceDescriptionIndex;
}
