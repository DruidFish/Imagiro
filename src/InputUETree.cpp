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
#include <cmath>

const string EVENT_NUMBER_COLUMN_NAME = "EventNumber";
const string EVENT_WEIGHT_COLUMN_NAME = "EventWeight";

//Default constructor - useless
InputUETree::InputUETree()
{
}

//Constructor taking arguments pointing to a particular Ntuple in a root file
InputUETree::InputUETree( string FilePath, string NtuplePath, string Description, unsigned int DescriptionIndex,
	       string CutName, double CutMinimum, double CutMaximum )
{
	sourceDescription = Description;
	sourceDescriptionIndex = DescriptionIndex;

	//Get the Ntuple from the file
	inputFile = new TFile( FilePath.c_str(), "READ" );
	if ( !inputFile->IsOpen() )
	{
		cerr << "Failed to open file " << FilePath << endl;
		exit(1);
	}
	wrappedNtuple = ( TTree* )inputFile->Get( NtuplePath.c_str() );
	if ( !wrappedNtuple )
	{
		cerr << "Failed to retrieve " << NtuplePath << " from " << FilePath << endl;
		exit(1);
	}

	//Loop over all columns/branches to find their names
	TIter branchIterator( wrappedNtuple->GetListOfLeaves() );
	TLeaf * nextLeaf;
	bool eventNumberFound = false;
	bool eventWeightFound = false;
	vector<string> otherColumnNames;
	while ( ( nextLeaf = ( TLeaf* )branchIterator() ) )
	{
		string leafName = nextLeaf->GetName();
		string leafType = nextLeaf->GetTypeName();

		//Find the special columns for EventNumber and EventWeight
		if ( leafType == "UInt_t" && leafName == EVENT_NUMBER_COLUMN_NAME )
		{
			wrappedNtuple->SetBranchAddress( EVENT_NUMBER_COLUMN_NAME.c_str(), &currentEventNumber, &eventNumberBranch );
			eventNumberFound = true;
		}
		else if ( leafType == "Double_t" )
		{
			//Hack at the moment to remove other formats

			if ( leafName == EVENT_WEIGHT_COLUMN_NAME )
			{
				//Event weight
				wrappedNtuple->SetBranchAddress( EVENT_WEIGHT_COLUMN_NAME.c_str(), &currentEventWeight, &eventWeightBranch );
				eventWeightFound = true;
			}
			else
			{
				//All other readable columns
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
	bool calculateCut = !( CutMaximum == 0.0 && CutMinimum == 0.0 && CutName == "NoCut" );
	currentRowNumber = 0;
	currentEventNumber = 0;
	for ( unsigned long rowIndex = 0; rowIndex < wrappedNtuple->GetEntries(); rowIndex++ )
	{
		//Load the row
		wrappedNtuple->GetEvent( rowIndex );

		//Perform a cut if required
		bool storeEvent = true;
		if ( calculateCut )
		{
			double cutValue = GetValue( CutName );
			storeEvent = ( cutValue >= CutMinimum && cutValue < CutMaximum );
		}

		//Make maps to the event
		if ( storeEvent )
		{
			eventNumberToRow[ currentEventNumber ] = rowIndex;
			validRows.push_back( rowIndex );
		}
	}

	//Settle on a valid row
	if ( validRows.size() > 0 )
	{
		 wrappedNtuple->GetEvent( validRows[ 0 ] );
		 currentRowNumber = 0;
	}
}

//Destructor
InputUETree::~InputUETree()
{
	//STL
	eventNumberToRow.clear();
	validRows.clear();
	columnNameToIndex.clear();
	branches.clear();
	currentValues.clear();

	//IO
	delete wrappedNtuple;
	inputFile->Close();
	delete inputFile;
}

//Change the Ntuple row being examined
bool InputUETree::ReadRow( unsigned long RowIndex )
{
	//Check if we're already there
	if ( currentRowNumber == RowIndex )
	{
		return true;
	}
	else
	{
		//Check if the row index is in range
		if ( RowIndex >= validRows.size() )
		{
			return false;
		}
		else
		{
			//Load the corresponding row from the file
			currentRowNumber = RowIndex;
			wrappedNtuple->GetEvent( validRows[ RowIndex ] );
			return true;
		}
	}
}

bool InputUETree::ReadNextRow()
{
	if ( currentRowNumber < validRows.size() - 1 )
	{
		//Load the next row
		currentRowNumber++;
		wrappedNtuple->GetEvent( validRows[ currentRowNumber ] );
		return true;
	}
	else
	{
		//No more rows
		return false;
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
			wrappedNtuple->GetEvent( eventIterator->second );

			//Find out which valid row that corresponds to
			for ( unsigned int rowIndex = 0; rowIndex < validRows.size(); rowIndex++ )
			{
				if ( eventIterator->second == validRows[ rowIndex ] )
				{
					currentRowNumber = rowIndex;
					break;
				}
			}
			return true;
		}
	}
}
bool InputUETree::ReadEvent( UInt_t EventNumber, unsigned int FileIndex )
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
		double storedValue = currentValues[ columnIterator->second ];

		if ( isnan( storedValue ) )
		{
			//cout << "Input NaN converted to 0" << endl;
			storedValue = 0.0;
		}

		return storedValue;
	}
}

//Get the number of rows
unsigned long InputUETree::NumberOfRows()
{
	return validRows.size();
}
unsigned long InputUETree::CurrentRow()
{
	return currentRowNumber;
}
unsigned int InputUETree::CurrentFile()
{
        return 0;
}

//Get the description of the source
string * InputUETree::Description()
{
	return &sourceDescription;
}

//Get the index of the source
unsigned int InputUETree::DescriptionIndex()
{
	return sourceDescriptionIndex;
}
