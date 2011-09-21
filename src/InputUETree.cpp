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
InputUETree::InputUETree( string FilePath, string NtuplePath, string Description, unsigned int DescriptionIndex, ObservableList * RelevanceChecker,
		string CutName, double CutMinimum, double CutMaximum )
{
	sourceDescription = Description;
	sourceDescriptionIndex = DescriptionIndex;
	totalRows = 0;
	currentRowNumber = 0;
	currentEventNumber = 0;
	bool calculateCut = !( CutMaximum == 0.0 && CutMinimum == 0.0 && CutName == "NoCut" );

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
	bool cutColumnFound = false;
	vector< string > otherColumnNames, vectorNames;
	unsigned int cutColumnIndex = 0;
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
			//Check the name
			if ( leafName == EVENT_WEIGHT_COLUMN_NAME )
			{
				//Event weight
				wrappedNtuple->SetBranchAddress( EVENT_WEIGHT_COLUMN_NAME.c_str(), &currentEventWeight, &eventWeightBranch );
				eventWeightFound = true;
			}
			else if ( calculateCut && leafName == CutName )
			{
				//The column with the cut variable
				cutColumnFound = true;
				cutColumnIndex = otherColumnNames.size();
				otherColumnNames.push_back( leafName );
			}
			else if ( RelevanceChecker->IsInList( leafName ) )
			{
				//All other relevant columns
				otherColumnNames.push_back( leafName );
			}
		}
		else if ( leafType == "vector<double>" )
		{

			if ( RelevanceChecker->IsInList( leafName ) )
			{
				//All relevant vectors
				vectorNames.push_back( leafName );
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

	//Set up the column with the cut variable (if there is one)
	//You mustn't edit the lengths of these two vectors once you start using SetBranchAddress with them
	currentValues = vector< double >( otherColumnNames.size(), 0 );
	valueBranches = vector< TBranch* >( otherColumnNames.size(), 0 );
	if ( calculateCut )
	{
		if ( cutColumnFound )
		{
			columnNameToIndex[ CutName ] = cutColumnIndex;
			wrappedNtuple->SetBranchAddress( CutName.c_str(), &( currentValues[ cutColumnIndex ] ), &( valueBranches[ cutColumnIndex ] ) );
		}
		else
		{
			cerr << "Ntuple contains no column called \"" << CutName << "\"" << endl;
			exit(1);
		}
	}

	//Loop over all rows to map event number to row index
	for ( unsigned long rowIndex = 0; rowIndex < wrappedNtuple->GetEntries(); rowIndex++ )
	{
		//Load the row
		//wrappedNtuple->GetEvent( rowIndex );
		eventNumberBranch->GetEvent( rowIndex );

		//Perform a cut if required
		bool storeEvent = true;
		if ( calculateCut )
		{
			valueBranches[ cutColumnIndex ]->GetEvent( rowIndex );
			double cutValue = currentValues[ cutColumnIndex ];//GetValue( CutName );
			storeEvent = ( cutValue >= CutMinimum && cutValue <= CutMaximum );
		}

		//Make maps to the event
		if ( storeEvent )
		{
			currentRowNumber = totalRows;
			eventNumberToExternalRow[ currentEventNumber ] = totalRows;
			externalRowToInternalRow[ totalRows ] = rowIndex;
			totalRows++;
		}
	}

	//Set up access to the other columns
	for ( unsigned int columnIndex = 0; columnIndex < otherColumnNames.size(); columnIndex++ )
	{
		if ( !calculateCut || columnIndex != cutColumnIndex )
		{
			string leafName = otherColumnNames[columnIndex];

			//Map the column name
			columnNameToIndex[ leafName ] = columnIndex;

			//Do the nasty root Ntuple access thing
			wrappedNtuple->SetBranchAddress( leafName.c_str(), &( currentValues[ columnIndex ] ), &( valueBranches[ columnIndex ] ) );
		}
	}
	currentVectors = vector< vector< double >* >( vectorNames.size(), 0 );
	vectorBranches = vector< TBranch* >( vectorNames.size(), 0 );
	for ( unsigned int columnIndex = 0; columnIndex < vectorNames.size(); columnIndex++ )
	{
		string leafName = vectorNames[columnIndex];

		//Map the column name
		vectorNameToIndex[ leafName ] = columnIndex;

		//Do the nasty root Ntuple access thing
		wrappedNtuple->SetBranchAddress( leafName.c_str(), &( currentVectors[ columnIndex ] ), &( vectorBranches[ columnIndex ] ) );
	}

	//Load a valid row
	if ( totalRows > 0 )
	{
		wrappedNtuple->GetEvent( externalRowToInternalRow[ 0 ] );
		currentRowNumber = 0;
	}
}

//Destructor
InputUETree::~InputUETree()
{
	//STL
	eventNumberToExternalRow.clear();
	externalRowToInternalRow.clear();
	columnNameToIndex.clear();
	valueBranches.clear();
	currentValues.clear();
	vectorBranches.clear();
	vectorNameToIndex.clear();
	currentVectors.clear();

	//IO
	delete wrappedNtuple;
	inputFile->Close();
	delete inputFile;
}

//Change the Ntuple row being examined
bool InputUETree::ReadRow( unsigned long RowIndex, unsigned int FileIndex )
{
	//Stupidity check
	if ( FileIndex != 0 )
	{
		cerr << "Asking for file " << FileIndex << " in InputUETree which only holds file 0" << endl;
		exit(1);
	}

	//Check if we're already there
	if ( currentRowNumber == RowIndex )
	{
		return true;
	}
	else
	{
		//Check if the row index is in range
		if ( RowIndex >= totalRows )
		{
			return false;
		}
		else
		{
			//Load the corresponding row from the file
			currentRowNumber = RowIndex;
			unsigned long loadRow = externalRowToInternalRow[ currentRowNumber ];

			//Slow
			//wrappedNtuple->GetEvent( loadRow );

			//Fast
			eventNumberBranch->GetEvent( loadRow );
			eventWeightBranch->GetEvent( loadRow );
			for ( unsigned int branchIndex = 0; branchIndex < valueBranches.size(); branchIndex++ )
			{
				valueBranches[ branchIndex ]->GetEvent( loadRow );
			}
			for ( unsigned int branchIndex = 0; branchIndex < vectorBranches.size(); branchIndex++ )
			{
				vectorBranches[ branchIndex ]->GetEvent( loadRow );
			}

			return true;
		}
	}
}

bool InputUETree::ReadEvent( UInt_t EventNumber, unsigned int FileIndex )
{
	//Stupidity check
	if ( FileIndex != 0 )
	{
		cerr << "Asking for file " << FileIndex << " in InputUETree which only holds file 0" << endl;
		exit(1);
	}

	//Check if we're already there
	if ( currentEventNumber == EventNumber )
	{
		return true;
	}
	else
	{
		//Look for an event with this number
		eventIterator = eventNumberToExternalRow.find( EventNumber );

		//Check if the event exists
		if ( eventIterator == eventNumberToExternalRow.end() )
		{
			return false;
		}
		else
		{
			//Load the corresponding row from the file
			currentRowNumber = eventIterator->second;
			unsigned long loadRow = externalRowToInternalRow[ currentRowNumber ];

			//Slow
			//wrappedNtuple->GetEvent( loadRow );

			//Fast
			eventNumberBranch->GetEvent( loadRow );
			eventWeightBranch->GetEvent( loadRow );
			for ( unsigned int branchIndex = 0; branchIndex < valueBranches.size(); branchIndex++ )
			{
				valueBranches[ branchIndex ]->GetEvent( loadRow );
			}
			for ( unsigned int branchIndex = 0; branchIndex < vectorBranches.size(); branchIndex++ )
			{
				vectorBranches[ branchIndex ]->GetEvent( loadRow );
			}

			return true;
		}
	}
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
vector< double > * InputUETree::GetVector( string VectorName )
{
	//Look for the vector
	columnIterator = vectorNameToIndex.find( VectorName );

	//Check if the column exists
	if ( columnIterator == vectorNameToIndex.end() )
	{
		//Not found
		cerr << "Vector named \"" << VectorName << "\" not found in Ntuple" << endl;
		exit(1);
	}
	else
	{
		//Found - return pointer to vector
		return currentVectors[ columnIterator->second ];
	}
}

//Get the number of rows and files
unsigned long InputUETree::NumberOfRows()
{
	return totalRows;
}
unsigned long InputUETree::CurrentRow()
{
	return currentRowNumber;
}
unsigned int InputUETree::NumberOfFiles()
{
	return 1;
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
