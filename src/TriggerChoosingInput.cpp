/**
  @class TriggerChoosingInput

  UE ANALYSIS SPECIFIC HARD CODE
  An interface to pick events from different files, according to LeadJetPt

  @author Benjamin M Wynne bwynne@cern.ch
  @date 06-05-2011
 */

#include "TriggerChoosingInput.h"
#include "InputUETree.h"
#include "TLeaf.h"
#include <iostream>
#include <cstdlib>
#include <cmath>

const string EVENT_NUMBER_COLUMN_NAME = "EventNumber";
const string EVENT_WEIGHT_COLUMN_NAME = "EventWeight";
const string LEAD_JET_PT_COLUMN_NAME = "LeadJetPt";
const string REPLACE_FOR_TRIGGER_IN_PATH = "TriggerName";
const char * TRIGGER_NAMES[] = { "L1_MBTS_1", "L1_J5", "L1_J15", "L1_J30", "L1_J55", "L1_J75", "L1_J95" };
const double TRIGGER_LOWER_BOUNDS[] = { 20000.0, 60000.0, 110000.0, 160000.0, 210000.0, 260000.0, 310000.0, 1E10 };

//Default constructor - useless
TriggerChoosingInput::TriggerChoosingInput()
{
}

//Constructor taking arguments pointing to a particular Ntuple in a root file
TriggerChoosingInput::TriggerChoosingInput( string FilePath, string NtuplePath, string Description, unsigned int DescriptionIndex )
{
	sourceDescription = Description;
	sourceDescriptionIndex = DescriptionIndex;
	totalRows = 0;

	//Check we can actually do the find-and-replace on the input path name
	if ( FilePath.find( REPLACE_FOR_TRIGGER_IN_PATH ) == string::npos )
	{
		//Can't find and replace - only one input file
		IFileInput * inputFile = new InputUETree( FilePath, NtuplePath, Description, DescriptionIndex );
		triggerInputs.push_back( inputFile );

		//Book-keeping
		totalRows = inputFile->NumberOfRows();
		currentRowNumber = inputFile->CurrentRow();
		currentFileNumber = 0;
		rowInCurrentFile = inputFile->CurrentRow();
	}
	else
	{
		//Open the input file corresponding to each trigger name
		vector< string > allTriggers( TRIGGER_NAMES, TRIGGER_NAMES + sizeof( TRIGGER_NAMES ) / sizeof( char* ) );
		vector< double > triggerLowerBounds( TRIGGER_LOWER_BOUNDS, TRIGGER_LOWER_BOUNDS + sizeof( TRIGGER_LOWER_BOUNDS ) / sizeof( double ) );
		currentRowNumber = 0;
		rowInCurrentFile = 0;
		for ( unsigned int triggerIndex = 0; triggerIndex < allTriggers.size(); triggerIndex++ )
		{
			//Find and replace the wildcard to make the input file name
			string inputTriggerFilePath = ReplaceString( FilePath, REPLACE_FOR_TRIGGER_IN_PATH, allTriggers[ triggerIndex ] );

			//Load the file for this trigger
			IFileInput * inputTriggerFile = new InputUETree( inputTriggerFilePath, NtuplePath, Description, DescriptionIndex,
					LEAD_JET_PT_COLUMN_NAME, triggerLowerBounds[ triggerIndex ], triggerLowerBounds[ triggerIndex + 1 ] );
			triggerInputs.push_back( inputTriggerFile );

			//Book-keeping
			totalRows += inputTriggerFile->NumberOfRows();
			currentFileNumber = triggerIndex;
			rowInCurrentFile = inputTriggerFile->CurrentRow();
			if ( triggerIndex == allTriggers.size() - 1 )
			{
				currentRowNumber += inputTriggerFile->CurrentRow();
			}
			else
			{
				currentRowNumber += inputTriggerFile->NumberOfRows();
			}
		}
	}
}

//Destructor
TriggerChoosingInput::~TriggerChoosingInput()
{
	for ( unsigned int triggerIndex = 0; triggerIndex < triggerInputs.size(); triggerIndex++ )
	{
		delete triggerInputs[ triggerIndex ];
	}
}

//Change the Ntuple row being examined
bool TriggerChoosingInput::ReadRow( unsigned long RowIndex )
{
	//Check if we are already there
	if ( RowIndex == currentRowNumber )
	{
		return true;
	}
	//Check if the row index is in range
	else if ( RowIndex >= totalRows )
	{
		return false;
	}
	else
	{
		//Load the corresponding row from the file
		currentRowNumber = RowIndex;
		for ( unsigned int triggerIndex = 0; triggerIndex < triggerInputs.size(); triggerIndex++ )
		{
			if ( RowIndex >= triggerInputs[ triggerIndex ]->NumberOfRows() )
			{
				RowIndex -= triggerInputs[ triggerIndex ]->NumberOfRows();
			}
			else
			{
				currentFileNumber = triggerIndex;
				rowInCurrentFile = RowIndex;
				return triggerInputs[ triggerIndex ]->ReadRow( RowIndex );
			}
		}
	}

	//Stop -Wall complaining
	return false;
}

bool TriggerChoosingInput::ReadNextRow()
{
	//Check if there are any more rows
	if ( currentRowNumber < totalRows - 1 )
	{
		//Check if the next row is within the current file
		if ( rowInCurrentFile < triggerInputs[ currentFileNumber ]->NumberOfRows() - 1 )
		{
			//Read the next row
			currentRowNumber++;
			rowInCurrentFile++;
			return triggerInputs[ currentFileNumber ]->ReadRow( rowInCurrentFile );
		}
		else
		{
			//Check if there are more files
			if ( currentFileNumber < triggerInputs.size() - 1 )
			{
				//Read the first row of the next file
				currentFileNumber++;
				currentRowNumber++;
				rowInCurrentFile = 0;
				return triggerInputs[ currentFileNumber ]->ReadRow( rowInCurrentFile );
			}
			else
			{
				//No more files - shouldn't get here
				return false;
			}
		}
	}
	else
	{
		return false;
	}
}

bool TriggerChoosingInput::ReadEvent( UInt_t EventNumber )
{
	//Check for this event in each input
	unsigned long newRowNumber = 0;
	for ( unsigned int triggerIndex = 0; triggerIndex < triggerInputs.size(); triggerIndex++ )
	{
		if ( triggerInputs[ triggerIndex ]->ReadEvent( EventNumber ) )
		{
			//Found the event in this file
			rowInCurrentFile = triggerInputs[ triggerIndex ]->CurrentRow();
			currentRowNumber = newRowNumber + triggerInputs[ triggerIndex ]->CurrentRow();
			currentFileNumber = triggerIndex;
			return true;
		}
		else
		{
			newRowNumber += triggerInputs[ triggerIndex ]->NumberOfRows();
		}
	}

	//If you get this far, the event number was not found
	return false;
}
bool TriggerChoosingInput::ReadEvent( UInt_t EventNumber, unsigned int FileIndex )
{
	return ReadEvent( EventNumber );
}

//Get the standard event number and weight information
UInt_t TriggerChoosingInput::EventNumber()
{
	return triggerInputs[ currentFileNumber ]->EventNumber();
}
double TriggerChoosingInput::EventWeight()
{
	return triggerInputs[ currentFileNumber ]->EventWeight();
}

//Get any other column value by name
double TriggerChoosingInput::GetValue( string VariableName )
{
	return triggerInputs[ currentFileNumber ]->GetValue( VariableName );
}

//Get the number of rows
unsigned long TriggerChoosingInput::NumberOfRows()
{
	return totalRows;
}
unsigned long TriggerChoosingInput::CurrentRow()
{
	return currentRowNumber;
}
unsigned int TriggerChoosingInput::CurrentFile()
{
	return currentFileNumber;
}

//Get the description of the source
string * TriggerChoosingInput::Description()
{
	return &sourceDescription;
}

//Get the index of the source
unsigned int TriggerChoosingInput::DescriptionIndex()
{
	return sourceDescriptionIndex;
}

//Replace any instances of a particular character in a string
string TriggerChoosingInput::ReplaceString( string Input, string FindString, string ReplaceString )
{
	string output;

	//Search the input character by character
	for ( unsigned int characterIndex = 0; characterIndex < Input.size(); characterIndex++ )
	{
		//Find out if the full search string is present
		bool isInstance = true;
		for ( unsigned int testIndex = 0; testIndex < FindString.size(); testIndex++ )
		{
			if ( Input[ characterIndex + testIndex ] != FindString[testIndex] )
			{
				isInstance = false;
				break;
			}
		}

		//Replace or push back
		if (isInstance)
		{
			for ( unsigned int replaceIndex = 0; replaceIndex < ReplaceString.size(); replaceIndex++ )
			{
				output.push_back( ReplaceString[replaceIndex] );
			}

			characterIndex += FindString.size() - 1;
		}
		else
		{
			output.push_back( Input[characterIndex] );
		}
	}

	return output;
}

