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
const char * TRIGGER_NAMES[] = { "L1_J5", "L1_J15", "L1_J30", "L1_J55", "L1_J75", "L1_J95" };
const double TRIGGER_LOWER_BOUNDS[] = { 60000.0, 110000.0, 160000.0, 210000.0, 260000.0, 310000.0 };

//Default constructor - useless
TriggerChoosingInput::TriggerChoosingInput()
{
}

//Constructor taking arguments pointing to a particular Ntuple in a root file
TriggerChoosingInput::TriggerChoosingInput( string FilePath, string NtuplePath, string Description, unsigned int DescriptionIndex )
{
	sourceDescription = Description;
	sourceDescriptionIndex = DescriptionIndex;

	//Open the input file corresponding to each trigger name
	vector< string > allTriggers( TRIGGER_NAMES, TRIGGER_NAMES + sizeof( TRIGGER_NAMES ) / sizeof( char* ) );
	vector< double > triggerLowerBounds( TRIGGER_LOWER_BOUNDS, TRIGGER_LOWER_BOUNDS + sizeof( TRIGGER_LOWER_BOUNDS ) / sizeof( double ) );
	for ( unsigned int triggerIndex = 0; triggerIndex < allTriggers.size(); triggerIndex++ )
	{
		//Load the file for this trigger
		currentFileNumber = triggerIndex;
		string inputTriggerFilePath = ReplaceString( FilePath, REPLACE_FOR_TRIGGER_IN_PATH, allTriggers[ triggerIndex ] );
		IFileInput * inputTriggerFile = new InputUETree( inputTriggerFilePath, NtuplePath, Description, DescriptionIndex );
		triggerInputs.push_back( inputTriggerFile );

		//Select all events from the file with the correct lead jet pT
		for ( unsigned int rowIndex = 0; rowIndex < inputTriggerFile->NumberOfRows(); rowIndex++ )
		{
			//Load the event
			inputTriggerFile->ReadRow( rowIndex );
			double leadJetPt = inputTriggerFile->GetValue( LEAD_JET_PT_COLUMN_NAME );

			//Reject the event if the pT is not in range for this trigger
			if ( triggerIndex == allTriggers.size() - 1 )
			{
				//No upper bound on highest trigger
				if ( leadJetPt < triggerLowerBounds[ triggerIndex ] )
				{
					continue;
				}
			}
			else
			{
				if ( leadJetPt < triggerLowerBounds[ triggerIndex ] || leadJetPt > triggerLowerBounds[ triggerIndex + 1 ] )
				{
					continue;
				}
			}

			//Make a mapping to the event
			currentRowNumber = eventNumbers.size();
			currentEventNumber = inputTriggerFile->EventNumber();
			eventNumbers.push_back( currentEventNumber );
			eventNumberToFile[ currentEventNumber ] = currentFileNumber;
			eventNumberToRow[ currentEventNumber ] = currentRowNumber;
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
	//Check if we're already there
	if ( currentRowNumber == RowIndex )
	{
		return true;
	}
	else
	{
		//Check if the row index is in range
		if ( RowIndex < 0 || RowIndex >= eventNumbers.size() )
		{
			return false;
		}
		else
		{
			//Load the corresponding row from the file
			currentRowNumber = RowIndex;
			currentEventNumber = eventNumbers[ RowIndex ];
			currentFileNumber = eventNumberToFile[ currentEventNumber ];
			triggerInputs[ currentFileNumber ]->ReadEvent( currentEventNumber );
			return true;
		}
	}
}
bool TriggerChoosingInput::ReadEvent( UInt_t EventNumber )
{
	//Check if we're already there
	if ( currentEventNumber == EventNumber )
	{
		return true;
	}
	else
	{
		//Look for an event with this number
		fileIterator = eventNumberToFile.find( EventNumber );

		//Check if the event exists
		if ( fileIterator == eventNumberToFile.end() )
		{
			return false;
		}
		else
		{
			//Load the corresponding row from the file
			currentRowNumber = eventNumberToRow[ EventNumber ];
			currentEventNumber = EventNumber;
			currentFileNumber = eventNumberToFile[ EventNumber ];
			triggerInputs[ currentFileNumber ]->ReadEvent( EventNumber );
			return true;
		}
	}
}
bool TriggerChoosingInput::ReadEvent( UInt_t EventNumber, unsigned int FileIndex )
{
	return ReadEvent( EventNumber );
}

//Get the standard event number and weight information
UInt_t TriggerChoosingInput::EventNumber()
{
	return currentEventNumber;
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
	return eventNumbers.size();
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

