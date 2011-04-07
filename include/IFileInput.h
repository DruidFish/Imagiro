/**
  @interface IFileInput

  A general interface for any class that loads data from Root files

  @author Benjamin M Wynne bwynne@cern.ch
  @date 06-04-2011
 */


#ifndef I_FILE_INPUT_H
#define I_FILE_INPUT_H

#include <string>
#include "Rtypes.h"

using namespace std;

class IFileInput
{
	public:
		//Access a particular event, return false if the event is not found
		virtual bool ReadRow( long RowIndex ) = 0;
		virtual bool ReadEvent( UInt_t EventNumber ) = 0;
		virtual bool ReadEvent( UInt_t EventNumber, int FileIndex ) = 0;

		//Get the standard event number and weight information
		virtual UInt_t EventNumber() = 0;
		virtual double EventWeight() = 0;

		//Get any other column value by name
		virtual double GetValue( string VariableName ) = 0;

		//Get the number of rows
		virtual long NumberOfRows() = 0;
		virtual long CurrentRow() = 0;
		virtual long CurrentFile() = 0;

		//Get the description of the source
		virtual string * Description() = 0;
		virtual int DescriptionIndex() = 0;
};

#endif
