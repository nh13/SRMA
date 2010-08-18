#ifndef BERROR_H_
#define BERROR_H_

#define BREAK_LINE "************************************************************\n"

/*! @enum
  @abstract              the type of action to be taken
  @field  Exit            exit the program
  @field  Warn            print a warning 
  @field  LastActionType  dummy action type
*/
enum {Exit, Warn, LastActionType};

/*! @enum
  @abstract                   the type of error
  @field  CommandLineArgument  improper command line argument
  @field  OutOfRange           value was out of range
  @field  ReallocMemory        memory re-allocation failure
  @field  MallocMemory         memory allocation failure
  @field  OpenFileError        could not open a file
  @field  ReadFileError        could not read from a file
  @field  WriteFileError       could not write from a file
  @field  EndOfFile            reached the end-of-file prematurely
  @field  ThreadError          error starting/joining threads
  @field  LastErrorType        dummy error type 
*/
enum {
	OutOfRange, 
	CommandLineArgument,
	ReallocMemory,
	MallocMemory,
	OpenFileError,
	ReadFileError,
	WriteFileError,
	EndOfFile,
	ThreadError,
	LastErrorType,
};

/*! @function
  @abstract              process an error based on the given action
  @param  function_name  the function name reporting the error
  @param  variable_name  the variable name or value associated with the error
  @param  action_type    the action to be taken
  @param  error_type     the error type 
*/
void srma_error(const char *function_name, const char *variable_name, int action_type, int error_type);

#endif
