/*
 *   C++ command line argument parser
 *
 *   Copyright (C) 2005 by
 *   Michael Hanke        michael.hanke@gmail.com
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 */
#ifndef __ARGVPARSER_H
#define __ARGVPARSER_H

#include <string>
#include <vector>
#include <map>
#include <list>

namespace CommandLineProcessing
{

/** Provides parsing and storage of POSIX-like command line arguments (argc, argv).
* To use this class for CLI-option parsing, first define a set of valid options using
* the defineOption() method. An option can have several attributes: it can be required
* itself (its missing is considered as an error) and/or it can require a value.
* Options with optional values can be realized by defining the option to not need a
* value and use the syntax '--option=value' or '-o=value' on the command line.
* Every option can have different alternative labels; see defineOptionAlternative().
* A usage description (see usageDescription()) can be generated from the set of defined options (also see the
* addErrorCode(), setIntroductoryDescription() and setHelpOption() methods).
* \note The implemented parsing algorithm requires that all options have to be given at
* the beginning of the command line. Everything on the commandline after the first
* non-option (or option value) is considered as an argument.
* \attention Short option labels (single letter options) must not be digits
* (the option string itself not a possible value).
* Valid syntaxes are:
* \li program --long-option value -sdfgh -u=5 -i 7 --last=val arg1 arg2 arg3
* Here is a small code example:
* \code
* #include
* ArgvParser cmd;  // the command line parser
*
* // init
* cmd.setIntroductoryDescription("This is foo written by bar.");
*
* //define error codes
* cmd.addErrorCode(0, "Success");
* cmd.addErrorCode(1, "Error");
*
* cmd.setHelpOption("h", "help", "Print this help");
*
* cmd.defineOption("version", ArgvParser::NoOptionAttribute, "Be verbose");
* cmd.defineOptionAlternative("verbose","v");
*
* cmd.defineOption("foo", ArgvParser::OptionRequiresValue, "Fooishness. Default value: 0");
*
* // finally parse and handle return codes (display help etc...)
* int result = cmd.parse(argc, argv);
*
* if (result != ArgvParser::NoParserError)
*   cout << cmd.parseErrorDescription(results);
*   exit(1);
* \endcode
*
* \author Michael Hanke
*/
class ArgvParser
{
public:
    typedef int OptionAttributes;
    typedef int ParserResults;
    typedef std::map<std::string, unsigned int> String2KeyMap;
    typedef std::map<unsigned int, OptionAttributes> Key2AttributeMap;
    typedef std::map<unsigned int, std::string> Key2StringMap;
    typedef std::vector<std::string> ArgumentContainer;

    ArgvParser();
    ~ArgvParser();

    /** Attributes for options. */
    enum
    {
        NoOptionAttribute = 0x00,
        OptionRequiresValue = 0x01,
        OptionRequired = 0x02
    };
    /** Return values of the parser. */
    enum
    {
        NoParserError = 0x00,
        ParserUnknownOption = 0x01,
        ParserMissingValue = 0x02,
        ParserOptionAfterArgument = 0x04,
        ParserMalformedMultipleShortOption = 0x08,
        ParserRequiredOptionMissing = 0x16,
        ParserHelpRequested = 0x32
    };

    /** Defines an option with optional attributes (required, ...) and an
    * additional (also optional) description. The description becomes part of the
    * generated usage help that can be requested by calling the usageDescription()
    * method.
    * \return Returns FALSE if there already is an option with this name
    * OR if a short option string (length == 1) is a digit. In that case no
    * action is peformed.
    */
    bool defineOption(const std::string& _name,
                      const std::string& _description = std::string(),
                      OptionAttributes _attributes = NoOptionAttribute);
    /** Define an alternative name for an option that was previously defined by
    * defineOption().
    * \return Returns FALSE if there already is an option with the alternative
    * name or no option with the original name OR if a short option string
    * (length == 1) is a digit. In that case no action is performed.
    */
    bool defineOptionAlternative(const std::string& _original,
                                 const std::string& _alternative);
    /** Returns whether _name is a defined option. */
    bool isDefinedOption(const std::string& _name) const;
    /** Returns whether _name is an option that was found while parsing
    * the command line arguments with the parse() method. In other word: This
    * method returns true if the string is an option AND it was given on the
    * parsed command line.
    */
    bool foundOption(const std::string& _name) const;
    /** Define a help option. If this option is found a special error code is
    * returned by the parse method.
    * \attention If this method is called twice without an intermediate call
    * to the reset() method the previously set help option will remain a valid
    * option but is not detected as the special help option and will therefore
    * not cause the parse() method to return the special help error code.
    * \return Returns FALSE if there already is an option defined that equals
    * the short or long name.
    */
    bool setHelpOption(const std::string& _longname = "h",
                       const std::string& _shortname = "help",
                       const std::string& _descr = "");
    /** Returns the number of read arguments. Arguments are efined as beeing
    * neither options nor option values and are specified at the end of the
    * command line after all options and their values. */
    unsigned int arguments() const;
    /** Returns the Nth argument. See arguments().
    * \return Argument string or an empty string if there was no argument of
    * that id.
    */
    std::string argument(unsigned int _number) const;
    /** Get the complete argument vector. The order of the arguments in the
    * vector is the same as on the commandline.
    */
    const std::vector<std::string>& allArguments() const;
    /** Add an error code and its description to the command line parser.
    * This will do nothing more than adding an entry to the usage description.
    */
    void addErrorCode(int _code, const std::string& _descr = "");
    /** Set some string as a general description, that will be printed before
    * the list of available options.
    */
    void setIntroductoryDescription(const std::string& _descr);
    /** Parse the command line arguments for all known options and arguments.
    * \return Error code with parsing result.
    * \retval NoParserError Everything went fine.
    * \retval ParserUnknownOption Unknown option was found.
    * \retval ParserMissingValue A value to a given option is missing.
    * \retval ParserOptionAfterArgument Option after an argument detected. All
    * options have to given before the first argument.
    * \retval ParserMalformedMultipleShortOption Malformed short option string.
    * \retval ParserRequiredOptionMissing Required option is missing.
    * \retval ParserHelpRequested Help option detected.
    */
    ParserResults parse(int _argc, char ** _argv);
    /** Return the value of an option.
    * \return Value of a commandline options given by the name of the option or
    * an empty string if there was no such option or the option required no
    * value.
    */
    std::string optionValue(const std::string& _option) const;
    /** Reset the parser. Call this function if you want to parse another set of
    * command line arguments with the same parser object.
    */
    void reset();
    /** Returns the name of the option that was responsible for a parser error.
      * An empty string is returned if no error occured at all.
      */
    const std::string& errorOption() const;
    /** This method can be used to evaluate parser error codes and generate a
    * human-readable description. In case of a help request error code the
    * usage description as returned by usageDescription() is printed.
    */
    std::string parseErrorDescription(ParserResults _error_code) const;
    /** Returns a string with the usage descriptions for all options. The
     * description string is formated to fit into a terminal of width _width.*/
    std::string usageDescription(unsigned int _width = 80) const;

private:
    /** Returns the key of a defined option with name _name or -1 if such option
     * is not defined. */
    int optionKey( const std::string& _name ) const;
    /** Returns a list of option names that are all alternative names associated
     * with a single key value.
     */
    std::list<std::string> getAllOptionAlternatives(unsigned int _key) const;

    /** The current maximum key value for an option. */
    unsigned int max_key;
    /** Map option names to a numeric key. */
    String2KeyMap option2key;

    /** Map option key to option attributes. */
    Key2AttributeMap option2attribute;

    /** Map option key to option description. */
    Key2StringMap option2descr;

    /** Map option key to option value. */
    Key2StringMap option2value;

    /** Map error code to its description. */
    std::map<int, std::string> errorcode2descr;

    /** Vector of command line arguments. */
    ArgumentContainer argument_container;

    /** General description to be returned as first part of the generated help page. */
    std::string intro_description;

    /** Holds the key for the help option. */
    unsigned int help_option;

    /** Holds the name of the option that was responsible for a parser error.
    */
    std::string error_option;
}; // class ArgvParser


// Auxillary functions

/** Returns whether the given string is a valid (correct syntax) option string.
 * It has to fullfill the following criteria:
 *  1. minimum length is 2 characters
 *  2. Start with '-'
 *  3. if if minimal length -> must not be '--'
 *  4. first short option character must not be a digit (to distinguish negative numbers)
 */
bool isValidOptionString(const std::string& _string);

/** Returns whether the given string is a valid (correct syntax) long option string.
 * It has to fullfill the following criteria:
 *  1. minimum length is 4 characters
 *  2. Start with '--'
 */
bool isValidLongOptionString(const std::string& _string);

/** Splits option and value string if they are given in the form 'option=value'.
* \return Returns TRUE if a value was found.
*/
bool splitOptionAndValue(const std::string& _string, std::string& _option,
                         std::string& _value);

/** String tokenizer using standard C++ functions. Taken from here:
 * http://gcc.gnu.org/onlinedocs/libstdc++/21_strings/howto.html#3
 * Splits the string _in by _delimiters and store the tokens in _container.
 */
template <typename Container>
void splitString(Container& _container, const std::string& _in,
                 const char* const _delimiters = " \t\n")
{
    const std::string::size_type len = _in.length();
    std::string::size_type i = 0;

    while ( i < len )
    {
        // eat leading whitespace
        i = _in.find_first_not_of (_delimiters, i);
        if (i == std::string::npos)
            return;   // nothing left but white space

        // find the end of the token
        std::string::size_type j = _in.find_first_of (_delimiters, i);

        // push token
        if (j == std::string::npos)
        {
            _container.push_back (_in.substr(i));
            return;
        }
        else
            _container.push_back (_in.substr(i, j-i));

        // set up for next loop
        i = j + 1;
    }
}

/** Returns true if the character is a digit (what else?). */
bool isDigit(const char& _char);

/** Build a vector of integers from a string of the form:
* '1,3-5,14,25-20'. This string will be expanded to a list of positive
* integers with the following elements: 1,3,4,5,14,25,24,23,22,21,20.
* All of the expanded elements will be added to the provided list.
* \return Returns FALSE if there was any syntax error in the given string
* In that case the function stops at the point where the error occured.
* Only elements processed up to that point will be added to the expanded
* list.
* \attention This function can only handle unsigned integers!
*/
bool expandRangeStringToUInt(const std::string& _string,
                             std::vector<unsigned int>& _expanded);
/** Returns a copy of _str with whitespace removed from front and back. */
std::string trimmedString(const std::string& _str);

/** Formats a string of an arbitrary length to fit a terminal of width
* _width and to be indented by _indent columns.
*/
std::string formatString(const std::string& _string,
                         unsigned int _width,
                         unsigned int _indent = 0);

}
; // namespace CommandLineProcessing

#endif // __CMDLINEPARSER_H
