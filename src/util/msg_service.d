/**
 * msg_service.d
 * Home to some commonly used messages in assert statements
 * and error messages.
 *
 * Author: Rowan J. Gollan
 */

module util.msg_service;

import std.string;

pure string brokenPreCondition(string variable,
			       size_t lineNo = __LINE__,
			       string fileName = __FILE__)
{
    return format("Pre-condition contract broken for %s on line %d in file %s\n",
		  variable, lineNo, fileName);
}

pure string brokenPostCondition(string variable,
				size_t lineNo = __LINE__,
				string fileName = __FILE__)
{
    return format("Post-condition contract broken for %s on line %d in file %s\n",
		  variable, lineNo, fileName);
}

string failedUnitTest(size_t lineNo = __LINE__,
		      string fileName = __FILE__)
{
    return format("Unit test failure on line %d in file %s\n",
		  lineNo, fileName);
}
