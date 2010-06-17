/*
 * Sqt2Pin.h
 *
 *  Created on: Jun 16, 2010
 *      Author: lukask
 */

#ifndef SQT2PIN_H_
#define SQT2PIN_H_

#include "SqtReader.h"


class Sqt2Pin {
public:
	Sqt2Pin();
	virtual ~Sqt2Pin();
	string greeter();
	bool parseOptions(int argc, char **argv);
	int run();

protected:
	ParseOptions parseOptions;
	string tokyoCabinetTmpFN;
};

int main(int argc, char **argv);

#endif /* SQT2PIN_H_ */
