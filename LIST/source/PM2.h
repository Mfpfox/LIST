/*
 * STM.h
 *
 *  Created on: May 7, 2018
 *      Author: nmalhis
 */

#ifndef PM2_H_
#define PM2_H_

#include <string>
#include <string.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <vector>
#include <set>
#include <algorithm>
#include <unistd.h> // sleep
#include <math.h> // log2

#include "const.h"
#include "Protein.h"

namespace std {

class PM2 {
public:
	Protein *protein;
	float center[32];
	float amplitude[32];

	PM2();
	void load_trained(string stmfile);
	float score_PM1(int ps);
	float score_PM2(int ps);

	virtual ~PM2();
};

} /* namespace std */

#endif /* PM2_H_ */
