/*
 * Position.h
 *
 *  Created on: Feb 6, 2019
 *      Author: nmalhis
 */

#ifndef POSITION_H_
#define POSITION_H_

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h> /* srand, rand */
#include <vector>
#include <set>
#include <algorithm>
#include <unistd.h> // sleep
#include <math.h> // log
#include <cmath>
#include <time.h>       /* time */
#include "const.h"

namespace std {

class Position {
public:
	int 	loci;
	short 	refi;
	float 	conservation;

	Position();
	Position(const Position &pos);
	Position(const string &line);

	virtual ~Position();
};

} /* namespace std */

#endif /* POSITION_H_ */
