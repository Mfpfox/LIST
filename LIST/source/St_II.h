/*
 * II.h
 *
 *  Created on: Apr 12, 2018
 *      Author: nmalhis
 */

#ifndef ST_II_H_
#define ST_II_H_

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <vector>
#include <set>
#include <algorithm>
#include <unistd.h> // sleep
#include <math.h> // log2

namespace std {

class St_II {
public:
	int ikey;
	int i2;
	// int lnum;
	St_II();
	St_II(const St_II &ii);

	bool operator<(const St_II &ii) const {
		return ikey < ii.ikey;
	}
	bool operator==(const St_II &ii) const {
		return ikey == ii.ikey;
	}
	virtual ~St_II();
};

} /* namespace std */

#endif /* ST_II_H_ */
