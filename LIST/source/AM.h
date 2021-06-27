/*
 * AM.h
 *
 *  Created on: May 5, 2018
 *      Author: nmalhis
 */

#ifndef AM_H_
#define AM_H_

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
#include "Protein.h"

namespace std {

class AM {
	float mx[2][20][20];
	float pr[20][20];
	float _sm2[20][20];
	float _sum[2];
public:
	AM();
	void load_trained(string rap_file);
	float score(short r, short a);

	virtual ~AM();

	static short _get_aa_idx(char aa) {
		switch (aa) {
		case 'A':
			return 0;
		case 'R':
			return 1;
		case 'N':
			return 2;
		case 'D':
			return 3;
		case 'C':
			return 4;
		case 'Q':
			return 5;
		case 'E':
			return 6;
		case 'G':
			return 7;
		case 'H':
			return 8;
		case 'I':
			return 9;
		case 'L':
			return 10;
		case 'K':
			return 11;
		case 'M':
			return 12;
		case 'F':
			return 13;
		case 'P':
			return 14;
		case 'S':
			return 15;
		case 'T':
			return 16;
		case 'W':
			return 17;
		case 'Y':
			return 18;
		case 'V':
			return 19;
		}
		return 20;
	}

	static char _get_aa(short aai) {
		switch (aai) {
		case 0:
			return 'A';
		case 1:
			return 'R';
		case 2:
			return 'N';
		case 3:
			return 'D';
		case 4:
			return 'C';
		case 5:
			return 'Q';
		case 6:
			return 'E';
		case 7:
			return 'G';
		case 8:
			return 'H';
		case 9:
			return 'I';
		case 10:
			return 'L';
		case 11:
			return 'K';
		case 12:
			return 'M';
		case 13:
			return 'F';
		case 14:
			return 'P';
		case 15:
			return 'S';
		case 16:
			return 'T';
		case 17:
			return 'W';
		case 18:
			return 'Y';
		case 19:
			return 'V';
		}
		return 'U';
	}

};

} /* namespace std */

#endif /* AM_H_ */
