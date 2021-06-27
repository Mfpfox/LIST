/*
 * Protein.h
 *
 *  Created on: May 5, 2018
 *      Author: nmalhis
 */

#ifndef PROTEIN_H_
#define PROTEIN_H_

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
#include "Position.h"
#include "STP.h"
// #include "St_SI.h"

namespace std {

class Protein {
public:

	string stp_path;

	string seqFile;
	string seq;

	vector <Position> positionsVec;
	vector <STP> stpVec;

	Protein();
	void clear();
	int load_seq(string _ac);
	int load_stp(string inSTP);

	virtual ~Protein();

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
			// default:
			// cerr << "Not an AA:\t" << aa << endl;
		}
		return 20;
	}

};

} /* namespace std */

#endif /* PROTEIN_H_ */
