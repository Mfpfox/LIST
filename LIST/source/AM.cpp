/*
 * AM.cpp
 *
 *  Created on: May 5, 2018
 *      Author: nmalhis
 */

#include "AM.h"

namespace std {

AM::AM() {
	_sum[0] = _sum[1] = 0;
	for (short ri = 0; ri < 20; ri++) {
		for (short ai = 0; ai < 20; ai++) {
			mx[0][ri][ai] = mx[1][ri][ai] = 0;
			pr[ri][ai] = 0;
		}
	}
}

void AM::load_trained(string rap_file) {
	float medians[] = { 5.1999, 6.05004, 4.45313, 5.45221, 3.87966, 4.09218,
			5.18468, 5.13501, 4.35445, 4.48237, 4.45219, 4.29586, 4.38185,
			4.3656, 5.01726, 4.15496, 4.37554, 3.06151, 4.81796, 4.99094 };

	for (short ri = 0; ri < 20; ri++) {
		for (short ai = 0; ai < 20; ai++) {
			_sm2[ri][ai] = 0.00815;
			pr[ri][ai] = medians[ri] / 10; // median Ts set ExAC scores
		}
	}
	string ifile(_ROOT_DATA);
	ifile.append(_LEARNED);
	ifile.append(rap_file);
	// cerr << "AM::load_trained....... : " << ifile << endl;
	ifstream fin;
	fin.open(ifile.c_str());
	if (!fin.is_open()) {
		cerr << "AM::load Error: can not open " << ifile << endl;
		return;
	}
	string junk;
	string line;
	short ri;
	short ai;
	getline(fin, line);
	while (!fin.eof()) {
		stringstream ss(line);
		ss >> ri;
		ss >> ai;
		if (ri >= 0 && ri < 20 && ai >= 0 && ai < 20) {
			ss >> _sm2[ri][ai];
			ss >> pr[ri][ai];
		}
		getline(fin, line);
	}
	fin.close();
}

float AM::score(short ri, short ai) {
	// float ret = 0;
	if (ri >= 0 && ri < 20 && ai >= 0 && ai < 20) {
		return pr[ri][ai];
	}
	return -1;
}

AM::~AM() {
	// TODO Auto-generated destructor stub
}

} /* namespace std */
