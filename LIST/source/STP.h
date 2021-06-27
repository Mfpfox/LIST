/*
 * STP.h
 *
 *  Created on: Apr 25, 2018
 *      Author: nmalhis
 *
 *  This class was created (Apr. 25 2018) to  be used by:
 *  Alignment::void _2_smsa_to_stp(string in_path, string out_path, string list_file, int radios, int idCount);
 */

#ifndef STP_H_
#define STP_H_

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

namespace std {

class STP {
	void _split_st3(string &st3, short &a_st, short &a_i, int &a_c);
	void _split_st2(string &st2, short &r_i, int &r_c);
public:
	int loci;
	char raa;
	// char msk;
	short ri;
	short ai_max;
	int depth;
	short a_staxa[20];
	short a_lIdentity[20];
	short a_sIdentity[20];
	int a_count[20];
	short stp_identity[32];
	int stp_count[32];

	static float pam_score[20];

	STP();
	STP(const STP &stp);
	STP(const string &line);

	float compute_PAM_scores();
	int get_alleles_cover();
	string to_string();

	virtual ~STP();
};

} /* namespace std */

#endif /* STP_H_ */
