/*
 * STP.cpp
 *
 *  Created on: Apr 25, 2018
 *      Author: nmalhis
 */

#include "STP.h"

namespace std {

STP::STP() {
	loci = 0;
	raa = 0;
	// msk = 0;
	ri = 0;
	ai_max = 0;
	depth = 0;
	for (short i = 0; i < 20; i++) {
		a_staxa[i] = 0;
		a_lIdentity[i] = 0;
		a_sIdentity[i] = 0;
		a_count[i] = 0;
	}
	for (short sti = 0; sti < 32; sti++) {
		stp_identity[sti] = 0;
		stp_count[sti] = 0;
	}
}

STP::STP(const STP &stp) {
	raa = stp.raa;
	ri = stp.ri;
	loci = stp.loci;
	ai_max = stp.ai_max;
	depth = stp.depth;
	for (short i = 0; i < 20; i++) {
		a_staxa[i] = stp.a_staxa[i];
		a_lIdentity[i] = stp.a_lIdentity[i];
		a_sIdentity[i] = stp.a_sIdentity[i];
		a_count[i] = stp.a_count[i];
	}
	for (short sti = 0; sti < 32; sti++) {
		stp_identity[sti] = stp.stp_identity[sti];
		stp_count[sti] = stp.stp_count[sti];
	}
}

STP::STP(const string& line) {
	string tmps;
	raa = 0;
	ri = 20;
	loci = 0;
	ai_max = 0;
	depth = 0;
	stringstream ss(line);
	ss >> loci;
// #ifndef _V2
//		ss >> tmps;
// #endif
	ss >> raa;
	ss >> ri;
	ss >> ai_max;
	ss >> depth;
	/*
#ifndef _V2
		for (short i = 0; i < 3; i++) {
			ss >> tmps;
		}
#endif//*/
	for (short i = 0; i < 20; i++) {
		ss >> tmps;
		_split_st3(tmps, a_staxa[i], a_lIdentity[i], a_count[i]);
	}
	for (short sti = 0; sti < 32; sti++) {
		ss >> tmps;
		_split_st2(tmps, stp_identity[sti], stp_count[sti]);
	}
}

void STP::_split_st3(string& st3, short &a_st, short &a_i, int &a_c) {
	a_st = 0;
	a_i = 0;
	a_c = 0;
	short cut1 = st3.find(':');
	short cut2 = st3.find(':', cut1 + 1);
	stringstream ss1(st3.substr(0, cut1));
	ss1 >> a_st;
	cut1++;
	stringstream ss2(st3.substr(cut1, (cut2 - cut1)));
	ss2 >> a_i;
	cut2++;
	stringstream ss3(st3.substr(cut2));
	ss3 >> a_c;
}

void STP::_split_st2(string& st2, short & r_i, int& r_c) {
	r_i = 0;
	r_c = 0;
	short cut1 = st2.find(':');
	stringstream ss1(st2.substr(0, cut1));
	ss1 >> r_i;
	cut1++;
	stringstream ss2(st2.substr(cut1));
	ss2 >> r_c;
}

float STP::compute_PAM_scores() {
	float ex20;
	//*
	for (short ai = 0; ai < 20; ai++) {
		if (a_lIdentity[ai] >= _PAM_IDEN_MIN) {
			pam_score[ai] = 1 - float(a_staxa[ai]) / 31;
		} else {
			pam_score[ai] = 1;
		}
	}
	ex20 = 0;
	for (short ai = 0; ai < 20; ai++) {
		if (ai != ri) {
			ex20 += (pam_score[ai] / 19);
		}
	}
	return ex20;
}

int STP::get_alleles_cover() {
	int ret = 0;
	for (short i = 0; i < 20; i++) {
		if (i != ri) {
			ret += a_count[i];
		}
	}
	return ret;
}

string STP::to_string() {
	stringstream ss;
	ss << loci << "\t" << raa << "\t" << ri << "\t" << ai_max << "\t" << depth;
	for (short i = 0; i < 20; i++) {
		ss << "\t" << a_staxa[i] << ":" << a_lIdentity[i] << ":" << a_count[i];
	}
	for (short sti = 0; sti < 32; sti++) {
		ss << "\t" << stp_identity[sti] << ":" << stp_count[sti];
	}
	return ss.str();
}

STP::~STP() {
	// TODO Auto-generated destructor stub
}

} /* namespace std */
