/*
 * SMSA.cpp
 *
 *  Created on: Oct 25, 2018
 *      Author: nmalhis
 */

#include "SMSA.h"

namespace std {


SMSA::SMSA() {
  SI = 0;
  taxa_id = 0;
  ST = 0;
  matches = 0;
  identity = 0;
  mm = 0;
  valid = false;
  sz = 0;
  first = 0;
  last = 0;
}

SMSA::SMSA(const SMSA &smsa) {
  SI = smsa.SI;
  taxa_id = smsa.taxa_id;
  ST = smsa.ST;
  matches = smsa.matches;
  identity = smsa.identity;
  mm = smsa.mm;
  valid = smsa.valid;
  sz = smsa.sz;
  first = smsa.first;
  last = smsa.last;
  ac.assign(smsa.ac);
  seq.assign(smsa.seq);
}

SMSA::SMSA(const string &line) {
  seq.assign("");// loaded next
  ac.assign(""); // loaded next
  valid = false; // loaded later
  ST = 0;        // loaded next
  matches = 0;
  identity = 0;
  mm = 0;
  taxa_id = 0;
  stringstream ss(line);
  ss >> ac;
  ss >> taxa_id;
  ss >> matches;
  ss >> mm;
  ss >> seq;
  for (first = 0; first < seq.size(); first++) {
    if (seq.at(first) != '.')
	break;
  }
  for (last = seq.size() - 1; last >= 0; last--) {
    if (seq.at(last) != '.')
	break;
  }
  sz = matches + mm;
  identity = float(matches) / sz;
  SI = short(identity * 100);
}

SMSA::SMSA(const string &line, const short &typ) {
  // ================================= not updated: matches and mm
  seq.assign("");
  ac.assign("");
  SI = 0;
  valid = false;
  sz = 0;
  ST = 0;
  matches = 0;
  identity = 0;
  mm = 0;
  taxa_id = 0;
  first = 0;
  last = 0;

  stringstream ss(line);
  ss >> ac;
  ss >> seq;
  ss >> matches;
  ss >> mm;
  ss >> seq;
  sz = matches + mm;
  identity = float(matches) / sz;
}

void std::SMSA::set(const SMSA& smsa) {
  SI = smsa.SI;
  taxa_id = smsa.taxa_id;
  ST = smsa.ST;
  matches = smsa.matches;
  identity = smsa.identity;
  mm = smsa.mm;
  valid = smsa.valid;
  sz = smsa.sz;
  first = smsa.first;
  last = smsa.last;
  ac.assign(smsa.ac);
  seq.assign(smsa.seq);
}

void std::SMSA::set(const string& line) {
  valid = false; // loaded later
  stringstream ss(line);
  ss >> ac;
  ss >> ST;
  matches = 0;
  mm = 0;
  ss >> seq;
  for (first = 0; first < seq.size(); first++) {
    if (seq.at(first) != '.')
	break;
  }
  for (last = seq.size() - 1; last >= 0; last--) {
    if (seq.at(last) != '.')
	break;
  }
}

void SMSA::clear() {
  SI = 0;
  taxa_id = 0;
  ST = 0;
  matches = 0;
  identity = 0;
  mm = 0;
  valid = false;
  sz = 0;
  first = 0;
  last = 0;
  ac.assign("");
  seq.assign("");
}

SMSA::~SMSA() {
	// TODO Auto-generated destructor stub
}

} /* namespace std */
