/*
 * STM.cpp
 *
 *  Created on: May 7, 2018
 *      Author: nmalhis
 */

#include "PM2.h"

namespace std {
  
  PM2::PM2() {
    protein = NULL;
    // pm2_iden_min = 0;
  }
  
  float PM2::score_PM2(int ps) {
    float ret = 0;
    float inst;
    for (short sti = 0; sti < 32; sti++) {
      if (protein->stpVec.at(ps).stp_identity[sti] > _PM2_IDEN_MIn
	  || (protein->stpVec.at(ps).stp_identity[sti] == _PM2_IDEN_MIn
	      && protein->stpVec.at(ps).stp_count[sti] >= 2)) {
	inst = log10(1 + protein->stpVec.at(ps).stp_identity[sti]);
      } else {
	inst = 0;
      }
      ret += (inst - center[sti]) * amplitude[sti];
    }
    return ret;
  }
  
  void PM2::load_trained(string stmfile) {
    string line;
    short xx;
    string ifile(_ROOT_DATA);
    ifile.append(_LEARNED);
    ifile.append(stmfile);
    ifstream fin;
    fin.open(ifile.c_str());
    if (!fin.is_open()) {
      cerr << "Error in STM::load_trained, can not open file: " << ifile
	   << endl;
      exit(0);
    }
    for (short sti = 0; sti < 32; sti++) {
      getline(fin, line);
      stringstream ss(line);
      ss >> xx;
      if (xx != sti) {
	cerr << "Error in STM::load_trained:\t" << sti << "\t" << xx
	     << endl;
	exit(0);
      }
      ss >> center[sti];
      ss >> amplitude[sti];
    }
    fin.close();
    
  }
  
  float PM2::score_PM1(int ps) {
    return protein->stpVec.at(ps).compute_PAM_scores();
  }
  
  PM2::~PM2() {
    // TODO Auto-generated destructor stub
  }
  
} /* namespace std */
