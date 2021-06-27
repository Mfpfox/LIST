/*
 * Almn.h
 *
 *  Created on: Aug 26, 2017
 *      Author: nmalhis
 */

#ifndef ALMN_H_
#define ALMN_H_

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
#include "TaxonomyTree.h"

namespace std {

  class Almn {
    
  public:
    struct SI {
    public:
      string st;
      int ii;
      SI(){
	ii = 0;
      }
      SI(const SI &si){
	st.assign(si.st);
	ii = si.ii;
      }
      
      bool operator<(const SI &si) const {
	return st.compare(si.st) < 0;
      }
      bool operator==(const SI &si) const {
	return st.compare(si.st) == 0;
      }
      virtual ~SI();
    };
    
    struct SIIv {
      string ac;
      short sTaxa;
      vector<string> msa_vec;
      vector<int> msa_e_vec;
    };
    
    int _length;
    int _qStart;
    string _query;
    int _qsstart;
    int _qssize;
    unsigned _qsz0;
    
    int lines_used;
    vector<SIIv> _siiv_vec;
    set<SI> _ac_set;
    
    TaxonomyTree TTree;
   
    short _process_MSA(string in_file);
    short _process_line(string &line);
    short _alpha_count(string st);
    short _save_msaSeq(string out_file);
    void process_MSA(string in_file, string out_file);
    void clear();
    
    Almn();
    virtual ~Almn();
  };
  
} /* namespace std */

#endif /* ALMN_H_ */
