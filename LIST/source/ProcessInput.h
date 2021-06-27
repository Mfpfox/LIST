/*
 * ProcessInput.h
 *
 *  Created on: Oct 25, 2018
 *      Author: nmalhis
 */

#ifndef PROCESSINPUT_H_
#define PROCESSINPUT_H_

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
#include "SMSA.h"
#include "STP.h"
#include "TaxonomyTree.h"

namespace std {

class ProcessInput {
  short _gap_edg_size;
  short _edge_score_min;
  
  short _sec_match_min;
  float _sec_iden_min;
  
  string querySectionSeq;
  int querySection_Start;
  int querySection_End; // after end
  
  string matchSectionSeq;
  
  string subjectSectionSeq;
  int sbjctSection_Start;
  int sbjctSection_End; // after end
  int subject_Size;
  
  vector<SMSA> smsaVec;
  vector<STP> stpVec;
  TaxonomyTree TTree;
  
  string _smsa_file;
  string _stp_file;
  ofstream _3_smsaout;
  ifstream _fin0;
 public:
  string thrds0;
  string thrds;
  SMSA _query;
  SMSA _subject;
  void align(string inSeq, string msaFile);
  void msa_to_smsa(string inSeq, string msaID, int dbCount);
  void smsa_to_stp(string inFile, string outFile);
  void smsa_to_stp_O(string inFile, string outFile);
  void _process_alignment();
  void _load_querySeq(string inSeq);
  void _new_subject(string &line);
  void _new_section(string &line);
  void _in_section(string &line);
  void _process_section();
  void _load_smsaVec();
  void _compu_mmm(SMSA &qSeq, SMSA &sSeq);
  void _process_smsaVec();
  void _process_smsaVec_O();
  int _get_LI(int loci, int pidx);

  void clean();
  ProcessInput();
  virtual ~ProcessInput();
  
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
};
 
} /* namespace std */

#endif /* PROCESSINPUT_H_ */
