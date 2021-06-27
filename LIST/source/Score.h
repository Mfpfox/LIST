/*
 * Score.h
 *
 *  Created on: May 5, 2018
 *      Author: nmalhis
 */

#ifndef SCORE_H_
#define SCORE_H_

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

#include "AM.h"
#include "STP.h"
#include "Protein.h"
#include "PM2.h"

namespace std {
  
  class Score {
    float minPVM[2];
    float maxPVM[2];
    float minPM1[2];
    float maxPM1[2];
    float minPM2[2];
    float maxPM2[2];
    float minAM[2];
    float maxAM[2];
    float _map[3][10001];
    float cover_Mr[4][803];
    float cover_ST20[4][803];
    float cover___ST[4][803];
    vector <float> map;
  public:
    AM am;
    PM2 pm2;
    Protein protein;
    
    Score();
    void _process_seq(string inSeq, string inSTP, string outFile);
    float map_it(float fv);
    void load__norm(string mnmx_1_file, string mnmx_2_file, string depth_file);
    void load_map(string map_file);
    void load_rescale(int ix);
    void clear();
    float Bayes(float f1, float f2, float weight1, float weight2);
    float _Bayes(float wf1, float wf2);
    float rescale_n(int idx, float fv);
    
    virtual ~Score();
    
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

#endif /* SCORE_H_ */
