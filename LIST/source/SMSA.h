/*
 * SMSA.h
 *
 *  Created on: Oct 25, 2018
 *      Author: nmalhis
 */

#ifndef SMSA_H_
#define SMSA_H_

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

namespace std {

class SMSA {
public:
   string head;
   string seq;
   string ac;
   int first;
   int last;
   int sz;
   int SI;
   int taxa_id;
   int mm; // mismatches
   short ST;
   bool valid;
   int matches;
   float identity;

   SMSA();
   SMSA(const SMSA &smsa);
   SMSA(const string &line);
   SMSA(const string &line, const short &typ);

   void set(const SMSA &smsa);
   void set(const string &line);
   void clear();
	virtual ~SMSA();
};

} /* namespace std */

#endif /* SMSA_H_ */
