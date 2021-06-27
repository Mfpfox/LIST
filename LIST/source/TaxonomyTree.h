/*
 * Taxa.h
 *
 *  Created on: Apr 13, 2018
 *      Author: nmalhis
 */

#ifndef TAXONOMYTREE_H_
#define TAXONOMYTREE_H_

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

#include "St_II.h"
#include "const.h"

namespace std {

class TaxonomyTree {
public:
	vector<int> taxaVec;
	set<St_II> lineageSet;

	TaxonomyTree();
	virtual ~TaxonomyTree();

	void fillLineage(int targetTaxa);
	int get_shared_taxa(int taxa_num);
	int get_parent(int taxa_num);

};

} /* namespace std */

#endif /* TAXONOMYTREE_H_ */
