/*
 * Taxa.cpp
 *
 *  Created on: Apr 13, 2018
 *      Author: nmalhis
 */

#include "TaxonomyTree.h"

namespace std {
  
  TaxonomyTree::TaxonomyTree() {
    string treeFile(_ROOT_DATA);
    treeFile.append("Taxa/tID_ptID.txt");
    ifstream fin_tree;
    St_II ii;
    string line;
    int tx;
    int pr;
    fin_tree.open(treeFile.c_str());
    if (!fin_tree.is_open()) {
      cerr << "Error in TaxonomyTree::TaxonomyTree(), can't open " << treeFile << endl;
      exit(1);
    }
    getline(fin_tree, line);
    while (!fin_tree.eof()) {
      if (line.size() > 2) {
	stringstream ss(line);
	ss >> tx;
	ss >> pr;
	if (tx == 0 || pr == 0 || tx == pr) {
	  // cerr << line << endl;
	}
	if (tx < 3000000) {
	  while ((int) taxaVec.size() < tx + 1)
	    taxaVec.push_back(0);
	  taxaVec.at(tx) = pr;
	} else {
	  cerr << "Bad line, Taxa is large:\t" << line << endl;
	}
      }
      getline(fin_tree, line);
    }
    fin_tree.close();
    // cerr << "taxaVec.size():\t" << taxaVec.size() << endl;
  }
  
  TaxonomyTree::~TaxonomyTree() {
    // TODO Auto-generated destructor stub
  }
  
  void TaxonomyTree::fillLineage(int targetTaxa) {
    St_II ii;
    ii.ikey = targetTaxa;
    ii.i2 = 31;
    lineageSet.insert(ii);
    while (ii.ikey > 1) {
      ii.ikey = taxaVec.at(ii.ikey);
      ii.i2--;
      if (ii.i2 > 0) {
	lineageSet.insert(ii);
      }
    }
  }
  
  int TaxonomyTree::get_shared_taxa(int taxa_num) {
    if (taxa_num >= (int) taxaVec.size() || taxa_num <= 0) {
      cerr << "Error in Taxa::get_shared_taxa: " << taxa_num << endl;
      return 0;
    }
    int ikeylast;
    St_II ii;
    set<St_II>::iterator it;
    ii.ikey = taxa_num;
    while ((it = lineageSet.find(ii)) == lineageSet.end()) {
      ikeylast = ii.ikey;
      ii.ikey = taxaVec.at(ii.ikey);
      if (ii.ikey == 1 || ikeylast == ii.ikey)
	return 0;
    }
    if (it != lineageSet.end()) {
      return it->i2;
    }
    return 0;
  }
  
  int TaxonomyTree::get_parent(int taxa_num) {
    if (taxa_num >= (int) taxaVec.size() || taxa_num <= 0) {
      cerr << "Error in Taxa::get_shared_taxa: " << taxa_num << endl;
      return 0;
    }
    return taxaVec.at(taxa_num);
  }
  
} /* namespace std */
