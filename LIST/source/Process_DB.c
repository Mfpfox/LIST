//============================================================================
// Name        : Process_DB.cpp
// Author      : Nawar Malhis
// Version     :
// Copyright   : Your copyright notice
//============================================================================
// /data2/P2/bak_msa/
// /olddisk/data/database/P2/bak_smsa/


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

using namespace std;

void process_uniprot(string spFile, string trFile, string dbPath);

int main(int argc, char* argv[]) {
  cerr << endl;
  cerr << " Processing UniProt databases for LIST 0.1 beta (2019)" << endl;
  cerr << " By Nawar Malhis, the University 0f British Columbia";
  cerr << " Michael Smith Laboratories" << endl;
  if (argc == 3) {
    cerr << "\nProcessing database" << endl;
    process_uniprot(argv[1], argv[2], "DB");
  } else {
    cerr << "./process_db UniProt/uniprot_sprot.fasta UniProt/uniprot_trembl.fasta " << endl;
  }
  return 0;
}

void process_uniprot(string spFile, string trFile, string dbPath) {
  string line, tmps;
  string mkdir("mkdir ");
  mkdir.append(dbPath);
  int ret1 = system(mkdir.c_str());

  
  vector<string> inDB;
  // vector<string> outDB;
  ofstream oVec;
  ofstream fout;
  inDB.push_back(spFile);
  inDB.push_back(trFile);
  string UP2(dbPath);
  UP2.append("/up2.fasta");
  oVec.open(UP2.c_str());
  bool skp;
  for (unsigned fl = 0; fl < inDB.size(); fl++) {
    skp = false;
    ifstream dbin;
    dbin.open(inDB.at(fl).c_str());
    if (!dbin.is_open()) {
      cerr << "Error in process_uniprot, can't open " << inDB.at(fl)
	   << endl;
      exit(1);
    }
    int cnt = 0;
    cerr << "Processing " << inDB.at(fl) << ": ";
    getline(dbin, line);
    while (!dbin.eof()) {
      if (line.size() > 0) {
	if (line.find("virus") != string::npos) {
	  skp = true;
	} else {
	  if (line.at(0) == '>') {
	    cnt++;
	    if (cnt % 100000 == 0)
	      cerr << '.';
	    skp = false;
	    stringstream ss(line);
	    ss >> tmps;
	    oVec << tmps;
	    int p1 = line.find("OX=");
	    int p2 = line.find(' ', p1 + 2);
	    int sz = p2 - p1;
	    oVec << " " << line.substr(p1 + 3, sz - 3) << endl;
	  } else {
	    if (!skp)
	      oVec << line << endl;
	  }
	}
      }
      getline(dbin, line);
    }
    dbin.close();
    cerr << endl;
  }
  oVec.close();
  
  
  cerr << "Formating database " << endl;
  string cmd("makeblastdb -in ");
  cmd.append(UP2);
  cmd.append(" -parse_seqids -dbtype prot -title DB -out ");
  cmd.append(dbPath);
  cmd.append("/DB");
  cerr << cmd << endl;
  int ret3 = system(cmd.c_str());
}
