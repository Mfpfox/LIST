/*
 * Protein.cpp
 *
 *  Created on: May 5, 2018
 *      Author: nmalhis
 */

#include "Protein.h"

namespace std {

  Protein::Protein() {
    // stp_path.assign(_ROOT_DATA); // "/data2/LIST_2_data/"
    /*
      stp_path.assign("/data2/LIST_2_data/");
      stp_path.append(_DATA);
      stp_path.append(_STP);//*/
  }
  
  int Protein::load_seq(string inSeq) {
    positionsVec.clear();
    stpVec.clear();
    // cerr << inSeq << endl;
    string seq_file(inSeq);
    // seq_file.append(inSeq);
    // seq_file.append(".fasta");
    ifstream seqin;
    seqin.open(seq_file.c_str());
    if (!seqin.is_open()) {
      cerr << "Can not open seq_file: " << seq_file << endl;
      return -1;
    }
    seq.assign(" ");
    string line;
    getline(seqin, line);
    stringstream ss(line);
    ss >> line;
    seqFile.assign(line.substr(1,15));
    getline(seqin, line);
    while (!seqin.eof()) {
      if (line.size() > 0) {
	seq.append(line);
      }
      getline(seqin, line);
    }
    seqin.close();
    // cerr << seq << endl;
    return seq.size() - 1;
  }
  
  int Protein::load_stp(string inSTP) {
    ifstream stpin;
    // string stp_file(_WD);
    // stp_file.append("input.stp");
    stpin.open(inSTP.c_str());
    if (!stpin.is_open()) {
      cerr << "can not open STP_File: " << inSTP << endl;
      return -1;
    } else {
      stpVec.clear();
      stpVec.push_back(STP());
      string line;
      getline(stpin, line);
      while (!stpin.eof()) {
	if (line.size() > 0) {
	  stpVec.push_back(STP(line));
	  if (stpVec.back().raa != 0
	      && seq.at(stpVec.size() - 1) != stpVec.back().raa) {
	  }
	  if (isalpha(stpVec.back().raa)
	      && stpVec.back().raa != seq.at(stpVec.back().loci)) {
	    cerr << "Error in Protein::load_stp, stp.raa not matched with seq " << seqFile
		 << "(" << stpVec.back().loci << ") " << stpVec.back().raa << " vs "
		 << seq.at(stpVec.back().loci) << endl;
	    return -1;
	  }
	}
	getline(stpin, line);
      }
      stpin.close();
      return stpVec.size() - 1;
    }
  }
  
  void Protein::clear() {
    seqFile.assign("");
    seq.assign("");
    
    positionsVec.clear();
    stpVec.clear();
  }
  
  Protein::~Protein() {
    // TODO Auto-generated destructor stub
  }
  
} /* namespace std */

