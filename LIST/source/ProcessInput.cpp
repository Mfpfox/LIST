/*
 * ProcessInput.cpp
 *
 *  Created on: Oct 25, 2018
 *      Author: nmalhis
 */

#include "ProcessInput.h"

namespace std {
  
  void ProcessInput::align(string inSeq, string msaFile) {
    int ret;
    string line;
    string cmd("blastp -outfmt 4 -evalue 0.01 -gapopen 11 -gapextend 2 ");
    cmd.append("-num_descriptions 150000  -num_alignments 150000 -num_threads ");
    cmd.append(thrds);
    cmd.append(" -query ");
    cmd.append(inSeq);
    cmd.append(" -db ");
    cmd.append(_DB_PATH);
    cmd.append(" -out ");
    cmd.append(msaFile);
    // cerr << endl << cmd << endl;
    ret = system(cmd.c_str());
  }

  void ProcessInput::smsa_to_stp(string inFile, string outFile) {
    _smsa_file.assign(inFile);
    _stp_file.assign(outFile);
    _load_smsaVec();
    _process_smsaVec();
    return;
  }
  
  void ProcessInput::smsa_to_stp_O(string inFile, string outFile) {
    _smsa_file.assign(inFile);
    _stp_file.assign(outFile);
    _load_smsaVec();
    // cerr << "_query\t" << _query.ac << "\t" << _query.sz << "\t" << _query.seq << endl;
    _process_smsaVec_O();
    return;
  }
  
  void ProcessInput::_load_smsaVec() {
    string smsa_line;
    smsaVec.clear();
    ifstream fin;
    fin.open(_smsa_file.c_str());
    if (!fin.is_open()) {
      cerr << "Error in Alignment::_2_smsa_to_stp, can not open "
	   << _smsa_file << endl;
      return;
    }
    getline(fin, smsa_line);
    stringstream ss(smsa_line);
    ss >> _query.ac;
    _query.ac.assign(_query.ac.substr(1));
    ss >> _query.sz;
    getline(fin, smsa_line);
    _query.seq.assign(smsa_line);
    getline(fin, smsa_line);
    while (!fin.eof()) {
      _subject.set(smsa_line);
      _compu_mmm(_query, _subject);
      // filter sequences longer than 20000 to identity greater than 0.45
      if (_subject.seq.size() < 20000 || _subject.identity > 0.45) {
	smsaVec.push_back(_subject);
      }
      getline(fin, smsa_line);
    }
    fin.close();
  }

  void ProcessInput::_compu_mmm(SMSA &qSeq, SMSA &sSeq) {
    sSeq.matches = 0;
    sSeq.mm = 0;
    sSeq.sz = 0;
    if (qSeq.seq.size() != sSeq.seq.size()) {
      cerr << "Error in ProcessInput::_compu_mmm(SMSA &qSeq, SMSA &sSeq)" << endl;
      exit(1);
    }
    for (unsigned i = 0; i < sSeq.seq.size(); i++) {
      if (toupper(sSeq.seq.at(i)) == qSeq.seq.at(i)) {
	sSeq.matches++;
      } else if (isdigit(sSeq.seq.at(i))) {
	sSeq.mm++;
      }
    }
    sSeq.sz = sSeq.matches + sSeq.mm;
    sSeq.identity = float(sSeq.matches) / sSeq.sz;
    sSeq.SI = short(float(100 * sSeq.matches) / sSeq.sz);
  }
  
  void ProcessInput::_process_smsaVec() {
    int liCount[20][32][STP_RADIOS * 2 + 2];
    stpVec.clear();
    /*
    if (_query.taxa_id == 0) {
      cerr << "ERROR in ProcessInput::_process_smsaVec()\t_query.taxa_id == 0"
	   << endl;
      exit(0);
    }
    // TTree.fillLineage(_query.taxa_id); // !!!!!!!!!!!!!!!!!!
    // _query.ST = _query.taxa_id;
    for (unsigned i = 0; i < smsaVec.size(); i++) {
      smsaVec.at(i).ST = smsaVec.at(i).taxa_id;
    }
    //*/
    for (unsigned pos = 0; pos < _query.seq.size(); pos++) { // =============================== POS
      for (short iAi = 0; iAi < 20; iAi++) {
	for (short iSt = 0; iSt < 32; iSt++) {
	  for (short iLi = 0; iLi < STP_RADIOS * 2 + 2; iLi++) {
	    liCount[iAi][iSt][iLi] = 0;
	  }
	}
      }
      STP stp;
      short refi = _get_aa_idx(_query.seq.at(pos));
      if (refi == 20) {
	stpVec.push_back(stp);
	continue;
      }
      for (unsigned i = 0; i < smsaVec.size(); i++) { // depth
	short ai = _get_aa_idx(smsaVec.at(i).seq.at(pos));
	short LI = _get_LI(pos, i);
	short sti = smsaVec.at(i).ST;
	if (ai < 20 && ai >= 0) {
	  if (sti > 0)
	    liCount[ai][sti][LI]++;
	  // start filter
	  if (smsaVec.at(i).identity >= STP_IDEN_MIN
	      && (smsaVec.at(i).matches >= STP_MATCH_MIN
		  || (smsaVec.size() <= 25
		      && smsaVec.at(i).matches
		      >= (smsaVec.size() * 75) / 100))) {
	    if (LI > STP_RADIOS) // ------(5)--- new added (+1)
	      stp.depth++;
	    
	    if (LI >= STP_RADIOS) { // + _4_rpmo ------------ STP
	      // stp.a_count[ai]++;
	      if (refi != ai) {
		stp.stp_count[sti]++;
		if (stp.stp_identity[sti] < LI) {
		  stp.stp_identity[sti] = LI;
		}
	      }
	    }
	    // ----------------------------------------------- PVM
	    if (LI > stp.a_lIdentity[ai]
		|| (LI == stp.a_lIdentity[ai]
		    && smsaVec.at(i).SI > stp.a_sIdentity[ai])
		|| (LI == stp.a_lIdentity[ai]
		    && smsaVec.at(i).SI == stp.a_sIdentity[ai]
		    && smsaVec.at(i).ST > stp.a_staxa[ai])) {
	      // --------------------------------------------------------
	      stp.a_lIdentity[ai] = LI;
	      stp.a_sIdentity[ai] = smsaVec.at(i).SI;
	      stp.a_staxa[ai] = smsaVec.at(i).ST;
	      // --------------------------------------------------------
	    }
	  }
	}
      } // for i=0
      
      for (short iAi = 0; iAi < 20; iAi++) {
	for (short iLi = 0; iLi < STP_RADIOS * 2 + 2; iLi++) {
	  for (short iSt = 30; iSt >= 0; iSt--) {
	    liCount[iAi][iSt][iLi] += liCount[iAi][iSt + 1][iLi];
	  }
	}
      }
      
      for (short iAi = 0; iAi < 20; iAi++) {
	short iLi = stp.a_lIdentity[iAi];
	short iST = stp.a_staxa[iAi];
	
	if (iLi > 2 && iST > 0) {
	  stp.a_count[iAi] = liCount[iAi][iST][iLi];
	  if (stp.depth > _HIGH_DEPTH) {
	    short x = iST;
	    while (liCount[iAi][x][iLi] < 2 && x > 0) {
	      x--;
	      stp.a_count[iAi] = liCount[iAi][x][iLi];
	      stp.a_staxa[iAi] = x;
	    }
	  }
	} else {
	  stp.a_count[iAi] = 0;
	}
      }
      
      short ai_max = 0;
      for (short ai = 0; ai < 20; ai++) {
	if (stp.a_lIdentity[ai] > stp.a_lIdentity[ai_max]
	    || (stp.a_lIdentity[ai] == stp.a_lIdentity[ai_max]
		&& stp.a_staxa[ai] > stp.a_staxa[ai_max])) {
	  ai_max = ai;
	}
	if (stp.a_lIdentity[ai_max] > 0) {
	  stp.ai_max = ai_max;
	} else {
	  stp.ai_max = -1;
	}
      }
      stpVec.push_back(stp);
    }
    
    // ------------------------------------------------------------ Output
    ofstream fout;
    fout.open(_stp_file.c_str());
    for (unsigned pos = 0; pos < stpVec.size(); pos++) {
      stpVec.at(pos).loci = pos + 1;
      stpVec.at(pos).raa = _query.seq.at(pos);
      stpVec.at(pos).ri = _get_aa_idx(_query.seq.at(pos));
      fout << stpVec.at(pos).to_string() << endl;
    }
    fout.close();
  }
  
  void ProcessInput::_process_smsaVec_O() {
    stpVec.clear();
    // cerr << "_query.seq.size():\t" << _query.seq.size() << endl;
    for (unsigned pos = 0; pos < _query.seq.size(); pos++) { // =============================== POS
      STP stp;
      short refi = _get_aa_idx(_query.seq.at(pos));
      if (refi == 20) {
	stpVec.push_back(stp);
	continue;
      }
      for (unsigned i = 0; i < smsaVec.size(); i++) { // depth
	short ai = _get_aa_idx(smsaVec.at(i).seq.at(pos));
	short LI = _get_LI(pos, i);
	short sti = smsaVec.at(i).ST;
	if (ai < 20 && ai >= 0) {
	  /*
	  cerr << pos << "\t" << i << "\t" << sti << "\t" << LI << endl;
	  if (i % 10 == 9)
	    sleep(1);
	  //*/
	  if (sti > 0) {
	    stp.a_count[ai]++;
	  }
	  // start filter
	  if (smsaVec.at(i).identity >= STP_IDEN_MIN
	      && (smsaVec.at(i).matches >= STP_MATCH_MIN
		  || (smsaVec.size() <= 25
		      && smsaVec.at(i).matches
		      >= (smsaVec.size() * 75) / 100))) {
	    if (LI > STP_RADIOS) // ------(5)--- new added (+1)
	      stp.depth++;
	    
	    // if (LI >= STP_RADIOS) { // + _4_rpmo ------------ STP
	      if (refi != ai) {
		stp.stp_count[sti]++;
		if (stp.stp_identity[sti] < LI) {
		  stp.stp_identity[sti] = LI;
		}
	      }
	      // }
	    // ----------------------------------------------- PVM
	    if (LI > stp.a_lIdentity[ai]
		|| (LI == stp.a_lIdentity[ai]
		    && smsaVec.at(i).SI > stp.a_sIdentity[ai])
		|| (LI == stp.a_lIdentity[ai]
		    && smsaVec.at(i).SI == stp.a_sIdentity[ai]
		    && sti > stp.a_staxa[ai])) {
	      // --------------------------------------------------------
	      stp.a_lIdentity[ai] = LI;
	      stp.a_sIdentity[ai] = smsaVec.at(i).SI;
	      stp.a_staxa[ai] = sti;
	      // --------------------------------------------------------
	    }
	  }
	}
      } // for i=0
      
      short ai_max = 0;
      for (short ai = 0; ai < 20; ai++) {
	if (stp.a_lIdentity[ai] > stp.a_lIdentity[ai_max]
	    || (stp.a_lIdentity[ai] == stp.a_lIdentity[ai_max]
		&& stp.a_staxa[ai] > stp.a_staxa[ai_max])) {
	  ai_max = ai;
	}
	if (stp.a_lIdentity[ai_max] > 0) {
	  stp.ai_max = ai_max;
	} else {
	  stp.ai_max = -1;
	}
      }
      stpVec.push_back(stp);
    }
    
    // ------------------------------------------------------------ Output
    ofstream fout;
    fout.open(_stp_file.c_str());
    for (unsigned pos = 0; pos < stpVec.size(); pos++) {
      stpVec.at(pos).loci = pos + 1;
      stpVec.at(pos).raa = _query.seq.at(pos);
      stpVec.at(pos).ri = _get_aa_idx(_query.seq.at(pos));
      fout << stpVec.at(pos).to_string() << endl;
    }
    fout.close();
  }
  
  int ProcessInput::_get_LI(int loci, int pidx) {
    short ret = 0;
    int sst;
    int sed;
    if (_query.seq.size() <= (2 * STP_RADIOS + 1)) {
      sst = 0;
      sed = _query.seq.size();
    } else {
      sst = loci - STP_RADIOS;
      if (sst < 0) {
	sst = 0;
      }
      sed = loci + STP_RADIOS + 1; // OK
      if (sed > smsaVec.at(pidx).seq.size()) {
	sed = smsaVec.at(pidx).seq.size();
      }
    }
    for (int i = sst; i < sed; i++) {
      if (i != loci
	  && _query.seq.at(i) == toupper(smsaVec.at(pidx).seq.at(i))) {
	ret++;
      }
    }
    if (loci == _query.first || loci == _query.last) {
      ret += ret / 2;
    }
    return ret;
  }

  void ProcessInput::clean() {
    smsaVec.clear();
    stpVec.clear();
    _edge_score_min = EDGE_SCORE_MIN_;
    _gap_edg_size = GAP_EDG_SIZE_;
    _sec_match_min = SEC_MATCH_MIN_;
    _sec_iden_min = SEC_IDEN_MIN_;
    querySectionSeq.assign("");
    matchSectionSeq.assign("");
    subjectSectionSeq.assign("");
    _smsa_file.assign("");
    _stp_file.assign("");
    querySection_Start = 0;
    querySection_End = 0;
    sbjctSection_Start = 0;
    sbjctSection_End = 0;
    subject_Size = 0;
  }
  
  void ProcessInput::msa_to_smsa(string inSeq, string msaID, int dbCount) {
    _load_querySeq(inSeq);
    _query.taxa_id = 9606;
    
    _smsa_file.assign(msaID);
    _smsa_file.append(".smsa");
    _3_smsaout.open(_smsa_file.c_str());
    _3_smsaout << _query.head << '\t' << _query.taxa_id << '\t'
	       << _query.seq.size();
    _3_smsaout << "\t0\t" << _query.seq << endl;
    
    for (int i = 0; i < dbCount; i++) {
      // cerr << ".";
      string str_part;
      stringstream ss;
      ss << i;
      ss >> str_part;
      string inFile_msa(msaID);
      inFile_msa.append("_");
      inFile_msa.append(str_part);
      inFile_msa.append(".msa");
      _fin0.open(inFile_msa.c_str());
      _process_alignment();
      _fin0.close();
    }
    _3_smsaout.close();
    return;
  }
  
  void ProcessInput::_process_alignment() {
    string line;
    int qSize;
    getline(_fin0, line);
    while (!_fin0.eof()) {
      if (line.find("Sequences producing significant alignments:")
	  < line.size()) {
	getline(_fin0, line); // blank
	break;
      } else if (line.find("Query=") < line.size()) {
	while (!(line.find("Length=") < line.size()))
	  getline(_fin0, line);
	stringstream ss;
	ss << line.substr(7);
	ss >> qSize;
      }
      getline(_fin0, line);
    }
    getline(_fin0, line);
    while (!_fin0.eof() && line.size() > 5) {
      getline(_fin0, line); // Skip the initial list
    }
    getline(_fin0, line);
    int blnk_count = 0;
    querySection_Start = -1;
    _subject.ac.assign("");
    while (!_fin0.eof()) {
      if (line.size() > 0) {
	//*
	if (line.at(0) == '>') {
	  if (blnk_count == 2) {
	    _process_section();
	  }
	  blnk_count = 0;
	  if (_subject.ac.size() > 3) {
	    _subject.clear();
	  }
	  _new_subject(line);
	} else if (line.find("Identities") < line.size()) {
	  _new_section(line);
	} else if (blnk_count == 2) {
	  _process_section();
	} else if (line.find("Query") == 0) {
	  _in_section(line);
	}
	//*/
	blnk_count = 0;
      } else {
	blnk_count++;
      }
      getline(_fin0, line);
    }
  }
  
  void ProcessInput::_new_subject(string &line) {
    stringstream ss1(line.substr(1));
    ss1 >> _subject.ac;
    ss1 >> _subject.taxa_id;
    while (!(line.find("Length=") < line.size()))
      getline(_fin0, line);
    stringstream ss2;
    ss2 << line.substr(7);
    ss2 >> subject_Size;
  }
  
  void ProcessInput::_new_section(string& line) {
    querySection_Start = -1;
    querySection_End = -1;
    sbjctSection_Start = -1;
    sbjctSection_End = -1;
    querySectionSeq.assign("");
    subjectSectionSeq.assign("");
    matchSectionSeq.assign("");
  }

  void ProcessInput::_in_section(string& line) {
    string tmps;
    int matchi;
    stringstream ss0(line);
    ss0 >> tmps;
    if (querySection_Start < 0)
      ss0 >> querySection_Start;
    else
      ss0 >> tmps;
    ss0 >> tmps;
    matchi = line.find(tmps, 5);
    querySectionSeq.append(tmps);
    ss0 >> querySection_End;
    
    getline(_fin0, line);
    matchSectionSeq.append(line.substr(matchi));
    
    getline(_fin0, line);
    stringstream ss2(line);
    ss2 >> tmps;
    if (sbjctSection_Start < 0)
      ss2 >> sbjctSection_Start;
    else
      ss2 >> tmps;
    ss2 >> tmps;
    subjectSectionSeq.append(tmps);
    ss2 >> sbjctSection_End;
    sbjctSection_End++;
  }
  
  void ProcessInput::_process_section() {
    if (querySection_Start >= 0) {
      if ((querySectionSeq.size() != matchSectionSeq.size())
	  || (querySectionSeq.size() != subjectSectionSeq.size())) {
	cerr << "Size Error: subAC: " << _subject.ac << "\t"
	     << sbjctSection_Start;
	cerr << "\t" << sbjctSection_End << endl;
	cerr << querySectionSeq << "| query section size: "
	     << querySectionSeq.size() << endl;
	cerr << matchSectionSeq << "| match section size: "
	     << matchSectionSeq.size() << endl;
	cerr << subjectSectionSeq << "| subject section size: "
	     << subjectSectionSeq.size() << endl;
	exit(0);
      }
      int gapps_inQuery = 0;
      string smsa_seq(string(_query.seq.size(), '.'));
      string match_seq(string(_query.seq.size(), '.'));
      vector<short> mask(_query.seq.size(), 0);
      int local = 0;
      for (int i = querySection_Start - 1;
	   (i - gapps_inQuery) < (int) _query.seq.size(); i++) {
	if (local >= (int) subjectSectionSeq.size()) {
	  break;
	}
	if (subjectSectionSeq.at(local) == '-') {
	  mask.at(i - gapps_inQuery) = 1;
	} else if (querySectionSeq.at(local) == '-') {
	  mask.at(i - gapps_inQuery) = 2;
	  gapps_inQuery++;
	} else {
	  smsa_seq.at(i - gapps_inQuery) = subjectSectionSeq.at(local);
	  match_seq.at(i - gapps_inQuery) = matchSectionSeq.at(local);
	}
	local++;
      }
      
      // ------------------------------------ _gap_edg_size (tolower)
      int lc_count = 0;
      for (unsigned k = 0; k < mask.size(); k++) {
	if (mask.at(k) > 0)
	  lc_count = _gap_edg_size;
	if (lc_count > 0 && isalpha(smsa_seq.at(k))) {
	  smsa_seq.at(k) = tolower(smsa_seq.at(k));
	  lc_count--;
	}
      }
      lc_count = 0;
      for (int k = mask.size() - 1; k >= 0; k--) {
	if (lc_count > 0 && isalpha(smsa_seq.at(k))) {
	  smsa_seq.at(k) = tolower(smsa_seq.at(k));
	  lc_count--;
	}
	if (mask.at(k) > 0)
	  lc_count = _gap_edg_size;
      }
      // ----------------------------------------- _edge_score_min
      if (querySection_End > (int) mask.size()) {
	cerr << "Fatal ERROD, ";
	cerr << "querySection_End:\t" << querySection_End << " : "
	     << mask.size() << endl;
	exit(1);
      }
      int pass = 0;
      int edgeStart = querySection_Start - 1;
      int edgeEnd = querySection_End;
      int score = 0;
      for (int k = querySection_Start - 1; k < querySection_End - 5; k++) {
	if (match_seq.at(k) == '+') {
	  score += 1;
	} else if (match_seq.at(k) == toupper(smsa_seq.at(k))) {
	  if (_query.seq.at(k) == smsa_seq.at(k)) {
	    score += 2;
	  } else {
	    score += 1;
	  }
	} else {
	  if (mask.at(k) > 0) {
	    score = -1;
	  } else {
	    score = 0;
	  }
	  edgeStart = k + 1;
	}
	if (score >= _edge_score_min) {
	  pass++;
	  break;
	}
      }
      score = 0;
      for (int k = querySection_End - 1; k > edgeStart + 5; k--) {
	if (match_seq.at(k) == '+') {
	  score += 1;
	} else if (match_seq.at(k) == toupper(smsa_seq.at(k))) {
	  if (_query.seq.at(k) == smsa_seq.at(k)) {
	    score += 2;
	  } else {
	    score += 1;
	  }
	} else {
	  if (mask.at(k) > 0) {
	    score = -1;
	  } else {
	    score = 0;
	  }
	  edgeEnd = k - 1;
	}
	if (score >= _edge_score_min) {
	  pass++;
	  break;
	}
      }
      
      if (pass == 2) {
	int k = querySection_Start - 1;
	for (; k < edgeStart; k++) {
	  smsa_seq.at(k) = '.';
	}
	k = edgeEnd + 1;
	for (; k < querySection_End; k++) {
	  smsa_seq.at(k) = '.';
	}
      }
      // --------------------------------------------- Identity
      int matches = 0;
      int mm = 0;
      int dots = 0;
      for (unsigned i = 0; i < smsa_seq.size(); i++) {
	if (smsa_seq.at(i) != '.') {
	  if (_query.seq.at(i) == toupper(smsa_seq.at(i))) {
	    matches++;
	  } else {
	    mm++;
	  }
	} else {
	  dots++;
	}
      }
      float sec_iden = float(matches) / float(matches + mm);
      if (sec_iden > _sec_iden_min && matches > _sec_match_min) {
	_3_smsaout << _subject.ac << '\t' << _subject.taxa_id;
	_3_smsaout << '\t' << matches << "\t" << mm << "\t" << smsa_seq
		   << endl;
      }
    } // if querySection_Start >= 0
    querySection_Start = -1;
    querySectionSeq.assign("");
	matchSectionSeq.assign("");
	subjectSectionSeq.assign("");
	return;
  }
  
  void ProcessInput::_load_querySeq(string inSeq) {
    //*
    string inSeqFile(inSeq);
    // inSeqFile.append(".fasta");
    
    string line;
    _fin0.open(inSeqFile.c_str());
    getline(_fin0, line);
    stringstream ss(line);
    ss >> line;
    _query.head.assign(line);
    _query.seq.assign("");
    getline(_fin0, line);
    while (!_fin0.eof()) {
      if (line.size() > 0) {
	if (line.at(0) != '>') {
	  _query.seq.append(line);
	}
      }
      getline(_fin0, line);
    }
    _fin0.close();
    _query.sz = _query.seq.size();
  }
  
  ProcessInput::ProcessInput() {
    _edge_score_min = EDGE_SCORE_MIN_;
    _gap_edg_size = GAP_EDG_SIZE_;
    _sec_match_min = SEC_MATCH_MIN_;
    _sec_iden_min = SEC_IDEN_MIN_;    
  }
  
  ProcessInput::~ProcessInput() {
    // TODO Auto-generated destructor stub
  }
  
} /* namespace std */
