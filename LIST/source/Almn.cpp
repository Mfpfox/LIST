/*
 * Almn.cpp
 *
 *  Created on: Aug 26, 2017
 *      Author: nmalhis
 */

#include "Almn.h"

namespace std {
  
  Almn::SI::~SI() {
  }
  
  Almn::Almn() {
    _length = 0;
    _qStart = 0;
    _query.assign("");
    _qsstart = 0;
    _qssize = 0;
    _qsz0 = 0;
    lines_used = 0;
  }
  
  Almn::~Almn() {
    // TODO Auto-generated destructor stub
  }
  
  void Almn::clear(){
    _length = 0;
    _qStart = 0;
    _query.assign("");
    _qsstart = 0;
    _qssize = 0;
    _qsz0 = 0;
    lines_used = 0;
    _siiv_vec.clear();
    _ac_set.clear();
  }

  void Almn::process_MSA(string in_file, string out_file) {
    short ret;
    ret = _process_MSA(in_file);
    if (ret < 0) {
    } else {
      _save_msaSeq(out_file);
    }
  }
  
  short Almn::_process_MSA(string in_file) {
    int count = 0;
    short ret = -10;
    set<SI>::iterator it;
    SI si;
    int ql_start;
    string q_line;
    int ql_end = 0;
    bool set_full = false;
    TTree.fillLineage(9606);
    _length = 0;
    _qStart = -1;
    lines_used = 0;
    _query.assign("");
    ifstream fin;
    fin.open(in_file.c_str());
    if (!fin.is_open()) {
      return -1;
    }
    _siiv_vec.clear();
    _ac_set.clear();
    string line;
    string tmps;
    getline(fin, line);
    while (!fin.eof()) {
      if (_length == 0) {
	if (line.find("Length=") == 0) {
	  stringstream ss(line.substr(7));
	  ss >> _length;
	  getline(fin, line);
	  getline(fin, line);
	  getline(fin, line);
	}
      } else {
	if (!set_full) {
	  if (line.size() > 5) {
	    if (line.find("Lambda") == 0) {
	      fin.close();
	      return -6;
	    }
	    int tx;
	    stringstream ss(line);
	    ss >> tmps;
	    ss >> tx;
	    for (unsigned i2 = 0; i2 < tmps.size(); i2++) {
	      if (tmps.at(i2) == '|')
		tmps.at(i2) = ' ';
	    }
	    stringstream ss2(tmps);
	    ss2 >> tmps;
	    ss2 >> si.st;
	    si.ii = _siiv_vec.size();
	    _ac_set.insert(si);
	    SIIv siiv;
	    _siiv_vec.push_back(siiv);
	    _siiv_vec.back().ac.assign(si.st);
	    _siiv_vec.back().sTaxa = TTree.get_shared_taxa(tx);
	  } else if (_ac_set.size() > 0) {
	    // cerr << "_ac_set.size():\t" << _ac_set.size() << endl;
	    set_full = true;
	    ret = -2; //
	  }
	} else { // ------------------------------------  set is full -------------
	  if (line.size() > 5) {
	    if (line.find("Query") == 0) {
	      _qsz0 = _query.size();
	      stringstream ss(line);
	      ss >> tmps;
	      ss >> ql_start;
	      if (_qStart < 0)
		_qStart = ql_start - 1;
	      if (!ql_start == ql_end + 1) {
		fin.close();
		return -3; // bad query sequence (based on number)
	      }
	      ss >> q_line;
	      ss >> ql_end;
	      _qsstart = line.find(q_line);
	      _qssize = q_line.size();
	      _query.append(q_line);
	      long memused = _query.size();
	      memused = memused * lines_used;
	      if (memused > 4000000000)
		return -7;
	      for (unsigned i = 0; i < _siiv_vec.size(); i++) {
		for (unsigned j = 0;
		     j < _siiv_vec.at(i).msa_e_vec.size(); j++) {
		  unsigned sz =
		    _siiv_vec.at(i).msa_vec.at(j).size();
		  if (sz < _qsz0) {
		    _siiv_vec.at(i).msa_vec.at(j).append(
							 string(_qsz0 - sz, '.'));
		  }
		}
	      }
	    } else if (line.find("Lambda") == 0) {
	      ret = 0; // loading is completed
	      break;
	    } else {
	      count += _process_line(line);
	    }
	  }
	}
      }
      getline(fin, line);
    }
    fin.close();
    return ret;
  }
  
  short Almn::_process_line(string& line) {
    short ret = 0;
    int sl_start = -1;
    int sl_end = -1;
    string sl_line;
    string tmps;
    set<SI>::iterator it;
    SI si;
    stringstream ss(line);
    ss >> si.st;
    ss >> sl_start;
    ss >> tmps;
    ss >> sl_end;
    if (sl_end < 1)
      sl_end = 1;
    it = _ac_set.find(si);
    if (it != _ac_set.end()) {
      int ii = it->ii;
      if (ii < (int) _siiv_vec.size()) {
	sl_line.assign(line.substr(_qsstart, _qssize));
	short countalpha = _alpha_count(sl_line);
	if (countalpha > 0) {
	  // ==============================================
	  int jj = -1;
	  for (unsigned j = 0; j < _siiv_vec.at(ii).msa_vec.size(); j++) {
	    if (_siiv_vec.at(ii).msa_e_vec.at(j) == sl_start - 1) {
	      if (_siiv_vec.at(ii).msa_vec.at(j).size()
		  + sl_line.size() == _query.size()) {
		jj = j;
		break;
	      }
	    }
	  }
	  if (jj == -1) {
	    jj = _siiv_vec.at(ii).msa_vec.size();
	    lines_used++;
	    _siiv_vec.at(ii).msa_vec.push_back(string(_qsz0, '.'));
	    _siiv_vec.at(ii).msa_e_vec.push_back(sl_end);
	    ret = 1;
	  } else {
	    _siiv_vec.at(ii).msa_e_vec.at(jj) = sl_end;
	  }
	  _siiv_vec.at(ii).msa_vec.at(jj).append(sl_line);
	} else {
	}
      } else {
	// cerr << "Alignments::_process_line: 3" << endl;
      }
    } else {
      // cerr << "Alignments::_process_line: 4" << endl;
      // cerr << line << endl;
    }
    return ret;
  }
  
  short Almn::_alpha_count(string st) {
    short ret = 0;
    for (unsigned i = 0; i < st.size(); i++) {
      if (isalpha(st.at(i)))
	ret++;
    }
    return ret;
  }

  short Almn::_save_msaSeq(string out_file) {
    short ret = 0;
    for (unsigned i = 0; i < _siiv_vec.size(); i++) {
      for (unsigned j = 0; j < _siiv_vec.at(i).msa_vec.size(); j++) {
	for (unsigned k = 0; k < _siiv_vec.at(i).msa_vec.at(j).size();
	     k++) {
	  if (_siiv_vec.at(i).msa_vec.at(j).at(k) == ' ')
	    _siiv_vec.at(i).msa_vec.at(j).at(k) = '.';
	}
	if (_siiv_vec.at(i).msa_vec.at(j).size() < _query.size()) {
	  _siiv_vec.at(i).msa_vec.at(j).append(
					       string(
						      _query.size()
						      - _siiv_vec.at(i).msa_vec.at(j).size(),
						      '.'));
	} else if (_siiv_vec.at(i).msa_vec.at(j).size() > _query.size()) {
	  // cerr << endl << _query << endl;
	  // cerr << _siiv_vec.at(i).ac << endl;
	  // cerr << _siiv_vec.at(i).msa_vec.at(j) << endl;
	  // cerr << endl;
	  // sleep(1);
	}
      }
    }
    for (unsigned i = 0; i < _query.size(); i++) {
      if (_query.at(i) == '-')
	_query.at(i) = '_';
    }
    if (_query.at(0) == '_' || _query.at(_query.size() - 1) == '_')
      return -14; // bad query sequence
    // -------------------------------------
    ofstream fout;
    fout.open(out_file.c_str());
    if (!fout.is_open())
      return -18;
    string spQuery(_length, '.');
    int spi = _qStart;
    for (unsigned i = 0; i < _query.size(); i++) {
      if (isalpha(_query.at(i))) {
	if (spi == _length) {
	  cerr << "Alignments::_process_MSA: 5.0" << endl;
	  exit(0);
	}
	spQuery.at(spi) = toupper(_query.at(i));
	spi++;
      }
    }
    if (spi > _length) {
      return -11;
    } else if ((spi - _qStart) < _length) {
      ret = _length - spi + _qStart;
    }
    int cut = out_file.find("_1/") + 3;
    cut = out_file.find_last_of("/\\");
    string tmp_ac(out_file.substr(cut+1));
    fout << ">" << tmp_ac.substr(0, tmp_ac.size() - 5) << "\t" << spQuery.size() << endl;
    fout << spQuery << endl;
    // ----------------------------------------------------------
    
    for (unsigned i = 0; i < _siiv_vec.size(); i++) {
      for (unsigned j = 0; j < _siiv_vec.at(i).msa_vec.size(); j++) {
	
	if (_siiv_vec.at(i).msa_vec.at(j).size() == _query.size()) {
	  float identity = 0;
	  string seq(_length, '.');
	  int si = _qStart;
	  unsigned faa = 1000000;
	  unsigned laa = 0;
	  short skip_f = 0;
	  for (unsigned k = 0; k < _query.size(); k++) {
	    if (isalpha(_siiv_vec.at(i).msa_vec.at(j).at(k))) {
	      if (faa == 1000000)
		faa = k;
	      laa = k;
	    }
	  }
	  // --------------------------------------------
	  for (unsigned k = 0; k < _query.size(); k++) {
	    // .....................................................
	    if (isalpha(_query.at(k))) {
	      if (isalpha(_siiv_vec.at(i).msa_vec.at(j).at(k))) {
		if (skip_f == 0) {
		  seq.at(si) = toupper(
				       _siiv_vec.at(i).msa_vec.at(j).at(k));
		} else {
		  seq.at(si) = tolower(
				       _siiv_vec.at(i).msa_vec.at(j).at(k));
		  skip_f--;
		}
	      } else {
		if (k > faa && k < laa) {
		  skip_f = 2;
		  int b2si = si;
		  b2si -= 2;
		  if (b2si < 0)
		    b2si = 0;
		  for (int bk = b2si; bk < si; bk++) {
		    identity--;
		    seq.at(bk) = tolower(seq.at(bk));
		  }
		}
	      }
	      si++;
	    } else {
	      if (isalpha(_siiv_vec.at(i).msa_vec.at(j).at(k))) {
		if (k > faa && k < laa) {
		  skip_f = 2;
		  int b2si = si;
		  b2si -= 2;
		  if (b2si < 0)
		    b2si = 0;
		  for (int bk = b2si; bk < si; bk++) {
		    identity--;
		    seq.at(bk) = tolower(seq.at(bk));
		  }
		}
	      } else {
	      }
	    }
	  }
	  fout << _siiv_vec.at(i).ac
	       << string(12 - _siiv_vec.at(i).ac.size(), ' ');
	  fout << _siiv_vec.at(i).sTaxa << "\t" << seq << endl;
	  // .....................................................
	  
	  // --------------------------------------------
	} else {
	  cerr << endl << "Error\t" << i << "\t" << j << "\t"
	       << _siiv_vec.at(i).msa_vec.at(j).size() << "\t"
	       << _query.size() << endl;
	  return -15;
	}
      }
    }
    //*/
    fout.close();
    
    return ret;
  }
  
} /* namespace std */
