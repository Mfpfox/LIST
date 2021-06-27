/*
 * Score.cpp
 *
 *  Created on: May 5, 2018
 *      Author: nmalhis
 */

#include "Score.h"

namespace std {

float STP::pam_score[20];

Score::Score() {
	am.load_trained("AM.txt");
	pm2.load_trained("PM2.txt");
}

  void Score::_process_seq(string inSeq, string inSTP, string outFile) {
    float _pm2, _pm1, conservation;
    float _am[20], _pvm[20];
    pm2.protein = &protein;
    ofstream fout;
    clear();
    float LIST_Score[20];
    float msaScore14;
    float msaScore15;
    protein.clear();
    if (protein.load_seq(inSeq) > 0) {
      if (protein.load_stp(inSTP) > 0) {
	fout.open(outFile.c_str());
	{
	  // ------------------------------------------------------------------
	  fout << "AC\tPos\tRef\tDepth\tConservation";
	  fout << "\tA\tR\tN\tD\tC\tQ\tE\tG\tH\tI\tL\tK\tM\tF\tP\tS\tT";
	  fout << "\tW\tY\tV";
	  fout << endl;
	  for (unsigned ps = 1; ps < protein.stpVec.size(); ps++) {
	    int _cvr = protein.stpVec.at(ps).depth;
	    if (_cvr > _L2SIZE + 1)
	      _cvr = _L2SIZE + 1;
	    short refi = protein.stpVec.at(ps).ri;
	    if (ps > 1) {
	      if (refi >= 0 && refi < 20) {
		if (protein.stpVec.at(ps).depth >= _MIN_DEPTH) {
		  _pm1 =
		    10
		    * protein.stpVec.at(ps).compute_PAM_scores();
		  _pm2 = 10 * pm2.score_PM2(ps);
		} else {
		  _pm2 = 2.45883;
		  _pm1 = 8.52292;
		}
		// cerr << ps << "\t" << _pm1 << "\t" << _pm2;
		_pm2 = (_pm2 - minPM2[0]) / maxPM2[0]; // rst_count
		_pm2 = _pm2 * 0.6 + 0.2;
		_pm1 = (_pm1 - minPM1[0]) / maxPM1[0]; // rst_count
		_pm1 = _pm1 * 0.6 + 0.2;
		
		_pm2 = _pm2 - cover_Mr[0][_cvr];
		_pm2 = (_pm2 - minPM2[1]) / maxPM2[1];
		_pm2 = _pm2 * 0.6 + 0.2;
		
		_pm1 = _pm1 - cover_ST20[0][_cvr];
		_pm1 = _pm1 * cover_ST20[2][_cvr]
		  / (cover_ST20[1][_cvr]);
		
		_pm1 = (_pm1 - minPM1[1]) / maxPM1[1];
		_pm1 = _pm1 * 0.6 + 0.2;
		
		for (short ai = 0; ai < 20; ai++) {
		  _pvm[ai] = _am[ai] = -1;
		}
		for (short ai = 0; ai < 20; ai++) {
		  if (ai != refi) {
		    if (protein.stpVec.at(ps).depth
			>= _MIN_DEPTH) {
		      _pvm[ai] = STP::pam_score[ai] * 10;
		    } else {
		      _pvm[ai] = 3.22581;
		    }
		    // cerr << "\t" << _pvm[ai];
		    _pvm[ai] = (_pvm[ai] - minPVM[0])
		      / maxPVM[0];
		    _pvm[ai] = _pvm[ai] * 0.6 + 0.2;
		    
		    _pvm[ai] = _pvm[ai] - cover___ST[0][_cvr];
		    _pvm[ai] = _pvm[ai] * cover___ST[2][_cvr]
		      / (cover___ST[1][_cvr]);
		    
		    _pvm[ai] = (_pvm[ai] - minPVM[1])
		      / maxPVM[1];
		    _pvm[ai] = _pvm[ai] * 0.6 + 0.2;
		    
		    _am[ai] = am.score(refi, ai) * 10;
		    if (_am[ai] < 0) {
		      _am[ai] = 4.91544;
		    }
		    _am[ai] = (_am[ai] - minAM[0]) / maxAM[0];
		    _am[ai] = (_am[ai] * 0.6 + 0.2) * 0.6 + 0.4;
		  } else {
		    // cerr << "\t0";
		  }
		}
		// cerr << endl;
		// cerr << ps << "\t" << _pm1 << "\t" << _pm2;
		for (short ai = 0; ai < 20; ai++) {
		  // cerr << "\t" << _pvm[ai];
		}
		// cerr << endl;
		// must add conservation and idp
		fout << protein.seqFile << "\t" << ps;
		fout << "\t" << _get_aa(refi);
		fout << "\t" << protein.stpVec.at(ps).depth;
		float msa14;
		float msa15;
		float msa16;
		msaScore14 = Bayes(_pm2, _pm1, 1, 1);
		msa14 = rescale_n(14, msaScore14);
		for (short ai = 0; ai < 20; ai++) {
		  if (ai != refi) {
		    msaScore15 = Bayes(_pvm[ai], msa14, 0.7, 1);
		    msa15 = rescale_n(15, msaScore15);
		    
		    LIST_Score[ai] = Bayes(_am[ai], msa15, 0.3,
					   1);
		    msa16 = rescale_n(16, LIST_Score[ai]);
		    LIST_Score[ai] = msa16;
		  } else {
		    LIST_Score[ai] = 0;
		  }
		}
		conservation = 0;
		for (short ai = 0; ai < 20; ai++) {
		  if (ai != refi) {
		    conservation += LIST_Score[ai];
		  }
		}
		fout << "\t" << (conservation / 19);
		for (short ai = 0; ai < 20; ai++) {
		  fout << "\t" << LIST_Score[ai];
		}
		fout << endl;
	      } else { // refi range
		fout << protein.seqFile << "\t" << ps;
		fout << "\t" << _get_aa(refi);
		fout << "\t" << protein.stpVec.at(ps).depth;
		for (short ai = 0; ai < 21; ai++) {
		  fout << "\t0"; // << LIST_Score[ai];
		}
		fout << endl;
	      }
	    } else { // first aa
	      fout << protein.seqFile << "\t" << ps;
	      fout << "\t" << _get_aa(refi);
	      fout << "\t" << protein.stpVec.at(ps).depth;
	      fout << "\t0.785";
	      for (short ai = 0; ai < 20; ai++) {
		if (ai != refi) {
		  fout << "\t0.785";
		} else {
		  fout << "\t0";
		}
	      }
	      fout << endl;
	    }
	  }
	  fout.close();
	}
      }
    }
    fout.close();
    
  }
  
void Score::load_rescale(int ix) {
	if (ix >= _RESCALE_BASE && ix < _RESCALE_BASE + 3) {
		ifstream fin;
		stringstream ss0;
		string line, tmps;
		string mx1(_ROOT_DATA);
		mx1.append(_LEARNED);
		mx1.append("rescale_");
		ss0 << ix;
		ss0 >> tmps;
		mx1.append(tmps);
		mx1.append(".txt");

		fin.open(mx1.c_str());
		for (int j = 0; j < 10001; j++) {
			getline(fin, line);
			stringstream ss(line);
			ss >> tmps;
			ss >> _map[ix - _RESCALE_BASE][j];
		}
		fin.close();
	} else {
		cerr << "Error in Score::load_rescale: " << ix << endl;
	}
}

void Score::load__norm(string mnmx_1_file, string mnmx_2_file,
		string depth_file) {
	load_rescale(14);
	load_rescale(15);
	load_rescale(16);
	string mx1(_ROOT_DATA);
	mx1.append(_LEARNED);
	string mx2(mx1);
	string dp(mx1);
	mx1.append(mnmx_1_file);
	mx2.append(mnmx_2_file);
	dp.append(depth_file);
	ifstream fin;
	string line;
	string tmps;
	fin.open(mx1.c_str());
	getline(fin, line);
	while (!fin.eof()) {
		if (line.find("\tAM") != string::npos) {
			stringstream ss(line);
			ss >> tmps;
			ss >> tmps;
			ss >> minAM[0];
			ss >> maxAM[0];
			maxAM[0] -= minAM[0];
			// cerr << minAM[0] << "\t" << maxAM[0] << endl;
		} else if (line.find("PM2") != string::npos) {
			stringstream ss(line);
			ss >> tmps;
			ss >> tmps;
			ss >> minPM2[0];
			ss >> maxPM2[0];
			maxPM2[0] -= minPM2[0];
		} else if (line.find("PM1") != string::npos) {
			stringstream ss(line);
			ss >> tmps;
			ss >> tmps;
			ss >> minPM1[0];
			ss >> maxPM1[0];
			maxPM1[0] -= minPM1[0];
		} else if (line.find("PAM") != string::npos) {
			stringstream ss(line);
			ss >> tmps;
			ss >> tmps;
			ss >> minPVM[0];
			ss >> maxPVM[0];
			maxPVM[0] -= minPVM[0];
		}
		getline(fin, line);
	}
	fin.close();

	fin.open(mx2.c_str());
	getline(fin, line);
	while (!fin.eof()) {
		if (line.find("\tAM") != string::npos) {
			stringstream ss(line);
			ss >> tmps;
			ss >> tmps;
			ss >> minAM[1];
			ss >> maxAM[1];
			maxAM[1] -= minAM[1];
			// cerr << minAM[1] << "\t" << maxAM[1] << endl;
		} else if (line.find("PM2") != string::npos) {
			stringstream ss(line);
			ss >> tmps;
			ss >> tmps;
			ss >> minPM2[1];
			ss >> maxPM2[1];
			maxPM2[1] -= minPM2[1];
		} else if (line.find("PM1") != string::npos) {
			stringstream ss(line);
			ss >> tmps;
			ss >> tmps;
			ss >> minPM1[1];
			ss >> maxPM1[1];
			maxPM1[1] -= minPM1[1];
		} else if (line.find("PAM") != string::npos) {
			stringstream ss(line);
			ss >> tmps;
			ss >> tmps;
			ss >> minPVM[1];
			ss >> maxPVM[1];
			maxPVM[1] -= minPVM[1];
		}
		getline(fin, line);
	}
	fin.close();
	short cv;
	short dt = 0;
	fin.open(dp.c_str());
	getline(fin, line);
	while (!fin.eof()) {
		if (line.size() < 7)
			dt = 0;
		if (line.find("PM1") != string::npos) {
			dt = 1;
		} else if (line.find("PAM") != string::npos) {
			dt = 2;
		} else if (line.find("PM2") != string::npos) {
			dt = 3;
		} else {
			if (dt > 0) {
				stringstream ss(line);
				ss >> cv;
				if (dt == 3) {
					ss >> cover_Mr[0][cv];
					ss >> cover_Mr[1][cv];
					ss >> cover_Mr[2][cv];
					ss >> cover_Mr[3][cv];
				} else if (dt == 2) {
					ss >> cover___ST[0][cv];
					ss >> cover___ST[1][cv];
					ss >> cover___ST[2][cv];
					ss >> cover___ST[3][cv];
				} else if (dt == 1) {
					ss >> cover_ST20[0][cv];
					ss >> cover_ST20[1][cv];
					ss >> cover_ST20[2][cv];
					ss >> cover_ST20[3][cv];
				}
				if (cv == 701)
					dt = 0;
			}
		}
		getline(fin, line);
	}
	fin.close();
}

void Score::load_map(string map_f) {
	ifstream fin;
	string line;
	string tmps;
	string map_file(_ROOT_DATA);
	float value;
	map.clear();
	map_file.append(_LEARNED);
	map_file.append(map_f);
	// cerr << map_file << endl;
	fin.open(map_file.c_str());
	getline(fin, line);
	while (!fin.eof()) {
		if (line.size() > 2) {
			stringstream ss(line);
			ss >> tmps;
			ss >> value;
			map.push_back(value);
		}
		getline(fin, line);
	}
	fin.close();
}

float Score::map_it(float fv) {
	int iv = int(fv * (_DistributionSIZE - 1) + 0.49);
	if (iv < 0)
		iv = 0;
	if (iv >= _DistributionSIZE)
		iv = _DistributionSIZE - 1;
	return map.at(iv);
}

void Score::clear() {
	protein.clear();
}

float Score::_Bayes(float wf1, float wf2) {
	if (wf1 == 0 || wf2 == 0)
		return 0;
	return ((wf1 * wf2) / ((wf1 * wf2) + ((1 - wf1) * (1 - wf2))));
}

float Score::Bayes(float f1, float f2, float weight1, float weight2) {
	float wf1 = f1 * weight1 + (1 - weight1) / 2;
	float wf2 = f2 * weight2 + (1 - weight2) / 2;
	return _Bayes(wf1, wf2);
}

float Score::rescale_n(int idx, float fv) {
	if (idx >= _RESCALE_BASE && idx < _RESCALE_BASE + 3) {
		int iv = int(fv * 10000 + 0.49);
		if (iv < 0)
			iv = 0;
		if (iv > 10000)
			iv = 10000;
		return _map[idx - _RESCALE_BASE][iv];
	} else {
		cerr << "Error in Score::rescale_n:\t" << idx << endl;
		exit(0);
	}
}

Score::~Score() {
	// TODO Auto-generated destructor stub
}

} /* namespace std */

