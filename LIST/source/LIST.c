//============================================================================
// Name        : LIST
// Author      : Nawar Malhis
// Version     : 2.0
//============================================================================

#include <iostream>
#include "const.h"
#include "Score.h"
#include "ProcessInput.h"
#include "Almn.h"

using namespace std;

void score_server_seq(string seq_file, string stp_file, string list_file);

int main(int argc, char *argv[]) {
  int pcount = argc;
  int thrd   = 4;
  cerr << endl;
  cerr << " LIST 0.1 beta (2019),";
  cerr << " Position conservation matrix for protein sequences" << endl;
  cerr << " By Nawar Malhis, the University 0f British Columbia Michael Smith Laboratories";
  cerr << endl;
  cerr << endl;
  vector <string> inV;
  vector <string> outV;
  if (pcount == 3) {
    string arg2(argv[2]);
    stringstream ss2(arg2);
    ss2 >> thrd;
    pcount = 2;
    if (thrd < 1) {	
      cerr << " Invalid thread counts " << argv[2] << endl;
      exit(1);
    }
  }
  if (pcount == 2) {
    string inPath(argv[1]);
    if (inPath.at(inPath.size() - 1) != '/')
      inPath.append("/");
    string rmd("rmdir tmp/");
    rmd.append(inPath);
    string cmd("mkdir tmp/");
    cmd.append(inPath);
    cmd.append(" 2> tmp/junk");
    int ret0 = system(cmd.c_str());
    cmd.assign("ls ");
    cmd.append(inPath);
    cmd.append(" > tmp/");
    cmd.append(inPath);
    cmd.append("__LIST_dir.tmp");
    int ret1 = system(cmd.c_str());
    string line, oname, iname;
    ifstream fin;
    string flist("tmp/");
    flist.append(inPath);
    flist.append("__LIST_dir.tmp");
    fin.open(flist.c_str());
    getline(fin, line);
    while(!fin.eof()) {
      if (line.find(".fasta") == line.size() - 6) {
	iname.assign(inPath);
	iname.append(line);
	oname.assign(inPath);
	oname.append(line.substr(0,line.size() - 6));
	oname.append(".list");
	inV.push_back(iname);
	outV.push_back(oname);
      }
      getline(fin, line);
    }
    fin.close();

    cmd.assign("rm ");
    cmd.append(flist);
    ret1 = system(cmd.c_str());
    string cln("rm ");
    cln.append("tmp/");
    cln.append(inPath);
    cln.append("*");
    ProcessInput procIn;
    stringstream ss;
    ss << thrd;
    ss >> procIn.thrds;
    procIn.thrds.append(" ");
    cerr << " Processing " << inV.size() << " fasta sequences in <" << inPath;
    cerr << "> using " << thrd << " threads" << endl << endl;
    vector <string> msaV;
    vector <string> smsaV;
    vector <string> stpV;
    for (unsigned i = 0; i < inV.size(); i++) {
      string id("tmp/");
      id.append(inV.at(i).substr(0,inV.at(i).size() - 6));
      string msaFile(id);
      msaFile.append(".msa");
      msaV.push_back(msaFile);
      string smsaFile(id);
      smsaFile.append(".smsa");
      smsaV.push_back(smsaFile);
      string stpFile(id);
      stpFile.append(".stp");
      stpV.push_back(stpFile);
    }

    Almn almn;
    for (unsigned i = 0; i < inV.size(); i++) {
      cerr << i+1 << "\t" << inV.at(i) << "\t" << outV.at(i) << endl;
      almn.clear();
      procIn.align(inV.at(i), msaV.at(i));
      almn.process_MSA(msaV.at(i), smsaV.at(i));
      procIn.clean();
      procIn.smsa_to_stp_O(smsaV.at(i), stpV.at(i));
      score_server_seq(inV.at(i), stpV.at(i), outV.at(i));
      ret1 = system(cln.c_str());
    }
    //*/
    ret1 = system(rmd.c_str());
  } else {
    cerr << "format: ./list 'data_path' [#seq  #threads]" << endl;
  }
  return 0;
}

void score_server_seq(string seq_file, string stp_file, string list_file) {
  Score s;
  s.load__norm("norm_1_mn_mx.txt", "norm_2_mn_mx.txt", "cover_tr.txt");
  s._process_seq(seq_file, stp_file, list_file);
}

