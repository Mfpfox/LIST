/*
 * Position.cpp
 *
 *  Created on: Feb 6, 2019
 *      Author: nmalhis
 */

#include "Position.h"

namespace std {

Position::Position() {
	loci = 0;
	refi = 0;
	conservation = 0;
}

Position::Position(const Position& pos) {
	loci = pos.loci;
	refi = pos.refi;
	conservation = pos.conservation;
}

Position::Position(const string& line) {
	loci = -1;
	refi = 0;
	conservation = 0;
	stringstream ss(line);
	ss >> loci;
	ss >> refi;
}

Position::~Position() {
	// TODO Auto-generated destructor stub
}

} /* namespace std */
