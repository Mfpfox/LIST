/*
 * const.h
 *
 *  Created on: Apr 30, 2018
 *      Author: nmalhis
 */
#ifndef CONST_H_
#define CONST_H_

// #define	_WD		"" // "/home/nmalhis/WD/"
#define _ROOT_DATA	"" // "/home/nmalhis/LIST_Server_data/"
#define _LEARNED	"learned_files/"
#define _DB_PATH	"DB/DB" 
// #define _DB_PATH	"st01TREMBL/ST01T"

#define EDGE_SCORE_MIN_	3
#define GAP_EDG_SIZE_	2
#define	SEC_MATCH_MIN_	20
#define	SEC_IDEN_MIN_ 	0.35

#define	STP_RADIOS	4 // 5
#define	STP_MATCH_MIN	20
#define	STP_IDEN_MIN	0.5
#define	_HIGH_DEPTH	100000

#define _RESCALE_BASE	14

#define _PM2_IDEN_MIn	5 // 6
#define _PAM_IDEN_MIN	4 // 4 // STP::compute_PAM_scores()
#define _MIN_DEPTH	3

#define _L2SIZE	700
#define _DistributionSIZE 10001

#endif /* CONST_H_ */
