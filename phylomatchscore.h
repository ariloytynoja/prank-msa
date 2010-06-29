/***************************************************************************
 *   Copyright (C) 2005 by Ari Loytynoja                                   *
 *   ari@sink                                                              *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef PHYLOMATCHSCORE_H
#define PHYLOMATCHSCORE_H

#include "sequence.h"
#include "ancestralsequence.h"
#include "terminalsequence.h"
#include "dbmatrix.h"

/**
 * Computes match scores for the current column.
 */

class PhyloMatchScore
{
	DbMatrix* fM;
	DbMatrix* bM;

	DbMatrix* flM;

	DbMatrix* idX;
	DbMatrix* idY;

	DbMatrix* match;
	DbMatrix* gap;

	Sequence* s1;
	Sequence* s2;
	AncestralSequence* a1;
	AncestralSequence* a2;
	TerminalSequence* t1;
	TerminalSequence* t2;

	int sl1;
	int sl2;
	int sfl1;
	int sfl2;

	int sAlpha;
	int nState;

	double small;
	double t;
	int k,m,n;
	double nullM1,nullM2;
	double matchBr1,matchBr2;

	// pointers to current functions
	void (PhyloMatchScore::*fwdp)(int,int);
	void (PhyloMatchScore::*bwdp)(int,int);
	void (PhyloMatchScore::*fullFwdp)(int,int);
	void (PhyloMatchScore::*fullBwdp)(int,int);

	// for two matrices
	void fwdMM(int j,int i);
	void bwdMM(int j,int i);
	void fullFwdMM(int j,int i);
	void fullBwdMM(int j,int i);

	// for two sequences
	void fwdSS(int j,int i);
	void bwdSS(int j,int i);
	void fullFwdSS(int j,int i);
	void fullBwdSS(int j,int i);

	// for a sequence and a matrix
	void fwdSM(int j,int i);
	void bwdSM(int j,int i);
	void fullFwdSM(int j,int i);
	void fullBwdSM(int j,int i);

	void fwdMS(int j,int i);
	void bwdMS(int j,int i);
	void fullFwdMS(int j,int i);
	void fullBwdMS(int j,int i);

	void logFwdMM(int j,int i);
	void logBwdMM(int j,int i);
	void logFullFwdMM(int j,int i);
	void logFullBwdMM(int j,int i);

	void logFwdSS(int j,int i);
	void logBwdSS(int j,int i);
	void logFullFwdSS(int j,int i);
	void logFullBwdSS(int j,int i);

	void logFwdSM(int j,int i);
	void logBwdSM(int j,int i);
	void logFullFwdSM(int j,int i);
	void logFullBwdSM(int j,int i);

	void logFwdMS(int j,int i);
	void logBwdMS(int j,int i);
	void logFullFwdMS(int j,int i);
	void logFullBwdMS(int j,int i);

	void computeSSMatrix();

public:
	~PhyloMatchScore();
	PhyloMatchScore(Sequence* seq1,Sequence* seq2);

	void computeFwd(int j,int i);
	void computeBwd(int j,int i);
	void computeFullFwd(int j,int i);
	void computeFullBwd(int j,int i);

	double fwdM(int k) { return fM->g(k); } // probability over all characters at the parent
	double bwdM(int k) { return bM->g(k); }

	double indelX(int k) { return idX->g(k); }
	double indelY(int k) { return idY->g(k); }

	double fullM(int k) { return flM->g(k); }

};

#endif
