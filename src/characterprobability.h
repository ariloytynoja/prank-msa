/***************************************************************************
 *   Copyright (C) 2005 by Ari Loytynoja   *
 *   ari@ebi.ac.uk   *
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
#ifndef CHARACTERPROBABILITY_H
#define CHARACTERPROBABILITY_H

#define PRINT(STR, VAR) std::cout<< STR " = "<< VAR << std::endl

/**
  A new way to compute the structure state probability; now phylogeny is taken into account.
 */

#include "sequence.h"
#include "ancestralsequence.h"
#include "terminalsequence.h"
#include "site.h"

class CharacterProbability{
    float fwdScore;
    float bwdScore;

	int nState;
	int sAlpha;

	int li;
	int ri;

	double small,sum,sum1,sum2;

	int skipMatch;

	AncestralSequence* a1;
	AncestralSequence* a2;
	TerminalSequence* t1;
	TerminalSequence* t2;

	void logScoresSS();
	void logScoresSM();
	void logScoresMS();
	void logScoresMM();
	void scoresSS();
	void scoresSM();
	void scoresMS();
	void scoresMM();

	Site *cSite;
	Site *pSite;
	Site *sSite;

    int j,k,l,m,n;
public:
    CharacterProbability(Sequence* sq1,Sequence* sq2);
    ~CharacterProbability();

    float getFwdScore() { return fwdScore; }
    float getBwdScore() { return bwdScore; }
};

#endif
