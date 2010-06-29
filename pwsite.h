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
#ifndef PWSITE_H
#define PWSITE_H

#include "flmatrix.h"
#include "intmatrix.h"
#include <iostream>

class PwSite {
private:
	static IntMatrix *lSite; // index of neighbours
	static IntMatrix *rSite;

	static IntMatrix *cIndex1; // character index in seq1
	static IntMatrix *lIndex1; // character index left of this in seq1
	static IntMatrix *rIndex1; // character index right of this in seq1

	static IntMatrix *cIndex2; // character index in seq2
	static IntMatrix *lIndex2; // character index left of this in seq2
	static IntMatrix *rIndex2; // character index right of this in seq2

	static IntMatrix *currMS; // match state used

	static IntMatrix *vfX; // for Viterbi path with a linear algorithm; ending probs for adjacent fragments
	static IntMatrix *vfY; // these three for forward start site
	static IntMatrix *vfM;
	static IntMatrix *vbX; // for Viterbi path with a linear algorithm; starting probs for adjacent fragments
	static IntMatrix *vbY; // these three for backward start site
	static IntMatrix *vbM;

	static int aSize;
	static int count;
    int in;
public:
	PwSite();
	PwSite(int i);
	~PwSite();

	void setMatrices(int longest,int slongest);
	void deleteMatrices();

	void next() { in = rSite->g(in); }
	void prev() { in = lSite->g(in); }

	void setNeighbours(PwSite *ls, PwSite *rs) {
		lSite->s(ls->getIndex(),in);
		rSite->s(rs->getIndex(),in);
		ls->setRSite(in);
		rs->setLSite(in);
	}

	void addNewSite() { in = count; count++; }
	void deleteLast() { count--; }

	void resetCounter() { count =2; }

	void setASize(int i){ aSize = i; }

	void setIndex(int n) { in = n; }
	int getIndex() { return in; }

	void index(int n) { in = n; }
	int index() { return in; }

	void setLSite(int i) { this->lSite->s(i,in); }
	int getLSite() { return this->lSite->g(in); }
	void setRSite(int i) { this->rSite->s(i,in); }
	int getRSite() { return this->rSite->g(in); }

	void cInd1(int i) { cIndex1->s(i,in);}
	void lInd1(int i) { lIndex1->s(i,in);}
	void rInd1(int i) { rIndex1->s(i,in);}

	void cInd2(int i) { cIndex2->s(i,in);}
	void lInd2(int i) { lIndex2->s(i,in);}
	void rInd2(int i) { rIndex2->s(i,in);}

	void currMatchState(int i) { currMS->s(i,in);}

	void vitfX(int i) { vfX->s(i,in);}
	void vitfY(int i) { vfY->s(i,in);}
	void vitfM(int i) { vfM->s(i,in);}
	void vitbX(int i) { vbX->s(i,in);}
	void vitbY(int i) { vbY->s(i,in);}
	void vitbM(int i) { vbM->s(i,in);}

	int cInd1() { return cIndex1->g(in);}
	int lInd1() { return lIndex1->g(in);}
	int rInd1() { return rIndex1->g(in);}

	int cInd2() { return cIndex2->g(in);}
	int lInd2() { return lIndex2->g(in);}
	int rInd2() { return rIndex2->g(in);}

	int currMatchState() { return currMS->g(in);}

	int vitfX() { return vfX->g(in);}
	int vitfY() { return vfY->g(in);}
	int vitfM() { return vfM->g(in);}
	int vitbX() { return vbX->g(in);}
	int vitbY() { return vbY->g(in);}
	int vitbM() { return vbM->g(in);}

};


#endif
