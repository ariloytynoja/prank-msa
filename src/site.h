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
#ifndef SITE_H
#define SITE_H

#define PRINT(STR, VAR) std::cout<< STR " = "<< VAR << std::endl

#include "intmatrix.h"
#include "flmatrix.h"
#include "dbmatrix.h"
#include "boolmatrix.h"

class Site
{
private:
    static BoolMatrix *anc;
    static BoolMatrix *nus;

    static IntMatrix *lSite; // index of neighbours
    static IntMatrix *rSite;

    static IntMatrix *cIndex1; // character index in seq1
    static IntMatrix *nIndex1; // character null index in seq1
    static IntMatrix *lIndex1; // character index left of this in seq1
    static IntMatrix *rIndex1; // character index right of this in seq1

    static IntMatrix *cIndex2; // character index in seq2
    static IntMatrix *nIndex2; // character null index in seq2
    static IntMatrix *lIndex2; // character index left of this in seq2
    static IntMatrix *rIndex2; // character index right of this in seq2

    static IntMatrix *currMS; // match state used
    static IntMatrix *currSS; // struct state used

    static IntMatrix *permIns; // permanent insertion

    static FlMatrix *vf; // for Viterbi path with a linear algorithm; ending probs for adjacent fragments
    static IntMatrix *vfM;
    static IntMatrix *vfS;
    static FlMatrix *vb;
    static IntMatrix *vbM;
    static IntMatrix *vbS;

    static FlMatrix* ffM;
    static FlMatrix* ffX;
    static FlMatrix* ffY;
    static FlMatrix* fbM;
    static FlMatrix* fbX;
    static FlMatrix* fbY;

    static DbMatrix* mcp;
    static FlMatrix* stp;
    static FlMatrix* pop;

    static int aSize;
    static int nState;

    static int count;
    int in;
public:
    Site();
    Site(int i);
    ~Site();

    void setMatrices(int longest,int slongest);
    void deleteMatrices();

    void next()
    {
        in = rSite->g(in);
    }
    void prev()
    {
        in = lSite->g(in);
    }

    void setNeighbours(Site *ls, Site *rs)
    {
        lSite->s(ls->getIndex(),in);
        rSite->s(rs->getIndex(),in);
        ls->setRSite(in);
        rs->setLSite(in);
    }

    void addNewSite()
    {
        in = count;
        count++;
    }
    void deleteLast()
    {
        count--;
    }

    void resetCounter()
    {
        count = 2;
    }
    int getLength()
    {
        return count;
    }

    void setIndex(int n)
    {
        in = n;
    }
    int getIndex()
    {
        return in;
    }

    void index(int n)
    {
        in = n;
    }
    int index()
    {
        return in;
    }

    void setLSite(int i)
    {
        this->lSite->s(i,in);
    }
    int getLSite()
    {
        return this->lSite->g(in);
    }
    void setRSite(int i)
    {
        this->rSite->s(i,in);
    }
    int getRSite()
    {
        return this->rSite->g(in);
    }

    void isAnchor(bool i)
    {
        anc->s(i,in);
    }
    bool isAnchor()
    {
        return anc->g(in);
    }
    void nullSite(bool i)
    {
        nus->s(i,in);
    }
    bool nullSite()
    {
        return nus->g(in);
    }

    void cInd1(int i)
    {
        cIndex1->s(i,in);
    }
    void nInd1(int i)
    {
        nIndex1->s(i,in);
    }
    void lInd1(int i)
    {
        lIndex1->s(i,in);
    }
    void rInd1(int i)
    {
        rIndex1->s(i,in);
    }

    void cInd2(int i)
    {
        cIndex2->s(i,in);
    }
    void nInd2(int i)
    {
        nIndex2->s(i,in);
    }
    void lInd2(int i)
    {
        lIndex2->s(i,in);
    }
    void rInd2(int i)
    {
        rIndex2->s(i,in);
    }

    void currMatchState(int i)
    {
        currMS->s(i,in);
    }
    void currModelState(int i)
    {
        currSS->s(i,in);
    }

    void permInsertion(int i)
    {
        permIns->s(i,in);
    }

    int cInd1()
    {
        return cIndex1->g(in);
    }
    int nInd1()
    {
        return nIndex1->g(in);
    }
    int lInd1()
    {
        return lIndex1->g(in);
    }
    int rInd1()
    {
        return rIndex1->g(in);
    }

    int cInd2()
    {
        return cIndex2->g(in);
    }
    int nInd2()
    {
        return nIndex2->g(in);
    }
    int lInd2()
    {
        return lIndex2->g(in);
    }
    int rInd2()
    {
        return rIndex2->g(in);
    }

    int currMatchState()
    {
        return currMS->g(in);
    }
    int currModelState()
    {
        return currSS->g(in);
    }

    int permInsertion()
    {
        return permIns->g(in);
    }

    void vitf(float i)
    {
        vf->s(i,in);
    }
    void vitfM(int i)
    {
        vfM->s(i,in);
    }
    void vitfS(int i)
    {
        vfS->s(i,in);
    }
    void vitb(float i)
    {
        vb->s(i,in);
    }
    void vitbM(int i)
    {
        vbM->s(i,in);
    }
    void vitbS(int i)
    {
        vbS->s(i,in);
    }

    float vitf()
    {
        return vf->g(in);
    }
    int vitfM()
    {
        return vfM->g(in);
    }
    int vitfS()
    {
        return vfS->g(in);
    }
    float vitb()
    {
        return vb->g(in);
    }
    int vitbM()
    {
        return vbM->g(in);
    }
    int vitbS()
    {
        return vbS->g(in);
    }

    void fullFwdX(float i,int k)
    {
        ffX->s(i,k,in);
    }
    void fullFwdY(float i,int k)
    {
        ffY->s(i,k,in);
    }
    void fullFwdM(float i,int k)
    {
        ffM->s(i,k,in);
    }
    void fullBwdX(float i,int k)
    {
        fbX->s(i,k,in);
    }
    void fullBwdY(float i,int k)
    {
        fbY->s(i,k,in);
    }
    void fullBwdM(float i,int k)
    {
        fbM->s(i,k,in);
    }

    float fullFwdX(int k)
    {
        return ffX->g(k,in);
    }
    float fullFwdY(int k)
    {
        return ffY->g(k,in);
    }
    float fullFwdM(int k)
    {
        return ffM->g(k,in);
    }
    float fullBwdX(int k)
    {
        return fbX->g(k,in);
    }
    float fullBwdY(int k)
    {
        return fbY->g(k,in);
    }
    float fullBwdM(int k)
    {
        return fbM->g(k,in);
    }

    void mlCharProb(double i,int k,int j)
    {
        mcp->s(i,k,j,in);
    }
    void stateProb(float i,int k)
    {
        stp->s(i,k,in);
    }
    void postProb(float i)
    {
        pop->s(i,in);
    }

    double mlCharProb(int k,int j)
    {
        return mcp->g(k,j,in);
    }
    float stateProb(int k)
    {
        return stp->g(k,in);
    }
    float postProb()
    {
        return pop->g(in);
    }

    void setASize(int i)
    {
        aSize = i;
    }
    void setNState(int i)
    {
        nState = i;
    }

};


#endif
