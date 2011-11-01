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
#ifndef HIRSCHBERG_H
#define HIRSCHBERG_H

#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include "sequence.h"
#include "phylomatchscore.h"
#include "flmatrix.h"
#include "dbmatrix.h"
#include "intmatrix.h"

class Hirschberg
{
protected:
    static int count;
    static int sAlpha;
    static int nState;

    Sequence* seq1;
    Sequence* seq2;
    PhyloMatchScore* msr;

    int sl1;
    int sl2;
    int mLen;
    int mLen2;
    int mSize;
    int maxIndex;

    int prevFwd;
    int currFwd;
    int currBwd;
    int nextBwd;

    static int alignmentNumber;
    static int matrixSize;

    static FlMatrix* fwdvX;
    static FlMatrix* fwdvY;
    static FlMatrix* fwdvM; // starting/ending values for Viterbi
    static FlMatrix* bwdvX;
    static FlMatrix* bwdvY;
    static FlMatrix* bwdvM;
    static FlMatrix* fwdxX;
    static FlMatrix* fwdxM;                  // starting/ending values for Viterbi (skip-X)
    static FlMatrix* bwdxX;
    static FlMatrix* bwdxM;
    static FlMatrix* fwdwX;
    static FlMatrix* fwdwM;                  // starting/ending values for Viterbi (skip-child-X)
    static FlMatrix* bwdwX;
    static FlMatrix* bwdwM;
    static FlMatrix* fwdyY;
    static FlMatrix* fwdyM;                  // starting/ending values for Viterbi (skip-Y)
    static FlMatrix* bwdyY;
    static FlMatrix* bwdyM;
    static FlMatrix* fwdzY;
    static FlMatrix* fwdzM;                  // starting/ending values for Viterbi (skip-child-Y)
    static FlMatrix* bwdzY;
    static FlMatrix* bwdzM;

    // matrices for the two rows kept in memory
    static DbMatrix* fVM1;
    static DbMatrix* fVX1;
    static DbMatrix* fVY1;
    static DbMatrix* fXM1;
    static DbMatrix* fXX1;
    static DbMatrix* fWM1;
    static DbMatrix* fWX1;
    static DbMatrix* fYM1;
    static DbMatrix* fYY1;
    static DbMatrix* fZM1;
    static DbMatrix* fZY1;

    static DbMatrix* fVM2;
    static DbMatrix* fVX2;
    static DbMatrix* fVY2;
    static DbMatrix* fXM2;
    static DbMatrix* fXX2;
    static DbMatrix* fWM2;
    static DbMatrix* fWX2;
    static DbMatrix* fYM2;
    static DbMatrix* fYY2;
    static DbMatrix* fZM2;
    static DbMatrix* fZY2;

    static DbMatrix* bVM1;
    static DbMatrix* bVX1;
    static DbMatrix* bVY1;
    static DbMatrix* bXM1;
    static DbMatrix* bXX1;
    static DbMatrix* bWM1;
    static DbMatrix* bWX1;
    static DbMatrix* bYM1;
    static DbMatrix* bYY1;
    static DbMatrix* bZM1;
    static DbMatrix* bZY1;

    static DbMatrix* bVM2;
    static DbMatrix* bVX2;
    static DbMatrix* bVY2;
    static DbMatrix* bXM2;
    static DbMatrix* bXX2;
    static DbMatrix* bWM2;
    static DbMatrix* bWX2;
    static DbMatrix* bYM2;
    static DbMatrix* bYY2;
    static DbMatrix* bZM2;
    static DbMatrix* bZY2;

    // matrices for pointers; just forward
    static IntMatrix* ptVM;
    static IntMatrix* ptVX;
    static IntMatrix* ptVY;
    static IntMatrix* ptXM;
    static IntMatrix* ptXX;
    static IntMatrix* ptWM;
    static IntMatrix* ptWX;
    static IntMatrix* ptYM;
    static IntMatrix* ptYY;
    static IntMatrix* ptZM;
    static IntMatrix* ptZY;

    // Temp variables
    int sX,sY,sM,sxX,sxM,swX,swM,syY,syM,szY,szM;    // state
    double mX,mY,mM,mxX,mxM,mwX,mwM,myY,myM,mzY,mzM; // max
    double cX,cY,cM,cxX,cxM,cwX,cwM,cyY,cyM,czY,czM; // current

    // pointers to the two rows
    DbMatrix* cfVX;
    DbMatrix* cfVY;
    DbMatrix* cfVM;
    DbMatrix* cfXX;
    DbMatrix* cfXM;
    DbMatrix* cfWX;
    DbMatrix* cfWM;
    DbMatrix* cfYY;
    DbMatrix* cfYM;
    DbMatrix* cfZY;
    DbMatrix* cfZM;

    DbMatrix* cbVX;
    DbMatrix* cbVY;
    DbMatrix* cbVM;
    DbMatrix* cbXX;
    DbMatrix* cbXM;
    DbMatrix* cbWX;
    DbMatrix* cbWM;
    DbMatrix* cbYY;
    DbMatrix* cbYM;
    DbMatrix* cbZY;
    DbMatrix* cbZM;

    DbMatrix* pVX;
    DbMatrix* pVY;
    DbMatrix* pVM;
    DbMatrix* pXX;
    DbMatrix* pXM;
    DbMatrix* pWX;
    DbMatrix* pWM;
    DbMatrix* pYY;
    DbMatrix* pYM;
    DbMatrix* pZY;
    DbMatrix* pZM;

    // tmp pointers needed during the re-pointing
    DbMatrix* tmpVX;
    DbMatrix* tmpVY;
    DbMatrix* tmpVM;
    DbMatrix* tmpXX;
    DbMatrix* tmpXM;
    DbMatrix* tmpWX;
    DbMatrix* tmpWM;
    DbMatrix* tmpYY;
    DbMatrix* tmpYM;
    DbMatrix* tmpZY;
    DbMatrix* tmpZM;

    //
    double maxFullScore;


    Site *beg;
    Site *end;
    Site* newsite;

    int nanch;
    IntMatrix* anchors;

    int totalSites;
    int countSites;
    int i,j;
    int k;
public:
    double small;
    Hirschberg();
    ~Hirschberg();

    void defineBegin();
    void defineSite(int i);
    void defineESite(int l,int r);
    void defineEnd();

    void getMidSite(int s1,int e1,int s2,int e2);
    void alignSeqs(Sequence* s1,Sequence* s2,PhyloMatchScore* pms);
    void divideSeq();

    double getMaxScore()
    {
        return maxFullScore;
    }

    bool rndBool();
    int rndInt(int i);
    double max(double a,double b);
    double max(double a,double b,double c);

    void printMatrix(std::string n,int i,DbMatrix* m);
    void printMatrix(std::string n,int i,IntMatrix* m);

    void initialiseIndex(int *min,int *max);
    void initialiseMatrices(int size);
    void cleanUp();
};

#endif

#ifndef STRUCTCELL_H
#define STRUCTCELL_H

struct Cell
{
    int prev;
    int curr;
    int k;
};
#endif
