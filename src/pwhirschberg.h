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
#ifndef PWHIRSCHBERG_H
#define PWHIRSCHBERG_H

#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include "pwsite.h"
#include "flmatrix.h"
#include "intmatrix.h"

class PwHirschberg
{
protected:
    static int count;
    int sAlpha;
    std::string alpha;

    int deltaX1,deltaX2,epsilonX,deltaY1,deltaY2,epsilonY;
    IntMatrix* substScores;

    std::string* seq1;
    std::string* seq2;

    int sl1;
    int sl2;
    int mLen;
    int mSize;
    int maxIndex;
    int small;

    int fwdvX;
    int fwdvY;
    int fwdvM; // starting/ending values for Viterbi
    int bwdvX;
    int bwdvY;
    int bwdvM;

    // matrices for the two rows kept in memory
    IntMatrix* fVM1;
    IntMatrix* fVX1;
    IntMatrix* fVY1;
    IntMatrix* fVM2;
    IntMatrix* fVX2;
    IntMatrix* fVY2;

    IntMatrix* bVM1;
    IntMatrix* bVX1;
    IntMatrix* bVY1;
    IntMatrix* bVM2;
    IntMatrix* bVX2;
    IntMatrix* bVY2;

    // matrices for pointers; just forward
    IntMatrix* ptVM;
    IntMatrix* ptVX;
    IntMatrix* ptVY;

    // Temp variables
    int sX,sY,sM;    // state
    int mX,mY,mM; // max
    int cX,cY,cM; // current

    // pointers to the two rows
    IntMatrix* cfVX;
    IntMatrix* cfVY;
    IntMatrix* cfVM;
    IntMatrix* cbVX;
    IntMatrix* cbVY;
    IntMatrix* cbVM;
    IntMatrix* pVX;
    IntMatrix* pVY;
    IntMatrix* pVM;

    // tmp pointers needed during the re-pointing
    IntMatrix* tmpVX;
    IntMatrix* tmpVY;
    IntMatrix* tmpVM;

    int maxFullScore;
    int matchScore;

    PwSite *beg;
    PwSite *end;
    PwSite* pwsite;

    int totalSites;
    int countSites;
    int i,j,k;
    static int depth;
public:

    PwHirschberg(int length);
    ~PwHirschberg();
    void setSequences(std::string* s1,std::string* s2);
    void setModel(IntMatrix* scores,int delta, int epsilon);

    void defineBegin();
    void defineESite(int l,int r);
    void defineEnd();

    void getMidSite(int s1,int e1,int s2,int e2);
    void alignSeqs();
    void getAnchors();
    void divideSeq();

    double getMaxScore()
    {
        return maxFullScore;
    }

    bool rndBool();
    int rndInt(int i);
    int max(int a,int b);
    int max(int a,int b,int c);

    void printMatrix(std::string n,int i,IntMatrix* m);

    void cleanUp();

    void computeFwd(int j,int i);
    void computeBwd(int j,int i);

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
