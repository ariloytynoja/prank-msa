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
#ifndef READALIGNMENT_H
#define READALIGNMENT_H

#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include "sequence.h"
#include "treenode.h"
#include "phylomatchscore.h"
#include "flmatrix.h"
#include "dbmatrix.h"
#include "intmatrix.h"

class ReadAlignment
{
    static int count;
    static int sAlpha;
    static int nState;

    Sequence* seq1;
    Sequence* seq2;
    PhyloMatchScore* msr;
    TreeNode *tnode;

    int sl1;
    int sl2;
    int maxIndex;

    static int matrixSize;

    static FlMatrix* vX;
    static FlMatrix* vY;
    static FlMatrix* vM;
    static FlMatrix* xX;
    static FlMatrix* xM;
    static FlMatrix* wX;
    static FlMatrix* wM;
    static FlMatrix* yY;
    static FlMatrix* yM;
    static FlMatrix* zY;
    static FlMatrix* zM;


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

    Site *beg;
    Site *end;
    Site* newsite;

    //
    double maxFullScore;

    int totalSites;
    int countSites;
    int i,j,k;
    double small;

    int random_seed;
public:
    ReadAlignment();
    ~ReadAlignment();

    bool readSeqs(Sequence* s1,Sequence* s2,PhyloMatchScore* pms,TreeNode* tnode,std::vector<int>* path);

    void defineBegin();
    void defineEnd();

    void initialiseMatrices(int size);
    void cleanUp();
    double getMaxScore()
    {
        return maxFullScore;
    }

    bool rndBool();
    int rndInt(int i);
    double max(double a,double b);
    double max(double a,double b,double c);

    void setRandomSeed(int i) {random_seed = i; srand(random_seed);}
};

#endif
