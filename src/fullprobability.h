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
#ifndef FULLPROBABILITY_H
#define FULLPROBABILITY_H

#include <string>
#include "sequence.h"
#include "phylomatchscore.h"
#include "site.h"
#include "dbmatrix.h"

/**
 * Forward and backward loops to compute full probability within a band.
 */

class FullProbability{

   Sequence* seq1;
   Sequence* seq2;

   int sAlpha;
   int nState;
   int width;

   IntMatrix* minBIndex;
   IntMatrix* maxBIndex;
   IntMatrix* diffIndex;

   double small;
   PhyloMatchScore *msr;

   double maxFwdScore;
   double maxBwdScore;

   // initialise the indeces for banding
   void initialiseIndex(Site* sites);

   int xgap;
   int i,j,k;
public:
    FullProbability();

    ~FullProbability();
   FullProbability(Sequence* s1,Sequence* s2,PhyloMatchScore *msr);

   void alignSeqs();
   void alignBand();

   double getMaxFwdScore() { return maxFwdScore; }
   double getMaxBwdScore() { return maxBwdScore; }

   void printMatrix(std::string s,int i,DbMatrix* n);
};

#endif
