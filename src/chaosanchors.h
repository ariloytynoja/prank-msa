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
#ifndef CHAOSANCHORS_H
#define CHAOSANCHORS_H

#include <string>
#include "intmatrix.h"

/**
 * anchoring with chaos. need more work!
 */

class ChaosAnchors
{

    std::string* left;
    std::string* right;

    IntMatrix* anchps;
    int nanch;

    IntMatrix* subst;
    int gOpen;
    int gExt;

    std::string alpha;
    int small;

    bool rndBool();
    int rndInt(int i);
    int max(int a,int b);
    int max(int a,int b,int c);

    int i,j;
public:
    ChaosAnchors(std::string* s1,std::string* s2);
    ~ChaosAnchors();
    IntMatrix* getAnchors(int *na);
    void alignRegions(int *c1,int *c2,int beg1,int end1,int beg2,int end2) ;
    void reverseAlignRegions(int *c1,int *c2,int beg1,int end1,int beg2,int end2) ;
    int matchScore(int site1,int site2);
};

#endif
