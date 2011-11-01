/***************************************************************************
 *   Copyright (C) 2005 by Ari Loytynoja                                   *
 *   ari@ebi.ac.uk                                                         *
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

#ifndef SEQUENCE_H
#define SEQUENCE_H

/**
 * Sequence described as character probabilities.
 */

#include <string>
#include <iostream>
#include "site.h"
#include "flmatrix.h"
#include "dbmatrix.h"
#include "intmatrix.h"

extern bool PRIORS;

class Sequence
{

protected:
    std::string charseq;
    std::string gappedseq;

    int seqLength;    // length of sequence
    int realLength;   // length when insertion-sites are skipped
    bool terminal;    // is/not terminal
    int sAlpha;

    bool hasPrior;
    int i,j,k;
public:
    ~Sequence();
    Sequence();

    bool isTerminal()
    {
        return terminal;
    }
    int length()
    {
        return seqLength;
    }
    int lengthF()
    {
        return realLength;
    }
    int gappedLength()
    {
        return gappedseq.length();
    }

    bool prealignedGapAt(int i)
    {
        return gappedseq.at(i)=='-';
    }

    int charAt(int )
    {
        return -1;
    }

    double mlCharProbAt(int j,int i,int k);
    double mlCharProbAtF(int j,int i,int k);

    int getLIndex(int i)
    {
        return i;
    }
    int getRIndex(int i)
    {
        return i;
    }

    std::string* getMLsequence()
    {
        return &charseq;
    }
    std::string* getGappedSeq()
    {
        return &gappedseq;
    }

    virtual bool isGap(int i) = 0;
    virtual bool isXGap(int i) = 0;
    virtual bool isYGap(int i) = 0;
    virtual bool isChildGap(int i) = 0;
    virtual bool hasNeighborGaps(int i) = 0;
    virtual bool isInsertion(int i) = 0;
    virtual bool isPermInsertion(int i) = 0;
    virtual void setPermInsertion(int i) = 0;

    virtual bool fwdGapStarts(int i) = 0;
    virtual bool fwdGapContinues(int i) = 0;
    virtual bool fwdGapEnds(int i) = 0;

    virtual bool fwdGapStartsNext(int i) = 0;
    virtual bool fwdGapContinuesNext(int i) = 0;
    virtual bool fwdGapEndsNext(int i) = 0;

    virtual bool bwdGapStarts(int i) = 0;
    virtual bool bwdGapContinues(int i) = 0;
    virtual bool bwdGapEnds(int i) = 0;
//
    virtual bool fwdChildGapStarts(int i) = 0;
    virtual bool fwdChildGapContinues(int i) = 0;
    virtual bool fwdChildGapEnds(int i) = 0;

    virtual bool fwdChildGapStartsNext(int i) = 0;
    virtual bool fwdChildGapContinuesNext(int i) = 0;
    virtual bool fwdChildGapEndsNext(int i) = 0;

    virtual bool bwdChildGapStarts(int i) = 0;
    virtual bool bwdChildGapContinues(int i) = 0;
    virtual bool bwdChildGapEnds(int i) = 0;

    void cleanSpace() {};
    void writeSequence(std::string ) {}

};

#endif
