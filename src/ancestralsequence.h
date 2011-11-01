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
#ifndef ANCESTRALSEQUENCE_H
#define ANCESTRALSEQUENCE_H

#include <sequence.h>

class AncestralSequence : public Sequence
{
    FlMatrix* seqmat;     // sequence data in profile format
    FlMatrix* logseqmat;     // sequence data in profile format
    DbMatrix* mlCharProb; // character probabilities given the tree below
    FlMatrix* postProb;   // alignment posterior probabilities
    FlMatrix* stateProb;  // state posterior probabilities
    FlMatrix* priorProb;

    IntMatrix* lcIndex;       // index of corresponding characters in child seqs
    IntMatrix* rcIndex;

    IntMatrix* xGapSite;      // is this gapsite
    IntMatrix* yGapSite;      // is this gapsite
    IntMatrix* childGapSite;  // is this gapsite in one of the children
    IntMatrix* insertionSite; //
    IntMatrix* permInsertionSite; //
    IntMatrix* realIndex;     // index for non-insertion sites

public:
    AncestralSequence();
    ~AncestralSequence();

    float postProbAt(int i)
    {
        return postProb->g(i);
    }
    float stateProbAt(int k,int i)
    {
        return stateProb->g(k,i);
    }

    int getLIndex(int i)
    {
        return lcIndex->g(i);
    }
    int getRIndex(int i)
    {
        return rcIndex->g(i);
    }

    bool isGap(int i)
    {
        if (i>=0 && i<length())
        {
            return xGapSite->g(i)==1 || yGapSite->g(i)==1;
        }
        else
        {
            return false;
        }
    }
    bool isXGap(int i)
    {
        if (i>=0 && i<length())
        {
            return xGapSite->g(i)==1;
        }
        else
        {
            return false;
        }
    }
    bool isYGap(int i)
    {
        if (i>=0 && i<length())
        {
            return yGapSite->g(i)==1;
        }
        else
        {
            return false;
        }
    }

    bool fwdGapStarts(int i)
    {
        return ( (!isXGap(i-2) && isXGap(i-1)) || (!isYGap(i-2) && isYGap(i-1)) );
    }
    bool fwdGapContinues(int i)
    {
        return ( (isXGap(i-2) && isXGap(i-1)) || (isYGap(i-2) && isYGap(i-1)) || (isXGap(i-2) && isYGap(i-1)) || (isYGap(i-2) && isXGap(i-1)) );
    }
    bool fwdGapEnds(int i)
    {
        return ( (isXGap(i-2) && !isXGap(i-1)) || (isYGap(i-2) && !isYGap(i-1)) );
    }

    bool fwdGapStartsNext(int i)
    {
        return ( (!isXGap(i-1) && isXGap(i)) || (!isYGap(i-1) && isYGap(i)) );
    }
    bool fwdGapContinuesNext(int i)
    {
        return ( (isXGap(i-1) && isXGap(i)) || (isYGap(i-1) && isYGap(i)) || (isXGap(i-1) && isYGap(i)) || (isYGap(i-1) && isXGap(i)) );
    }
    bool fwdGapEndsNext(int i)
    {
        return ( (isXGap(i-1) && !isXGap(i)) || (isYGap(i-1) && !isYGap(i)) );
    }

    bool bwdGapStarts(int i)
    {
        return ( (isXGap(i-1) && !isXGap(i)) || (isYGap(i-1) && !isYGap(i)) );
    }
    bool bwdGapContinues(int i)
    {
        return ( (isXGap(i-1) && isXGap(i)) || (isYGap(i-1) && isYGap(i)) || (isXGap(i-1) && isYGap(i)) || (isYGap(i-1) && isXGap(i)) );
    }
    bool bwdGapEnds(int i)
    {
        return ( (!isXGap(i-1) && isXGap(i)) || (!isYGap(i-1) && isYGap(i)) );
    }
    //
    bool fwdChildGapStarts(int i)
    {
        return ( !isChildGap(i-2) && isChildGap(i-1) && ( isXGap(i-1) || isYGap(i-1) ) );
    }
    bool fwdChildGapContinues(int i)
    {
        return ( isChildGap(i-2) && isChildGap(i-1) );
    }
    bool fwdChildGapEnds(int i)
    {
        return ( isChildGap(i-2) && !isChildGap(i-1) && ( isXGap(i-2) || isYGap(i-2) ) );
    }

    bool fwdChildGapStartsNext(int i)
    {
        return ( !isChildGap(i-1) && isChildGap(i) && ( isXGap(i) || isYGap(i) ) );
    }
    bool fwdChildGapContinuesNext(int i)
    {
        return ( isChildGap(i-1) && isChildGap(i) );
    }
    bool fwdChildGapEndsNext(int i)
    {
        return ( isChildGap(i-1) && !isChildGap(i) && ( isXGap(i-1) || isYGap(i-1) ) );
    }

    bool bwdChildGapStarts(int i)
    {
        return ( isChildGap(i-1) && !isChildGap(i) && ( isXGap(i-1) || isYGap(i-1) ) );
    }
    bool bwdChildGapContinues(int i)
    {
        return ( isChildGap(i-1) && isChildGap(i) );
    }
    bool bwdChildGapEnds(int i)
    {
        return ( !isChildGap(i-1) && isChildGap(i) && ( isXGap(i) || isYGap(i) ) );
    }

    void setChildGaps(Sequence *l,Sequence *r);
    void setRealIndex(bool left);

    double mlCharProbAt(int j,int i,int k);
    double mlCharProbAtF(int j,int i,int k);


    void cleanSpace();

    bool isInsertion(int i)
    {
        if (i>=0 && i<length())
        {
            return insertionSite->g(i)==1;
        }
        else
        {
            return false;
        }
    };
    bool isPermInsertion(int i)
    {
        if (i>=0 && i<length())
        {
            return permInsertionSite->g(i)==1;
        }
        else
        {
            return false;
        }
    }
    void setPermInsertion(int i)
    {
        if (i>=0 && i<length())
        {
            permInsertionSite->s(1,i);
        }
    }

    bool isChildGap(int i)
    {
        if (i>=0 && i<length() && childGapSite->g(i)==1)
        {
            return true;
        }
        else
        {
            return false;
        }
    };

    bool hasNeighborGaps(int i)
    {
        if ( isGap(i-2) || isGap(i-1) || isGap(i) || isGap(i+1) || isGap(i+2) || isChildGap(i-2) || isChildGap(i-1) || isChildGap(i) || isChildGap(i+1) || isChildGap(i+2) )
        {
            return true;
        }
        else
        {
            return false;
        }
    };

    void writeSequence(std::string name);

    void setGappedSeq(std::string *s)
    {
        gappedseq = *s;
    }
};

#endif
