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
#ifndef TERMINALSEQUENCE_H
#define TERMINALSEQUENCE_H

#include <string>
#include <iostream>
#include <sequence.h>
#include <string.h>

class TerminalSequence : public Sequence
{
    IntMatrix* seqvec;    // shortcut for unambiguous terminal seqs: profile is useless!

    bool isXGap(int )
    {
        return false;
    }
    bool isYGap(int )
    {
        return false;
    }

public:

    TerminalSequence(std::string* s);
    ~TerminalSequence();

    bool isGap(int )
    {
        return false;
    }
    int charAt(int i)
    {
        return seqvec->g(i);    // index unambiguous character at site i
    }

    bool isChildGap(int )
    {
        return false;
    }
    bool hasNeighborGaps(int )
    {
        return false;
    }
    bool isInsertion(int )
    {
        return false;
    }
    bool isPermInsertion(int )
    {
        return false;
    }
    void setPermInsertion(int ) {}

    bool fwdGapStarts(int )
    {
        return false;
    }
    bool fwdGapContinues(int )
    {
        return false;
    }
    bool fwdGapEnds(int )
    {
        return false;
    }

    bool fwdGapStartsNext(int )
    {
        return false;
    }
    bool fwdGapContinuesNext(int )
    {
        return false;
    }
    bool fwdGapEndsNext(int )
    {
        return false;
    }

    bool bwdGapStarts(int )
    {
        return false;
    }
    bool bwdGapContinues(int )
    {
        return false;
    }
    bool bwdGapEnds(int )
    {
        return false;
    }
    //
    bool fwdChildGapStarts(int )
    {
        return false;
    }
    bool fwdChildGapContinues(int )
    {
        return false;
    }
    bool fwdChildGapEnds(int )
    {
        return false;
    }

    bool fwdChildGapStartsNext(int )
    {
        return false;
    }
    bool fwdChildGapContinuesNext(int )
    {
        return false;
    }
    bool fwdChildGapEndsNext(int )
    {
        return false;
    }

    bool bwdChildGapStarts(int )
    {
        return false;
    }
    bool bwdChildGapContinues(int )
    {
        return false;
    }
    bool bwdChildGapEnds(int )
    {
        return false;
    }

    void removeGaps(std::string *si)
    {
        std::string s = "";
        for (int i=0; i<(int)si->length(); i++)
        {
            char c = si->at(i);
            if (c!='-')
            {
                s+=c;
            }
        }
        *si = s;
    }

};

#endif
