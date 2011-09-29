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

#include <iostream>
#include <cmath>
#include "config.h"
#include "phylomatchscore.h"

using namespace std;

PhyloMatchScore::~PhyloMatchScore()
{
// 	cout<<"delete PhyloMatchScore"<<endl;
    delete flM;
    delete fM;
    delete bM;
    delete idX;
    delete idY;
	if (s1->isTerminal() && s2->isTerminal()) {
		delete match;
		delete gap;
	}
}

PhyloMatchScore::PhyloMatchScore(Sequence* seq1,Sequence* seq2)
{
//        cout<<"create PhyloMatchScore"<<endl;
    s1 = seq1;
    s2 = seq2;

    sl1 = s1->length(); // length for Viterbi
    sl2 = s2->length();

    sfl1 = s1->lengthF(); // length for full probability (possibly shorter sequences)
    sfl2 = s2->lengthF();

    sAlpha = hmm->getASize();
    nState = hmm->getNStates();

    fM = new DbMatrix(nState,"fM");
    bM = new DbMatrix(nState,"bM");

    flM = new DbMatrix(nState,"flM");

    idX = new DbMatrix(nState,"idX");
    idY = new DbMatrix(nState,"idY");

    small = -HUGE_VAL;


    if (LOGVALUES) {
        if (s1->isTerminal() && s2->isTerminal()) {
            t1 = static_cast<TerminalSequence*>(s1);
            t2 = static_cast<TerminalSequence*>(s2);
            fwdp = &PhyloMatchScore::logFwdSS;
            bwdp = &PhyloMatchScore::logBwdSS;
            fullFwdp = &PhyloMatchScore::logFwdSS;
            fullBwdp = &PhyloMatchScore::logBwdSS;
            computeSSMatrix();
        }
        else if (s1->isTerminal() && !s2->isTerminal()) {
            t1 = static_cast<TerminalSequence*>(s1);
            a2 = static_cast<AncestralSequence*>(s2);
            fwdp = &PhyloMatchScore::logFwdSM;
            bwdp = &PhyloMatchScore::logBwdSM;
            fullFwdp = &PhyloMatchScore::logFullFwdSM;
            fullBwdp = &PhyloMatchScore::logFullBwdSM;
        }
        else if (!s1->isTerminal() && s2->isTerminal()) {
            a1 = static_cast<AncestralSequence*>(s1);
            t2 = static_cast<TerminalSequence*>(s2);
            fwdp = &PhyloMatchScore::logFwdMS;
            bwdp = &PhyloMatchScore::logBwdMS;
            fullFwdp = &PhyloMatchScore::logFullFwdMS;
            fullBwdp = &PhyloMatchScore::logFullBwdMS;
        }
        else if (!s1->isTerminal() && !s2->isTerminal()) {
            a1 = static_cast<AncestralSequence*>(s1);
            a2 = static_cast<AncestralSequence*>(s2);
            fwdp = &PhyloMatchScore::logFwdMM;
            bwdp = &PhyloMatchScore::logBwdMM;
            fullFwdp = &PhyloMatchScore::logFullFwdMM;
            fullBwdp = &PhyloMatchScore::logFullBwdMM;
        }
    } else {
        if (s1->isTerminal() && s2->isTerminal()) {
            t1 = static_cast<TerminalSequence*>(s1);
            t2 = static_cast<TerminalSequence*>(s2);
            fwdp = &PhyloMatchScore::fwdSS;
            bwdp = &PhyloMatchScore::bwdSS;
            fullFwdp = &PhyloMatchScore::fwdSS;
            fullBwdp = &PhyloMatchScore::bwdSS;
            computeSSMatrix();
        }
        else if (s1->isTerminal() && !s2->isTerminal()) {
            t1 = static_cast<TerminalSequence*>(s1);
            a2 = static_cast<AncestralSequence*>(s2);
            fwdp = &PhyloMatchScore::fwdSM;
            bwdp = &PhyloMatchScore::bwdSM;
            fullFwdp = &PhyloMatchScore::fullFwdSM;
            fullBwdp = &PhyloMatchScore::fullBwdSM;
        }
        else if (!s1->isTerminal() && s2->isTerminal()) {
            a1 = static_cast<AncestralSequence*>(s1);
            t2 = static_cast<TerminalSequence*>(s2);
            fwdp = &PhyloMatchScore::fwdMS;
            bwdp = &PhyloMatchScore::bwdMS;
            fullFwdp = &PhyloMatchScore::fullFwdMS;
            fullBwdp = &PhyloMatchScore::fullBwdMS;
        }
        else if (!s1->isTerminal() && !s2->isTerminal()) {
            a1 = static_cast<AncestralSequence*>(s1);
            a2 = static_cast<AncestralSequence*>(s2);
            fwdp = &PhyloMatchScore::fwdMM;
            bwdp = &PhyloMatchScore::bwdMM;
            fullFwdp = &PhyloMatchScore::fullFwdMM;
            fullBwdp = &PhyloMatchScore::fullBwdMM;
        }
    }
}

void PhyloMatchScore::computeFwd(int j,int i)
{
    (this->*fwdp)(j,i);
}

void PhyloMatchScore::computeBwd(int j,int i)
{
    (this->*bwdp)(j,i);
}

void PhyloMatchScore::computeFullFwd(int j,int i)
{
    (this->*fullFwdp)(j,i);
}

void PhyloMatchScore::computeFullBwd(int j,int i)
{
    (this->*fullBwdp)(j,i);
}

/////////

void PhyloMatchScore::fwdMM(int j,int i)
{
    FOR(k,nState) {

        fM->s(0,k);
        idX->s(0,k);
        idY->s(0,k);

        nullM1=nullM2=0.0;

        FOR(n,sAlpha) {

            matchBr1=matchBr2=0.0;

            FOR(m,sAlpha) {
                // match
                if (j>0 && i>0) {
                    matchBr1 += hmm->charSbProbL(k,n,m)*a1->mlCharProbAt(m,j-1,k);
                    matchBr2 += hmm->charSbProbR(k,n,m)*a2->mlCharProbAt(m,i-1,k);
                }

                // x-gap
                if (j>0) {
                    if (n==0) {
                        idX->a(hmm->charBgFreq(k,m)*a1->mlCharProbAt(m,j-1,k),k);
                        nullM1 += hmm->nullBgFreq(m)*a1->mlCharProbAt(m,j-1,k);
                    }
                }

                // y-gap
                if (i>0) {
                    if (n==0) {
                        idY->a(hmm->charBgFreq(k,m)*a2->mlCharProbAt(m,i-1,k),k);
                        nullM2 += hmm->nullBgFreq(m)*a2->mlCharProbAt(m,i-1,k);
                    }
                }
            }

            t = hmm->charBgFreq(k,n)*matchBr1*matchBr2/(nullM1*nullM2);
            fM->a( t, k);
        }

        fM->clog( k );

        idX->d( nullM1, k );
        idX->clog( k );

        idY->d( nullM2, k );
        idY->clog( k );
    }

    return;
}

void PhyloMatchScore::fwdSS(int j,int i)
{

    FOR(k,nState) {
        if (j>0 && i>0) {
            fM->s( match->g( t1->charAt(j-1), t2->charAt(i-1), k ), k );
            flM->s( match->g( t1->charAt(j-1), t2->charAt(i-1), k ), k );
        }
        else {
            fM->s( 0, k );
            flM->s( 0, k );
        }
        if (j>0)
            idX->s( gap->g( t1->charAt(j-1), k ), k );
        else
            idX->s(0,k);
        if (i>0)
            idY->s( gap->g( t2->charAt(i-1), k ), k );
        else
            idY->s(0,k);
    }

    return;

    FOR(k,nState) {

        fM->s(0,k);
        idX->s(0,k);
        idY->s(0,k);

        nullM1=nullM2=0.0;

        FOR(n,sAlpha) {

            matchBr1=matchBr2=0.0;

            // match
            if (j>0 && i>0) {
                matchBr1 += hmm->charSbProbL(k,n,t1->charAt(j-1));
                matchBr2 += hmm->charSbProbR(k,n,t2->charAt(i-1));
            }

            // x-gap
            if (j>0) {
                if (n==0) {
                    idX->a(hmm->charBgFreq(k,t1->charAt(j-1)),k);
                    nullM1 += hmm->nullBgFreq(t1->charAt(j-1));
                }
            }

            // y-gap
            if (i>0) {
                if (n==0) {
                    idY->a(hmm->charBgFreq(k,t2->charAt(i-1)),k);
                    nullM2 += hmm->nullBgFreq(t2->charAt(i-1));
                }
            }

            t = hmm->charBgFreq(k,n)*matchBr1*matchBr2/(nullM1*nullM2);

            fM->a( t, k );
        }

        fM->clog( k );

        idX->d( nullM1, k );
        idX->clog( k );

        idY->d( nullM2, k );
        idY->clog( k );
    }

    return;
}

void PhyloMatchScore::fwdSM(int j,int i)
{
    FOR(k,nState) {

        fM->s(0,k);
        idX->s(0,k);
        idY->s(0,k);

        nullM1=nullM2=0.0;

        FOR(n,sAlpha) {

            matchBr1=matchBr2=0.0;

            FOR(m,sAlpha) {
                // match
                if (j>0 && i>0) {
                    matchBr2 += hmm->charSbProbR(k,n,m)*a2->mlCharProbAt(m,i-1,k);
                }

                // y-gap
                if (i>0) {
                    if (n==0) {
                        idY->a(hmm->charBgFreq(k,m)*a2->mlCharProbAt(m,i-1,k),k);
                        nullM2 += hmm->nullBgFreq(m)*a2->mlCharProbAt(m,i-1,k);
                    }
                }
            }

            // match
            if (j>0 && i>0) {
                matchBr1 += hmm->charSbProbL(k,n,t1->charAt(j-1));
            }

            // x-gap
            if (j>0) {
                if (n==0) {
                    idX->a(hmm->charBgFreq(k,t1->charAt(j-1)),k);
                    nullM1 += hmm->nullBgFreq(t1->charAt(j-1));
                }
            }

            t = hmm->charBgFreq(k,n)*matchBr1*matchBr2/(nullM1*nullM2);

            fM->a( t, k );
        }

        fM->clog( k );

        idX->d( nullM1, k );
        idX->clog( k );

        idY->d( nullM2, k );
        idY->clog( k );
    }

    return;
}


void PhyloMatchScore::fwdMS(int j,int i)
{
    FOR(k,nState) {

        fM->s(0,k);
        idX->s(0,k);
        idY->s(0,k);

        nullM1=nullM2=0.0;

        FOR(n,sAlpha) {

            matchBr1=matchBr2=0.0;

            FOR(m,sAlpha) {
                // match
                if (j>0 && i>0) {
                    matchBr1 += hmm->charSbProbL(k,n,m)*a1->mlCharProbAt(m,j-1,k);
                }

                // x-gap
                if (j>0) {
                    if (n==0) {
                        idX->a(hmm->charBgFreq(k,m)*a1->mlCharProbAt(m,j-1,k),k);
                        nullM1 += hmm->nullBgFreq(m)*a1->mlCharProbAt(m,j-1,k);
                    }
                }
            }

            // match
            if (j>0 && i>0) {
                matchBr2 += hmm->charSbProbR(k,n,t2->charAt(i-1));
            }

            // y-gap
            if (i>0) {
                if (n==0) {
                    idY->a(hmm->charBgFreq(k,t2->charAt(i-1)),k);
                    nullM2 += hmm->nullBgFreq(t2->charAt(i-1));
                }
            }

            t = hmm->charBgFreq(k,n)*matchBr1*matchBr2/(nullM1*nullM2);

            fM->a( t, k );
        }

        fM->clog( k );

        idX->d( nullM1, k );
        idX->clog( k );

        idY->d( nullM2, k );
        idY->clog( k );
    }

    return;
}


void PhyloMatchScore::bwdMM(int j,int i)
{
    FOR(k,nState) {

        bM->s(0,k);
        idX->s(0,k);
        idY->s(0,k);

        nullM1=nullM2=0.0;

        FOR(n,sAlpha) {

            matchBr1=matchBr2=0.0;

            FOR(m,sAlpha) {
                // match
                if (j<sl1 && i<sl2) {
                    matchBr1 += hmm->charSbProbL(k,n,m)*a1->mlCharProbAt(m,j,k);
                    matchBr2 += hmm->charSbProbR(k,n,m)*a2->mlCharProbAt(m,i,k);
                }

                // x-gap
                if (j<sl1) {
                    if (n==0) {
                        idX->a(hmm->charBgFreq(k,m)*a1->mlCharProbAt(m,j,k),k);
                        nullM1 += hmm->nullBgFreq(m)*a1->mlCharProbAt(m,j,k);
                    }
                }

                // y-gap
                if (i<sl2) {
                    if (n==0) {
                        idY->a(hmm->charBgFreq(k,m)*a2->mlCharProbAt(m,i,k),k);
                        nullM2 += hmm->nullBgFreq(m)*a2->mlCharProbAt(m,i,k);
                    }
                }
            }

            t = hmm->charBgFreq(k,n)*matchBr1*matchBr2/(nullM1*nullM2);
            bM->a( t, k );
        }

        bM->clog( k );

        idX->d( nullM1, k );
        idX->clog( k );

        idY->d( nullM2, k );
        idY->clog( k );

    }
    return;
}

void PhyloMatchScore::bwdSS(int j,int i)
{
    FOR(k,nState) {
        if (j<sl1 && i<sl2) {
            bM->s( match->g( t1->charAt(j), t2->charAt(i), k ), k );
            flM->s( match->g( t1->charAt(j), t2->charAt(i), k ), k );
        }
        else {
            bM->s( 0, k );
            flM->s( 0, k );
        }
        if (j<sl1)
            idX->s( gap->g( t1->charAt(j), k ), k );
        else
            idX->s(0,k);
        if (i<sl2)
            idY->s( gap->g( t2->charAt(i), k ), k );
        else
            idY->s(0,k);
    }

    return;


    FOR(k,nState) {

        bM->s(0,k);
        idX->s(0,k);
        idY->s(0,k);

        nullM1=nullM2=0.0;

        FOR(n,sAlpha) {

            matchBr1=matchBr2=0.0;

            // match
            if (j<sl1 && i<sl2) {
                matchBr1 += hmm->charSbProbL(k,n,t1->charAt(j));
                matchBr2 += hmm->charSbProbR(k,n,t2->charAt(i));
            }

            // x-gap
            if (j<sl1) {
                if (n==0) {
                    idX->a(hmm->charBgFreq(k,t1->charAt(j)),k);
                    nullM1 += hmm->nullBgFreq(t1->charAt(j));
                }
            }

            // y-gap
            if (i<sl2) {
                if (n==0) {
                    idY->a(hmm->charBgFreq(k,t2->charAt(i)),k);
                    nullM2 += hmm->nullBgFreq(t2->charAt(i));
                }
            }

            t = hmm->charBgFreq(k,n)*matchBr1*matchBr2/(nullM1*nullM2);
            bM->a( t, k );
        }

        bM->clog( k );

        idX->d( nullM1, k );
        idX->clog( k );

        idY->d( nullM2, k );
        idY->clog( k );

    }
    return;
}

void PhyloMatchScore::bwdSM(int j,int i)
{
    FOR(k,nState) {

        bM->s(0,k);
        idX->s(0,k);
        idY->s(0,k);

        nullM1=nullM2=0.0;

        FOR(n,sAlpha) {

            matchBr1=matchBr2=0.0;

            FOR(m,sAlpha) {
                // match
                if (j<sl1 && i<sl2) {
                    matchBr2 += hmm->charSbProbR(k,n,m)*a2->mlCharProbAt(m,i,k);
                }

                // y-gap
                if (i<sl2) {
                    if (n==0) {
                        idY->a(hmm->charBgFreq(k,m)*a2->mlCharProbAt(m,i,k),k);
                        nullM2 += hmm->nullBgFreq(m)*a2->mlCharProbAt(m,i,k);
                    }
                }
            }

            // match
            if (j<sl1 && i<sl2) {
                matchBr1 += hmm->charSbProbL(k,n,t1->charAt(j));
            }

            // x-gap
            if (j<sl1) {
                if (n==0) {
                    idX->a(hmm->charBgFreq(k,t1->charAt(j)),k);
                    nullM1 += hmm->nullBgFreq(t1->charAt(j));
                }
            }

            t = hmm->charBgFreq(k,n)*matchBr1*matchBr2/(nullM1*nullM2);
            bM->a( t, k );
        }

        bM->clog( k );

        idX->d( nullM1, k );
        idX->clog( k );

        idY->d( nullM2, k );
        idY->clog( k );

    }
    return;
}

void PhyloMatchScore::bwdMS(int j,int i)
{
    FOR(k,nState) {

        bM->s(0,k);
        idX->s(0,k);
        idY->s(0,k);

        nullM1=nullM2=0.0;

        FOR(n,sAlpha) {

            matchBr1=matchBr2=0.0;

            FOR(m,sAlpha) {
                // match
                if (j<sl1 && i<sl2) {
                    matchBr1 += hmm->charSbProbL(k,n,m)*a1->mlCharProbAt(m,j,k);
                }

                // x-gap
                if (j<sl1) {
                    if (n==0) {
                        idX->a(hmm->charBgFreq(k,m)*a1->mlCharProbAt(m,j,k),k);
                        nullM1 += hmm->nullBgFreq(m)*a1->mlCharProbAt(m,j,k);
                    }
                }
            }

            // match
            if (j<sl1 && i<sl2) {
                matchBr2 += hmm->charSbProbR(k,n,t2->charAt(i));
            }

            // y-gap
            if (i<sl2) {
                if (n==0) {
                    idY->a(hmm->charBgFreq(k,t2->charAt(i)),k);
                    nullM2 += hmm->nullBgFreq(t2->charAt(i));
                }
            }

            t = hmm->charBgFreq(k,n)*matchBr1*matchBr2/(nullM1*nullM2);
            bM->a( t, k );
        }

        bM->clog( k );

        idX->d( nullM1, k );
        idX->clog( k );

        idY->d( nullM2, k );
        idY->clog( k );

    }
    return;
}

// full probability
void PhyloMatchScore::fullFwdMM(int j,int i)
{

    FOR(k,nState) {

        flM->s(0,k);

        idX->s(0,k);
        idY->s(0,k);

        nullM1 = nullM2 = 0.0;

        FOR(n,sAlpha) {

            matchBr1=matchBr2=0.0;

            FOR(m,sAlpha) {
                // match
                if (j>0 && i>0) {
                    matchBr1 += hmm->charSbProbL(k,n,m)*a1->mlCharProbAtF(m,j-1,k);
                    matchBr2 += hmm->charSbProbR(k,n,m)*a2->mlCharProbAtF(m,i-1,k);
                }

                // x-gap
                if (j>0) {
                    if (n==0) {
                        idX->a(hmm->charBgFreq(k,m)*a1->mlCharProbAtF(m,j-1,k),k);
                        nullM1 += hmm->nullBgFreq(m)*a1->mlCharProbAtF(m,j-1,k); // added ",k"
                    }
                }

                // y-gap
                if (i>0) {
                    if (n==0) {
                        idY->a(hmm->charBgFreq(k,m)*a2->mlCharProbAtF(m,i-1,k),k);
                        nullM2 += hmm->nullBgFreq(m)*a2->mlCharProbAtF(m,i-1,k); // added ",k"
                    }
                }
            }

            flM->a( hmm->charBgFreq(k,n)*matchBr1*matchBr2/(nullM1*nullM2), k );
        }

        if (flM->g(k)!=0)
            flM->clog( k );
        else
            flM->s( -HUGE_VAL, k );

        if (idX->g(k)!=0) {
            idX->s( idX->g( k )/nullM1, k );
            idX->clog( k );
        } else {
            idX->s( -HUGE_VAL, k );
        }

        if (idY->g(k)!=0) {
            idY->s( idY->g( k )/nullM2, k );
            idY->clog( k );
        } else {
            idY->s( -HUGE_VAL, k );
        }
    }

    return;
}

void PhyloMatchScore::fullFwdSS(int j,int i)
{
    FOR(k,nState) {

        flM->s(0,k);

        idX->s(0,k);
        idY->s(0,k);

        nullM1 = nullM2 = 0.0;

        FOR(n,sAlpha) {

            matchBr1=matchBr2=0.0;

            // match
            if (j>0 && i>0) {
                matchBr1 += hmm->charSbProbL(k,n,t1->charAt(j-1));
                matchBr2 += hmm->charSbProbR(k,n,t2->charAt(i-1));
            }

            // x-gap
            if (j>0) {
                if (n==0) {
                    idX->a(hmm->charBgFreq(k,t1->charAt(j-1)),k);
                    nullM1 += hmm->nullBgFreq(t1->charAt(j-1));
                }
            }

            // y-gap
            if (i>0) {
                if (n==0) {
                    idY->a(hmm->charBgFreq(k,t2->charAt(i-1)),k);
                    nullM2 += hmm->nullBgFreq(t2->charAt(i-1));
                }
            }

            flM->a( hmm->charBgFreq(k,n)*matchBr1*matchBr2/(nullM1*nullM2), k );
        }

        if (flM->g(k)!=0)
            flM->clog( k );
        else
            flM->s( -HUGE_VAL, k );

        if (idX->g(k)!=0) {
            idX->s( idX->g( k )/nullM1, k );
            idX->clog( k );
        } else {
            idX->s( -HUGE_VAL, k );
        }

        if (idY->g(k)!=0) {
            idY->s( idY->g( k )/nullM2, k );
            idY->clog( k );
        } else {
            idY->s( -HUGE_VAL, k );
        }
    }

    return;
}

void PhyloMatchScore::fullFwdSM(int j,int i)
{

    FOR(k,nState) {

        flM->s(0,k);

        idX->s(0,k);
        idY->s(0,k);

        nullM1 = nullM2 = 0.0;

        FOR(n,sAlpha) {

            matchBr1=matchBr2=0.0;

            FOR(m,sAlpha) {
                // match
                if (j>0 && i>0) {
                    matchBr2 += hmm->charSbProbR(k,n,m)*a2->mlCharProbAtF(m,i-1,k);
                }

                // y-gap
                if (i>0) {
                    if (n==0) {
                        idY->a(hmm->charBgFreq(k,m)*a2->mlCharProbAtF(m,i-1,k),k);
                        nullM2 += hmm->nullBgFreq(m)*a2->mlCharProbAtF(m,i-1,k); // added ",k"
                    }
                }
            }

            // match
            if (j>0 && i>0) {
                matchBr1 += hmm->charSbProbL(k,n,t1->charAt(j-1));
            }

            // x-gap
            if (j>0) {
                if (n==0) {
                    idX->a(hmm->charBgFreq(k,t1->charAt(j-1)),k);
                    nullM1 += hmm->nullBgFreq(t1->charAt(j-1));
                }
            }

            flM->a( hmm->charBgFreq(k,n)*matchBr1*matchBr2/(nullM1*nullM2), k );
        }

        if (flM->g(k)!=0)
            flM->clog( k );
        else
            flM->s( -HUGE_VAL, k );

        if (idX->g(k)!=0) {
            idX->s( idX->g( k )/nullM1, k );
            idX->clog( k );
        } else {
            idX->s( -HUGE_VAL, k );
        }

        if (idY->g(k)!=0) {
            idY->s( idY->g( k )/nullM2, k );
            idY->clog( k );
        } else {
            idY->s( -HUGE_VAL, k );
        }
    }

    return;
}

void PhyloMatchScore::fullFwdMS(int j,int i)
{

    FOR(k,nState) {

        flM->s(0,k);

        idX->s(0,k);
        idY->s(0,k);

        nullM1 = nullM2 = 0.0;

        FOR(n,sAlpha) {

            matchBr1=matchBr2=0.0;

            FOR(m,sAlpha) {
                // match
                if (j>0 && i>0) {
                    matchBr1 += hmm->charSbProbL(k,n,m)*a1->mlCharProbAtF(m,j-1,k);
                }

                // x-gap
                if (j>0) {
                    if (n==0) {
                        idX->a(hmm->charBgFreq(k,m)*a1->mlCharProbAtF(m,j-1,k),k);
                        nullM1 += hmm->nullBgFreq(m)*a1->mlCharProbAtF(m,j-1,k); // added ",k"
                    }
                }
            }

            // match
            if (j>0 && i>0) {
                matchBr2 += hmm->charSbProbR(k,n,t2->charAt(i-1));
            }

            // y-gap
            if (i>0) {
                if (n==0) {
                    idY->a(hmm->charBgFreq(k,t2->charAt(i-1)),k);
                    nullM2 += hmm->nullBgFreq(t2->charAt(i-1));
                }
            }

            flM->a( hmm->charBgFreq(k,n)*matchBr1*matchBr2/(nullM1*nullM2), k );
        }

        if (flM->g(k)!=0)
            flM->clog( k );
        else
            flM->s( -HUGE_VAL, k );

        if (idX->g(k)!=0) {
            idX->s( idX->g( k )/nullM1, k );
            idX->clog( k );
        } else {
            idX->s( -HUGE_VAL, k );
        }

        if (idY->g(k)!=0) {
            idY->s( idY->g( k )/nullM2, k );
            idY->clog( k );
        } else {
            idY->s( -HUGE_VAL, k );
        }
    }

    return;
}

void PhyloMatchScore::fullBwdMM(int j,int i) {

    FOR(k,nState) {

        flM->s(0,k);

        idX->s(0,k);
        idY->s(0,k);

        nullM1 = nullM2 = 0.0;

        FOR(n,sAlpha) {

            matchBr1=matchBr2=0.0;

            FOR(m,sAlpha) {
                // match
                if (j<sfl1 && i<sfl2) {
                    matchBr1 += hmm->charSbProbL(k,n,m)*a1->mlCharProbAtF(m,j,k);
                    matchBr2 += hmm->charSbProbR(k,n,m)*a2->mlCharProbAtF(m,i,k);
                }

                // x-gap
                if (j<sfl1) {
                    if (n==0) {
                        idX->a(hmm->charBgFreq(k,m)*a1->mlCharProbAtF(m,j,k),k);
                        nullM1 += hmm->nullBgFreq(m)*a1->mlCharProbAtF(m,j,k); // added ",k"
                    }
                }

                // y-gap
                if (i<sfl2) {
                    if (n==0) {
                        idY->a(hmm->charBgFreq(k,m)*a2->mlCharProbAtF(m,i,k),k);
                        nullM2 += hmm->nullBgFreq(m)*a2->mlCharProbAtF(m,i,k); // added ",k"
                    }
                }
            }

            flM->a( hmm->charBgFreq(k,n)*matchBr1*matchBr2/(nullM1*nullM2), k );

        }

        if (flM->g(k)!=0)
            flM->clog( k );
        else
            flM->s( -HUGE_VAL, k );

        if (idX->g(k)!=0) {
            idX->s( idX->g( k )/nullM1, k );
            idX->clog( k );
        } else {
            idX->s( -HUGE_VAL, k );
        }

        if (idY->g(k)!=0) {
            idY->s( idY->g( k )/nullM2, k );
            idY->clog( k );
        } else {
            idY->s( -HUGE_VAL, k );
        }

    }

    return;
}

void PhyloMatchScore::fullBwdSS(int j,int i) {

    FOR(k,nState) {

        flM->s(0,k);

        idX->s(0,k);
        idY->s(0,k);

        nullM1 = nullM2 = 0.0;

        FOR(n,sAlpha) {

            matchBr1=matchBr2=0.0;

            // match
            if (j<sfl1 && i<sfl2) {
                matchBr1 += hmm->charSbProbL(k,n,t1->charAt(j));
                matchBr2 += hmm->charSbProbR(k,n,t2->charAt(i));
            }

            // x-gap
            if (j<sfl1) {
                if (n==0) {
                    idX->a(hmm->charBgFreq(k,t1->charAt(j)),k);
                    nullM1 += hmm->nullBgFreq(t1->charAt(j));
                }
            }

            // y-gap
            if (i<sfl2) {
                if (n==0) {
                    idY->a(hmm->charBgFreq(k,t2->charAt(i)),k);
                    nullM2 += hmm->nullBgFreq(t2->charAt(i));
                }
            }

            flM->a( hmm->charBgFreq(k,n)*matchBr1*matchBr2/(nullM1*nullM2), k );

        }

        if (flM->g(k)!=0)
            flM->clog( k );
        else
            flM->s( -HUGE_VAL, k );

        if (idX->g(k)!=0) {
            idX->s( idX->g( k )/nullM1, k );
            idX->clog( k );
        } else {
            idX->s( -HUGE_VAL, k );
        }

        if (idY->g(k)!=0) {
            idY->s( idY->g( k )/nullM2, k );
            idY->clog( k );
        } else {
            idY->s( -HUGE_VAL, k );
        }

    }

    return;
}

void PhyloMatchScore::fullBwdSM(int j,int i) {

    FOR(k,nState) {

        flM->s(0,k);

        idX->s(0,k);
        idY->s(0,k);

        nullM1 = nullM2 = 0.0;

        FOR(n,sAlpha) {

            matchBr1=matchBr2=0.0;

            FOR(m,sAlpha) {
                // match
                if (j<sfl1 && i<sfl2) {
                    matchBr2 += hmm->charSbProbR(k,n,m)*a2->mlCharProbAtF(m,i,k);
                }

                // y-gap
                if (i<sfl2) {
                    if (n==0) {
                        idY->a(hmm->charBgFreq(k,m)*a2->mlCharProbAtF(m,i,k),k);
                        nullM2 += hmm->nullBgFreq(m)*a2->mlCharProbAtF(m,i,k); // added ",k"
                    }
                }
            }

            // match
            if (j<sfl1 && i<sfl2) {
                matchBr1 += hmm->charSbProbL(k,n,t1->charAt(j));
            }

            // x-gap
            if (j<sfl1) {
                if (n==0) {
                    idX->a(hmm->charBgFreq(k,t1->charAt(j)),k);
                    nullM1 += hmm->nullBgFreq(t1->charAt(j));
                }
            }

            flM->a( hmm->charBgFreq(k,n)*matchBr1*matchBr2/(nullM1*nullM2), k );

        }

        if (flM->g(k)!=0)
            flM->clog( k );
        else
            flM->s( -HUGE_VAL, k );

        if (idX->g(k)!=0) {
            idX->s( idX->g( k )/nullM1, k );
            idX->clog( k );
        } else {
            idX->s( -HUGE_VAL, k );
        }

        if (idY->g(k)!=0) {
            idY->s( idY->g( k )/nullM2, k );
            idY->clog( k );
        } else {
            idY->s( -HUGE_VAL, k );
        }

    }

    return;
}

void PhyloMatchScore::fullBwdMS(int j,int i) {

    FOR(k,nState) {

        flM->s(0,k);

        idX->s(0,k);
        idY->s(0,k);

        nullM1 = nullM2 = 0.0;

        FOR(n,sAlpha) {

            matchBr1=matchBr2=0.0;

            FOR(m,sAlpha) {
                // match
                if (j<sfl1 && i<sfl2) {
                    matchBr1 += hmm->charSbProbL(k,n,m)*a1->mlCharProbAtF(m,j,k);
                }

                // x-gap
                if (j<sfl1) {
                    if (n==0) {
                        idX->a(hmm->charBgFreq(k,m)*a1->mlCharProbAtF(m,j,k),k);
                        nullM1 += hmm->nullBgFreq(m)*a1->mlCharProbAtF(m,j,k); // added ",k"
                    }
                }
            }

            // match
            if (j<sfl1 && i<sfl2) {
                matchBr2 += hmm->charSbProbR(k,n,t2->charAt(i));
            }

            // y-gap
            if (i<sfl2) {
                if (n==0) {
                    idY->a(hmm->charBgFreq(k,t2->charAt(i)),k);
                    nullM2 += hmm->nullBgFreq(t2->charAt(i));
                }
            }

            flM->a( hmm->charBgFreq(k,n)*matchBr1*matchBr2/(nullM1*nullM2), k );

        }

        if (flM->g(k)!=0)
            flM->clog( k );
        else
            flM->s( -HUGE_VAL, k );

        if (idX->g(k)!=0) {
            idX->s( idX->g( k )/nullM1, k );
            idX->clog( k );
        } else {
            idX->s( -HUGE_VAL, k );
        }

        if (idY->g(k)!=0) {
            idY->s( idY->g( k )/nullM2, k );
            idY->clog( k );
        } else {
            idY->s( -HUGE_VAL, k );
        }

    }

    return;
}


void PhyloMatchScore::logFwdMM(int j,int i)
{
    FOR(k,nState) {

        fM->s(small,k);
        idX->s(small,k);
        idY->s(small,k);

        nullM1=nullM2=small;

        FOR(n,sAlpha) {

            matchBr1=matchBr2=small;

            FOR(m,sAlpha) {
                // match
                if (j>0 && i>0) {
                    matchBr1 = sumLogs(matchBr1,hmm->logCharSbProbL(k,n,m)+a1->mlCharProbAt(m,j-1,k));
                    matchBr2 = sumLogs(matchBr2,hmm->logCharSbProbR(k,n,m)+a2->mlCharProbAt(m,i-1,k));
                }

                // x-gap
                if (j>0) {
                    if (n==0) {
                        idX->alog( hmm->logCharBgFreq(k,m)+a1->mlCharProbAt(m,j-1,k), k );
                        nullM1 = sumLogs(nullM1,hmm->logNullBgFreq(m)+a1->mlCharProbAt(m,j-1,k));
                    }
                }

                // y-gap
                if (i>0) {
                    if (n==0) {
                        idY->alog( hmm->logCharBgFreq(k,m)+a2->mlCharProbAt(m,i-1,k), k );
                        nullM2 = sumLogs(nullM2,hmm->logNullBgFreq(m)+a2->mlCharProbAt(m,i-1,k));
                    }
                }
            }

            t = hmm->logCharBgFreq(k,n)+matchBr1+matchBr2-(nullM1+nullM2);
            fM->alog( t, k );

        }

        idX->a( -1*nullM1, k );
        idY->a( -1*nullM2, k );

    }
    return;
}

void PhyloMatchScore::logFwdSS(int j,int i)
{
    FOR(k,nState) {
        if (j>0 && i>0) {
            fM->s( match->g( t1->charAt(j-1), t2->charAt(i-1), k ), k );
            flM->s( fM->g( k ), k );
        }
        else {
            fM->s( small, k );
            flM->s( small, k );
        }
        if (j>0)
            idX->s( gap->g( t1->charAt(j-1), k ), k );
        else
            idX->s(small,k);
        if (i>0)
            idY->s( gap->g( t2->charAt(i-1), k ), k );
        else
            idY->s(small,k);
    }

    return;

    FOR(k,nState) {

        fM->s(small,k);
        idX->s(small,k);
        idY->s(small,k);

        nullM1=nullM2=small;

        FOR(n,sAlpha) {

            matchBr1=matchBr2=small;

            // match
            if (j>0 && i>0) {
                matchBr1 = sumLogs(matchBr1,hmm->logCharSbProbL(k,n,t1->charAt(j-1)));
                matchBr2 = sumLogs(matchBr2,hmm->logCharSbProbR(k,n,t2->charAt(i-1)));
            }

            // x-gap
            if (j>0) {
                if (n==0) {
                    idX->alog( hmm->logCharBgFreq(k,t1->charAt(j-1)), k );
                    nullM1 = sumLogs(nullM1,hmm->logNullBgFreq(t1->charAt(j-1)));
                }
            }

            // y-gap
            if (i>0) {
                if (n==0) {
                    idY->alog( hmm->logCharBgFreq(k,t2->charAt(i-1)), k );
                    nullM2 = sumLogs(nullM2,hmm->logNullBgFreq(t2->charAt(i-1)));
                }
            }

            t = hmm->logCharBgFreq(k,n)+matchBr1+matchBr2-(nullM1+nullM2);
            fM->alog( t, k );

        }

        idX->a( -1*nullM1, k );
        idY->a( -1*nullM2, k );

    }
    return;
}

void PhyloMatchScore::logFwdSM(int j,int i)
{
    FOR(k,nState) {

        fM->s(small,k);
        idX->s(small,k);
        idY->s(small,k);

        nullM1=nullM2=small;

        FOR(n,sAlpha) {

            matchBr1=matchBr2=small;

            FOR(m,sAlpha) {
                // match
                if (j>0 && i>0) {
                    matchBr2 = sumLogs(matchBr2,hmm->logCharSbProbR(k,n,m)+a2->mlCharProbAt(m,i-1,k));
                }

                // y-gap
                if (i>0) {
                    if (n==0) {
                        idY->alog( hmm->logCharBgFreq(k,m)+a2->mlCharProbAt(m,i-1,k), k );
                        nullM2 = sumLogs(nullM2,hmm->logNullBgFreq(m)+a2->mlCharProbAt(m,i-1,k));
                    }
                }
            }

            // match
            if (j>0 && i>0) {
                matchBr1 = sumLogs(matchBr1,hmm->logCharSbProbL(k,n,t1->charAt(j-1)));
            }

            // x-gap
            if (j>0) {
                if (n==0) {
                    idX->alog( hmm->logCharBgFreq(k,t1->charAt(j-1)), k );
                    nullM1 = sumLogs(nullM1,hmm->logNullBgFreq(t1->charAt(j-1)));
                }
            }

            t = hmm->logCharBgFreq(k,n)+matchBr1+matchBr2-(nullM1+nullM2);
            fM->alog( t, k );

        }

        idX->a( -1*nullM1, k );
        idY->a( -1*nullM2, k );

    }
    return;
}


void PhyloMatchScore::logFwdMS(int j,int i)
{
    FOR(k,nState) {

        fM->s(small,k);
        idX->s(small,k);
        idY->s(small,k);

        nullM1=nullM2=small;

        FOR(n,sAlpha) {

            matchBr1=matchBr2=small;

            FOR(m,sAlpha) {
                // match
                if (j>0 && i>0) {
                    matchBr1 = sumLogs(matchBr1,hmm->logCharSbProbL(k,n,m)+a1->mlCharProbAt(m,j-1,k));
                }

                // x-gap
                if (j>0) {
                    if (n==0) {
                        idX->alog( hmm->logCharBgFreq(k,m)+a1->mlCharProbAt(m,j-1,k), k );
                        nullM1 = sumLogs(nullM1,hmm->logNullBgFreq(m)+a1->mlCharProbAt(m,j-1,k));
                    }
                }
            }

            // match
            if (j>0 && i>0) {
                matchBr2 = sumLogs(matchBr2,hmm->logCharSbProbR(k,n,t2->charAt(i-1)));
            }

            // y-gap
            if (i>0) {
                if (n==0) {
                    idY->alog( hmm->logCharBgFreq(k,t2->charAt(i-1)), k );
                    nullM2 = sumLogs(nullM2,hmm->logNullBgFreq(t2->charAt(i-1)));
                }
            }

            t = hmm->logCharBgFreq(k,n)+matchBr1+matchBr2-(nullM1+nullM2);
            fM->alog( t, k );

        }

        idX->a( -1*nullM1, k );
        idY->a( -1*nullM2, k );

    }
    return;
}


void PhyloMatchScore::logBwdMM(int j,int i)
{
    FOR(k,nState) {

        bM->s(small,k);
        idX->s(small,k);
        idY->s(small,k);

        nullM1=nullM2=small;

        FOR(n,sAlpha) {

            matchBr1=matchBr2=small;

            FOR(m,sAlpha) {
                // match
                if (j<sl1 && i<sl2) {
                    matchBr1 = sumLogs(matchBr1,hmm->logCharSbProbL(k,n,m)+a1->mlCharProbAt(m,j,k));
                    matchBr2 = sumLogs(matchBr2,hmm->logCharSbProbR(k,n,m)+a2->mlCharProbAt(m,i,k));
                }

                // x-gap
                if (j<sl1) {
                    if (n==0) {
                        idX->alog( hmm->logCharBgFreq(k,m)+a1->mlCharProbAt(m,j,k), k );
                        nullM1 = sumLogs(nullM1,hmm->logNullBgFreq(m)+a1->mlCharProbAt(m,j,k));
                    }
                }

                // y-gap
                if (i<sl2) {
                    if (n==0) {
                        idY->alog( hmm->logCharBgFreq(k,m)+a2->mlCharProbAt(m,i,k), k );
                        nullM2 = sumLogs(nullM2,hmm->logNullBgFreq(m)+a2->mlCharProbAt(m,i,k));
                    }
                }
            }

            t = hmm->logCharBgFreq(k,n)+matchBr1+matchBr2-(nullM1+nullM2);
            bM->alog( t, k );

        }

        idX->a( -1*nullM1, k );
        idY->a( -1*nullM2, k );

    }
    return;
}

void PhyloMatchScore::logBwdSS(int j,int i)
{
    FOR(k,nState) {
        if (j<sl1 && i<sl2) {
            bM->s( match->g( t1->charAt(j), t2->charAt(i), k ), k );
            flM->s( bM->g( k ), k );
        }
        else {
            bM->s( small, k );
            flM->s( small, k );
        }
        if (j<sl1)
            idX->s( gap->g( t1->charAt(j), k ), k );
        else
            idX->s(small,k);
        if (i<sl2)
            idY->s( gap->g( t2->charAt(i), k ), k );
        else
            idY->s(small,k);
    }

    return;

    FOR(k,nState) {

        bM->s(small,k);
        idX->s(small,k);
        idY->s(small,k);

        nullM1=nullM2=small;

        FOR(n,sAlpha) {

            matchBr1=matchBr2=small;

            // match
            if (j<sl1 && i<sl2) {
                matchBr1 = sumLogs(matchBr1,hmm->logCharSbProbL(k,n,t1->charAt(j)));
                matchBr2 = sumLogs(matchBr2,hmm->logCharSbProbR(k,n,t2->charAt(i)));
            }

            // x-gap
            if (j<sl1) {
                if (n==0) {
                    idX->alog( hmm->logCharBgFreq(k,t1->charAt(j)), k );
                    nullM1 = sumLogs(nullM1,hmm->logNullBgFreq(t1->charAt(j)));
                }
            }

            // y-gap
            if (i<sl2) {
                if (n==0) {
                    idY->alog( hmm->logCharBgFreq(k,t2->charAt(i)), k );
                    nullM2 = sumLogs(nullM2,hmm->logNullBgFreq(t2->charAt(i)));
                }
            }

            t = hmm->logCharBgFreq(k,n)+matchBr1+matchBr2-(nullM1+nullM2);
            bM->alog( t, k );

        }

        idX->a( -1*nullM1, k );
        idY->a( -1*nullM2, k );

    }
    return;
}
void PhyloMatchScore::logBwdSM(int j,int i)
{
    FOR(k,nState) {

        bM->s(small,k);

        idX->s(small,k);
        idY->s(small,k);

        nullM1=nullM2=small;

        FOR(n,sAlpha) {

            matchBr1=matchBr2=small;

            FOR(m,sAlpha) {
                // match
                if (j<sl1 && i<sl2) {
                    matchBr2 = sumLogs(matchBr2,hmm->logCharSbProbR(k,n,m)+a2->mlCharProbAt(m,i,k));
                }

                // y-gap
                if (i<sl2) {
                    if (n==0) {
                        idY->alog( hmm->logCharBgFreq(k,m)+a2->mlCharProbAt(m,i,k), k );
                        nullM2 = sumLogs(nullM2,hmm->logNullBgFreq(m)+a2->mlCharProbAt(m,i,k));
                    }
                }
            }

            // match
            if (j<sl1 && i<sl2) {
                matchBr1 = sumLogs(matchBr1,hmm->logCharSbProbL(k,n,t1->charAt(j)));
            }

            // x-gap
            if (j<sl1) {
                if (n==0) {
                    idX->alog( hmm->logCharBgFreq(k,t1->charAt(j)), k );
                    nullM1 = sumLogs(nullM1,hmm->logNullBgFreq(t1->charAt(j)));
                }
            }

            t = hmm->logCharBgFreq(k,n)+matchBr1+matchBr2-(nullM1+nullM2);
            bM->alog( t, k );

        }

        idX->a( -1*nullM1, k );
        idY->a( -1*nullM2, k );

    }
    return;
}
void PhyloMatchScore::logBwdMS(int j,int i)
{
    FOR(k,nState) {

        bM->s(small,k);

        idX->s(small,k);
        idY->s(small,k);

        nullM1=nullM2=small;

        FOR(n,sAlpha) {

            matchBr1=matchBr2=small;

            FOR(m,sAlpha) {
                // match
                if (j<sl1 && i<sl2) {
                    matchBr1 = sumLogs(matchBr1,hmm->logCharSbProbL(k,n,m)+a1->mlCharProbAt(m,j,k));
                }

                // x-gap
                if (j<sl1) {
                    if (n==0) {
                        idX->alog( hmm->logCharBgFreq(k,m)+a1->mlCharProbAt(m,j,k), k );
                        nullM1 = sumLogs(nullM1,hmm->logNullBgFreq(m)+a1->mlCharProbAt(m,j,k));
                    }
                }
            }

            // match
            if (j<sl1 && i<sl2) {
                matchBr2 = sumLogs(matchBr2,hmm->logCharSbProbR(k,n,t2->charAt(i)));
            }

            // y-gap
            if (i<sl2) {
                if (n==0) {
                    idY->alog( hmm->logCharBgFreq(k,t2->charAt(i)), k );
                    nullM2 = sumLogs(nullM2,hmm->logNullBgFreq(t2->charAt(i)));
                }
            }

            t = hmm->logCharBgFreq(k,n)+matchBr1+matchBr2-(nullM1+nullM2);
            bM->alog( t, k );

        }

        idX->a( -1*nullM1, k );
        idY->a( -1*nullM2, k );

    }
    return;
}

void PhyloMatchScore::logFullFwdMM(int j,int i)
{
    FOR(k,nState) {

        flM->s(small,k);

        idX->s(small,k);
        idY->s(small,k);

        nullM1 = nullM2 = small;

        FOR(n,sAlpha) {

            matchBr1=matchBr2=small;

            FOR(m,sAlpha) {
                // match
                if (j>0 && i>0) {
                    matchBr1 = sumLogs(matchBr1,hmm->logCharSbProbL(k,n,m)+a1->mlCharProbAtF(m,j-1,k));
                    matchBr2 = sumLogs(matchBr2,hmm->logCharSbProbR(k,n,m)+a2->mlCharProbAtF(m,i-1,k));
                }

                // x-gap
                if (j>0) {
                    if (n==0) {
                        idX->alog( hmm->logCharBgFreq(k,m)+a1->mlCharProbAtF(m,j-1,k), k );
                        nullM1 = sumLogs(nullM1,hmm->logNullBgFreq(m)+a1->mlCharProbAtF(m,j-1,k));
                    }
                }

                // y-gap
                if (i>0) {
                    if (n==0) {
                        idY->alog( hmm->logCharBgFreq(k,m)+a2->mlCharProbAtF(m,i-1,k), k );
                        nullM2 = sumLogs(nullM2,hmm->logNullBgFreq(m)+a2->mlCharProbAtF(m,i-1,k));
                    }
                }
            }

            flM->alog( hmm->logCharBgFreq(k,n)+matchBr1+matchBr2-(nullM1+nullM2), k );
        }

        idX->a( -1*nullM1, k );
        idY->a( -1*nullM2, k );

    }

    return;
}

void PhyloMatchScore::logFullFwdSS(int j,int i)
{
    FOR(k,nState) {

        flM->s(small,k);

        idX->s(small,k);
        idY->s(small,k);

        nullM1 = nullM2 = small;

        FOR(n,sAlpha) {

            matchBr1=matchBr2=small;

            // match
            if (j>0 && i>0) {
                matchBr1 = sumLogs(matchBr1,hmm->logCharSbProbL(k,n,t1->charAt(j-1)));
                matchBr2 = sumLogs(matchBr2,hmm->logCharSbProbR(k,n,t2->charAt(i-1)));
            }

            // x-gap
            if (j>0) {
                if (n==0) {
                    idX->alog( hmm->logCharBgFreq(k,t1->charAt(j-1)), k );
                    nullM1 = sumLogs(nullM1,hmm->logNullBgFreq(t1->charAt(j-1)));
                }
            }

            // y-gap
            if (i>0) {
                if (n==0) {
                    idY->alog( hmm->logCharBgFreq(k,t2->charAt(i-1)), k );
                    nullM2 = sumLogs(nullM2,hmm->logNullBgFreq(t2->charAt(i-1)));
                }
            }

            flM->alog( hmm->logCharBgFreq(k,n)+matchBr1+matchBr2-(nullM1+nullM2), k );
        }

        idX->a( -1*nullM1, k );
        idY->a( -1*nullM2, k );

    }

    return;
}

void PhyloMatchScore::logFullFwdSM(int j,int i)
{
    FOR(k,nState) {

        flM->s(small,k);

        idX->s(small,k);
        idY->s(small,k);

        nullM1 = nullM2 = small;

        FOR(n,sAlpha) {

            matchBr1=matchBr2=small;


            FOR(m,sAlpha) {
                // match
                if (j>0 && i>0) {
                    matchBr2 = sumLogs(matchBr2,hmm->logCharSbProbR(k,n,m)+a2->mlCharProbAtF(m,i-1,k));
                }

                // y-gap
                if (i>0) {
                    if (n==0) {
                        idY->alog( hmm->logCharBgFreq(k,m)+a2->mlCharProbAtF(m,i-1,k), k );
                        nullM2 = sumLogs(nullM2,hmm->logNullBgFreq(m)+a2->mlCharProbAtF(m,i-1,k));
                    }
                }
            }

            // match
            if (j>0 && i>0) {
                matchBr1 = sumLogs(matchBr1,hmm->logCharSbProbL(k,n,t1->charAt(j-1)));
            }

            // x-gap
            if (j>0) {
                if (n==0) {
                    idX->alog( hmm->logCharBgFreq(k,t1->charAt(j-1)), k );
                    nullM1 = sumLogs(nullM1,hmm->logNullBgFreq(t1->charAt(j-1)));
                }
            }

            flM->alog( hmm->logCharBgFreq(k,n)+matchBr1+matchBr2-(nullM1+nullM2), k );
        }

        idX->a( -1*nullM1, k );
        idY->a( -1*nullM2, k );

    }

    return;
}

void PhyloMatchScore::logFullFwdMS(int j,int i)
{
    FOR(k,nState) {

        flM->s(small,k);

        idX->s(small,k);
        idY->s(small,k);

        nullM1 = nullM2 = small;

        FOR(n,sAlpha) {

            matchBr1=matchBr2=small;


            FOR(m,sAlpha) {
                // match
                if (j>0 && i>0) {
                    matchBr1 = sumLogs(matchBr1,hmm->logCharSbProbL(k,n,m)+a1->mlCharProbAtF(m,j-1,k));
                }

                // x-gap
                if (j>0) {
                    if (n==0) {
                        idX->alog( hmm->logCharBgFreq(k,m)+a1->mlCharProbAtF(m,j-1,k), k );
                        nullM1 = sumLogs(nullM1,hmm->logNullBgFreq(m)+a1->mlCharProbAtF(m,j-1,k));
                    }
                }
            }

            // match
            if (j>0 && i>0) {
                matchBr2 = sumLogs(matchBr2,hmm->logCharSbProbR(k,n,t2->charAt(i-1)));
            }

            // y-gap
            if (i>0) {
                if (n==0) {
                    idY->alog( hmm->logCharBgFreq(k,t2->charAt(i-1)), k );
                    nullM2 = sumLogs(nullM2,hmm->logNullBgFreq(t2->charAt(i-1)));
                }
            }

            flM->alog( hmm->logCharBgFreq(k,n)+matchBr1+matchBr2-(nullM1+nullM2), k );
        }

        idX->a( -1*nullM1, k );
        idY->a( -1*nullM2, k );

    }

    return;
}


void PhyloMatchScore::logFullBwdMM(int j,int i)
{
    FOR(k,nState) {

        flM->s(small,k);

        idX->s(small,k);
        idY->s(small,k);

        nullM1 = nullM2 = small;

        FOR(n,sAlpha) {

            matchBr1=matchBr2=small;

            FOR(m,sAlpha) {
                // match
                if (j<sfl1 && i<sfl2) {
                    matchBr1 = sumLogs(matchBr1,hmm->logCharSbProbL(k,n,m)+a1->mlCharProbAtF(m,j,k));
                    matchBr2 = sumLogs(matchBr2,hmm->logCharSbProbR(k,n,m)+a2->mlCharProbAtF(m,i,k));
                }

                // x-gap
                if (j<sfl1) {
                    if (n==0) {
                        idX->alog( hmm->logCharBgFreq(k,m)+a1->mlCharProbAtF(m,j,k), k );
                        nullM1 = sumLogs(nullM1,hmm->logNullBgFreq(m)+a1->mlCharProbAtF(m,j,k));
                    }
                }

                // y-gap
                if (i<sfl2) {
                    if (n==0) {
                        idY->alog( hmm->logCharBgFreq(k,m)+a2->mlCharProbAtF(m,i,k), k );
                        nullM2 = sumLogs(nullM2,hmm->logNullBgFreq(m)+a2->mlCharProbAtF(m,i,k));
                    }
                }
            }

            flM->alog( hmm->logCharBgFreq(k,n)+matchBr1+matchBr2-(nullM1+nullM2), k );

        }

        idX->a( -1*nullM1, k );
        idY->a( -1*nullM2, k );

    }

    return;
}

void PhyloMatchScore::logFullBwdSS(int j,int i)
{
    FOR(k,nState) {

        flM->s(small,k);

        idX->s(small,k);
        idY->s(small,k);

        nullM1 = nullM2 = small;

        FOR(n,sAlpha) {

            matchBr1=matchBr2=small;

            // match
            if (j<sfl1 && i<sfl2) {
                matchBr1 = sumLogs(matchBr1,hmm->logCharSbProbL(k,n,t1->charAt(j)));
                matchBr2 = sumLogs(matchBr2,hmm->logCharSbProbR(k,n,t2->charAt(i)));
            }

            // x-gap
            if (j<sfl1) {
                if (n==0) {
                    idX->alog( hmm->logCharBgFreq(k,t1->charAt(j)), k );
                    nullM1 = sumLogs(nullM1,hmm->logNullBgFreq(t1->charAt(j)));
                }
            }

            // y-gap
            if (i<sfl2) {
                if (n==0) {
                    idY->alog( hmm->logCharBgFreq(k,t2->charAt(i)), k );
                    nullM2 = sumLogs(nullM2,hmm->logNullBgFreq(t2->charAt(i)));
                }
            }

            flM->alog( hmm->logCharBgFreq(k,n)+matchBr1+matchBr2-(nullM1+nullM2), k );

        }

        idX->a( -1*nullM1, k );
        idY->a( -1*nullM2, k );

        /*		cout<<flM->g(k)<<" "<<idX->g(k)<<" "<<idY->g(k)<<endl;*/
    }

    return;
}

void PhyloMatchScore::logFullBwdSM(int j,int i)
{
    FOR(k,nState) {

        flM->s(small,k);

        idX->s(small,k);
        idY->s(small,k);

        nullM1 = nullM2 = small;

        FOR(n,sAlpha) {

            matchBr1=matchBr2=small;


            FOR(m,sAlpha) {
                // match
                if (j<sfl1 && i<sfl2) {
                    matchBr2 = sumLogs(matchBr2,hmm->logCharSbProbR(k,n,m)+a2->mlCharProbAtF(m,i,k));
                }

                // y-gap
                if (i<sfl2) {
                    if (n==0) {
                        idY->alog( hmm->logCharBgFreq(k,m)+a2->mlCharProbAtF(m,i,k), k );
                        nullM2 = sumLogs(nullM2,hmm->logNullBgFreq(m)+a2->mlCharProbAtF(m,i,k));
                    }
                }
            }

            // match
            if (j<sfl1 && i<sfl2) {
                matchBr1 = sumLogs(matchBr1,hmm->logCharSbProbL(k,n,t1->charAt(j)));
            }

            // x-gap
            if (j<sfl1) {
                if (n==0) {
                    idX->alog( hmm->logCharBgFreq(k,t1->charAt(j)), k );
                    nullM1 = sumLogs(nullM1,hmm->logNullBgFreq(t1->charAt(j)));
                }
            }

            flM->alog( hmm->logCharBgFreq(k,n)+matchBr1+matchBr2-(nullM1+nullM2), k );

        }

        idX->a( -1*nullM1, k );
        idY->a( -1*nullM2, k );

    }

    return;
}

void PhyloMatchScore::logFullBwdMS(int j,int i)
{
    FOR(k,nState) {

        flM->s(small,k);

        idX->s(small,k);
        idY->s(small,k);

        nullM1 = nullM2 = small;

        FOR(n,sAlpha) {

            matchBr1=matchBr2=small;

            FOR(m,sAlpha) {
                // match
                if (j<sfl1 && i<sfl2) {
                    matchBr1 = sumLogs(matchBr1,hmm->logCharSbProbL(k,n,m)+a1->mlCharProbAtF(m,j,k));
                }

                // x-gap
                if (j<sfl1) {
                    if (n==0) {
                        idX->alog( hmm->logCharBgFreq(k,m)+a1->mlCharProbAtF(m,j,k), k );
                        nullM1 = sumLogs(nullM1,hmm->logNullBgFreq(m)+a1->mlCharProbAtF(m,j,k));
                    }
                }
            }

            // match
            if (j<sfl1 && i<sfl2) {
                matchBr2 = sumLogs(matchBr2,hmm->logCharSbProbR(k,n,t2->charAt(i)));
            }

            // y-gap
            if (i<sfl2) {
                if (n==0) {
                    idY->alog( hmm->logCharBgFreq(k,t2->charAt(i)), k );
                    nullM2 = sumLogs(nullM2,hmm->logNullBgFreq(t2->charAt(i)));
                }
            }

            flM->alog( hmm->logCharBgFreq(k,n)+matchBr1+matchBr2-(nullM1+nullM2), k );

        }

        idX->a( -1*nullM1, k );
        idY->a( -1*nullM2, k );

    }

    return;
}



void PhyloMatchScore::computeSSMatrix()
{
    int fas = hmm->getFullASize();

    match = new DbMatrix(fas,fas,nState);
    gap = new DbMatrix(fas,nState);
    if (LOGVALUES) {
        match->initialise(small);
        gap->initialise(small);
    }
    else {
        match->initialise(0);
        gap->initialise(0);
    }

    int i,j;

	if (LOGVALUES) {
		FOR(k,nState) {
			FOR(i,fas) {
				FOR(j,fas) {

					nullM1 = hmm->logNullBgFreq(j);
					nullM2 = hmm->logNullBgFreq(i);

					t = small;

					FOR(n,sAlpha) {
						matchBr1 = hmm->logCharSbProbL(k,n,j);
						matchBr2 = hmm->logCharSbProbR(k,n,i);

                        t = sumLogs(t, hmm->logCharBgFreq(k,n)+matchBr1+matchBr2-(nullM1+nullM2));

                    }

                    match->s( t, i, j, k );
                }

				gap->s( hmm->logCharBgFreq(k,i)-hmm->logNullBgFreq(i), i, k );

            }
        }
    }
    else {

        FOR(k,nState) {
            FOR(i,fas) {
                FOR(j,fas) {

					nullM1 = hmm->nullBgFreq(j);
					nullM2 = hmm->nullBgFreq(i);

					t = 0.0;

                    FOR(n,sAlpha) {

						matchBr1 = hmm->charSbProbL(k,n,j);
						matchBr2 = hmm->charSbProbR(k,n,i);

						t += hmm->charBgFreq(k,n)*matchBr1*matchBr2/(nullM1*nullM2);

                    }

                    match->s( log(t), i, j, k );

                }

				gap->s( hmm->charBgFreq(k,i)/hmm->nullBgFreq(i), i, k );
                gap->clog( i, k );

            }
        }
    }

    if (0){
      cout<<fas<<endl;
        cout<<"match"<<endl;
        match->print();
        cout<<endl;
        cout<<"gap"<<endl;
        gap->print();
        cout<<endl;
    }
}
