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

#include <cstdlib>
#include <iostream>
#include <cmath>
#include "config.h"
#include "characterprobability.h"

using namespace std;

CharacterProbability::~CharacterProbability()
{
}

CharacterProbability::CharacterProbability(Sequence* sq1,Sequence* sq2)
{
    nState = hmm->getNStates();
    sAlpha = hmm->getASize();

	sSite = new Site();

	cSite = new Site();
	cSite->index(1);
	cSite->prev();

	cSite->index(0);
	cSite->next();

    li=0;
    ri=0;

    small = -HUGE_VAL;

    skipMatch = -1;

    if (LOGVALUES) {
        if (sq1->isTerminal() && sq2->isTerminal()) {
            t1 = static_cast<TerminalSequence*>(sq1);
            t2 = static_cast<TerminalSequence*>(sq2);
            logScoresSS();
        }
        else if (sq1->isTerminal() && !sq2->isTerminal()) {
            t1 = static_cast<TerminalSequence*>(sq1);
            a2 = static_cast<AncestralSequence*>(sq2);
            logScoresSM();
        }
        else if (!sq1->isTerminal() && sq2->isTerminal()) {
            a1 = static_cast<AncestralSequence*>(sq1);
            t2 = static_cast<TerminalSequence*>(sq2);
            logScoresMS();
        }
        else if (!sq1->isTerminal() && !sq2->isTerminal()) {
            a1 = static_cast<AncestralSequence*>(sq1);
            a2 = static_cast<AncestralSequence*>(sq2);
            logScoresMM();
        }
        else {
            cout<<"CharacterProbability(Sequence* sq1,Sequence* sq2)"<<endl;
            exit(-1);
        }
    }
    else {
        if (sq1->isTerminal() && sq2->isTerminal()) {
            t1 = static_cast<TerminalSequence*>(sq1);
            t2 = static_cast<TerminalSequence*>(sq2);
            scoresSS();
        }
        else if (sq1->isTerminal() && !sq2->isTerminal()) {
            t1 = static_cast<TerminalSequence*>(sq1);
            a2 = static_cast<AncestralSequence*>(sq2);
            scoresSM();
        }
        else if (!sq1->isTerminal() && sq2->isTerminal()) {
            a1 = static_cast<AncestralSequence*>(sq1);
            t2 = static_cast<TerminalSequence*>(sq2);
            scoresMS();
        }
        else if (!sq1->isTerminal() && !sq2->isTerminal()) {
            a1 = static_cast<AncestralSequence*>(sq1);
            a2 = static_cast<AncestralSequence*>(sq2);
            scoresMM();
        }
        else {
            cout<<"CharacterProbability(Sequence* sq1,Sequence* sq2)"<<endl;
            exit(-1);
        }
    }


    // Compute a kind of fwd/bwd probability
    // Save fwd values in the list and fix the probabilities on the way back
    //
    DbMatrix* score = new DbMatrix(nState,"score");
    DbMatrix* vec1 = new DbMatrix(nState,"vec1");
    DbMatrix* vec2 = new DbMatrix(nState,"vec2");

    DbMatrix* prev;
    DbMatrix* curr;
    DbMatrix* temp;

    prev = vec1;
    curr = vec2;

    double t;

	cSite->index(0);
	cSite->next();

	pSite = new Site();
	pSite->index(0);

    int moveFrom;
    int moveTo = cSite->currMatchState();

    int lastMove;

    FOR(k,nState) {
        if (moveTo==0 || moveTo==3 || moveTo==9){
            t = hmm->structBgFreq(k)+ hmm->probWX(k);
            pSite->stateProb( t, k );
            prev->s( t, k );
        } else if (moveTo==1 || moveTo==7 || moveTo==13){
            t = hmm->structBgFreq(k)+ hmm->probWY(k);
            pSite->stateProb( t, k );
            prev->s( t, k );
        } else if (moveTo==2 || moveTo==5 || moveTo==8 || moveTo==11 || moveTo==14){
            t = hmm->structBgFreq(k)+ hmm->probWM(k);
            pSite->stateProb( t, k );
            prev->s( t, k );
        }
    }


    for (;cSite->index()!=1; cSite->next(),pSite->next()) {

        moveFrom = pSite->currMatchState();
        moveTo = cSite->currMatchState();

        if (moveFrom<0||moveFrom>14)
            moveFrom = 2;

        FOR(l,nState) {

            double score;

            if (LOGVALUES) {

                score=small;
                FOR(j,sAlpha) {
                    score = sumLogs(score,hmm->logCharBgFreq(l,j)+cSite->mlCharProb(l,j));
                }

            } else {

                score=0;
                FOR(j,sAlpha) {
                    score += hmm->charBgFreq(l,j)*cSite->mlCharProb(l,j);
                }
                score = log(score);
            }

            if(score<-100)
                score = -100;

            if ((FOREVER || FOREVER_OLD) && moveTo>2)
            {
                if (LOGVALUES)
                    score = 0;
                else
                    score = 1;
            }
            sum=-HUGE_VAL;

            FOR(k,nState) {

                if (moveFrom==0 || moveFrom==3 || moveFrom==5 || moveFrom==9 || moveFrom==11) {

                    if (moveTo==0 || moveTo==3 || moveTo==5 || moveTo==9 || moveTo==11) {
                        sum = sumLogs(sum,prev->g(k) + hmm->probXX(k,l) + score);
                    } else if (moveTo==1 || moveTo==7 || moveTo==8 || moveTo==13 || moveTo==14) {
                        sum = sumLogs(sum,prev->g(k) + hmm->probXY(k,l) + score);
                    } else if (moveTo==2) {
                        sum = sumLogs(sum,prev->g(k) + hmm->probXM(k,l) + score);
                    } else {
                        cout<<"CharacterProbability: impossible pointer "<<moveFrom<<" "<<moveTo<<endl;
                    }

                } else if (moveFrom==1 || moveFrom==7 || moveFrom==8 || moveFrom==13 || moveFrom==14) {

                    if (moveTo==0 || moveTo==3 || moveTo==5 || moveTo==9 || moveTo==11) {
                        sum = sumLogs(sum,prev->g(k) + hmm->probYX(k,l) + score);
                    } else if (moveTo==1 || moveTo==7 || moveTo==8 || moveTo==13 || moveTo==14) {
                        sum = sumLogs(sum,prev->g(k) + hmm->probYY(k,l) + score);
                    } else if (moveTo==2) {
                        sum = sumLogs(sum,prev->g(k) + hmm->probYM(k,l) + score);
                    } else {
                        cout<<"CharacterProbability: impossible pointer "<<moveFrom<<" "<<moveTo<<endl;
                    }

                } else if (moveFrom==2) {

                    if (moveTo==0 || moveTo==3 || moveTo==5 || moveTo==9 || moveTo==11) {
                        sum = sumLogs(sum,prev->g(k) + hmm->probMX(k,l) + score);
                    } else if (moveTo==1 || moveTo==7 || moveTo==8 || moveTo==13 || moveTo==14) {
                        sum = sumLogs(sum,prev->g(k) + hmm->probMY(k,l) + score);
                    } else if (moveTo==2) {
                        sum = sumLogs(sum,prev->g(k) + hmm->probMM(k,l) + score);
                    } else {
                        cout<<"CharacterProbability: impossible pointer "<<moveFrom<<" "<<moveTo<<endl;
                    }
                }
            }
            cSite->stateProb(sum,l);
            curr->s(sum,l);

        }
        temp = prev;
        prev = curr;
        curr = temp;


        lastMove = moveTo;
    }


    double full = -HUGE_VAL;
    FOR(k,nState) {
        if (lastMove==0 || lastMove==3 || lastMove==5 || lastMove==9 || lastMove==11)
            full = sumLogs(full,prev->g(k)+hmm->probXW(k));
        else if (lastMove==1 || lastMove==7 || lastMove==8 || lastMove==13 || lastMove==14)
            full = sumLogs(full,prev->g(k)+hmm->probYW(k));
        else if (lastMove==2)
            full = sumLogs(full,prev->g(k)+hmm->probMW(k));
    }

    pSite->index(1); // previous site
    pSite->prev();
    cSite->index(1); // current site
    cSite->prev();

     cSite->prev();

    prev = vec1;
    curr = vec2;

    moveTo = pSite->currMatchState();
    moveFrom = cSite->currMatchState();

    double all=0;
    FOR(k,nState) {
        if (moveTo==0 || moveTo==3 || moveTo==5 || moveTo==9 || moveTo==11){
            prev->s( hmm->probXW(k), k );
        } else if (moveTo==1 || moveTo==7 || moveTo==8 || moveTo==13 || moveTo==14) {
            prev->s( hmm->probYW(k), k );
        } else if (moveTo==2){
            prev->s( hmm->probMW(k), k );
        }
        t = exp(pSite->stateProb(k)+prev->g(k) - full);
        pSite->stateProb( t, k );

        all+=pSite->stateProb(k);
    }

    FOR(k,nState) {
        pSite->stateProb( pSite->stateProb(k)/all, k );
    }

    cSite->next();pSite->next();
    do {
        cSite->prev();pSite->prev();

        moveTo = pSite->currMatchState();
        moveFrom = cSite->currMatchState();

        if (moveFrom<0||moveFrom>14)
            moveFrom = 2;

        if (LOGVALUES) {
            FOR(l,nState) {

                if ((FOREVER || FOREVER_OLD) && moveTo>2) {
                    score->s(0,l);

                } else {

                    score->s(small,l);

                    FOR(j,sAlpha) {
                        score->alog( hmm->logCharBgFreq(l,j)+pSite->mlCharProb(l,j), l );
                    }
                }
            }
        } else {
            FOR(l,nState) {

                if ((FOREVER || FOREVER_OLD) && moveTo>2) {
                    score->s(1,l);

                } else {

                    score->s(0,l);

                    FOR(j,sAlpha) {
                        score->a( hmm->charBgFreq(l,j)*pSite->mlCharProb(l,j), l );
                    }
                    score->clog(l);
                }
            }
        }

        FOR(k,nState) {

            sum=-HUGE_VAL;

            FOR(l,nState) {

                if (moveFrom==0 || moveFrom==3 || moveFrom==5 || moveFrom==9 || moveFrom==11) {

                    if (moveTo==0 || moveTo==3 || moveTo==5 || moveTo==9 || moveTo==11) {
                        sum = sumLogs(sum,hmm->probXX(k,l) + score->g(l) + prev->g(l));
                    } else if (moveTo==1 || moveTo==7 || moveTo==8 || moveTo==13 || moveTo==14) {
                        sum = sumLogs(sum,hmm->probXY(k,l) + score->g(l) + prev->g(l));
                    } else if (moveTo==2) {
                        sum = sumLogs(sum,hmm->probXM(k,l) + score->g(l) + prev->g(l));
                    } else {
                        cout<<"CharacterProbabilityL: impossible pointer "<<moveFrom<<" "<<moveTo<<endl;
                    }

                } else if (moveFrom==1 || moveFrom==7 || moveFrom==8 || moveFrom==13 || moveFrom==14) {

                    if (moveTo==0 || moveTo==3 || moveTo==5 || moveTo==9 || moveTo==11) {
                        sum = sumLogs(sum,hmm->probYX(k,l) + score->g(l) + prev->g(l));
                    } else if (moveTo==1 || moveTo==7 || moveTo==8 || moveTo==13 || moveTo==14) {
                        sum = sumLogs(sum,hmm->probYY(k,l) + score->g(l) + prev->g(l));
                    } else if (moveTo==2) {
                        sum = sumLogs(sum,hmm->probYM(k,l) + score->g(l) + prev->g(l));
                    } else {
                        cout<<"CharacterProbabilityL: impossible pointer "<<moveFrom<<" "<<moveTo<<endl;
                    }

                } else if (moveFrom==2) {

                    if (moveTo==0 || moveTo==3 || moveTo==5 || moveTo==9 || moveTo==11) {
                        sum = sumLogs(sum,hmm->probMX(k,l) + score->g(l) + prev->g(l));
                    } else if (moveTo==1 || moveTo==7 || moveTo==8 || moveTo==13 || moveTo==14) {
                        sum = sumLogs(sum,hmm->probMY(k,l) + score->g(l) + prev->g(l));
                    } else if (moveTo==2) {
                        sum = sumLogs(sum,hmm->probMM(k,l) + score->g(l) + prev->g(l));
                    } else {
                        cout<<"CharacterProbabilityL: impossible pointer "<<moveFrom<<" "<<moveTo<<endl;
                    }
                }
            }
            curr->s(sum,k);
        }
        all=0;


        FOR(k,nState) {
            t = exp(cSite->stateProb(k)+curr->g(k) - full);
            cSite->stateProb( t, k );
            all+=t;
        }

        FOR(k,nState) {
			cSite->stateProb( cSite->stateProb(k)/all, k );
		}

        temp = prev;
        prev = curr;
        curr = temp;

        lastMove = moveTo;

	} while (cSite->index()!=0);

    fwdScore = full;

    full = -HUGE_VAL;

    FOR(k,nState) {
        if (lastMove==0 || lastMove==3 || lastMove==9)
            full = sumLogs(full,prev->g(k)+hmm->structBgFreq(k)+hmm->probWX(k));
        else if (lastMove==1 || lastMove==7 || lastMove==13)
            full = sumLogs(full,prev->g(k)+hmm->structBgFreq(k)+hmm->probWY(k));
        else if (lastMove==2 || moveTo==5 || moveTo==8 || moveTo==11 || moveTo==14)
            full = sumLogs(full,prev->g(k)+hmm->structBgFreq(k)+hmm->probWM(k));
        else
            cout<<"CharacterProbability: impossible last "<<lastMove<<endl;
    }

    bwdScore = full;

    delete score;
    delete vec1;
    delete vec2;

	delete pSite;
	delete cSite;
	delete sSite;

}

void CharacterProbability::logScoresSS()
{

    for (;cSite->index()!=1; cSite->next()) {

        cSite->permInsertion(0);

        // X-gap
        //
        if (cSite->currMatchState()==0){

            skipMatch = -1;

            FOR(k,nState) {
                FOR(n,sAlpha) {
                    sum = hmm->logCharSbProbL(k,n,t1->charAt(li));
                    cSite->mlCharProb( sum, k, n );
                }
            }

            li++;

            // Y-gap
            //
        } else if (cSite->currMatchState()==1){

            skipMatch = -1;

            FOR(k,nState) {
                FOR(n,sAlpha) {
                    sum = hmm->logCharSbProbR(k,n,t2->charAt(ri));
                    cSite->mlCharProb( sum, k, n );
                }
            }

            ri++;

            // match
            //
        } else if (cSite->currMatchState()==2){

            FOR(k,nState) {
                FOR(n,sAlpha) {
                    sum1 = hmm->logCharSbProbL(k,n,t1->charAt(li));
                    sum2 = hmm->logCharSbProbR(k,n,t2->charAt(ri));
                    cSite->mlCharProb( sum1+sum2, k, n );
                }
            }

            li++;
            ri++;

        } else {
			cout<<"CharacterProbability: impossible state "<<cSite->index()<<" "<<cSite->currMatchState()<<endl;
        }
    }
}

void CharacterProbability::logScoresSM()
{

    for (;cSite->index()!=1; cSite->next()) {

      if( a2->isPermInsertion(ri) )
        cSite->permInsertion(1);
      else
        cSite->permInsertion(0);

        // X-gap
        //
        if (cSite->currMatchState()==0){

            skipMatch = -1;

            FOR(k,nState) {
                FOR(n,sAlpha) {
                    sum = hmm->logCharSbProbL(k,n,t1->charAt(li));
                    cSite->mlCharProb( sum, k, n );
                }
            }

            li++;

            // Y-gap
            //
        } else if (cSite->currMatchState()==1){

            skipMatch = -1;

            FOR(k,nState) {
                FOR(n,sAlpha) {
                    sum=small;
                    FOR(m,sAlpha) {
                        sum = sumLogs(sum,hmm->logCharSbProbR(k,n,m)+a2->mlCharProbAt(m,ri,k));
                    }
                    cSite->mlCharProb( sum, k, n );
                }
            }

            ri++;

            // match
            //
        } else if (cSite->currMatchState()==2){

            FOR(k,nState) {
                FOR(n,sAlpha) {
                    sum1 = hmm->logCharSbProbL(k,n,t1->charAt(li));
                    sum2=small;
                    FOR(m,sAlpha) {
                        sum2 = sumLogs(sum2,hmm->logCharSbProbR(k,n,m)+a2->mlCharProbAt(m,ri,k));
                    }
                    cSite->mlCharProb( sum1+sum2, k, n );
                }
            }

            li++;
            ri++;

            // insertion skipped; copy old values
            //
        } else if (cSite->currMatchState()==7 || cSite->currMatchState()==8 || cSite->currMatchState()==13 || cSite->currMatchState()==14){

            if (FOREVER && ( cSite->currMatchState()==8 || cSite->currMatchState()==14 ) ){ /*e090626*/
                if (skipMatch<0) {
                    sSite->index(cSite->index());
                }
                skipMatch = 1;
            } else {
                skipMatch = -1;
            }

            if (FOREVER_OLD) {
                FOR(k,nState) {
                    FOR(n,sAlpha) {
                        cSite->mlCharProb( small, k, n );
                    }
                }
            } else {
                FOR(k,nState) {
                    FOR(n,sAlpha) {
                        cSite->mlCharProb( a2->mlCharProbAt(n,ri,k), k, n );
                    }
                }
            }

            ri++;

        } else {
            cout<<"CharacterProbability: impossible state "<<cSite->currMatchState()<<endl;
        }

        if (FOREVER && cSite->currMatchState()==2 && skipMatch>0) {
            if(TERMF || ( sSite->getLSite()!=0 && cSite->index()!=1 ))
            {
                for (;sSite->index()!=cSite->index(); sSite->next()) {
                    FOR(k,nState) {
                        FOR(n,sAlpha) {
                            sSite->mlCharProb( small, k, n );
                        }
                    }
                    sSite->permInsertion(1);
                }
            }
            skipMatch = -1;
        }

    }

    // 030310
    if (FOREVER && skipMatch>0) {
        if(TERMF || ( sSite->getLSite()!=0 && cSite->index()!=1 ))
        {
            for (;sSite->index()!=cSite->index(); sSite->next()) {
                FOR(k,nState) {
                    FOR(n,sAlpha) {
                        sSite->mlCharProb( small, k, n );
                    }
                }
                sSite->permInsertion(1);
            }
        }
        skipMatch = -1;
    }


}

void CharacterProbability::logScoresMS()
{

	for (;cSite->index()!=1; cSite->next()) {

      if( a1->isPermInsertion(li) )
        cSite->permInsertion(1);
      else
        cSite->permInsertion(0);

        // X-gap
        //
        if (cSite->currMatchState()==0){

            skipMatch = -1;

            FOR(k,nState) {
                FOR(n,sAlpha) {
                    sum=small;
                    FOR(m,sAlpha) {
                        sum = sumLogs(sum,hmm->logCharSbProbL(k,n,m)+a1->mlCharProbAt(m,li,k));
                    }
                    cSite->mlCharProb( sum, k, n );
                }
            }

            li++;

            // Y-gap
            //
        } else if (cSite->currMatchState()==1){

            skipMatch = -1;

            FOR(k,nState) {
                FOR(n,sAlpha) {
                    sum = hmm->logCharSbProbR(k,n,t2->charAt(ri));
                    cSite->mlCharProb( sum, k, n );
                }
            }

            ri++;

            // match
            //
        } else if (cSite->currMatchState()==2){

            FOR(k,nState) {
                FOR(n,sAlpha) {
                    sum1=small;
                    FOR(m,sAlpha) {
                        sum1 = sumLogs(sum1,hmm->logCharSbProbL(k,n,m)+a1->mlCharProbAt(m,li,k));
                    }
                    sum2 = hmm->logCharSbProbR(k,n,t2->charAt(ri));
                    cSite->mlCharProb( sum1+sum2, k, n );
                }
            }

            li++;
            ri++;

        } else if (cSite->currMatchState()==3 || cSite->currMatchState()==5 || cSite->currMatchState()==9 || cSite->currMatchState()==11){

            if (FOREVER && ( cSite->currMatchState()==5 || cSite->currMatchState()==11 ) ){  /*e090626*/
                if (skipMatch<0) {
                    sSite->index(cSite->index());
                }
                skipMatch = 1;
            } else {
                skipMatch = -1;
            }

            if (FOREVER_OLD) {
                FOR(k,nState) {
                    FOR(n,sAlpha) {
                        cSite->mlCharProb( small, k, n );
                    }
                }
            } else {
                FOR(k,nState) {
                    FOR(n,sAlpha) {
                        cSite->mlCharProb( a1->mlCharProbAt(n,li,k), k, n );
                    }
                }
            }

            li++;

        } else {
            cout<<"CharacterProbability: impossible state "<<cSite->currMatchState()<<endl;
        }

        if (FOREVER && cSite->currMatchState()==2 && skipMatch>0) {
            if(TERMF || ( sSite->getLSite()!=0 && cSite->index()!=1 ))
            {
                for (;sSite->index()!=cSite->index(); sSite->next()) {
                    FOR(k,nState) {
                        FOR(n,sAlpha) {
                            sSite->mlCharProb( small, k, n );
                        }
                    }
                    sSite->permInsertion(1);
                }
            }
            skipMatch = -1;
        }

    }

    //030310
    if (FOREVER && skipMatch>0) {
        if(TERMF || ( sSite->getLSite()!=0 && cSite->index()!=1 ))
        {
            for (;sSite->index()!=cSite->index(); sSite->next()) {
                FOR(k,nState) {
                    FOR(n,sAlpha) {
                        sSite->mlCharProb( small, k, n );
                    }
                }
                sSite->permInsertion(1);
            }
        }
        skipMatch = -1;
    }

}

void CharacterProbability::logScoresMM()
{

	// Compute the probability of alternative characters given the tree below
    //
    for (;cSite->index()!=1; cSite->next()) {

      if( a1->isPermInsertion(li) || a2->isPermInsertion(ri) )
        cSite->permInsertion(1);
      else
        cSite->permInsertion(0);

        // X-gap
        //
        if (cSite->currMatchState()==0){

            skipMatch = -1;

            FOR(k,nState) {
                FOR(n,sAlpha) {
                    sum=small;
                    FOR(m,sAlpha) {
                        sum = sumLogs(sum,hmm->logCharSbProbL(k,n,m)+a1->mlCharProbAt(m,li,k));
                    }
                    cSite->mlCharProb( sum, k, n );
                }
            }

            li++;

            // Y-gap
            //
        } else if (cSite->currMatchState()==1){

            skipMatch = -1;

            FOR(k,nState) {
                FOR(n,sAlpha) {
                    sum=small;
                    FOR(m,sAlpha) {
                        sum = sumLogs(sum,hmm->logCharSbProbR(k,n,m)+a2->mlCharProbAt(m,ri,k));
                    }
                    cSite->mlCharProb( sum, k, n );
                }
            }

            ri++;

            // match
            //
        } else if (cSite->currMatchState()==2){

            FOR(k,nState) {
                FOR(n,sAlpha) {
                    sum1=small;
                    sum2=small;
                    FOR(m,sAlpha) {
                        sum1 = sumLogs(sum1,hmm->logCharSbProbL(k,n,m)+a1->mlCharProbAt(m,li,k));
                        sum2 = sumLogs(sum2,hmm->logCharSbProbR(k,n,m)+a2->mlCharProbAt(m,ri,k));
                    }
                    cSite->mlCharProb( sum1+sum2, k, n );
                }
            }

            li++;
            ri++;

            // insertion skipped; copy old values
            //
        } else if (cSite->currMatchState()==3 || cSite->currMatchState()==5 || cSite->currMatchState()==9 || cSite->currMatchState()==11){

            if (FOREVER && ( cSite->currMatchState()==5 || cSite->currMatchState()==11 ) ){ /*e090626*/
                if (skipMatch<0) {
                    sSite->index(cSite->index());
                }
                skipMatch = 1;
            } else {
                skipMatch = -1;
            }

            if (FOREVER_OLD) {
                FOR(k,nState) {
                    FOR(n,sAlpha) {
                        cSite->mlCharProb( small, k, n );
                    }
                }
            } else {
                FOR(k,nState) {
                    FOR(n,sAlpha) {
                        cSite->mlCharProb( a1->mlCharProbAt(n,li,k), k, n );
                    }
                }
            }

            li++;

            // insertion skipped; copy old values
            //
        } else if (cSite->currMatchState()==7 || cSite->currMatchState()==8 || cSite->currMatchState()==13 || cSite->currMatchState()==14){

            if (FOREVER && ( cSite->currMatchState()==8 || cSite->currMatchState()==14 ) ){ /*e090626*/
                if (skipMatch<0) {
                    sSite->index(cSite->index());
                }
                skipMatch = 1;
            } else {
                skipMatch = -1;
            }

            if (FOREVER_OLD) {
                FOR(k,nState) {
                    FOR(n,sAlpha) {
                        cSite->mlCharProb( small, k, n );
                    }
                }
            } else {
                FOR(k,nState) {
                    FOR(n,sAlpha) {
                        cSite->mlCharProb( a2->mlCharProbAt(n,ri,k), k, n );
                    }
                }
            }

            ri++;

        } else {
            cout<<"CharacterProbability: impossible state "<<cSite->currMatchState()<<endl;
        }

        if (FOREVER && cSite->currMatchState()==2 && skipMatch>0) {
            if(TERMF || ( sSite->getLSite()!=0 && sSite->getRSite()!=1 ))
            {
                for (;sSite->index()!=cSite->index(); sSite->next()) {
                    FOR(k,nState) {
                        FOR(n,sAlpha) {
                            sSite->mlCharProb( small, k, n );
                        }
                    }
                    sSite->permInsertion(1);
                }
            }
            skipMatch = -1;
        }

    }

    //030310
    if (FOREVER && skipMatch>0) {
        if(TERMF || ( sSite->getLSite()!=0 && sSite->getRSite()!=1 ))
        {
            for (;sSite->index()!=cSite->index(); sSite->next()) {
                FOR(k,nState) {
                    FOR(n,sAlpha) {
                        sSite->mlCharProb( small, k, n );
                    }
                }
                sSite->permInsertion(1);
            }
        }
        skipMatch = -1;
    }

}

void CharacterProbability::scoresSS()
{

    for (;cSite->index()!=1; cSite->next()) {

        cSite->permInsertion(0);

        // X-gap
        //
        if (cSite->currMatchState()==0){

            skipMatch = -1;

            FOR(k,nState) {
                FOR(n,sAlpha) {
                    sum = hmm->charSbProbL(k,n,t1->charAt(li));
                    cSite->mlCharProb( sum, k, n );
                }
            }

            li++;

            // Y-gap
            //
        } else if (cSite->currMatchState()==1){

            skipMatch = -1;

            FOR(k,nState) {
                FOR(n,sAlpha) {
                    sum = hmm->charSbProbR(k,n,t2->charAt(ri));
                    cSite->mlCharProb( sum, k, n );
                }
            }

            ri++;

            // match
            //
        } else if (cSite->currMatchState()==2){

            FOR(k,nState) {
                FOR(n,sAlpha) {
                    sum1 = hmm->charSbProbL(k,n,t1->charAt(li));
                    sum2 = hmm->charSbProbR(k,n,t2->charAt(ri));
                    cSite->mlCharProb( sum1*sum2, k, n );
                }
            }

            li++;
            ri++;

        } else {
			cout<<"CharacterProbability: impossible state "<<cSite->index()<<" "<<cSite->currMatchState()<<endl;
        }

    }

}

void CharacterProbability::scoresSM()
{

	for (;cSite->index()!=1; cSite->next()) {

      if( a2->isPermInsertion(ri) )
        cSite->permInsertion(1);
      else
        cSite->permInsertion(0);

        // X-gap
        //
        if (cSite->currMatchState()==0){

            skipMatch = -1;

            FOR(k,nState) {
                FOR(n,sAlpha) {
                    sum = hmm->charSbProbL(k,n,t1->charAt(li));
                    cSite->mlCharProb( sum, k, n );
                }
            }

            li++;

            // Y-gap
            //
        } else if (cSite->currMatchState()==1){

            skipMatch = -1;

            FOR(k,nState) {
                FOR(n,sAlpha) {
                    sum=0;
                    FOR(m,sAlpha) {
                        sum += hmm->charSbProbR(k,n,m)*a2->mlCharProbAt(m,ri,k);
                    }
                    cSite->mlCharProb( sum, k, n );
                }
            }

            ri++;

            // match
            //
        } else if (cSite->currMatchState()==2){

            FOR(k,nState) {
                FOR(n,sAlpha) {
                    sum1 = hmm->charSbProbL(k,n,t1->charAt(li));
                    sum2=0;
                    FOR(m,sAlpha) {
                        sum2 += hmm->charSbProbR(k,n,m)*a2->mlCharProbAt(m,ri,k);
                    }
                    cSite->mlCharProb( sum1*sum2, k, n );
                }
            }

            li++;
            ri++;

            // insertion skipped; copy old values
            //
        } else if (cSite->currMatchState()==7 || cSite->currMatchState()==8 || cSite->currMatchState()==13 || cSite->currMatchState()==14){

            if (FOREVER && ( cSite->currMatchState()==8 || cSite->currMatchState()==14 ) ){ /*e090626*/
                if (skipMatch<0) {
                    sSite->index(cSite->index());
                }
                skipMatch = 1;
            } else {
                skipMatch = -1;
            }

            if (FOREVER_OLD) {
                FOR(k,nState) {
                    FOR(n,sAlpha) {
                        cSite->mlCharProb( 0, k, n );
                    }
                }
            } else {
                FOR(k,nState) {
                    FOR(n,sAlpha) {
                        cSite->mlCharProb( a2->mlCharProbAt(n,ri,k), k, n );
                    }
                }
            }

            ri++;

        } else {
            cout<<"CharacterProbability: impossible state "<<cSite->currMatchState()<<endl;
        }

        //030310
        if (FOREVER && cSite->currMatchState()==2 && skipMatch>0) {
            if(TERMF || ( sSite->getLSite()!=0 && sSite->getRSite()!=1 ))
            {
                for (;sSite->index()!=cSite->index(); sSite->next()) {
                    FOR(k,nState) {
                        FOR(n,sAlpha) {
                            sSite->mlCharProb( 0, k, n );
                        }
                    }
                    sSite->permInsertion(1);
                }
            }
            skipMatch = -1;
        }
    }

    //030310
    if (FOREVER && skipMatch>0) {
        if(TERMF || ( sSite->getLSite()!=0 && cSite->index()!=1 ))
        {
            for (;sSite->index()!=cSite->index(); sSite->next()) {
                FOR(k,nState) {
                    FOR(n,sAlpha) {
                        sSite->mlCharProb( 0, k, n );
                    }
                }
                sSite->permInsertion(1);
            }
        }
        skipMatch = -1;
    }

}

void CharacterProbability::scoresMS()
{

	for (;cSite->index()!=1; cSite->next()) {

      if( a1->isPermInsertion(li) )
        cSite->permInsertion(1);
      else
        cSite->permInsertion(0);

        // X-gap
        //
        if (cSite->currMatchState()==0){

            skipMatch = -1;

            FOR(k,nState) {
                FOR(n,sAlpha) {
                    sum=0;
                    FOR(m,sAlpha) {
                        sum += hmm->charSbProbL(k,n,m)*a1->mlCharProbAt(m,li,k);
                    }
                    cSite->mlCharProb( sum, k, n );
                }
            }

            li++;

            // Y-gap
            //
        } else if (cSite->currMatchState()==1){

            skipMatch = -1;

            FOR(k,nState) {
                FOR(n,sAlpha) {
                    sum = hmm->charSbProbR(k,n,t2->charAt(ri));
                    cSite->mlCharProb( sum, k, n );
                }
            }

            ri++;

            // match
            //
        } else if (cSite->currMatchState()==2){

            FOR(k,nState) {
                FOR(n,sAlpha) {
                    sum1=0;
                    sum2 = hmm->charSbProbR(k,n,t2->charAt(ri));
                    FOR(m,sAlpha) {
                        sum1 += hmm->charSbProbL(k,n,m)*a1->mlCharProbAt(m,li,k);
                    }
                    cSite->mlCharProb( sum1*sum2, k, n );
                }
            }

            li++;
            ri++;

            // insertion skipped; copy old values
            //
        } else if (cSite->currMatchState()==3 || cSite->currMatchState()==5 || cSite->currMatchState()==9 || cSite->currMatchState()==11){

            if (FOREVER && ( cSite->currMatchState()==5 || cSite->currMatchState()==11 ) ){ /*e090626*/
                if (skipMatch<0) {
                    sSite->index(cSite->index());
                }
                skipMatch = 1;
            } else {
                skipMatch = -1;
            }

            if (FOREVER_OLD) {
                FOR(k,nState) {
                    FOR(n,sAlpha) {
                        cSite->mlCharProb( 0, k, n );
                    }
                }
            } else {
                FOR(k,nState) {
                    FOR(n,sAlpha) {
                        cSite->mlCharProb( a1->mlCharProbAt(n,li,k), k, n );
                    }
                }
            }

            li++;

        } else {
            cout<<"CharacterProbability: impossible state "<<cSite->currMatchState()<<endl;
        }

        //030310
        if (FOREVER && cSite->currMatchState()==2 && skipMatch>0) {
            if(TERMF || ( sSite->getLSite()!=0 && sSite->getRSite()!=1 ))
            {
                for (;sSite->index()!=cSite->index(); sSite->next()) {
                    FOR(k,nState) {
                        FOR(n,sAlpha) {
                            sSite->mlCharProb( 0, k, n );
                        }
                    }
                    sSite->permInsertion(1);
                }
            }
            skipMatch = -1;
        }
    }

    //030310
    if (FOREVER && skipMatch>0) {
        if(TERMF || ( sSite->getLSite()!=0 && cSite->index()!=1 ))
        {
            for (;sSite->index()!=cSite->index(); sSite->next()) {
                FOR(k,nState) {
                    FOR(n,sAlpha) {
                        sSite->mlCharProb( 0, k, n );
                    }
                }
                sSite->permInsertion(1);
            }
        }
        skipMatch = -1;
    }

}

void CharacterProbability::scoresMM()
{

	for (;cSite->index()!=1; cSite->next()) {

        if( a1->isPermInsertion(li) || a2->isPermInsertion(ri) )
          cSite->permInsertion(1);
        else
          cSite->permInsertion(0);

        // X-gap
        //
        if (cSite->currMatchState()==0){

            skipMatch = -1;

            FOR(k,nState) {
                FOR(n,sAlpha) {
                    sum=0;
                    FOR(m,sAlpha) {
                        sum += hmm->charSbProbL(k,n,m)*a1->mlCharProbAt(m,li,k);
                    }
                    cSite->mlCharProb( sum, k, n );
                }
            }

            li++;

            // Y-gap
            //
        } else if (cSite->currMatchState()==1){

            skipMatch = -1;

            FOR(k,nState) {
                FOR(n,sAlpha) {
                    sum=0;
                    FOR(m,sAlpha) {
                        sum += hmm->charSbProbR(k,n,m)*a2->mlCharProbAt(m,ri,k);
                    }
                    cSite->mlCharProb( sum, k, n );
                }
            }

            ri++;

            // match
            //
        } else if (cSite->currMatchState()==2){

            FOR(k,nState) {
                FOR(n,sAlpha) {
                    sum1=0;
                    sum2=0;
                    FOR(m,sAlpha) {
                        sum1 += hmm->charSbProbL(k,n,m)*a1->mlCharProbAt(m,li,k);
                        sum2 += hmm->charSbProbR(k,n,m)*a2->mlCharProbAt(m,ri,k);
                    }
                    cSite->mlCharProb( sum1*sum2, k, n );
                }
            }

            li++;
            ri++;

            // insertion skipped; copy old values
            //
        } else if (cSite->currMatchState()==3 || cSite->currMatchState()==5 || cSite->currMatchState()==9 || cSite->currMatchState()==11){

            if (FOREVER && ( cSite->currMatchState()==5 || cSite->currMatchState()==11 ) ){ /*e090626*/
                if (skipMatch<0) {
                    sSite->index(cSite->index());
                }
                skipMatch = 1;
            } else {
                skipMatch = -1;
            }

            if (FOREVER_OLD) {
                FOR(k,nState) {
                    FOR(n,sAlpha) {
                        cSite->mlCharProb( 0, k, n );
                    }
                }
            } else {
                FOR(k,nState) {
                    FOR(n,sAlpha) {
                        cSite->mlCharProb( a1->mlCharProbAt(n,li,k), k, n );
                    }
                }
            }

            li++;

            // insertion skipped; copy old values
            //
        } else if (cSite->currMatchState()==7 || cSite->currMatchState()==8 || cSite->currMatchState()==13 || cSite->currMatchState()==14){

            if (FOREVER && ( cSite->currMatchState()==8 || cSite->currMatchState()==14 ) ){ /*e090626*/
                if (skipMatch<0) {
                    sSite->index(cSite->index());
                }
                skipMatch = 1;
            } else {
                skipMatch = -1;
            }

            if (FOREVER_OLD) {
                FOR(k,nState) {
                    FOR(n,sAlpha) {
                        cSite->mlCharProb( 0, k, n );
                    }
                }
            } else {
                FOR(k,nState) {
                    FOR(n,sAlpha) {
                        cSite->mlCharProb( a2->mlCharProbAt(n,ri,k), k, n );
                    }
                }
            }

            ri++;

        } else {
            cout<<"CharacterProbability: impossible state "<<cSite->currMatchState()<<endl;
        }

        //030310
        if (FOREVER && cSite->currMatchState()==2 && skipMatch>0) {
            if(TERMF || ( sSite->getLSite()!=0 && sSite->getRSite()!=1 ))
            {
                for (;sSite->index()!=cSite->index(); sSite->next()) {
                    FOR(k,nState) {
                        FOR(n,sAlpha) {
                            sSite->mlCharProb( 0, k, n );
                        }
                    }
                    sSite->permInsertion(1);
                }
            }
            skipMatch = -1;
        }
    }

    //030310
    if (FOREVER && skipMatch>0) {
        if(TERMF || ( sSite->getLSite()!=0 && cSite->index()!=1 ))
        {
            for (;sSite->index()!=cSite->index(); sSite->next()) {
                FOR(k,nState) {
                    FOR(n,sAlpha) {
                        sSite->mlCharProb( 0, k, n );
                    }
                }
                sSite->permInsertion(1);
            }
        }
        skipMatch = -1;
    }

}

