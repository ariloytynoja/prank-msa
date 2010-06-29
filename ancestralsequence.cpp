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
#include <string>
#include <iostream>
#include <fstream>
#include "ancestralsequence.h"
#include "config.h"

using namespace std;

AncestralSequence::~AncestralSequence()
{
    if (LOGVALUES)
        delete logseqmat;
    else
        delete seqmat;

    delete lcIndex;
    delete rcIndex;
    delete insertionSite;
    delete permInsertionSite;

    if (DOPOST && postProb!=0){
        delete postProb;
    }
    if (stateProb!=0){
        delete stateProb;
    }
    if (realIndex!=0) {
        delete realIndex;
    }


    if (mlCharProb!=0) {
        delete mlCharProb;
        delete childGapSite;
        delete xGapSite;
        delete yGapSite;
    }

}

// Define an internal sequence (matrix) from a list of alignment sites
//
AncestralSequence::AncestralSequence()
        : Sequence(){

    terminal = false;
    hasPrior = false;

    string alpha = hmm->getAlphabet();
    sAlpha = alpha.length();
    if (CODON)
        sAlpha = 61;

    int nState = hmm->getNStates();

	Site *sites = new Site();

    seqLength = realLength = sites->getLength()-2;
    charseq = "";

    lcIndex = new IntMatrix(seqLength,"lcIndex");
    rcIndex = new IntMatrix(seqLength,"rcIndex");

    mlCharProb = new DbMatrix(nState,sAlpha,seqLength,"mlCharProb");
    if (DOPOST)
        postProb = new FlMatrix(seqLength,"postProb");

    stateProb = new FlMatrix(nState,seqLength,"stateProb");
    realIndex = 0;

    if (LOGVALUES) {
        logseqmat = new FlMatrix(sAlpha,seqLength,"logseqmat");
        logseqmat->initialise(-HUGE_VAL);
    } else {
        seqmat = new FlMatrix(sAlpha,seqLength,"seqmat");
        seqmat->initialise(0);
    }

    xGapSite = new IntMatrix(seqLength,"xGapsite");
    yGapSite = new IntMatrix(seqLength,"yGapsite");
    xGapSite->initialise(0);
    yGapSite->initialise(0);

    childGapSite = new IntMatrix(seqLength,"childGapSite");     // if site was gap earlier
    insertionSite = new IntMatrix(seqLength,"insertionSite");
    permInsertionSite = new IntMatrix(seqLength,"permInsertionSite");
     permInsertionSite->initialise(0);

    int i=0;
    int li=0;
    int ri=0;

    DbMatrix* meanCharProb = new DbMatrix(sAlpha,"meanCharProb");
    double t;

    int maxYGap = 0;
    int maxXGap = 0;
    int thisYGap = 0;
    int thisXGap = 0;


	sites->index(0);
	sites->next();

	while(sites->index()!=1) {

        if (DOPOST)
            postProb->s( sites->postProb(), i );

        double sum;

        if (LOGVALUES) {
            meanCharProb->initialise(-HUGE_VAL);
            sum = -HUGE_VAL;

            FOR(k,nState) {

                FOR(j,sAlpha) {
                    t = sites->mlCharProb(k,j);
                    mlCharProb->s( t, k, j, i );
                    meanCharProb->alog( t, j );
                    sum = sumLogs(sum,t);
                }
                stateProb->s( sites->stateProb(k), k, i );
            }

            if (sum>-HUGE_VAL) {
                FOR(j,sAlpha) {
                    logseqmat->s( meanCharProb->g(j)-sum, j, i );
                }
            } else {
                FOR(j,sAlpha) {
                    logseqmat->s( -HUGE_VAL, j, i );
                }
            }

        } else {
            meanCharProb->initialise(0);
            sum = 0;

            FOR(k,nState) {

                FOR(j,sAlpha) {
                    t = sites->mlCharProb(k,j);
                    mlCharProb->s( t, k, j, i );
                    meanCharProb->a( t, j );
                    sum += t;
                }
                stateProb->s( sites->stateProb(k), k, i );
            }

            if (sum>0) {
                FOR(j,sAlpha) {
                    seqmat->s( meanCharProb->g(j)/sum, j, i );
                }
            } else {
                FOR(j,sAlpha) {
                    seqmat->s( 0, j, i );
                }
            }
        }


        if (sites->currMatchState()==0){
            lcIndex->s( li++, i );
            rcIndex->s( -1, i );
            xGapSite->s( 1, i );
            yGapSite->s( 0, i );
            insertionSite->s( 0, i );
            thisXGap++;
        } else if (sites->currMatchState()==1){
            lcIndex->s( -1, i );
            rcIndex->s( ri++, i );
            xGapSite->s( 0, i );
            yGapSite->s( 1, i );
            insertionSite->s( 0, i );
            thisYGap++;
        } else if (sites->currMatchState()==2){
            lcIndex->s( li++, i );
            rcIndex->s( ri++, i );
            xGapSite->s( 0, i );
            yGapSite->s( 0, i );
            insertionSite->s( 0, i );
            thisXGap = thisYGap = 0;
        } else if (sites->currMatchState()==3 || sites->currMatchState()==5 || sites->currMatchState()==9 || sites->currMatchState()==11){
            lcIndex->s( li++, i );
            rcIndex->s( -1, i );
            xGapSite->s( 1, i );
            yGapSite->s( 0, i );
            insertionSite->s( 1, i );
        } else if (sites->currMatchState()==7 || sites->currMatchState()==8 || sites->currMatchState()==13 || sites->currMatchState()==14){
            lcIndex->s( -1, i );
            rcIndex->s( ri++, i );
            xGapSite->s( 0, i );
            yGapSite->s( 1, i );
            insertionSite->s( 1, i );
        }

        if(sites->permInsertion()){
          permInsertionSite->s( 1, i );
        }

        if (thisXGap>maxXGap)
            maxXGap = thisXGap;
        if (thisYGap>maxYGap)
            maxYGap = thisYGap;

        /////// build charseq here and avoid doing that later

        if (insertionSite->g(i)){
            charseq += "-";
        } else {
            if (LOGVALUES) {
                float ms = -HUGE_VAL;

                int mi = -1;
                FOR(j,sAlpha) {
                    if (logseqmat->g(j,i) >= ms){
                        ms = logseqmat->g(j,i);
                        mi = j;
                    }
                }
                if (mi>=0)
                    charseq += alpha.at(mi);

            } else {
                float ms = 0;

                int mi = -1;
                FOR(j,sAlpha) {
                    if (seqmat->g(j,i) >= ms){
                        ms = seqmat->g(j,i);
                        mi = j;
                    }
                }
                if (mi>=0)
                    charseq += alpha.at(mi);
            }
        }

        ///////

        if (NOISE>1) {
            cout<<i<<"/"<<seqLength<<": ";
            if (LOGVALUES) {
                FOR(j,sAlpha) {
                    cout<<exp(logseqmat->g(j,i))<<" ";
                }
            } else {
                FOR(j,sAlpha) {
                    cout<<seqmat->g(j,i)<<" ";
                }
            }
            cout<<": "<<xGapSite->g(i)<<" "<<yGapSite->g(i)<<" "<<childGapSite->g(i)<<"; "<<sites->currMatchState()<<": ";
            if (DOPOST)
                cout<<postProb->g(i)<<": ";
            FOR(k,nState) {
                cout<<sites->stateProb(k)<<", ";
            }
            cout<<endl;
        }
        i++;
		sites->next();
	}
    delete meanCharProb;

    if (PATCHMISSING && (maxXGap > missingLimit || maxYGap > missingLimit)) {
        if (NOISE>0)
            cout<<"patching missing data: "<<maxXGap<<" "<<maxYGap<<endl;

        int lastYMatch = 0;
        int lastXMatch = 0;
        thisYGap = 0;
        thisXGap = 0;
        FOR(i,seqLength) {
            if (xGapSite->g( i ) == 1) {
                thisXGap++;
            } else {
                if (thisXGap>missingLimit) {
                    for (int j=lastXMatch;j<i-1;j++) {
                        if (insertionSite->g( j ) == 0) {
                            xGapSite->s( 0, j );
                        }
                    }
                    if (NOISE>0)
                        cout<<"patchX: "<<lastXMatch<<" "<<i-1<<endl;
                }
                thisXGap = 0;
                lastXMatch = i;
            }

            if (yGapSite->g( i ) == 1) {
                thisYGap++;
            } else {
                if (thisYGap>missingLimit) {
                    for (int j=lastYMatch;j<i-1;j++) {
                        if (insertionSite->g( j ) == 0) {
                            yGapSite->s( 0, j );
                        }
                    }
                    if (NOISE>0)
                        cout<<"patchY: "<<lastYMatch<<" "<<i-1<<endl;
                }
                thisYGap = 0;
                lastYMatch = i;
            }
        }
        if (thisXGap>missingLimit) {
            for (int j=lastXMatch;j<seqLength;j++) {
                if (insertionSite->g( j ) == 0) {
                    xGapSite->s( 0, j );
                }
            }
            if (NOISE>0)
                cout<<"patchX: "<<lastXMatch<<" "<<i-1<<endl;
        }

        if (thisYGap>missingLimit) {
            for (int j=lastYMatch;j<seqLength;j++) {
                if (insertionSite->g( j ) == 0) {
                    yGapSite->s( 0, j );
                }
            }
            if (NOISE>0)
                cout<<"patchY: "<<lastYMatch<<" "<<i-1<<endl;
        }
    }

	delete sites;
}

// Gaps in the child seqs
//
void AncestralSequence::setChildGaps(Sequence *l,Sequence *r)
{
    if (NOISE>1) {
        cout<<"Set child gaps:"<<endl;
    }
    FOR(i,seqLength) {
        if ((l->isGap(lcIndex->g(i)) && r->isGap(rcIndex->g(i))) ||
                (l->isGap(lcIndex->g(i)) && rcIndex->g(i)<0) ||
                (lcIndex->g(i)<0 && r->isGap(rcIndex->g(i))))
            childGapSite->s(1,i);
        else
            childGapSite->s(0,i);

        if (NOISE>1) {
            cout<<i<<"/"<<this->length()<<": ";
            cout<<xGapSite->g(i)<<" "<<yGapSite->g(i)<<" "<<childGapSite->g(i)<<endl;
        }

    }
}

// Sites may be marked as insertion sites after an alignment to an outgroup. The index
// of character sites needs to be updated, as the posterior probability computation skips
// the insertion sites.
//
void AncestralSequence::setRealIndex(bool left)
{
// cout<<"left: "<<left<<endl;
	Site *sites = new Site();
	sites->index(0);
	sites->next();

    int diffX=0;
    int diffY=0;

	while(sites->index()!=1) {
//       cout<<sites->index()<<"; "<<sites->nInd1()<<" "<<diffX<<"; "<<sites->nInd2()<<" "<<diffY<<endl;
		if (sites->currMatchState()==0 || sites->currMatchState()==1 || sites->currMatchState()==2){
			sites->nInd1( sites->nInd1() - diffX );
			sites->nInd2( sites->nInd2() - diffY );
        } else if (( sites->currMatchState()==3 || sites->currMatchState()==5 || sites->currMatchState()==9 || sites->currMatchState()==11 ) && left){
            diffX++;
        } else if (( sites->currMatchState()==7 || sites->currMatchState()==8 || sites->currMatchState()==13 || sites->currMatchState()==14 ) && !left){
            diffY++;
        }
		sites->next();
    }

	sites->index(0);
	sites->next();

	IntMatrix* tmpIndex = new IntMatrix(sites->getLength(),"tmpIndex");
    int i=0;
    int h=0;

	while(sites->index()!=1) {
//       cout<<sites->currMatchState()<<" "<<i<<" "<<h<<endl;
        if (sites->currMatchState()==0 && left){
            tmpIndex->s( i++, h++ );
        } else if (sites->currMatchState()==1 && !left){
            tmpIndex->s( i++, h++ );
        } else if (sites->currMatchState()==2){
            tmpIndex->s( i++, h++ );
        } else if (( sites->currMatchState()==3 || sites->currMatchState()==5 || sites->currMatchState()==9 || sites->currMatchState()==11 ) && left){
            i++;
        } else if (( sites->currMatchState()==7 || sites->currMatchState()==8 || sites->currMatchState()==13 || sites->currMatchState()==14 ) && !left){
            i++;
        }
		sites->next();
	}

    realLength = h;

    if(h>0){
      realIndex = new IntMatrix(h,"realIndex");  // index for non-insertion sites
      FOR(j,h) {
          realIndex->s( tmpIndex->g(j), j );
      }
    } else {
      realIndex = new IntMatrix(1,"realIndex");
    }
    delete tmpIndex;

    if (NOISE>1) {
        cout<<"new index: ";
        FOR(j,h) {
            cout<<realIndex->g(j)<<",";
        }
        cout<<endl;
    }

	delete sites;

}


void AncestralSequence::cleanSpace()
{
    delete mlCharProb;
    delete childGapSite;
    delete xGapSite;
    delete yGapSite;

    mlCharProb = 0;
}


// Get the probability of characters in different structures given the tree.
// Takes into account the phylogeny.
//
double AncestralSequence::mlCharProbAt(int j,int i,int k)
{
    return mlCharProb->g(k,j,i);
}

// Same for insertions skipped
//
double AncestralSequence::mlCharProbAtF(int j,int i,int k)
{
    realIndex->g(i);

    return mlCharProb->g(k,j,realIndex->g(i));
}

void AncestralSequence::writeSequence(string name)
{
    char str[10];
    ofstream output((name+".seq").c_str());
    FOR(i,seqLength) {
        if (LOGVALUES) {
            FOR(j,sAlpha) {
                sprintf(str,"%.4f ",exp(logseqmat->g(j,i) ) );
                output<<str;
            }
        } else {
            FOR(j,sAlpha) {
                sprintf(str,"%.4f ",seqmat->g(j,i) );
                output<<str;
            }
        }
    }
}
