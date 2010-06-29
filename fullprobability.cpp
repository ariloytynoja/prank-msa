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

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include "config.h"
#include "fullprobability.h"

using namespace std;

FullProbability::~FullProbability()
{
    if (!FULLFULL) {
        delete minBIndex;
        delete maxBIndex;
        delete diffIndex;
    }
}
FullProbability::FullProbability(Sequence* s1,Sequence* s2,PhyloMatchScore* pms)
{
    seq1 = s1;
    seq2 = s2;
    msr = pms;

	sAlpha = hmm->getASize();
    nState = hmm->getNStates();

    small = -HUGE_VAL;
}

void FullProbability::initialiseIndex(Site *sites)
{

    minBIndex = new IntMatrix(seq2->lengthF()+1,"minBIndex");
    maxBIndex = new IntMatrix(seq2->lengthF()+1,"maxBIndex");
    diffIndex = new IntMatrix(seq2->lengthF()+1,"diffIndex");

    minBIndex->initialise(0);
    maxBIndex->initialise(0);
    diffIndex->initialise(0);

    int mini = 0;
    int minj = 0;
    int maxi = 0;
    int maxj = 0;

	sites->index(0);
	sites->next();

	while (sites->nullSite()){
        sites->next();
    }
	while(sites->index()!=1 && !sites->nullSite()) {

        if (sites->currMatchState()!=1) {
            diffIndex->a( 1, mini );
        }


        if (sites->currMatchState()!=0) {
            mini++;
        }
        if (sites->currMatchState()!=1) {
            minj++;
        }

        maxBIndex->s( maxj, maxi );

        if (sites->currMatchState()!=0) {
            minBIndex->s( minj, mini );
        } else if (sites->currMatchState()==0 && mini>0) {
            minBIndex->s( minBIndex->g(mini-1), mini );
        } else if (sites->currMatchState()==0){
            minBIndex->s( -1, mini );
        }

        if (sites->currMatchState()!=0) {
            maxi++;
        }
        if (sites->currMatchState()!=1) {
            maxj++;
        }

        sites->next();

		while(sites->index()!=1 && sites->nullSite()) {
			sites->next();
        }
    }

    int sl1 = seq1->lengthF();
    int sl2 = seq2->lengthF();

    minBIndex->s(0,0);
    minBIndex->s( minBIndex->g(sl2-1), sl2 );
    maxBIndex->s( sl1, sl2 );
    diffIndex->s( sl1-minBIndex->g(sl2), sl2 );


    if (FULLBAND) {
        int sum = 0;
        int msum = 0;
        int sLen2 = seq2->lengthF();

        for (int i=0;i<FBW+5 && i<sLen2;i++){
            sum += diffIndex->g(i);
        }

        msum = sum;
        for (int i=FBW+5;i<sLen2;i++){
            sum += diffIndex->g(i);
            sum -= diffIndex->g(i-FBW-5);
            if (sum>msum)
                msum=sum;
        }

        width = 2*msum+20;
        width = max(width,2*FBW+20);
    }
}

void FullProbability::alignSeqs()
{
	Site *sites = new Site();
	sites->index(0);
	sites->next();

    if (!FULLFULL)
        initialiseIndex(sites);


	sites->index(0);

	while (sites->nullSite()){
		sites->next();
	}



    int mLen1 = seq1->lengthF()+1;                 // short cuts
    int mLen2 = seq2->lengthF()+1;

    if (NOISE>1)
        cout<<"seq1 length:"<<mLen1<<" seq2 length:"<<mLen2<<endl;

    DbMatrix* matM1 = new DbMatrix(nState,mLen1,"matM1");
    DbMatrix* matX1 = new DbMatrix(nState,mLen1,"matX1");
    DbMatrix* matY1 = new DbMatrix(nState,mLen1,"matY1");

    DbMatrix* matM2 = new DbMatrix(nState,mLen1,"matM2");
    DbMatrix* matX2 = new DbMatrix(nState,mLen1,"matX2");
    DbMatrix* matY2 = new DbMatrix(nState,mLen1,"matY2");

    DbMatrix* curM = matM1;
    DbMatrix* curX = matX1;
    DbMatrix* curY = matY1;

    DbMatrix* prevM = matM2;
    DbMatrix* prevX = matX2;
    DbMatrix* prevY = matY2;

    DbMatrix* tmpM;
    DbMatrix* tmpX;
    DbMatrix* tmpY;

    if (FULLFULL) {
        curX->initialise(small); curY->initialise(small); curM->initialise(small);
        prevX->initialise(small); prevY->initialise(small); prevM->initialise(small);
    } else {
        int si = min(mLen2-1,FBW);
        for (int j=0;j<maxBIndex->g(si)+10 && j<mLen1;j++) {
            FOR(k,nState) {
                curX->s(small,k,j); curY->s(small,k,j); curM->s(small,k,j);
                prevM->s(small,k,j); prevX->s(small,k,j); prevY->s(small,k,j);
            }
        }
    }

    // Temp variables
    //
    double cX,cY,cM; // current

    // Iterate through the matrix
    //
    FOR(i,mLen2) {
        FOR(j,mLen1) {

            if (i==0 && j==0) { // Corner: starting values

                FOR(k,nState) {
                    curX->s( hmm->structBgFreq(k)+ hmm->probWX(k), k, 0 );
                    curY->s( hmm->structBgFreq(k)+ hmm->probWY(k), k, 0 );
                    curM->s( hmm->structBgFreq(k)+ hmm->probWM(k), k, 0 );
                }
                continue;

            }

            // compute values if banding not used or values are within the band
            if (FULLFULL ||
                    ( j>minBIndex->g(i)-FBW-1 && j<maxBIndex->g(i)+FBW+1 ) ||
                    ( i-FBW>=0 && i+FBW>=mLen2 && j>minBIndex->g(i-FBW) ) ||
                    ( i-FBW<0 && i+FBW<mLen2 && j<maxBIndex->g(i+FBW) ) ||
                    ( i-FBW>=0 && i+FBW<mLen2 && j>minBIndex->g(i-FBW) && j<maxBIndex->g(i+FBW) ) ){ /*e090626*/

                msr->computeFullFwd(j,i);

                FOR(k,nState) {

                    if (i==0 && j>0) { // only X-gaps are possible

                        // move into X-matrix
                        //
                        cX=cY=cM=small;

                        int l = hmm->transIndY(k,0);
                        while (l>=0) {
                            cX = sumLogs(cX,curX->g(l,j-1) + hmm->probXX(l,k) + msr->indelX(k));
                            cX = sumLogs(cX,curY->g(l,j-1) + hmm->probYX(l,k) + msr->indelX(k));
                            cX = sumLogs(cX,curM->g(l,j-1) + hmm->probMX(l,k) + msr->indelX(k));

                            l = hmm->transIndY(k,l+1);
                        }

                        curX->s( cX, k, j );
                        curY->s( cY, k, j );
                        curM->s( cM, k, j );
                        continue;

                    } else if (i>0 && j==0) { // only Y-gaps are possible

                        // move into Y-matrix
                        //
                        cX=cY=cM=small;

                        int l = hmm->transIndY(k,0);
                        while (l>=0) {
                            cY = sumLogs(cY,prevX->g(l,j) + hmm->probXY(l,k) + msr->indelY(k));
                            cY = sumLogs(cY,prevY->g(l,j) + hmm->probYY(l,k) + msr->indelY(k));
                            cY = sumLogs(cY,prevM->g(l,j) + hmm->probMY(l,k) + msr->indelY(k));

                            l = hmm->transIndY(k,l+1);
                        }

                        curX->s( cX, k, j );
                        curY->s( cY, k, j );
                        curM->s( cM, k, j );
                        continue;

                    } else {  // so far, the moves have been exceptional; from now on they are "normal"

                        // all moves
                        //
                        cX=cY=cM=small;

                        int l = hmm->transIndY(k,0);
                        while (l>=0) {
                            cX = sumLogs(cX, curX->g(l,j-1) + hmm->probXX(l,k) + msr->indelX(k));
                            cX = sumLogs(cX, curY->g(l,j-1) + hmm->probYX(l,k) + msr->indelX(k));
                            cX = sumLogs(cX, curM->g(l,j-1) + hmm->probMX(l,k) + msr->indelX(k));

                            cY = sumLogs(cY, prevX->g(l,j) + hmm->probXY(l,k) + msr->indelY(k));
                            cY = sumLogs(cY, prevY->g(l,j) + hmm->probYY(l,k) + msr->indelY(k));
                            cY = sumLogs(cY, prevM->g(l,j) + hmm->probMY(l,k) + msr->indelY(k));

                            cM = sumLogs(cM, prevX->g(l,j-1) + hmm->probXM(l,k) + msr->fullM(k));
                            cM = sumLogs(cM, prevY->g(l,j-1) + hmm->probYM(l,k) + msr->fullM(k));
                            cM = sumLogs(cM, prevM->g(l,j-1) + hmm->probMM(l,k) + msr->fullM(k));

                            l = hmm->transIndY(k,l+1);
                        }

                        curX->s( cX, k, j );
                        curY->s( cY, k, j );
                        curM->s( cM, k, j );
                    }
                }

                // wipe out the old values and surround the band area with -inf's
            } else if (
                ( j>minBIndex->g(i)-FBW-2 && j<maxBIndex->g(i)+FBW+2 ) ||
                ( i-1>=0 && i+1<mLen2 && j>minBIndex->g(i-1)-FBW-2 && j<maxBIndex->g(i+1)+FBW+2 ) ||
                ( i-FBW-1>=0 && i+FBW+1>=mLen2 && j>minBIndex->g(i-FBW-1)-1 ) ||
                ( i-FBW-1<0 && i+FBW+1<mLen2 && j<maxBIndex->g(i+FBW+1)+1 ) ||
                ( i-FBW-1>=0 && i+FBW+1<mLen2 && j>minBIndex->g(i-FBW-1)-1 && j<maxBIndex->g(i+FBW+1)+1 ) ) { /*e090626*/

                FOR(k,nState) {
                    curX->s( small, k, j );
                    curY->s( small, k, j );
                    curM->s( small, k, j );
                }
            }
        }

		while(sites->index()!=1 && sites->nInd2()==i) {

            FOR(k,nState) {
                sites->fullFwdX( curX->g(k,sites->nInd1()), k );
				sites->fullFwdY( curY->g(k,sites->nInd1()), k );
                sites->fullFwdM( curM->g(k,sites->nInd1()), k );
            }

            sites->next();

			while (sites->index()!=1 && sites->nullSite()){
                sites->next();
            }
        }


        if (NOISE>2) {
            printMatrix("fx",i,curX);
            printMatrix("fy",i,curY);
            printMatrix("fm",i,curM);
        }

        tmpM = prevM;
        tmpX = prevX;
        tmpY = prevY;

        prevM = curM;
        prevX = curX;
        prevY = curY;

        curM = tmpM;
        curX = tmpX;
        curY = tmpY;

    }

    maxFwdScore = small;
    FOR(k,nState) {
        maxFwdScore = sumLogs(maxFwdScore,
                              sumLogs(prevX->g(k,mLen1-1)+hmm->probXW(k),
                                      sumLogs(prevY->g(k,mLen1-1)+hmm->probYW(k),
                                              prevM->g(k,mLen1-1)+hmm->probMW(k))));
    }

	sites->index(1);
	sites->prev();

    while (sites->nullSite()){
        sites->prev();
    }

    curM = matM1;
    curX = matX1;
    curY = matY1;

    prevM = matM2;
    prevX = matX2;
    prevY = matY2;

    if (FULLFULL) {

        RFOR(j,mLen1-1) {
            FOR(k,nState) {
                curX->s( small, k, j); curY->s( small, k, j); curM->s( small, k, j);
                prevM->s( small, k, j); prevX->s( small, k, j); prevY->s( small, k, j);
            }
        }

    } else {
        int si = min(mLen2-1,FBW);

        for (int j=mLen1-1;j>minBIndex->g(mLen2-si-1)-10 && j>=0;j--) {
            FOR(k,nState) {
                curX->s( small, k, j); curY->s( small, k, j); curM->s( small, k, j);
                prevM->s( small, k, j); prevX->s( small, k, j); prevY->s( small, k, j);
            }
        }

    }

    RFOR(i,mLen2-1) {
        RFOR(j,mLen1-1) {

            if (i==mLen2-1 && j==mLen1-1) { // Corner: starting values

                FOR(k,nState) {

                    curX->s( hmm->probXW(k), k, j );
                    curY->s( hmm->probYW(k), k, j );
                    curM->s( hmm->probMW(k), k, j );
                }
                continue;

            }

            // compute values if banding not used or values are within th eband
            if (FULLFULL ||
                    ( j>minBIndex->g(i)-FBW-1 && j<maxBIndex->g(i)+FBW+1 ) ||
                    ( i-FBW>=0 && i+FBW>=mLen2 && j>minBIndex->g(i-FBW) ) ||
                    ( i-FBW<0 && i+FBW<mLen2 && j<maxBIndex->g(i+FBW) ) ||
                    ( i-FBW>=0 && i+FBW<mLen2 && j>minBIndex->g(i-FBW) && j<maxBIndex->g(i+FBW) ) ) { /*e090626*/


                // Compute the substitution prices
                //
                msr->computeFullBwd(j,i);
                FOR(k,nState) {

                    if (i==mLen2-1 && j<mLen1-1) { // only x-gaps possible

                        // move into X-matrix
                        //
                        cX=cY=cM=small;

                        int l = hmm->transIndX(k,0);
                        while (l>=0) {
                            cX = sumLogs(cX, hmm->probXX(k,l) + msr->indelX(l) + curX->g(l,j+1) );
                            cY = sumLogs(cY, hmm->probYX(k,l) + msr->indelX(l) + curX->g(l,j+1) );
                            cM = sumLogs(cM, hmm->probMX(k,l) + msr->indelX(l) + curX->g(l,j+1) );

                            l = hmm->transIndX(k,l+1);
                        }

                        curX->s( cX, k, j );
                        curY->s( cY, k, j );
                        curM->s( cM, k, j );

                    } else if (i<mLen2-1 && j==mLen1-1) { // only y-gaps possible

                        // move into Y-matrix
                        //
                        cX=cY=cM=small;

                        int l = hmm->transIndX(k,0);
                        while (l>=0) {
                            cX = sumLogs(cX, hmm->probXY(k,l) + msr->indelY(l) + prevY->g(l,j));
                            cY = sumLogs(cY, hmm->probYY(k,l) + msr->indelY(l) + prevY->g(l,j));
                            cM = sumLogs(cM, hmm->probMY(k,l) + msr->indelY(l) + prevY->g(l,j));

                            l = hmm->transIndX(k,l+1);
                        }

                        curX->s( cX, k, j );
                        curY->s( cY, k, j );
                        curM->s( cM, k, j );

                    } else if (i<mLen2-1 && j<mLen1-1) { // everything possible

                        // all moves
                        //
                        cX=cY=cM=small;

                        int l = hmm->transIndX(k,0);
                        while (l>=0) {
                            cX = sumLogs(cX, hmm->probXX(k,l) + msr->indelX(l) + curX->g(l,j+1));
                            cX = sumLogs(cX, hmm->probXY(k,l) + msr->indelY(l) + prevY->g(l,j));
                            cX = sumLogs(cX, hmm->probXM(k,l) + msr->fullM(l) + prevM->g(l,j+1));

                            cY = sumLogs(cY, hmm->probYX(k,l) + msr->indelX(l) + curX->g(l,j+1));
                            cY = sumLogs(cY, hmm->probYY(k,l) + msr->indelY(l) + prevY->g(l,j));
                            cY = sumLogs(cY, hmm->probYM(k,l) + msr->fullM(l) + prevM->g(l,j+1));

                            cM = sumLogs(cM, hmm->probMX(k,l) + msr->indelX(l) + curX->g(l,j+1));
                            cM = sumLogs(cM, hmm->probMY(k,l) + msr->indelY(l) + prevY->g(l,j));
                            cM = sumLogs(cM, hmm->probMM(k,l) + msr->fullM(l) + prevM->g(l,j+1));

                            l = hmm->transIndX(k,l+1);
                        }

                        curX->s( cX, k, j );
                        curY->s( cY, k, j );
                        curM->s( cM, k, j );

                    } else {
                        cout<<"FullProbability::error"<<endl;
                        exit(1);
                    }
                }

                // wipe out the old values and surround the band area with -inf's
            } else if (
                ( j>minBIndex->g(i)-FBW-3 && j<maxBIndex->g(i)+FBW+3 ) ||
                ( i-1>=0 && i+1<mLen2 && j>minBIndex->g(i-1)-FBW-3 && j<maxBIndex->g(i+1)+FBW+3 ) ||
                ( i-FBW-1>=0 && i+FBW+1>=mLen2 && j>minBIndex->g(i-FBW-1)-1 ) ||
                ( i-FBW-1<0 && i+FBW+1<mLen2 && j<maxBIndex->g(i+FBW+1)+1 ) ||
                ( i-FBW-1>=0 && i+FBW+1<mLen2 && j>minBIndex->g(i-FBW-1)-1 && j<maxBIndex->g(i+FBW+1)+1 ) ) { /*e090626*/

                FOR(k,nState) {
                    curX->s( small, k, j );
                    curY->s( small, k, j );
                    curM->s( small, k, j );
                }

            }
        }

        while (sites->nInd2()==i && sites->index()!=0) {

            for (int k=0;k<nState;k++) {
                sites->fullBwdX( curX->g(k, sites->nInd1()), k );
                sites->fullBwdY( curY->g(k, sites->nInd1()), k );
                sites->fullBwdM( curM->g(k, sites->nInd1()), k );
            }

            sites->prev();

            while (sites->index()!=0 && sites->nullSite()){
                sites->prev();
            }
        }

		if (sites->nInd2()==i && sites->index()==0) {

			for (int k=0;k<nState;k++) {
				sites->fullBwdX( curX->g(k, sites->nInd1()), k );
				sites->fullBwdY( curY->g(k, sites->nInd1()), k );
				sites->fullBwdM( curM->g(k, sites->nInd1()), k );
			}
		}

        if (NOISE>2) {
            printMatrix("bx",i,curX);
            printMatrix("by",i,curY);
            printMatrix("bm",i,curM);
        }

        tmpM = prevM;
        tmpX = prevX;
        tmpY = prevY;

        prevM = curM;
        prevX = curX;
        prevY = curY;

        curM = tmpM;
        curX = tmpX;
        curY = tmpY;

    }

    maxBwdScore = small;
    for (int k=0;k<nState;k++) {
        maxBwdScore = sumLogs(maxBwdScore,
                              sumLogs(prevX->g(k,0) + hmm->structBgFreq(k) + hmm->probWX(k),
                                      sumLogs(prevY->g(k,0) + hmm->structBgFreq(k) + hmm->probWY(k),
                                              prevM->g(k,0) + hmm->structBgFreq(k) + hmm->probWM(k))));
    }


    delete matM1;
    delete matX1;
    delete matY1;

    delete matM2;
    delete matX2;
    delete matY2;

	delete sites;
}

void FullProbability::alignBand()
{

	Site *sites = new Site();
	sites->index(0);
	sites->next();

    initialiseIndex(sites);

    int mLen1 = seq1->lengthF()+1;                 // short cuts
    int mLen2 = seq2->lengthF()+1;

    if (NOISE>1)
        cout<<"seq1 length:"<<mLen1<<" seq2 length:"<<mLen2<<endl;

    DbMatrix* matM1 = new DbMatrix(nState,width,"matM1");
    DbMatrix* matX1 = new DbMatrix(nState,width,"matX1");
    DbMatrix* matY1 = new DbMatrix(nState,width,"matY1");

    DbMatrix* matM2 = new DbMatrix(nState,width,"matM2");
    DbMatrix* matX2 = new DbMatrix(nState,width,"matX2");
    DbMatrix* matY2 = new DbMatrix(nState,width,"matY2");

    DbMatrix* curM = matM1;
    DbMatrix* curX = matX1;
    DbMatrix* curY = matY1;

    DbMatrix* prevM = matM2;
    DbMatrix* prevX = matX2;
    DbMatrix* prevY = matY2;

    DbMatrix* tmpM;
    DbMatrix* tmpX;
    DbMatrix* tmpY;

    curX->initialise(small); curY->initialise(small); curM->initialise(small);
    prevM->initialise(small); prevX->initialise(small); prevY->initialise(small);


    // Temp variables
    //
    double cX,cY,cM; // current
    int i=0; int cj=0; int rj; int dif = 0;

    // Iterate through the matrix
    //
	while (sites->nullSite()){
		sites->next();
	}

    while (sites->index()!=1) {

        while (sites->nullSite()){
            sites->next();
			if (sites->index()!=1)
                break;
        }
        FOR(j,width) {

            rj = cj + j - width/2;

            if (i==0 && rj==0) { // Corner: starting values

                FOR(k,nState) {
                    curX->s( hmm->structBgFreq(k)+ hmm->probWX(k), k , j );
                    curY->s( hmm->structBgFreq(k)+ hmm->probWY(k), k , j );
                    curM->s( hmm->structBgFreq(k)+ hmm->probWM(k), k , j );
                }
                continue;

            }

            // compute values if banding not used or values are within the band
            if (rj>=0 && rj<mLen1 &&
                    ( ( rj>minBIndex->g(i)-FBW-1 && rj<maxBIndex->g(i)+FBW+1 ) ||
                      ( i-FBW>=0 && i+FBW>=mLen2 && rj>minBIndex->g(i-FBW) ) ||
                      ( i-FBW<0 && i+FBW<mLen2 && rj<maxBIndex->g(i+FBW) ) ||
                      ( i-FBW>=0 && i+FBW<mLen2 && rj>minBIndex->g(i-FBW) && rj<maxBIndex->g(i+FBW) ) ) ) { /*e090626*/

                msr->computeFullFwd(rj,i);

                FOR(k,nState) {

                    if (i==0 && rj>0) { // only X-gaps are possible

                        // move into X-matrix
                        //
                        cX=cY=cM=small;

                        int l = hmm->transIndY(k,0);
                        while (l>=0) {
                            cX = sumLogs(cX,curX->g(l,j-1) + hmm->probXX(l,k) + msr->indelX(k));
                            cX = sumLogs(cX,curY->g(l,j-1) + hmm->probYX(l,k) + msr->indelX(k));
                            cX = sumLogs(cX,curM->g(l,j-1) + hmm->probMX(l,k) + msr->indelX(k));

                            l = hmm->transIndY(k,l+1);
                        }

                        curX->s( cX, k, j );
                        curY->s( cY, k, j );
                        curM->s( cM, k, j );
                        continue;

                    } else if (i>0 && rj==0) { // only Y-gaps are possible

                        // move into Y-matrix
                        //
                        cX=cY=cM=small;

                        int l = hmm->transIndY(k,0);
                        while (l>=0) {
                            cY = sumLogs(cY,prevX->g(l,j+dif) + hmm->probXY(l,k) + msr->indelY(k));
                            cY = sumLogs(cY,prevY->g(l,j+dif) + hmm->probYY(l,k) + msr->indelY(k));
                            cY = sumLogs(cY,prevM->g(l,j+dif) + hmm->probMY(l,k) + msr->indelY(k));

                            l = hmm->transIndY(k,l+1);
                        }

                        curX->s( cX, k, j );
                        curY->s( cY, k, j );
                        curM->s( cM, k, j );
                        continue;

                    } else {  // so far, the moves have been exceptional; from now on they are "normal"

                        // all moves
                        //
                        cX=cY=cM=small;

                        int l = hmm->transIndY(k,0);
                        while (l>=0) {
                            cX = sumLogs(cX, curX->g(l,j-1) + hmm->probXX(l,k) + msr->indelX(k));
                            cX = sumLogs(cX, curY->g(l,j-1) + hmm->probYX(l,k) + msr->indelX(k));
                            cX = sumLogs(cX, curM->g(l,j-1) + hmm->probMX(l,k) + msr->indelX(k));

                            cY = sumLogs(cY, prevX->g(l,j+dif) + hmm->probXY(l,k) + msr->indelY(k));
                            cY = sumLogs(cY, prevY->g(l,j+dif) + hmm->probYY(l,k) + msr->indelY(k));
                            cY = sumLogs(cY, prevM->g(l,j+dif) + hmm->probMY(l,k) + msr->indelY(k));

                            cM = sumLogs(cM, prevX->g(l,j-1+dif) + hmm->probXM(l,k) + msr->fullM(k));
                            cM = sumLogs(cM, prevY->g(l,j-1+dif) + hmm->probYM(l,k) + msr->fullM(k));
                            cM = sumLogs(cM, prevM->g(l,j-1+dif) + hmm->probMM(l,k) + msr->fullM(k));

                            l = hmm->transIndY(k,l+1);
                        }

                        curX->s( cX, k, j );
                        curY->s( cY, k, j );
                        curM->s( cM, k, j );
                    }
                }

                // wipe out the old values and surround the band area with -inf's
            } else if (
                ( rj>minBIndex->g(i)-FBW-2 && rj<maxBIndex->g(i)+FBW+2 ) ||
                ( i-1>=0 && i+1<mLen2 && rj>minBIndex->g(i-1)-FBW-2 && rj<maxBIndex->g(i+1)+FBW+2 ) ||
                ( i-FBW-1>=0 && i+FBW+1>=mLen2 && rj>minBIndex->g(i-FBW-1)-1 ) ||
                ( i-FBW-1<0 && i+FBW+1<mLen2 && rj<maxBIndex->g(i+FBW+1)+1 ) ||
                ( i-FBW-1>=0 && i+FBW+1<mLen2 && rj>minBIndex->g(i-FBW-1)-1 && rj<maxBIndex->g(i+FBW+1)+1 ) ) { /*e090626*/

                FOR(k,nState) {
                    curX->s(small, k, j );
                    curY->s(small, k, j );
                    curM->s(small, k, j );
                }
            }
        }


        FOR(k,nState) {
            sites->fullFwdX( curX->g(k, width/2), k );
            sites->fullFwdY( curY->g(k, width/2), k );
            sites->fullFwdM( curM->g(k, width/2), k );
        }
        sites->next();

        while (sites->index()!=1 && sites->nullSite()){
            sites->next();
        }

        xgap=0;

		while (sites->index()!=1 && sites->currMatchState()==0){

            xgap++;

            FOR(k,nState) {
                sites->fullFwdX( curX->g(k, width/2 + xgap), k );
                sites->fullFwdY( curY->g(k, width/2 + xgap), k );
                sites->fullFwdM( curM->g(k, width/2 + xgap), k );
            }

            cj++;
			sites->next();
        }

		while (sites->index()!=1 && sites->nullSite()){
			sites->next();
        }

        if (NOISE>2) {
            printMatrix("fx",i,curX);
            printMatrix("fy",i,curY);
            printMatrix("fm",i,curM);
        }

        if (sites->currMatchState()==2)
            cj++;

        dif=diffIndex->g(i);

        i++;


        tmpM = prevM;
        tmpX = prevX;
        tmpY = prevY;

        prevM = curM;
        prevX = curX;
        prevY = curY;

        curM = tmpM;
        curX = tmpX;
        curY = tmpY;

    }

    maxFwdScore = small;
    FOR(k,nState) {
        maxFwdScore = sumLogs(maxFwdScore,
                              sumLogs(prevX->g(k, width/2 + xgap)+hmm->probXW(k),
                                      sumLogs(prevY->g(k, width/2 + xgap)+hmm->probYW(k),
                                              prevM->g(k, width/2 + xgap)+hmm->probMW(k))));
    }


    curM = matM1;
    curX = matX1;
    curY = matY1;

    prevM = matM2;
    prevX = matX2;
    prevY = matY2;

    curX->initialise(small); curY->initialise(small); curM->initialise(small);
    prevX->initialise(small); prevY->initialise(small); prevM->initialise(small);

    i=mLen2-1; cj=mLen1-1; dif = 0;


    // Iterate through the matrix
    //
    sites->index(1);
	sites->prev();

    while (sites->nullSite()){
		sites->prev();
	}

    while (sites->index()!=0) {

        while (sites->nullSite()){
            sites->prev();
        }

        RFOR(j,width-1) {

            rj = cj + j - width/2;

            if (i==mLen2-1 && rj==mLen1-1) { // Corner: starting values

                FOR(k,nState) {

                    curX->s( hmm->probXW(k), k, j );
                    curY->s( hmm->probYW(k), k, j );
                    curM->s( hmm->probMW(k), k, j );
                }
                continue;

            }

            // compute values if banding not used or values are within the band
            if (rj>=0 && rj<mLen1 &&
                    ( ( rj>minBIndex->g(i)-FBW-1 && rj<maxBIndex->g(i)+FBW+1 ) ||
                      ( i-FBW>=0 && i+FBW>=mLen2 && rj>minBIndex->g(i-FBW) ) ||
                      ( i-FBW<0 && i+FBW<mLen2 && rj<maxBIndex->g(i+FBW) ) ||
                      ( i-FBW>=0 && i+FBW<mLen2 && rj>minBIndex->g(i-FBW) && rj<maxBIndex->g(i+FBW) ) ) ) { /*e090626*/

                // Compute the substitution prices
                //
                msr->computeFullBwd(rj,i);

                FOR(k,nState) {

                    if (i==mLen2-1 && rj<mLen1-1) { // only x-gaps possible

                        // move into X-matrix
                        //
                        cX=cY=cM=small;

                        int l = hmm->transIndX(k,0);
                        while (l>=0) {
                            cX = sumLogs(cX, hmm->probXX(k,l) + msr->indelX(l) + curX->g(l,j+1));
                            cY = sumLogs(cY, hmm->probYX(k,l) + msr->indelX(l) + curX->g(l,j+1));
                            cM = sumLogs(cM, hmm->probMX(k,l) + msr->indelX(l) + curX->g(l,j+1));

                            l = hmm->transIndX(k,l+1);
                        }

                        curX->s(cX, k, j );
                        curY->s(cY, k, j );
                        curM->s(cM, k, j );

                    } else if (i<mLen2-1 && rj==mLen1-1) { // only y-gaps possible

                        // move into Y-matrix
                        //
                        cX=cY=cM=small;

                        int l = hmm->transIndX(k,0);
                        while (l>=0) {
                            cX = sumLogs(cX, hmm->probXY(k,l) + msr->indelY(l) + prevY->g(l,j-dif));
                            cY = sumLogs(cY, hmm->probYY(k,l) + msr->indelY(l) + prevY->g(l,j-dif));
                            cM = sumLogs(cM, hmm->probMY(k,l) + msr->indelY(l) + prevY->g(l,j-dif));

                            l = hmm->transIndX(k,l+1);
                        }

                        curX->s(cX, k, j );
                        curY->s(cY, k, j );
                        curM->s(cM, k, j );

                    } else if (i<mLen2-1 && rj<mLen1-1) { // everything possible

                        // all moves
                        //
                        cX=cY=cM=small;

                        int l = hmm->transIndX(k,0);
                        while (l>=0) {
                            cX = sumLogs(cX, hmm->probXX(k,l) + msr->indelX(l) + curX->g(l,j+1));
                            cX = sumLogs(cX, hmm->probXY(k,l) + msr->indelY(l) + prevY->g(l,j-dif));
                            cX = sumLogs(cX, hmm->probXM(k,l) + msr->fullM(l) + prevM->g(l,j+1-dif));

                            cY = sumLogs(cY, hmm->probYX(k,l) + msr->indelX(l) + curX->g(l,j+1));
                            cY = sumLogs(cY, hmm->probYY(k,l) + msr->indelY(l) + prevY->g(l,j-dif));
                            cY = sumLogs(cY, hmm->probYM(k,l) + msr->fullM(l) + prevM->g(l,j+1-dif));

                            cM = sumLogs(cM, hmm->probMX(k,l) + msr->indelX(l) + curX->g(l,j+1));
                            cM = sumLogs(cM, hmm->probMY(k,l) + msr->indelY(l) + prevY->g(l,j-dif));
                            cM = sumLogs(cM, hmm->probMM(k,l) + msr->fullM(l) + prevM->g(l,j+1-dif));

                            l = hmm->transIndX(k,l+1);
                        }

                        curX->s(cX, k, j );
                        curY->s(cY, k, j );
                        curM->s(cM, k, j );

                    } else {
                        cout<<"FullProbability::error"<<endl;
                        exit(1);
                    }
                }

                // wipe out the old values and surround the band area with -inf's
            } else if (
                ( rj>minBIndex->g(i)-FBW-2 && rj<maxBIndex->g(i)+FBW+2 ) ||
                ( i-1>=0 && i+1<mLen2 && rj>minBIndex->g(i-1)-FBW-2 && rj<maxBIndex->g(i+1)+FBW+2 ) ||
                ( i-FBW-1>=0 && i+FBW+1>=mLen2 && rj>minBIndex->g(i-FBW-1)-1 ) ||
                ( i-FBW-1<0 && i+FBW+1<mLen2 && rj<maxBIndex->g(i+FBW+1)+1 ) ||
                ( i-FBW-1>=0 && i+FBW+1<mLen2 && rj>minBIndex->g(i-FBW-1)-1 && rj<maxBIndex->g(i+FBW+1)+1 ) ) { /*e090626*/

                FOR(k,nState) {
                    curX->s( small, k, j );
                    curY->s( small, k, j );
                    curM->s( small, k, j );
                }

            }
        }

        FOR(k,nState) {
            sites->fullBwdX( curX->g(k, width/2), k );
            sites->fullBwdY( curY->g(k, width/2), k );
            sites->fullBwdM( curM->g(k, width/2), k );
        }

        if (sites->currMatchState()==2)
            cj--;

        sites->prev();

        while (sites->index()!=0 && sites->nullSite()){
            sites->prev();
        }

        xgap = 0;

		while (sites->index()!=0 && sites->currMatchState()==0){

            xgap++;

            FOR(k,nState) {
                sites->fullBwdX( curX->g(k, width/2 - xgap), k );
                sites->fullBwdY( curY->g(k, width/2 - xgap), k );
                sites->fullBwdM( curM->g(k, width/2 - xgap), k );
            }

            cj--;
			sites->prev();
		}

		while (sites->index()!=0 && sites->nullSite()){
            sites->prev();
        }

        if (NOISE>2) {
            printMatrix("bx",i,curX);
            printMatrix("by",i,curY);
            printMatrix("bm",i,curM);
        }

        if (i>=0)
            dif=diffIndex->g(i);

        i--;


        tmpM = prevM;
        tmpX = prevX;
        tmpY = prevY;

        prevM = curM;
        prevX = curX;
        prevY = curY;

        curM = tmpM;
        curX = tmpX;
        curY = tmpY;

    }


    maxBwdScore = small;
    FOR(k,nState) {
        maxBwdScore = sumLogs(maxBwdScore,
                              sumLogs(prevX->g(k, width/2-1-xgap) + hmm->structBgFreq(k) + hmm->probWX(k),
                                      sumLogs(prevY->g(k, width/2-1-xgap) + hmm->structBgFreq(k) + hmm->probWY(k),
                                              prevM->g(k, width/2-1-xgap) + hmm->structBgFreq(k) + hmm->probWM(k))));
    }


    delete matM1;
    delete matX1;
    delete matY1;

    delete matM2;
    delete matX2;
    delete matY2;

	delete sites;
}


void FullProbability::printMatrix(string n,int i,DbMatrix* m)
{
    cout<<n<<i<<": ";
    m->print();
}
