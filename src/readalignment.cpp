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
#include <fstream>
#include "readalignment.h"
#include "config.h"

using namespace std;

ReadAlignment::~ReadAlignment()
{
// 	delete beg;
// 	delete end;
// 	delete newsite;

}

ReadAlignment::ReadAlignment()
{
    count = 2;
    small = -HUGE_VAL;

}

void ReadAlignment::cleanUp()
{
    if (NOISE>1)
        cout<<"ReadAlignment::cleanUp()"<<endl;

    // clean up
    delete vX;
    delete vY;
    delete vM;
    delete xX;
    delete xM;
    delete yY;
    delete yM;
    delete wX;
    delete wM;
    delete zY;
    delete zM;


    delete ptVM;
    delete ptVX;
    delete ptVY;
    delete ptXM;
    delete ptXX;
    delete ptWM;
    delete ptWX;
    delete ptYM;
    delete ptYY;
    delete ptZM;
    delete ptZY;

}


void ReadAlignment::initialiseMatrices(int size)
{

    if (NOISE>1)
        cout<<"ReadAlignment::initialiseMatrices("<<size<<")"<<endl;

    matrixSize = size;

    sAlpha = hmm->getASize();
    nState = hmm->getNStates();

    // Initialize matrices
    //
    vX = new FlMatrix(nState,size,"vX");
    vY = new FlMatrix(nState,size,"vY");
    vM = new FlMatrix(nState,size,"vM");
    xX = new FlMatrix(nState,size,"xX");
    xM = new FlMatrix(nState,size,"xM");
    yY = new FlMatrix(nState,size,"yY");
    yM = new FlMatrix(nState,size,"yM");
    wX = new FlMatrix(nState,size,"wX");
    wM = new FlMatrix(nState,size,"wM");
    zY = new FlMatrix(nState,size,"zY");
    zM = new FlMatrix(nState,size,"zM");

    ptVM = new IntMatrix(nState,size,"ptVM");   // matrices for the backward pointers
    ptVX = new IntMatrix(nState,size,"ptVX");
    ptVY = new IntMatrix(nState,size,"ptVY");
    ptXM = new IntMatrix(nState,size,"ptXM");
    ptXX = new IntMatrix(nState,size,"ptXX");
    ptWM = new IntMatrix(nState,size,"ptWM");
    ptWX = new IntMatrix(nState,size,"ptWX");
    ptYM = new IntMatrix(nState,size,"ptYM");
    ptYY = new IntMatrix(nState,size,"ptYY");
    ptZM = new IntMatrix(nState,size,"ptZM");
    ptZY = new IntMatrix(nState,size,"ptZY");

}

int ReadAlignment::count = 2;
int ReadAlignment::nState;
int ReadAlignment::sAlpha;
int ReadAlignment::matrixSize;

FlMatrix* ReadAlignment::vX;
FlMatrix* ReadAlignment::vY;
FlMatrix* ReadAlignment::vM;
FlMatrix* ReadAlignment::xX;
FlMatrix* ReadAlignment::xM;
FlMatrix* ReadAlignment::wX;
FlMatrix* ReadAlignment::wM;
FlMatrix* ReadAlignment::yY;
FlMatrix* ReadAlignment::yM;
FlMatrix* ReadAlignment::zY;
FlMatrix* ReadAlignment::zM;

IntMatrix* ReadAlignment::ptVM;
IntMatrix* ReadAlignment::ptVX;
IntMatrix* ReadAlignment::ptVY;
IntMatrix* ReadAlignment::ptXM;
IntMatrix* ReadAlignment::ptXX;
IntMatrix* ReadAlignment::ptWM;
IntMatrix* ReadAlignment::ptWX;
IntMatrix* ReadAlignment::ptYM;
IntMatrix* ReadAlignment::ptYY;
IntMatrix* ReadAlignment::ptZM;
IntMatrix* ReadAlignment::ptZY;

bool ReadAlignment::readSeqs(Sequence* s1,Sequence* s2,PhyloMatchScore *pms,TreeNode* tn,vector<int>* path)
{

    seq1 = s1;
    seq2 = s2;

    sl1 = s1->length();
    sl2 = s2->length();

    totalSites = seq1->length()+seq2->length();
    countSites = 0;

    msr = pms;
    tnode = tn;


    FOR(k,nState)
    {
        if (NOTGAP)
        {
            vX->s(hmm->structBgFreq(k),k,0);
            vY->s(hmm->structBgFreq(k),k,0);
        }
        else
        {
            vX->s(hmm->structBgFreq(k)+hmm->probWX(k),k,0);
            vY->s(hmm->structBgFreq(k)+hmm->probWY(k),k,0);
        }
        vM->s(hmm->structBgFreq(k)+hmm->probWM(k),k,0);

    }

    vector<int>::iterator mi = path->begin();
    int move = *mi;

    int i=1;
    int j=1;
    int s = 1;
    unsigned int ii;
    for (;; s++)
    {

        if (SCREEN && totalSites>0 && countSites%reportLimit==0 && verbose == true)
        {
            FOR(ii,message.length())
            {
                cout<<'\b';
            }

            char prop[10];
            sprintf(prop,": %i",countSites*100/totalSites);
            message = currentNode+prop+"% computed                    ";

            cout<<message;
            cout.flush();
        }

        // Compute the substitution prices
        //
        msr->computeFwd( j, i );

        FOR(k,nState)
        {

            sX=sY=sM=sxX=sxM=syY=syM=swX=swM=szY=szM=-1;
            mX=mY=mM=mxX=mxM=myY=myM=mwX=mwM=mzY=mzM=small;
            cX=cY=cM=cxX=cxM=cyY=cyM=cwX=cwM=czY=czM=small;

            if (move==0)
            {

                if (seq1->fwdGapStarts( j ))   // flagged gap starts in seq1
                {

                    cxX = vX->g(k,s-1);
                    if (cxX > mxX)
                    {
                        mxX = cxX;
                        sxX = k*15+0;
                    }

                    cxM = vM->g(k,s-1);
                    if (cxM > mxM)
                    {
                        mxM = cxM;
                        sxM = k*15+2;
                    }

                    if (seq2->fwdGapEnds( i ) || i==sl2 && seq2->fwdGapContinues( i ))   // ..and another closes is seq2
                    {

                        cxM = yM->g(k,s-1);
                        if (cxM > mxM)
                        {
                            mxM = cxM;
                            sxM = k*15+8;
                        }
                    }
                    if (seq2->fwdChildGapEnds( i ) || i==sl2 && seq2->fwdChildGapContinues( i ))   // ..and another closes is seq2 child
                    {

                        cxM = zM->g(k,s-1);
                        if (cxM > mxM)
                        {
                            mxM = cxM;
                            sxM = k*15+14;
                        }
                    }
                }

                if (seq1->fwdGapContinues( j ))   // flagged gap continues in seq1
                {

                    cxX = xX->g(k,s-1);
                    if (cxX > mxX)
                    {
                        mxX = cxX;
                        sxX = k*15+3;
                    }

                    cxM = xM->g(k,s-1);
                    if (cxM > mxM)
                    {
                        mxM = cxM;
                        sxM = k*15+5;
                    }
                }


                if (seq1->fwdGapEnds( j ) || j==sl1 && seq1->fwdGapContinues( j ) )   // flagged gap ends in seq1
                {

                    int l = hmm->transIndY(k,0);
                    while (l>=0)
                    {

                        cX = max(xX->g(l,s-1) + hmm->probXX(l,k),
                                 small,
                                 xM->g(l,s-1) + hmm->probMX(l,k)) + msr->indelX(k);
                        if (cX > mX)
                        {
                            mX = cX;
                            sX = l*15+3+maxIndex;
                        }

                        l = hmm->transIndY(k,l+1);
                    }
                }

                if (seq1->fwdChildGapStarts( j ))   // flagged gap starts in seq1 child
                {

                    cwX = vX->g(k,s-1);
                    if (cwX > mwX)
                    {
                        mwX = cwX;
                        swX = k*15+0;
                    }

                    cwM = vM->g(k,s-1);
                    if (cwM > mwM)
                    {
                        mwM = cwM;
                        swM = k*15+2;
                    }

                    if (seq2->fwdGapEnds( i ) || i==sl2 && seq2->fwdGapContinues( i ))   // ..and another closes is seq2
                    {

                        cwM = yM->g(k,s-1);
                        if (cwM > mwM)
                        {
                            mwM = cwM;
                            swM = k*15+8;
                        }
                    }
                    if (seq2->fwdChildGapEnds( i ) || i==sl2 && seq2->fwdChildGapContinues( i ))   // ..and another closes is seq2 child
                    {

                        cwM = zM->g(k,s-1);
                        if (cwM > mwM)
                        {
                            mwM = cwM;
                            swM = k*15+14;
                        }
                    }
                }

                if (seq1->fwdChildGapContinues( j ))   // flagged gap continues in seq1
                {

                    cwX = wX->g(k,s-1);
                    if (cwX > mwX)
                    {
                        mwX = cwX;
                        swX = k*15+9;
                    }

                    cwM = wM->g(k,s-1);
                    if (cwM > mwM)
                    {
                        mwM = cwM;
                        swM = k*15+11;
                    }
                }

                if (seq1->fwdChildGapEnds( j )  || j==sl1 && seq1->fwdChildGapContinues( j ))   // flagged gap ends in seq1 child
                {

                    int l = hmm->transIndY(k,0);
                    while (l>=0)
                    {

                        cX = max(wX->g(l,s-1) + hmm->probXX(l,k),
                                 small,
                                 wM->g(l,s-1) + hmm->probMX(l,k)) + msr->indelX(k);
                        if (cX > mX)
                        {
                            mX = cX;
                            sX = l*15+9+maxIndex;
                        }

                        l = hmm->transIndY(k,l+1);
                    }
                }

                if (seq2->fwdGapEnds( i ) || i==sl2 && seq2->fwdGapContinues( i ) )   // flagged gap ends in seq2; X-gap goes right so earlier
                {

                    int l = hmm->transIndY(k,0);
                    while (l>=0)
                    {

                        cX = max(small,
                                 yY->g(l,s-1) + hmm->probYX(l,k),
                                 yM->g(l,s-1) + hmm->probMX(l,k)) + msr->indelX(k);
                        if (cX > mX)
                        {
                            mX = cX;
                            sX = l*15+6+maxIndex;
                        }

                        l = hmm->transIndY(k,l+1);
                    }
                }

                if (seq2->fwdChildGapEnds( i ) || i==sl2 && seq2->fwdChildGapContinues( i ))   // flagged gap ends in seq2 child; X-gap goes right so earlier
                {

                    int l = hmm->transIndY(k,0);
                    while (l>=0)
                    {

                        cX = max(small,
                                 zY->g(l,s-1) + hmm->probYX(l,k),
                                 zM->g(l,s-1) + hmm->probMX(l,k)) + msr->indelX(k);
                        if (cX > mX)
                        {
                            mX = cX;
                            sX = l*15+12+maxIndex;
                        }

                        l = hmm->transIndY(k,l+1);
                    }
                }

                int l = hmm->transIndY(k,0);
                while (l>=0)
                {

                    cX = max(vX->g(l,s-1) + hmm->probXX(l,k),
                             vY->g(l,s-1) + hmm->probYX(l,k),
                             vM->g(l,s-1) + hmm->probMX(l,k)) + msr->indelX(k);
                    if (cX > mX)
                    {
                        mX = cX;
                        sX = l*15+maxIndex;
                    }

                    l = hmm->transIndY(k,l+1);
                }

            }
            else if (move==1)
            {

                if (seq1->fwdGapEnds( j ) || j==sl1 && seq1->fwdGapContinues( j ))   // flagged gap ends in seq1; Y-gap goes down so earlier
                {

                    int l = hmm->transIndY(k,0);
                    while (l>=0)
                    {

                        cY = max(xX->g(l,s-1) + hmm->probXY(l,k),
                                 small,
                                 xM->g(l,s-1) + hmm->probMY(l,k)) + msr->indelY(k);
                        if (cY > mY)
                        {
                            mY = cY;
                            sY = l*15+3+maxIndex;
                        }

                        l = hmm->transIndY(k,l+1);
                    }
                }


                if (seq1->fwdChildGapEnds( j )  || j==sl1 && seq1->fwdChildGapContinues( j ))   // flagged gap ends in seq1 child; Y-gap goes down so earlier
                {

                    int l = hmm->transIndY(k,0);
                    while (l>=0)
                    {

                        cY = max(wX->g(l,s-1) + hmm->probXY(l,k),
                                 small,
                                 wM->g(l,s-1) + hmm->probMY(l,k)) + msr->indelY(k);
                        if (cY > mY)
                        {
                            mY = cY;
                            sY = l*15+9+maxIndex;
                        }

                        l = hmm->transIndY(k,l+1);
                    }
                }

                if (seq2->fwdGapStarts( i ))   // flagged gap starts in seq2
                {

                    cyY = vY->g(k,s-1);
                    if (cyY > myY)
                    {
                        myY = cyY;
                        syY = k*15+1;
                    }

                    cyM = vM->g(k,s-1);
                    if (cyM > myM)
                    {
                        myM = cyM;
                        syM = k*15+2;
                    }

                    if (seq1->fwdGapEnds( j )  || j==sl1 && seq1->fwdChildGapContinues( j ))   // .. and another closes in seq1
                    {

                        cyM = xM->g(k,s-1);
                        if (cyM > myM)
                        {
                            myM = cyM;
                            syM = k*15+5;
                        }
                    }
                    if (seq1->fwdChildGapEnds( j )  || j==sl1 && seq1->fwdChildGapContinues( j ))   // .. and another closes in seq1 child
                    {

                        cyM = wM->g(k,s-1);
                        if (cyM > myM)
                        {
                            myM = cyM;
                            syM = k*15+11;
                        }
                    }
                }

                if (seq2->fwdGapContinues( i ))   // flagged gap continues in seq2
                {

                    cyY = yY->g(k,s-1);
                    if (cyY > myY)
                    {
                        myY = cyY;
                        syY = k*15+7;
                    }

                    cyM = yM->g(k,s-1);
                    if (cyM > myM)
                    {
                        myM = cyM;
                        syM = k*15+8;
                    }
                }


                if (seq2->fwdGapEnds( i )  || i==sl2 && seq2->fwdGapContinues( i ) )   // flagged gap ends in seq2
                {

                    int l = hmm->transIndY(k,0);
                    while (l>=0)
                    {
                        cY = max(small,
                                 yY->g(l,s-1) + hmm->probYY(l,k),
                                 yM->g(l,s-1) + hmm->probMY(l,k)) + msr->indelY(k);
                        if (cY > mY)
                        {
                            mY = cY;
                            sY = l*15+6+maxIndex;
                        }

                        l = hmm->transIndY(k,l+1);
                    }
                }


                if (seq2->fwdChildGapStarts( i ))   // flagged gap starts in seq2 child
                {

                    czY = vY->g(k,s-1);
                    if (czY > mzY)
                    {
                        mzY = czY;
                        szY = k*15+1;
                    }

                    czM = vM->g(k,s-1);
                    if (czM > mzM)
                    {
                        mzM = czM;
                        szM = k*15+2;
                    }

                    if (seq1->fwdGapEnds( j )  || j==sl1 && seq1->fwdGapContinues( j ))   // .. and another closes in seq1
                    {

                        czM = xM->g(k,s-1);
                        if (czM > mzM)
                        {
                            mzM = czM;
                            szM = k*15+5;
                        }
                    }
                    if (seq1->fwdChildGapEnds( j ) || j==sl1 && seq1->fwdChildGapContinues( j ))   // .. and another closes in seq1 child
                    {

                        czM = wM->g(k,s-1);
                        if (czM > mzM)
                        {
                            mzM = czM;
                            szM = k*15+11;
                        }
                    }
                }

                if (seq2->fwdChildGapContinues( i ))   // flagged gap continues in seq2
                {

                    czY = zY->g(k,s-1);
                    if (czY > mzY)
                    {
                        mzY = czY;
                        szY = k*15+13;
                    }

                    czM = zM->g(k,s-1);
                    if (czM > mzM)
                    {
                        mzM = czM;
                        szM = k*15+14;
                    }
                }


                if (seq2->fwdChildGapEnds( i )  || i==sl2 && seq2->fwdChildGapContinues( i ))   // flagged gap ends in seq2 child
                {

                    int l = hmm->transIndY(k,0);
                    while (l>=0)
                    {
                        cY = max(small,
                                 zY->g(l,s-1) + hmm->probYY(l,k),
                                 zM->g(l,s-1) + hmm->probMY(l,k)) + msr->indelY(k);
                        if (cY > mY)
                        {
                            mY = cY;
                            sY = l*15+12+maxIndex;
                        }

                        l = hmm->transIndY(k,l+1);
                    }
                }

                int l = hmm->transIndY(k,0);
                while (l>=0)
                {

                    cY = max(vX->g(l,s-1) + hmm->probXY(l,k),
                             vY->g(l,s-1) + hmm->probYY(l,k),
                             vM->g(l,s-1) + hmm->probMY(l,k)) + msr->indelY(k);
                    if (cY > mY)
                    {
                        mY = cY;
                        sY = l*15+maxIndex;
                    }

                    l = hmm->transIndY(k,l+1);
                }

            }
            else if (move==2)
            {

                if (seq1->fwdGapEnds( j ) || j==sl1 && seq1->fwdGapContinues( j ))   // flagged gap ends in seq1
                {

                    int l = hmm->transIndY(k,0);
                    while (l>=0)
                    {

                        cM = max(xX->g(l,s-1) + hmm->probXM(l,k),
                                 small,
                                 xM->g(l,s-1) + hmm->probMM(l,k)) + msr->fwdM(k);
                        if (cM > mM)
                        {
                            mM = cM;
                            sM = l*15+3+maxIndex;
                        }

                        l = hmm->transIndY(k,l+1);
                    }
                }

                if (seq1->fwdChildGapEnds( j ) || j==sl1 && seq1->fwdChildGapContinues( j ))   // flagged gap ends in seq1 child
                {

                    int l = hmm->transIndY(k,0);
                    while (l>=0)
                    {

                        cM = max(wX->g(l,s-1) + hmm->probXM(l,k),
                                 small,
                                 wM->g(l,s-1) + hmm->probMM(l,k)) + msr->fwdM(k);
                        if (cM > mM)
                        {
                            mM = cM;
                            sM = l*15+9+maxIndex;
                        }

                        l = hmm->transIndY(k,l+1);
                    }
                }

                if (seq2->fwdGapEnds( i ) || i==sl2 && seq2->fwdGapContinues( i ))   // flagged gap ends in seq2
                {

                    int l = hmm->transIndY(k,0);
                    while (l>=0)
                    {

                        cM = max(small,
                                 yY->g(l,s-1) + hmm->probYM(l,k),
                                 yM->g(l,s-1) + hmm->probMM(l,k)) + msr->fwdM(k);
                        if (cM > mM)
                        {
                            mM = cM;
                            sM = l*15+6+maxIndex;
                        }

                        l = hmm->transIndY(k,l+1);
                    }
                }

                if (seq2->fwdChildGapEnds( i ) || i==sl2 && seq2->fwdChildGapContinues( i ))   // flagged gap ends in seq2 child
                {

                    int l = hmm->transIndY(k,0);
                    while (l>=0)
                    {
                        cM = max(small,
                                 zY->g(l,s-1) + hmm->probYM(l,k),
                                 zM->g(l,s-1) + hmm->probMM(l,k)) + msr->fwdM(k);
                        if (cM > mM)
                        {
                            mM = cM;
                            sM = l*15+12+maxIndex;
                        }

                        l = hmm->transIndY(k,l+1);
                    }
                }

                int l = hmm->transIndY(k,0);
                while (l>=0)
                {

                    cM = max(vX->g(l,s-1) + hmm->probXM(l,k),
                             vY->g(l,s-1) + hmm->probYM(l,k),
                             vM->g(l,s-1) + hmm->probMM(l,k)) + msr->fwdM(k);
                    if (cM > mM)
                    {
                        mM = cM;
                        sM = l*15+maxIndex;
                    }

                    l = hmm->transIndY(k,l+1);
                }

            }

            vX->s(mX,k,s);
            vY->s(mY,k,s);
            vM->s(mM,k,s);

            ptVX->s(sX,k,s);
            ptVY->s(sY,k,s);
            ptVM->s(sM,k,s);

            xX->s(mxX,k,s);
            xM->s(mxM,k,s);

            ptXX->s(sxX,k,s);
            ptXM->s(sxM,k,s);

            yY->s(myY,k,s);
            yM->s(myM,k,s);

            ptYY->s(syY,k,s);
            ptYM->s(syM,k,s);

            wX->s(mwX,k,s);
            wM->s(mwM,k,s);

            ptWX->s(swX,k,s);
            ptWM->s(swM,k,s);

            zY->s(mzY,k,s);
            zM->s(mzM,k,s);

            ptZY->s(szY,k,s);
            ptZM->s(szM,k,s);

//            cout<<s<<": ("<<j<<") "<<seq1->fwdGapStarts( j )<<" "<<seq1->fwdChildGapStarts( j )<<"; "<<seq1->fwdGapContinues( j )<<" "<<seq1->fwdChildGapContinues( j )<<"; "<<seq1->fwdGapEnds( j )<<" "<<seq1->fwdChildGapEnds( j )<<"; "<<seq1->fwdGapEndsNext( j )<<" "<<seq1->fwdChildGapEndsNext( j )<<" : ";
//            cout     <<"("<<i<<") "<<seq2->fwdGapStarts( i )<<" "<<seq2->fwdChildGapStarts( i )<<"; "<<seq2->fwdGapContinues( i )<<" "<<seq2->fwdChildGapContinues( i )<<"; "<<seq2->fwdGapEnds( i )<<" "<<seq2->fwdChildGapEnds( i )<<"; "<<seq2->fwdGapEndsNext( i )<<" "<<seq2->fwdChildGapEndsNext( i )<<"\n";
//            cout<<s<<": "<<move<<" : "<<mX<<" "<<mY<<" "<<mM<<"; "<<mxX<<" "<<mxM<<"; "<<myY<<" "<<myM<<"; "<<mwX<<" "<<mwM<<"; "<<mzY<<" "<<mzM<<"\n";
        }
        if (move==0)
        {
            if (j<sl1)
                j++;
        }
        else if (move==1)
        {
            if (i<sl2)
                i++;
        }
        else if (move==2)
        {
            if (j<sl1)
                j++;
            if (i<sl2)
                i++;
            countSites++;
        }
        else
        {
            cout<<"impossible pointer "<<move<<endl;
            exit(-1);
        }
        countSites++;

        mi++;
        if (mi!=path->end())
        {
            move = *mi;
        }
        else
        {
            break;
        }
// 		cout<<move<<endl;
    }
    s++;

    FOR(k,nState)
    {
        if (NOTGAP)
        {
            vX->s(vX->g(k,s-1)+hmm->structBgFreq(k),k,s); // no gap penalty for terminal gaps
            vY->s(vY->g(k,s-1)+hmm->structBgFreq(k),k,s);
        }
        else
        {
            vX->s(vX->g(k,s-1)+hmm->structBgFreq(k)+hmm->probXW(k),k,s);
            vY->s(vY->g(k,s-1)+hmm->structBgFreq(k)+hmm->probYW(k),k,s);
        }
        vM->s(vM->g(k,s-1)+hmm->structBgFreq(k)+hmm->probMW(k),k,s);
    }

    if (seq1->bwdGapStarts( sl1 ) || seq1->bwdGapContinues( sl1 ))
    {
        FOR(k,nState)
        {
            xX->s(xX->g(k,s-1)+hmm->structBgFreq(k),k,s);
            xM->s(xM->g(k,s-1)+hmm->structBgFreq(k)+hmm->probMW(k),k,s);
        }
    }
    else
    {
        FOR(k,nState)
        {
            xX->s(small,k,s);
            xM->s(small,k,s);
        }
    }
    if (seq2->bwdGapStarts( sl2 ) || seq2->bwdGapContinues( sl2 ))
    {
        FOR(k,nState)
        {
            yY->s(yY->g(k,s-1)+hmm->structBgFreq(k),k,s);
            yM->s(yM->g(k,s-1)+hmm->structBgFreq(k)+hmm->probMW(k),k,s);
        }
    }
    else
    {
        FOR(k,nState)
        {
            yY->s(small,k,s);
            yM->s(small,k,s);
        }
    }
    if (seq1->bwdChildGapStarts( sl1 ) || seq1->bwdChildGapContinues( sl1 ))
    {
        FOR(k,nState)
        {
            wX->s(wX->g(k,s-1)+hmm->structBgFreq(k),k,s);
            wM->s(wM->g(k,s-1)+hmm->structBgFreq(k)+hmm->probMW(k),k,s);
        }
    }
    else
    {
        FOR(k,nState)
        {
            wX->s(small,k,s);
            wM->s(small,k,s);
        }
    }
    if (seq2->bwdChildGapStarts( sl2 ) || seq2->bwdChildGapContinues( sl2 ))
    {
        FOR(k,nState)
        {
            zY->s(zY->g(k,s-1)+hmm->structBgFreq(k),k,s);
            zM->s(zM->g(k,s-1)+hmm->structBgFreq(k)+hmm->probMW(k),k,s);
        }
    }
    else
    {
        FOR(k,nState)
        {
            zY->s(small,k,s);
            zM->s(small,k,s);
        }
    }

    maxFullScore = small;
    int pointer = -1;
    FOR(k,nState)
    {
        if (vX->g(k,s)>maxFullScore)
        {
            maxFullScore=vX->g(k,s);
            pointer = k*15+0;
        }
        if (vY->g(k,s)>maxFullScore)
        {
            maxFullScore=vY->g(k,s);
            pointer = k*15+1;
        }
        if (vM->g(k,s)>maxFullScore)
        {
            maxFullScore=vM->g(k,s);
            pointer = k*15+2;
        }
        if (xX->g(k,s)>maxFullScore)
        {
            maxFullScore=xX->g(k,s);
            pointer = k*15+3;
        }
        if (xM->g(k,s)>maxFullScore)
        {
            maxFullScore=xM->g(k,s);
            pointer = k*15+5;
        }
        if (yY->g(k,s)>maxFullScore)
        {
            maxFullScore=yY->g(k,s);
            pointer = k*15+7;
        }
        if (yM->g(k,s)>maxFullScore)
        {
            maxFullScore=yM->g(k,s);
            pointer = k*15+8;
        }
        if (wX->g(k,s)>maxFullScore)
        {
            maxFullScore=wX->g(k,s);
            pointer = k*15+9;
        }
        if (wM->g(k,s)>maxFullScore)
        {
            maxFullScore=wM->g(k,s);
            pointer = k*15+11;
        }
        if (zY->g(k,s)>maxFullScore)
        {
            maxFullScore=zY->g(k,s);
            pointer = k*15+13;
        }
        if (zM->g(k,s)>maxFullScore)
        {
            maxFullScore=zM->g(k,s);
            pointer = k*15+14;
        }
    }

    countSites=0;

    beg = new Site(0);
    end = new Site(1);
    newsite = new Site();

    defineBegin();
    defineEnd();

    newsite->resetCounter();

    int proc = -1;
    int state = -1;
    j = sl1;
    i = sl2;
    s--;
    for (; s>0; s--)
    {
        proc = pointer/15;
        state = pointer%15;

        newsite->addNewSite();
        newsite->isAnchor(false);
        newsite->currModelState(proc);
        newsite->currMatchState(state);

//        cout<<"ns: "<<s<<"; "<<proc<<" "<<state<<"; ("<<i<<" "<<j<<"): "<<newsite->index()<<endl;
        if (state==0)
        {
            newsite->cInd1(j);
            newsite->cInd2(-1);
            newsite->nInd1(j);
            newsite->nInd2(i);
            newsite->nullSite(false);
            j--;
            pointer = ptVX->g(proc,s);
        }
        else if (state==1)
        {
            newsite->cInd1(-1);
            newsite->cInd2(i);
            newsite->nInd1(j);
            newsite->nInd2(i);
            newsite->nullSite(false);
            i--;
            pointer = ptVY->g(proc,s);
        }
        else if (state==2)
        {
            newsite->cInd1(j);
            newsite->cInd2(i);
            newsite->nInd1(j);
            newsite->nInd2(i);
            newsite->nullSite(false);
            i--;
            j--;
            pointer = ptVM->g(proc,s);
            countSites++;
        }
        else if (state==3)
        {
            newsite->cInd1(j);
            newsite->cInd2(-1);
            newsite->nInd1(-1);
            newsite->nInd2(-1);
            newsite->nullSite(true);
            j--;
            pointer = ptXX->g(proc,s);
        }
        else if (state==5)
        {
            newsite->cInd1(j);
            newsite->cInd2(-1);
            newsite->nInd1(-1);
            newsite->nInd2(-1);
            newsite->nullSite(true);
            j--;
            pointer = ptXM->g(proc,s);
        }
        else if (state==7)
        {
            newsite->cInd1(-1);
            newsite->cInd2(i);
            newsite->nInd1(-1);
            newsite->nInd2(-1);
            newsite->nullSite(true);
            i--;
            pointer = ptYY->g(proc,s);
        }
        else if (state==8)
        {
            newsite->cInd1(-1);
            newsite->cInd2(i);
            newsite->nInd1(-1);
            newsite->nInd2(-1);
            newsite->nullSite(true);
            i--;
            pointer = ptYM->g(proc,s);
        }
        else if (state==9)
        {
            newsite->cInd1(j);
            newsite->cInd2(-1);
            newsite->nInd1(-1);
            newsite->nInd2(-1);
            newsite->nullSite(true);
            j--;
            pointer = ptWX->g(proc,s);
        }
        else if (state==11)
        {
            newsite->cInd1(j);
            newsite->cInd2(-1);
            newsite->nInd1(-1);
            newsite->nInd2(-1);
            newsite->nullSite(true);
            j--;
            pointer = ptWM->g(proc,s);
        }
        else if (state==13)
        {
            newsite->cInd1(-1);
            newsite->cInd2(i);
            newsite->nInd1(-1);
            newsite->nInd2(-1);
            newsite->nullSite(true);
            i--;
            pointer = ptZY->g(proc,s);
        }
        else if (state==14)
        {
            newsite->cInd1(-1);
            newsite->cInd2(i);
            newsite->nInd1(-1);
            newsite->nInd2(-1);
            newsite->nullSite(true);
            i--;
            pointer = ptZM->g(proc,s);
        }
        else
        {
            delete beg;
            delete end;
            delete newsite;

            return false;

//            cout<<"something wrong"<<endl;
//            cout<<"ns: "<<s<<"; "<<proc<<" "<<state<<"; ("<<i<<" "<<j<<"): "<<newsite->index()<<endl;
//            exit(-1);
        }
        countSites++;



        newsite->setNeighbours(beg,end);
        end->prev();

        if (SCREEN && totalSites>0 && countSites%reportLimit==0 && verbose == true)
        {
            FOR(ii,message.length())
            {
                cout<<'\b';
            }

            char prop[10];
            sprintf(prop,": %i",countSites*100/totalSites);
            message = currentNode+prop+"% computed                    ";

            cout<<message;
            cout.flush();
        }

//  		cout<<pointer<<endl;
    }

    delete beg;
    delete end;
    delete newsite;

    return true;
}



void ReadAlignment::defineBegin()
{
    beg->index(0);
    beg->isAnchor(false);
    beg->currMatchState(-1);
    beg->currModelState(-1);
    beg->nullSite(false);

    beg->cInd1(0);
    beg->cInd2(0);
    beg->nInd1(0);
    beg->nInd2(0);
    beg->rInd1(0);
    beg->rInd2(-1); // before the start
    beg->lInd1(0);
    beg->lInd2(0);

    beg->vitf(small);
    beg->vitfM(-1);
    beg->vitfS(-1);

    beg->vitb(small);
    beg->vitbM(-1);
    beg->vitbS(-1);

}


void ReadAlignment::defineEnd()
{
    end->index(1);
    end->isAnchor(false);
    end->nullSite(true);

    end->cInd1(-1);
    end->nInd1(-1);
    end->cInd2(-1);
    end->nInd2(-1);
    end->rInd1(seq1->length());
    end->rInd2(seq2->length());
    end->lInd1(-1);
    end->lInd2(-1);

    end->vitf(small);
    end->vitfM(-1);
    end->vitfS(-1);

    end->vitb(small);
    end->vitbM(-1);
    end->vitbS(-1);

}


bool ReadAlignment::rndBool()
{
    if(REPRODUCIBLE)
        srand(random_seed);

    double p = (double)rand()/(double)RAND_MAX;
    if (p>0.5)
        return true;
    else
        return false;
}

int ReadAlignment::rndInt(int i)
{
    if(REPRODUCIBLE)
        srand(random_seed);

    return (int)(i*(rand()/(RAND_MAX+1.0)));
}


double ReadAlignment::max(double a,double b)
{
    if (a==small && b==small)
    {
        return a;
    }
    else if (a>b)
    {
        return a;
    }
    else if (a<b)
    {
        return b;
    }
    else
    {
        if (rndBool())
        {
            return a;
        }
        else
        {
            return b;
        }
    }
}

double ReadAlignment::max(double a,double b, double c)
{
    if (a==small && b==small && c==small)
    {
        maxIndex = 0;
        return a;
    }
    else if (a>b && a>c)
    {
        maxIndex = 0;
        return a;
    }
    else if (a<b && b>c)
    {
        maxIndex = 1;
        return b;
    }
    else if (a<c && b<c)
    {
        maxIndex = 2;
        return c;
    }
    else if (a>b && a==c)
    {
        if (rndBool())
        {
            maxIndex = 0;
            return a;
        }
        else
        {
            maxIndex = 2;
            return c;
        }
    }
    else if (a>c && a==b)
    {
        if (rndBool())
        {
            maxIndex = 0;
            return a;
        }
        else
        {
            maxIndex = 1;
            return b;
        }
    }
    else if (a<b && b==c)
    {
        if (rndBool())
        {
            maxIndex = 1;
            return b;
        }
        else
        {
            maxIndex = 2;
            return c;
        }
    }
    else
    {
        int i = rndInt(3);
        maxIndex = i;
        if (i==0 || i==3)
        {
            return a;
        }
        else if (i==1)
        {
            return b;
        }
        else if (i==2)
        {
            return c;
        }
        else
        {
            cout <<"ReadAlignment::random number error: i="<<i<<endl;
            exit(1);
        }
    }
}

