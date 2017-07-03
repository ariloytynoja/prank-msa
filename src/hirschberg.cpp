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
#include <sstream>
#include "config.h"
#include "exonerate_reads.h"
#include "hirschberg.h"

using namespace std;

extern string tmpNodeName;

Hirschberg::~Hirschberg()
{
    delete beg;
    delete end;
    delete newsite;

}

int Hirschberg::alignmentNumber = 0;

void Hirschberg::cleanUp()
{
    if (NOISE>1)
        cout<<"Hirschberg::cleanUp()"<<endl;

//     cout<<"cleanUp()"<<endl;

    // clean up

    delete fwdvX;
    delete fwdvY;
    delete fwdvM;
    delete bwdvX;
    delete bwdvY;
    delete bwdvM;

    delete fwdxX;
    delete fwdxM;
    delete bwdxX;
    delete bwdxM;

    delete fwdyY;
    delete fwdyM;
    delete bwdyY;
    delete bwdyM;

    delete fwdwX;
    delete fwdwM;
    delete bwdwX;
    delete bwdwM;

    delete fwdzY;
    delete fwdzM;
    delete bwdzY;
    delete bwdzM;

//
    delete fVM1;
    delete fVX1;
    delete fVY1;
    delete fXM1;
    delete fXX1;
    delete fWM1;
    delete fWX1;
    delete fYM1;
    delete fYY1;
    delete fZM1;
    delete fZY1;

    delete fVM2;
    delete fVX2;
    delete fVY2;
    delete fXM2;
    delete fXX2;
    delete fWM2;
    delete fWX2;
    delete fYM2;
    delete fYY2;
    delete fZM2;
    delete fZY2;

    delete bVM1;
    delete bVX1;
    delete bVY1;
    delete bXM1;
    delete bXX1;
    delete bWM1;
    delete bWX1;
    delete bYM1;
    delete bYY1;
    delete bZM1;
    delete bZY1;

    delete bVM2;
    delete bVX2;
    delete bVY2;
    delete bXM2;
    delete bXX2;
    delete bWM2;
    delete bWX2;
    delete bYM2;
    delete bYY2;
    delete bZM2;
    delete bZY2;

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


Hirschberg::Hirschberg()
{
    count = 2;
    small = -HUGE_VAL;

    beg = new Site(0);
    end = new Site(1);
    newsite = new Site();

    random_seed = (unsigned)time(NULL);
}

void Hirschberg::initialiseMatrices(int size)
{

    if (NOISE>1)
        cout<<"Hirschberg::initialiseMatrices("<<size<<")"<<endl;

    matrixSize = size;

    sAlpha = hmm->getASize();
    nState = hmm->getNStates();

    // Initialize matrices
    //
    fwdvX = new FlMatrix(nState,"fwdvX");
    fwdvY = new FlMatrix(nState,"fwdvY");
    fwdvM = new FlMatrix(nState,"fwdvM");
    bwdvX = new FlMatrix(nState,"bwdvX");
    bwdvY = new FlMatrix(nState,"bwdvY");
    bwdvM = new FlMatrix(nState,"bwdvM");

    fwdxX = new FlMatrix(nState,"fwdxX");
    fwdxM = new FlMatrix(nState,"fwdxM");
    bwdxX = new FlMatrix(nState,"bwdxX");
    bwdxM = new FlMatrix(nState,"bwdxM");

    fwdyY = new FlMatrix(nState,"fwdyY");
    fwdyM = new FlMatrix(nState,"fwdyM");
    bwdyY = new FlMatrix(nState,"bwdyY");
    bwdyM = new FlMatrix(nState,"bwdyM");

//
    fwdwX = new FlMatrix(nState,"fwdwX");
    fwdwM = new FlMatrix(nState,"fwdwM");
    bwdwX = new FlMatrix(nState,"bwdwX");
    bwdwM = new FlMatrix(nState,"bwdwM");

    fwdzY = new FlMatrix(nState,"fwdzY");
    fwdzM = new FlMatrix(nState,"fwdzM");
    bwdzY = new FlMatrix(nState,"bwdzY");
    bwdzM = new FlMatrix(nState,"bwdzM");

//
    fVX1 = new DbMatrix(nState,size,"fVX1");  // matrices for fwd Viterbi scores & skipped (X-gap, Y-gap) scores
    fVY1 = new DbMatrix(nState,size,"fVY1");
    fVM1 = new DbMatrix(nState,size,"fVM1");
    fXX1 = new DbMatrix(nState,size,"fXX1");
    fXM1 = new DbMatrix(nState,size,"fXM1");
    fWX1 = new DbMatrix(nState,size,"fWX1");
    fWM1 = new DbMatrix(nState,size,"fWM1");
    fYY1 = new DbMatrix(nState,size,"fYY1");
    fYM1 = new DbMatrix(nState,size,"fYM1");
    fZY1 = new DbMatrix(nState,size,"fZY1");
    fZM1 = new DbMatrix(nState,size,"fZM1");

    fVX1->initialise(small);
    fVY1->initialise(small);
    fVM1->initialise(small);
    fXX1->initialise(small);
    fXM1->initialise(small);
    fYY1->initialise(small);
    fYM1->initialise(small);
    fWX1->initialise(small);
    fWM1->initialise(small);
    fZY1->initialise(small);
    fZM1->initialise(small);

    fVX2 = new DbMatrix(nState,size,"fVX2");  // matrices for fwd Viterbi scores & skipped (X-gap, Y-gap) scores
    fVY2 = new DbMatrix(nState,size,"fVY2");
    fVM2 = new DbMatrix(nState,size,"fVM2");
    fXX2 = new DbMatrix(nState,size,"fXX2");
    fXM2 = new DbMatrix(nState,size,"fXM2");
    fWX2 = new DbMatrix(nState,size,"fWX2");
    fWM2 = new DbMatrix(nState,size,"fWM2");
    fYY2 = new DbMatrix(nState,size,"fYY2");
    fYM2 = new DbMatrix(nState,size,"fYM2");
    fZY2 = new DbMatrix(nState,size,"fZY2");
    fZM2 = new DbMatrix(nState,size,"fZM2");

    fVX2->initialise(small);
    fVY2->initialise(small);
    fVM2->initialise(small);
    fXX2->initialise(small);
    fXM2->initialise(small);
    fYY2->initialise(small);
    fYM2->initialise(small);
    fWX2->initialise(small);
    fWM2->initialise(small);
    fZY2->initialise(small);
    fZM2->initialise(small);

    bVM1 = new DbMatrix(nState,size,"bVM1");  // matrices for bwd Viterbi scores & skipped (X-gap, Y-gap) scores
    bVX1 = new DbMatrix(nState,size,"bVX1");
    bVY1 = new DbMatrix(nState,size,"bVY1");
    bXX1 = new DbMatrix(nState,size,"bXX1");
    bXM1 = new DbMatrix(nState,size,"bXM1");
    bWX1 = new DbMatrix(nState,size,"bWX1");
    bWM1 = new DbMatrix(nState,size,"bWM1");
    bYY1 = new DbMatrix(nState,size,"bYY1");
    bYM1 = new DbMatrix(nState,size,"bYM1");
    bZY1 = new DbMatrix(nState,size,"bZY1");
    bZM1 = new DbMatrix(nState,size,"bZM1");

    bVX1->initialise(small);
    bVY1->initialise(small);
    bVM1->initialise(small);
    bXX1->initialise(small);
    bXM1->initialise(small);
    bYY1->initialise(small);
    bYM1->initialise(small);
    bWX1->initialise(small);
    bWM1->initialise(small);
    bZY1->initialise(small);
    bZM1->initialise(small);

    bVM2 = new DbMatrix(nState,size,"bVM2");  // matrices for bwd Viterbi scores & skipped (X-gap, Y-gap) scores
    bVX2 = new DbMatrix(nState,size,"bVX2");
    bVY2 = new DbMatrix(nState,size,"bVY2");
    bXX2 = new DbMatrix(nState,size,"bXX2");
    bXM2 = new DbMatrix(nState,size,"bXM2");
    bWX2 = new DbMatrix(nState,size,"bWX2");
    bWM2 = new DbMatrix(nState,size,"bWM2");
    bYY2 = new DbMatrix(nState,size,"bYY2");
    bYM2 = new DbMatrix(nState,size,"bYM2");
    bZY2 = new DbMatrix(nState,size,"bZY2");
    bZM2 = new DbMatrix(nState,size,"bZM2");

    bVX2->initialise(small);
    bVY2->initialise(small);
    bVM2->initialise(small);
    bXX2->initialise(small);
    bXM2->initialise(small);
    bYY2->initialise(small);
    bYM2->initialise(small);
    bWX2->initialise(small);
    bWM2->initialise(small);
    bZY2->initialise(small);
    bZM2->initialise(small);

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

    ptVX->initialise(-1);
    ptVY->initialise(-1);
    ptVM->initialise(-1);
    ptXX->initialise(-1);
    ptXM->initialise(-1);
    ptYY->initialise(-1);
    ptYM->initialise(-1);
    ptWX->initialise(-1);
    ptWM->initialise(-1);
    ptZY->initialise(-1);
    ptZM->initialise(-1);

}

int Hirschberg::count = 2;
int Hirschberg::nState;
int Hirschberg::sAlpha;
int Hirschberg::matrixSize;

// Site* Hirschberg::beg;
// Site* Hirschberg::end;
// Site* Hirschberg::newsite;

FlMatrix* Hirschberg::fwdvX;
FlMatrix* Hirschberg::fwdvY;
FlMatrix* Hirschberg::fwdvM;
FlMatrix* Hirschberg::bwdvX;
FlMatrix* Hirschberg::bwdvY;
FlMatrix* Hirschberg::bwdvM;
FlMatrix* Hirschberg::fwdxX;
FlMatrix* Hirschberg::fwdxM;
FlMatrix* Hirschberg::bwdxX;
FlMatrix* Hirschberg::bwdxM;
FlMatrix* Hirschberg::fwdwX;
FlMatrix* Hirschberg::fwdwM;
FlMatrix* Hirschberg::bwdwX;
FlMatrix* Hirschberg::bwdwM;
FlMatrix* Hirschberg::fwdyY;
FlMatrix* Hirschberg::fwdyM;
FlMatrix* Hirschberg::bwdyY;
FlMatrix* Hirschberg::bwdyM;
FlMatrix* Hirschberg::fwdzY;
FlMatrix* Hirschberg::fwdzM;
FlMatrix* Hirschberg::bwdzY;
FlMatrix* Hirschberg::bwdzM;

DbMatrix* Hirschberg::fVM1;
DbMatrix* Hirschberg::fVX1;
DbMatrix* Hirschberg::fVY1;
DbMatrix* Hirschberg::fXM1;
DbMatrix* Hirschberg::fXX1;
DbMatrix* Hirschberg::fWM1;
DbMatrix* Hirschberg::fWX1;
DbMatrix* Hirschberg::fYM1;
DbMatrix* Hirschberg::fYY1;
DbMatrix* Hirschberg::fZM1;
DbMatrix* Hirschberg::fZY1;

DbMatrix* Hirschberg::fVM2;
DbMatrix* Hirschberg::fVX2;
DbMatrix* Hirschberg::fVY2;
DbMatrix* Hirschberg::fXM2;
DbMatrix* Hirschberg::fXX2;
DbMatrix* Hirschberg::fWM2;
DbMatrix* Hirschberg::fWX2;
DbMatrix* Hirschberg::fYM2;
DbMatrix* Hirschberg::fYY2;
DbMatrix* Hirschberg::fZM2;
DbMatrix* Hirschberg::fZY2;

DbMatrix* Hirschberg::bVM1;
DbMatrix* Hirschberg::bVX1;
DbMatrix* Hirschberg::bVY1;
DbMatrix* Hirschberg::bXM1;
DbMatrix* Hirschberg::bXX1;
DbMatrix* Hirschberg::bWM1;
DbMatrix* Hirschberg::bWX1;
DbMatrix* Hirschberg::bYM1;
DbMatrix* Hirschberg::bYY1;
DbMatrix* Hirschberg::bZM1;
DbMatrix* Hirschberg::bZY1;

DbMatrix* Hirschberg::bVM2;
DbMatrix* Hirschberg::bVX2;
DbMatrix* Hirschberg::bVY2;
DbMatrix* Hirschberg::bXM2;
DbMatrix* Hirschberg::bXX2;
DbMatrix* Hirschberg::bWM2;
DbMatrix* Hirschberg::bWX2;
DbMatrix* Hirschberg::bYM2;
DbMatrix* Hirschberg::bYY2;
DbMatrix* Hirschberg::bZM2;
DbMatrix* Hirschberg::bZY2;

IntMatrix* Hirschberg::ptVM;
IntMatrix* Hirschberg::ptVX;
IntMatrix* Hirschberg::ptVY;
IntMatrix* Hirschberg::ptXM;
IntMatrix* Hirschberg::ptXX;
IntMatrix* Hirschberg::ptWM;
IntMatrix* Hirschberg::ptWX;
IntMatrix* Hirschberg::ptYM;
IntMatrix* Hirschberg::ptYY;
IntMatrix* Hirschberg::ptZM;
IntMatrix* Hirschberg::ptZY;

void Hirschberg::alignSeqs(Sequence* s1,Sequence* s2,PhyloMatchScore *pms)
{

    alignmentNumber++;
    if (NOISE>0) cout<<"Alignment number: "<<alignmentNumber<<"."<<endl;

    seq1 = s1;
    seq2 = s2;

    sl1 = s1->length();
    sl2 = s2->length();

    totalSites = seq1->length()+seq2->length();
    countSites = 0;

    msr = pms;

    newsite->resetCounter();

    defineBegin();
    defineEnd();

    unsigned int ii = 0;
    if (SCREEN && totalSites>0 && verbose == true)
    {
        FOR( ii,message.length())
        {
            cout<<'\b';
        }

        char prop[10];
        sprintf(prop,": %i",countSites*100/totalSites);
        message = currentNode+prop+"% aligned                    ";

        cout<<message;
        cout.flush();
    }


    if (EXONERATE)
    {
        defineBegin();

        if (NOISE>0)
            cout<<"lengths: "<<sl1<<" and "<<sl2<<" ("<<anchSkipDist<<")"<<endl;

        if (sl1>anchSkipDist && sl2>anchSkipDist)
        {

            vector<hit> exonerate_hits;
            Exonerate_reads er;
            er.local_alignment(seq1->getMLsequence(),seq2->getMLsequence(),&exonerate_hits, true);

            vector<pair<int,int> > anchor_pairs;

            for (int i=0; i<exonerate_hits.size(); i++)
            {
                hit h = exonerate_hits.at(i);
                if (NOISE>1)
                    cout<<"e "<<h.query<<" "<<h.node<<" "<<h.score<<" "<<h.q_start<<" "<<h.q_end<<" "<<h.q_strand<<" "<<h.t_start<<" "<<h.t_end<<" "<<h.t_strand<<"\n";

                if(h.q_start < h.q_end && h.t_start < h.t_end)
                {
                    for (int j=0; j+h.q_start<h.q_end; j+=10)
                    {
                        if (j+h.q_start>5 && j+h.t_start>5)
                            anchor_pairs.push_back(make_pair(j+h.q_start,j+h.t_start));
                    }
                }
                else
                {
                    cout<<"\nAlignment anchoring indicates a reverse match: check the input data.\n";
                }
            }

            if (anchor_pairs.size()>0)
            {
                for (int i=0; i<anchor_pairs.size(); i++)
                {

                    if (NOISE>1)
                        cout<<" ex anchor "<<anchor_pairs.at(i).first<<","<<anchor_pairs.at(i).second<<" *"<<endl;

                    if ( SKIPGAPANCH && ( seq1->hasNeighborGaps(anchor_pairs.at(i).first) || seq2->hasNeighborGaps(anchor_pairs.at(i).second) ) )
                    {
                        if (NOISE>1)
                            cout<<"drop anchor "<<anchor_pairs.at(i).first<<","<<anchor_pairs.at(i).second<<" [gap]"<<endl;
                        continue;
                    }

                    defineESite(anchor_pairs.at(i).first,anchor_pairs.at(i).second);

                    if (NOISE>1)
                    {
                        cout<<" beg: "<<beg->index()<<" "<<beg->lInd1()<<" "<<beg->lInd2()<<" | ";
                        cout<<" anc: "<<end->index()<<" "<<end->lInd1()<<" "<<end->lInd2()<<endl;
                    }

                    if (end->lInd1() - beg->lInd1() > matrixSize)
                    {
                        cleanUp();
                        initialiseMatrices(end->lInd1() - beg->lInd1());
                    }

                    divideSeq();
                }
            }

            if (NOISE>1)
            {
                cout<<" beg: "<<beg->index()<<" "<<beg->lInd1()<<" "<<beg->lInd2()<<" | ";
                cout<<" end: "<<end->index()<<" "<<end->lInd1()<<" "<<end->lInd2()<<endl;
            }

            defineEnd();

            if (seq1->length()+1-beg->lInd1() > matrixSize)
            {
                if (nanch>0)
                    cleanUp();
                initialiseMatrices(seq1->length() + 1 - beg->lInd1());
            }



        }
        else
        {
            initialiseMatrices(sl1+1);

            defineEnd();
        }

        divideSeq();


        // plain alignment, no anchors
    }
    else
    {

// 		cout<<sl1<<" "<<sl2<<" "<<matrixSize<<endl;

        if (sl1+1>matrixSize)
        {
            cleanUp();
            initialiseMatrices((int)(((float)sl1+1)*initialMatrixSize));
        }

        if (sl2+1>matrixSize)
        {
            cleanUp();
            initialiseMatrices((int)(((float)sl2+1)*initialMatrixSize));
        }

        defineBegin();
        defineEnd();

        divideSeq();

    }

}


void Hirschberg::defineBegin()
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

void Hirschberg::defineSite(int idx)
{
    end->index(1);
    end->isAnchor(true);
    end->nullSite(true);

    end->cInd1(anchors->g(0,idx));
    end->nInd1(anchors->g(0,idx));
    end->lInd1(anchors->g(0,idx));
    end->cInd2(anchors->g(1,idx));
    end->nInd2(anchors->g(1,idx));
    end->lInd2(anchors->g(1,idx));
    end->rInd1(anchors->g(0,idx)-1);
    end->rInd2(anchors->g(1,idx)-1);

    end->vitf(small);
    end->vitfM(-1);
    end->vitfS(-1);

    end->vitb(small);
    end->vitbM(-1);
    end->vitbS(-1);
}

void Hirschberg::defineESite(int l,int r)
{
    end->index(1);
    end->isAnchor(true);
    end->nullSite(true);

    end->cInd1(l);
    end->nInd1(l);
    end->lInd1(l);
    end->cInd2(r);
    end->nInd2(r);
    end->lInd2(r);
    end->rInd1(l-1);
    end->rInd2(r-1);

    end->vitf(small);
    end->vitfM(-1);
    end->vitfS(-1);

    end->vitb(small);
    end->vitbM(-1);
    end->vitbS(-1);
}

void Hirschberg::defineEnd()
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

void Hirschberg::divideSeq()
{

    fwdvX->initialise(small);
    fwdvY->initialise(small);
    fwdvM->initialise(small);
    bwdvX->initialise(small);
    bwdvY->initialise(small);
    bwdvM->initialise(small);

    fwdxX->initialise(small);
    fwdxM->initialise(small);
    bwdxX->initialise(small);
    bwdxM->initialise(small);
    fwdyY->initialise(small);
    fwdyM->initialise(small);
    bwdyY->initialise(small);
    bwdyM->initialise(small);

    fwdwX->initialise(small);
    fwdwM->initialise(small);
    bwdwX->initialise(small);
    bwdwM->initialise(small);
    fwdzY->initialise(small);
    fwdzM->initialise(small);
    bwdzY->initialise(small);
    bwdzM->initialise(small);

    if (beg->index() == 0)
    {
        FOR(k,nState)
        {
            if (NOTGAP)
            {
                fwdvX->s(hmm->structBgFreq(k),k);
                fwdvY->s(hmm->structBgFreq(k),k);
            }
            else
            {
                fwdvX->s(hmm->structBgFreq(k)+hmm->probWX(k),k);
                fwdvY->s(hmm->structBgFreq(k)+hmm->probWY(k),k);
            }
            fwdvM->s(hmm->structBgFreq(k)+hmm->probWM(k),k);
        }
    }
    else
    {
        if (beg->vitfM() == 0)
            fwdvX->s(beg->vitf(),beg->vitfS());
        else if (beg->vitfM() == 1)
            fwdvY->s(beg->vitf(),beg->vitfS());
        else if (beg->vitfM() == 2)
            fwdvM->s(beg->vitf(),beg->vitfS());
        else if (beg->vitfM() == 3)
            fwdxX->s(beg->vitf(),beg->vitfS());
        else if (beg->vitfM() == 5)
            fwdxM->s(beg->vitf(),beg->vitfS());
        else if (beg->vitfM() == 7)
            fwdyY->s(beg->vitf(),beg->vitfS());
        else if (beg->vitfM() == 8)
            fwdyM->s(beg->vitf(),beg->vitfS());
        else if (beg->vitfM() == 9)
            fwdwX->s(beg->vitf(),beg->vitfS());
        else if (beg->vitfM() == 11)
            fwdwM->s(beg->vitf(),beg->vitfS());
        else if (beg->vitfM() == 13)
            fwdzY->s(beg->vitf(),beg->vitfS());
        else if (beg->vitfM() == 14)
            fwdzM->s(beg->vitf(),beg->vitfS());
        else
        {
            cout<<"hirschberg initialisation: impossible fwd state '"<<beg->vitfM()<<"'"<<endl;
        }
    }

    if (end->index() == 1 && !end->isAnchor())
    {
        FOR(k,nState)
        {
            if (NOTGAP)
            {
                bwdvX->s(hmm->structBgFreq(k),k); // no gap penalty for terminal gaps
                bwdvY->s(hmm->structBgFreq(k),k);
            }
            else
            {
                bwdvX->s(hmm->structBgFreq(k)+hmm->probXW(k),k);
                bwdvY->s(hmm->structBgFreq(k)+hmm->probYW(k),k);
            }
            bwdvM->s(hmm->structBgFreq(k)+hmm->probMW(k),k);
        }

        if (seq1->bwdGapStarts( sl1 ) || seq1->bwdGapContinues( sl1 ))
        {
            FOR(k,nState)
            {
                bwdxX->s(hmm->structBgFreq(k),k);
                bwdxM->s(hmm->structBgFreq(k)+hmm->probMW(k),k);
            }
        }
        if (seq2->bwdGapStarts( sl2 ) || seq2->bwdGapContinues( sl2 ))
        {
            FOR(k,nState)
            {
                bwdyY->s(hmm->structBgFreq(k),k);
                bwdyM->s(hmm->structBgFreq(k)+hmm->probMW(k),k);
            }
        }
        if (seq1->bwdChildGapStarts( sl1 ) || seq1->bwdChildGapContinues( sl1 ))
        {
            FOR(k,nState)
            {
                bwdwX->s(hmm->structBgFreq(k),k);
                bwdwM->s(hmm->structBgFreq(k)+hmm->probMW(k),k);
            }
        }
        if (seq2->bwdChildGapStarts( sl2 ) || seq2->bwdChildGapContinues( sl2 ))
        {
            FOR(k,nState)
            {
                bwdzY->s(hmm->structBgFreq(k),k);
                bwdzM->s(hmm->structBgFreq(k)+hmm->probMW(k),k);
            }
        }
    }
    else if (end->isAnchor())
    {
        bwdvX->initialise(0);
        bwdvY->initialise(0);
        bwdvM->initialise(0);

        if (seq1->bwdGapStarts( end->cInd1()-1 ) )
        {
            bwdxX->initialise(0);
            bwdxM->initialise(0);
        }
        if (seq2->bwdGapStarts( end->cInd2()-1 ) )
        {
            bwdyY->initialise(0);
            bwdyM->initialise(0);
        }
        if (seq1->bwdChildGapStarts( end->cInd1()-1 ) )
        {
            bwdwX->initialise(0);
            bwdwM->initialise(0);
        }
        if (seq2->bwdChildGapStarts( end->cInd2()-1 ) )
        {
            bwdzY->initialise(0);
            bwdzM->initialise(0);
        }
    }
    else
    {
        if (end->vitbM() == 0)
            bwdvX->s(end->vitb(),end->vitbS());
        else if (end->vitbM() == 1)
            bwdvY->s(end->vitb(),end->vitbS());
        else if (end->vitbM() == 2)
            bwdvM->s(end->vitb(), end->vitbS());
        else if (end->vitbM() == 3)
            bwdxX->s(end->vitb(),end->vitbS());
        else if (end->vitbM() == 5)
            bwdxM->s(end->vitb(),end->vitbS());
        else if (end->vitbM() == 7)
            bwdyY->s(end->vitb(),end->vitbS());
        else if (end->vitbM() == 8)
            bwdyM->s(end->vitb(),end->vitbS());
        else if (end->vitbM() == 9)
            bwdwX->s(end->vitb(),end->vitbS());
        else if (end->vitbM() == 11)
            bwdwM->s(end->vitb(),end->vitbS());
        else if (end->vitbM() == 13)
            bwdzY->s(end->vitb(),end->vitbS());
        else if (end->vitbM() == 14)
            bwdzM->s(end->vitb(),end->vitbS());
        else
        {
            cout<<"hirschberg initialisation: impossible bwd state '"<<end->vitbM()<<"'"<<endl;
        }
    }

    getMidSite(beg->lInd1(),end->rInd1(),beg->lInd2(),end->rInd2());

    if (newsite->index()%reportLimit==0)
    {
        if (SCREEN && verbose == true)
        {
            unsigned int ii;
            FOR(ii,message.length())
            {
                cout<<'\b';
            }

            char prop[10];
            sprintf(prop,": %i",countSites*100/totalSites);
            message = currentNode+prop+"% aligned                    ";

            cout<<message;
            cout.flush();
        }
        else if (NOISE>0)
        {
            cout<<currentNode<<": "<<countSites*100/totalSites<<"% aligned"<<endl;
        }
    }

    if (end->isAnchor() && (end->rInd1() - newsite->lInd1() < anchDropDist || end->rInd2() - newsite->lInd2() < anchDropDist))
    {

        if (NOISE>1)
            cout<<"new site "<<newsite->lInd1()<<","<<newsite->lInd2()<<"; drop anchor ("<<beg->lInd1()<<","<<end->rInd1()<<"; "<<beg->lInd2()<<","<<end->rInd2()<<")"<<endl;
        newsite->deleteLast();

    }
    else
    {

        newsite->setNeighbours(beg,end);
// 		cout<<newsite->index()<<" "<<newsite->cInd1()<<" "<<newsite->cInd2()<<" "<<newsite->getLSite()<<" "<<newsite->getRSite()<<endl;

        beg->next();

        // do right loop
        if (beg->index()!=0)
        {
            end->prev();
            beg->index(end->getLSite());

            // seqs still have chars on right
            if ( beg->lInd1() < end->rInd1() || beg->lInd2() < end->rInd2() )
            {
                divideSeq();
            }

            beg->index(end->index());
            end->next();
        }

        // do left loop
        if (beg->lInd1() < end->rInd1() || beg->lInd2() < end->rInd2() )
        {
            divideSeq();
        }

        beg->index(end->getLSite());
        end->index(beg->getRSite());

    }
}


void Hirschberg::getMidSite(int s1,int e1,int s2,int e2)
{

    int h = (s2+e2)/2+1;  // midpoint
    if (s2==e2)           // exception for zero-length seq2
        h = (s2+e2)/2;

    int s1Beg = s1;        // seq1 start site
    int s1Len = e1-s1;     // seq1 length

//    int s2Len = e2-s2;     // seq2 length

    int s2bBeg = s2;       // seq2_begin start site
    int s2bLen = h-s2;     // seq2_begin length
    if (s2==e2)
        s2bLen=0;

    int s2eBeg = h;        // seq2_end start site
    int s2eLen = e2-(h+1); // seq2_end length

    mLen = s1Len+1;        // short cuts
    mSize = nState*mLen;

    // Define pointers for current & previous row
    //
    cfVX = fVX1;
    cfVY = fVY1;
    cfVM = fVM1;
    cfXX = fXX1;
    cfXM = fXM1;
    cfYY = fYY1;
    cfYM = fYM1;
    cfWX = fWX1;
    cfWM = fWM1;
    cfZY = fZY1;
    cfZM = fZM1;

    pVX = fVX2;
    pVY = fVY2;
    pVM = fVM2;
    pXX = fXX2;
    pXM = fXM2;
    pYY = fYY2;
    pYM = fYM2;
    pWX = fWX2;
    pWM = fWM2;
    pZY = fZY2;
    pZM = fZM2;


    // A loop through the first half of seq2
    //
    FOR(i,s2bLen+1)
    {

        // A loop through seq1
        //
        FOR(j,mLen)
        {

            // Starting: set the corner values
            //
            if (i==0 && j==0 )
            {

                FOR(k,nState)
                {

                    // set starting values
                    cfVM->s(fwdvM->g(k),k,j);
                    cfVX->s(fwdvX->g(k),k,j);
                    cfVY->s(fwdvY->g(k),k,j);

                    cfXM->s(fwdxM->g(k),k,j);
                    cfXX->s(fwdxX->g(k),k,j);
                    cfYM->s(fwdyM->g(k),k,j);
                    cfYY->s(fwdyY->g(k),k,j);

                    cfWM->s(fwdwM->g(k),k,j);
                    cfWX->s(fwdwX->g(k),k,j);
                    cfZM->s(fwdzM->g(k),k,j);
                    cfZY->s(fwdzY->g(k),k,j);

                    ptVX->s(0,k,j);
                    ptVY->s(1,k,j);
                    ptVM->s(2,k,j);

                    ptXX->s(3,k,j);
                    ptXM->s(5,k,j);
                    ptYY->s(7,k,j);
                    ptYM->s(8,k,j);

                    ptWX->s(9,k,j);
                    ptWM->s(11,k,j);
                    ptZY->s(13,k,j);
                    ptZM->s(14,k,j);

                }
                continue;
            }


            // Compute the substitution prices
            //
            msr->computeFwd( s1Beg + j, s2bBeg + i );

            if (i==0 && j>0)   // only X-gaps are possible
            {


                FOR(k,nState)
                {

                    sX=sY=sM=sxX=sxM=syY=syM=swX=swM=szY=szM=-1;
                    mX=mY=mM=mxX=mxM=myY=myM=mwX=mwM=mzY=mzM=small;
                    cX=cY=cM=cxX=cxM=cyY=cyM=cwX=cwM=czY=czM=small;

                    if (seq1->fwdGapStarts( s1Beg + j ))   // flagged gap starts in seq1
                    {

                        cxX = cfVX->g(k,j-1);
                        if (cxX > mxX)
                        {
                            mxX = cxX;
                            sxX = k*15+0;
                        }

                        cxM = cfVM->g(k,j-1);
                        if (cxM > mxM)
                        {
                            mxM = cxM;
                            sxM = k*15+2;
                        }

                        if (seq2->fwdGapEndsNext( s2bBeg + i ))   // ..and another closes in seq2
                        {

                            cxM = cfYM->g(k,j-1);
                            if (cxM > mxM)
                            {
                                mxM = cxM;
                                sxM = k*15+8;
                            }
                        }
                        if (seq2->fwdChildGapEndsNext( s2bBeg + i ))   // ..and another closes in seq2 child
                        {

                            cxM = cfZM->g(k,j-1);
                            if (cxM > mxM)
                            {
                                mxM = cxM;
                                sxM = k*15+14;
                            }
                        }
                    }

                    if (seq1->fwdGapContinues( s1Beg + j ))   // flagged gap continues in seq1
                    {

                        cxX = cfXX->g(k,j-1);
                        if (cxX > mxX)
                        {
                            mxX = cxX;
                            sxX = k*15+3;
                        }

                        cxM = cfXM->g(k,j-1);
                        if (cxM > mxM)
                        {
                            mxM = cxM;
                            sxM = k*15+5;
                        }
                    }

                    if (seq1->fwdGapEnds( s1Beg + j ))     // flagged gap ends in seq1
                    {

                        int l = hmm->transIndY(k,0);
                        while (l>=0)
                        {

                            cX = max(cfXX->g(l,j-1) + hmm->probXX(l,k),
                                     small,
                                     cfXM->g(l,j-1) + hmm->probMX(l,k)) + msr->indelX(k);
                            if (cX > mX)
                            {
                                mX = cX;
                                sX = l*15+3+maxIndex;
                            }

                            l = hmm->transIndY(k,l+1);
                        }
                    }


                    if (seq1->fwdChildGapStarts( s1Beg + j ))   // flagged gap starts in seq1 child
                    {

                        cwX = cfVX->g(k,j-1);
                        if (cwX > mwX)
                        {
                            mwX = cwX;
                            swX = k*15+0;
                        }

                        cwM = cfVM->g(k,j-1);
                        if (cwM > mwM)
                        {
                            mwM = cwM;
                            swM = k*15+2;
                        }

                        if (seq2->fwdGapEndsNext( s2bBeg + i ))   // ..and another closes in seq2
                        {

                            cwM = cfYM->g(k,j-1);
                            if (cwM > mwM)
                            {
                                mwM = cwM;
                                swM = k*15+8;
                            }
                        }
                        if (seq2->fwdChildGapEndsNext( s2bBeg + i ))   // ..and another closes in seq2 child
                        {

                            cwM = cfZM->g(k,j-1);
                            if (cwM > mwM)
                            {
                                mwM = cwM;
                                swM = k*15+14;
                            }
                        }
                    }

                    if (seq1->fwdChildGapContinues( s1Beg + j ))   // flagged gap continues in seq1 child
                    {

                        cwX = cfWX->g(k,j-1);
                        if (cwX > mwX)
                        {
                            mwX = cwX;
                            swX = k*15+9;
                        }

                        cwM = cfWM->g(k,j-1);
                        if (cwM > mwM)
                        {
                            mwM = cwM;
                            swM = k*15+11;
                        }
                    }

                    if (seq1->fwdChildGapEnds( s1Beg + j ))     // flagged gap ends in seq1 child
                    {

                        int l = hmm->transIndY(k,0);
                        while (l>=0)
                        {

                            cX = max(cfWX->g(l,j-1) + hmm->probXX(l,k),
                                     small,
                                     cfWM->g(l,j-1) + hmm->probMX(l,k)) + msr->indelX(k);
                            if (cX > mX)
                            {
                                mX = cX;
                                sX = l*15+9+maxIndex;
                            }

                            l = hmm->transIndY(k,l+1);
                        }
                    }

                    if (seq2->fwdGapEndsNext( s2bBeg + i ))   // flagged gap ends in seq2
                    {

                        int l = hmm->transIndY(k,0);
                        while (l>=0)
                        {

                            cX = max(small,
                                     cfYY->g(l,j-1) + hmm->probYX(l,k),
                                     cfYM->g(l,j-1) + hmm->probMX(l,k)) + msr->indelX(k);
                            if (cX > mX)
                            {
                                mX = cX;
                                sX = l*15+6+maxIndex;
                            }

                            l = hmm->transIndY(k,l+1);
                        }
                    }

                    if (seq2->fwdChildGapEndsNext( s2bBeg + i ))   // flagged gap ends in seq2 child
                    {

                        int l = hmm->transIndY(k,0);
                        while (l>=0)
                        {

                            cX = max(small,
                                     cfZY->g(l,j-1) + hmm->probYX(l,k),
                                     cfZM->g(l,j-1) + hmm->probMX(l,k)) + msr->indelX(k);
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

                        cX = max(cfVX->g(l,j-1) + hmm->probXX(l,k),
                                 cfVY->g(l,j-1) + hmm->probYX(l,k),
                                 cfVM->g(l,j-1) + hmm->probMX(l,k)) + msr->indelX(k);
                        if (cX > mX)
                        {
                            mX = cX;
                            sX = l*15+maxIndex;
                        }
                        l = hmm->transIndY(k,l+1);
                    }

                    cfVX->s(mX,k,j);
                    cfVY->s(mY,k,j);
                    cfVM->s(mM,k,j);

                    ptVX->s(sX,k,j);
                    ptVY->s(sY,k,j);
                    ptVM->s(sM,k,j);

                    cfXX->s(mxX,k,j);
                    cfXM->s(mxM,k,j);

                    ptXX->s(sxX,k,j);
                    ptXM->s(sxM,k,j);

                    cfYY->s(myY,k,j);
                    cfYM->s(myM,k,j);

                    ptYY->s(syY,k,j);
                    ptYM->s(syM,k,j);

                    cfWX->s(mwX,k,j);
                    cfWM->s(mwM,k,j);

                    ptWX->s(swX,k,j);
                    ptWM->s(swM,k,j);

                    cfZY->s(mzY,k,j);
                    cfZM->s(mzM,k,j);

                    ptZY->s(szY,k,j);
                    ptZM->s(szM,k,j);
                }

            }
            else if (i>0 && j==0)   // only Y-gaps are possible
            {

                FOR(k,nState)
                {

                    sX=sY=sM=sxX=sxM=syY=syM=swX=swM=szY=szM=-1;
                    mX=mY=mM=mxX=mxM=myY=myM=mwX=mwM=mzY=mzM=small;
                    cX=cY=cM=cxX=cxM=cyY=cyM=cwX=cwM=czY=czM=small;

                    if (seq2->fwdGapStarts( s2bBeg + i ))      // flagged gap starts in seq2
                    {

                        cyY = pVY->g(k,j);
                        if (cyY > myY)
                        {
                            myY = cyY;
                            syY = k*15+1;
                        }

                        cyM = pVM->g(k,j);
                        if (cyM > myM)
                        {
                            myM = cyM;
                            syM = k*15+2;
                        }

                        if (seq1->fwdGapEndsNext( s1Beg + j ))   // .. and another closes in seq1
                        {

                            cyM = pXM->g(k,j);
                            if (cyM > myM)
                            {
                                myM = cyM;
                                syM = k*15+5;
                            }
                        }
                        if (seq1->fwdChildGapEndsNext( s1Beg + j ))   // .. and another closes in seq1 child
                        {

                            cyM = pWM->g(k,j);
                            if (cyM > myM)
                            {
                                myM = cyM;
                                syM = k*15+11;
                            }
                        }
                    }

                    if (seq2->fwdGapContinues( s2bBeg + i ))   // flagged gap continues in seq2
                    {

                        cyY = pYY->g(k,j);
                        if (cyY > myY)
                        {
                            myY = cyY;
                            syY = k*15+7;
                        }

                        cyM = pYM->g(k,j);
                        if (cyM > myM)
                        {
                            myM = cyM;
                            syM = k*15+8;
                        }
                    }

                    if (seq2->fwdGapEnds( s2bBeg + i ))   // flagged gap ends in seq2
                    {

                        int l = hmm->transIndY(k,0);
                        while (l>=0)
                        {

                            cY = max(small,
                                     pYY->g(l,j) + hmm->probYY(l,k),
                                     pYM->g(l,j) + hmm->probMY(l,k)) + msr->indelY(k);
                            if (cY > mY)
                            {
                                mY = cY;
                                sY = l*15+6+maxIndex;
                            }

                            l = hmm->transIndY(k,l+1);
                        }
                    }


                    if (seq2->fwdChildGapStarts( s2bBeg + i ))      // flagged gap starts in seq2 child
                    {

                        czY = pVY->g(k,j);
                        if (czY > mzY)
                        {
                            mzY = czY;
                            szY = k*15+1;
                        }

                        czM = pVM->g(k,j);
                        if (czM > mzM)
                        {
                            mzM = czM;
                            szM = k*15+2;
                        }

                        if (seq1->fwdGapEndsNext( s1Beg + j ))   // .. and another closes in seq1
                        {

                            czM = pXM->g(k,j);
                            if (czM > mzM)
                            {
                                mzM = czM;
                                szM = k*15+5;
                            }
                        }
                        if (seq1->fwdChildGapEndsNext( s1Beg + j ))   // .. and another closes in seq1 child
                        {

                            czM = pWM->g(k,j);
                            if (czM > mzM)
                            {
                                mzM = czM;
                                szM = k*15+11;
                            }
                        }
                    }

                    if (seq2->fwdChildGapContinues( s2bBeg + i ))   // flagged gap continues in seq2 child
                    {

                        czY = pZY->g(k,j);
                        if (czY > mzY)
                        {
                            mzY = czY;
                            szY = k*15+13;
                        }

                        czM = pZM->g(k,j);
                        if (czM > mzM)
                        {
                            mzM = czM;
                            szM = k*15+14;
                        }
                    }

                    if (seq2->fwdChildGapEnds( s2bBeg + i ))   // flagged gap ends in seq2 child
                    {

                        int l = hmm->transIndY(k,0);
                        while (l>=0)
                        {

                            cY = max(small,
                                     pZY->g(l,j) + hmm->probYY(l,k),
                                     pZM->g(l,j) + hmm->probMY(l,k)) + msr->indelY(k);
                            if (cY > mY)
                            {
                                mY = cY;
                                sY = l*15+12+maxIndex;
                            }

                            l = hmm->transIndY(k,l+1);
                        }
                    }

                    if (seq1->fwdGapEndsNext( s1Beg + j ))   // flagged gap ends in seq1
                    {

                        int l = hmm->transIndY(k,0);
                        while (l>=0)
                        {

                            cY = max(pXX->g(l,j) + hmm->probXY(l,k),
                                     small,
                                     pXM->g(l,j) + hmm->probMY(l,k)) + msr->indelY(k);
                            if (cY > mY)
                            {
                                mY = cY;
                                sY = l*15+3+maxIndex;
                            }

                            l = hmm->transIndY(k,l+1);
                        }
                    }

                    if (seq1->fwdChildGapEndsNext( s1Beg + j ))   // flagged gap ends in seq1
                    {

                        int l = hmm->transIndY(k,0);
                        while (l>=0)
                        {

                            cY = max(pWX->g(l,j) + hmm->probXY(l,k),
                                     small,
                                     pWM->g(l,j) + hmm->probMY(l,k)) + msr->indelY(k);
                            if (cY > mY)
                            {
                                mY = cY;
                                sY = l*15+9+maxIndex;
                            }

                            l = hmm->transIndY(k,l+1);
                        }
                    }

                    int l = hmm->transIndY(k,0);
                    while (l>=0)
                    {

                        cY = max(pVX->g(l,j) + hmm->probXY(l,k),
                                 pVY->g(l,j) + hmm->probYY(l,k),
                                 pVM->g(l,j) + hmm->probMY(l,k)) + msr->indelY(k);
                        if (cY > mY)
                        {
                            mY = cY;
                            sY = l*15+maxIndex;
                        }

                        l = hmm->transIndY(k,l+1);
                    }

                    cfVX->s(mX,k,j);
                    cfVY->s(mY,k,j);
                    cfVM->s(mM,k,j);

                    ptVX->s(sX,k,j);
                    ptVY->s(sY,k,j);
                    ptVM->s(sM,k,j);

                    cfXX->s(mxX,k,j);
                    cfXM->s(mxM,k,j);

                    ptXX->s(sxX,k,j);
                    ptXM->s(sxM,k,j);

                    cfYY->s(myY,k,j);
                    cfYM->s(myM,k,j);

                    ptYY->s(syY,k,j);
                    ptYM->s(syM,k,j);

                    cfWX->s(mwX,k,j);
                    cfWM->s(mwM,k,j);

                    ptWX->s(swX,k,j);
                    ptWM->s(swM,k,j);

                    cfZY->s(mzY,k,j);
                    cfZM->s(mzM,k,j);

                    ptZY->s(szY,k,j);
                    ptZM->s(szM,k,j);
                }

            }
            else    // so far, the moves have been exceptional; from now on they are "normal"
            {

                FOR(k,nState)
                {

                    sX=sY=sM=sxX=sxM=syY=syM=swX=swM=szY=szM=-1;
                    mX=mY=mM=mxX=mxM=myY=myM=mwX=mwM=mzY=mzM=small;
                    cX=cY=cM=cxX=cxM=cyY=cyM=cwX=cwM=czY=czM=small;

                    if (seq1->fwdGapStarts( s1Beg + j ))   // flagged gap starts in seq1
                    {

                        cxX = cfVX->g(k,j-1);
                        if (cxX > mxX)
                        {
                            mxX = cxX;
                            sxX = k*15+0;
                        }

                        cxM = cfVM->g(k,j-1);
                        if (cxM > mxM)
                        {
                            mxM = cxM;
                            sxM = k*15+2;
                        }

                        if (seq2->fwdGapEndsNext( s2bBeg + i ))   // ..and another closes is seq2
                        {

                            cxM = cfYM->g(k,j-1);
                            if (cxM > mxM)
                            {
                                mxM = cxM;
                                sxM = k*15+8;
                            }
                        }
                        if (seq2->fwdChildGapEndsNext( s2bBeg + i ))   // ..and another closes is seq2 child
                        {

                            cxM = cfZM->g(k,j-1);
                            if (cxM > mxM)
                            {
                                mxM = cxM;
                                sxM = k*15+14;
                            }
                        }
                    }

                    if (seq1->fwdGapContinues( s1Beg + j ))   // flagged gap continues in seq1
                    {

                        cxX = cfXX->g(k,j-1);
                        if (cxX > mxX)
                        {
                            mxX = cxX;
                            sxX = k*15+3;
                        }

                        cxM = cfXM->g(k,j-1);
                        if (cxM > mxM)
                        {
                            mxM = cxM;
                            sxM = k*15+5;
                        }
                    }

                    if (seq1->fwdGapEnds( s1Beg + j ))   // flagged gap ends in seq1
                    {

                        int l = hmm->transIndY(k,0);
                        while (l>=0)
                        {

                            cX = max(cfXX->g(l,j-1) + hmm->probXX(l,k),
                                     small,
                                     cfXM->g(l,j-1) + hmm->probMX(l,k)) + msr->indelX(k);
                            if (cX > mX)
                            {
                                mX = cX;
                                sX = l*15+3+maxIndex;
                            }

                            cM = max(pXX->g(l,j-1) + hmm->probXM(l,k),
                                     small,
                                     pXM->g(l,j-1) + hmm->probMM(l,k)) + msr->fwdM(k);
                            if (cM > mM)
                            {
                                mM = cM;
                                sM = l*15+3+maxIndex;
                            }

                            l = hmm->transIndY(k,l+1);
                        }
                    }

                    if (seq1->fwdGapEndsNext( s1Beg + j ))   // flagged gap ends in seq1; Y-gap goes down so earlier
                    {

                        int l = hmm->transIndY(k,0);
                        while (l>=0)
                        {

                            cY = max(pXX->g(l,j) + hmm->probXY(l,k),
                                     small,
                                     pXM->g(l,j) + hmm->probMY(l,k)) + msr->indelY(k);
                            if (cY > mY)
                            {
                                mY = cY;
                                sY = l*15+3+maxIndex;
                            }

                            l = hmm->transIndY(k,l+1);
                        }
                    }


                    if (seq1->fwdChildGapStarts( s1Beg + j ))   // flagged gap starts in seq1 child
                    {

                        cwX = cfVX->g(k,j-1);
                        if (cwX > mwX)
                        {
                            mwX = cwX;
                            swX = k*15+0;
                        }

                        cwM = cfVM->g(k,j-1);
                        if (cwM > mwM)
                        {
                            mwM = cwM;
                            swM = k*15+2;
                        }

                        if (seq2->fwdGapEndsNext( s2bBeg + i ))   // ..and another closes is seq2
                        {

                            cwM = cfYM->g(k,j-1);
                            if (cwM > mwM)
                            {
                                mwM = cwM;
                                swM = k*15+8;
                            }
                        }
                        if (seq2->fwdChildGapEndsNext( s2bBeg + i ))   // ..and another closes is seq2 child
                        {

                            cwM = cfZM->g(k,j-1);
                            if (cwM > mwM)
                            {
                                mwM = cwM;
                                swM = k*15+14;
                            }
                        }
                    }

                    if (seq1->fwdChildGapContinues( s1Beg + j ))   // flagged gap continues in seq1
                    {

                        cwX = cfWX->g(k,j-1);
                        if (cwX > mwX)
                        {
                            mwX = cwX;
                            swX = k*15+9;
                        }

                        cwM = cfWM->g(k,j-1);
                        if (cwM > mwM)
                        {
                            mwM = cwM;
                            swM = k*15+11;
                        }
                    }

                    if (seq1->fwdChildGapEnds( s1Beg + j ))   // flagged gap ends in seq1 child
                    {

                        int l = hmm->transIndY(k,0);
                        while (l>=0)
                        {

                            cX = max(cfWX->g(l,j-1) + hmm->probXX(l,k),
                                     small,
                                     cfWM->g(l,j-1) + hmm->probMX(l,k)) + msr->indelX(k);
                            if (cX > mX)
                            {
                                mX = cX;
                                sX = l*15+9+maxIndex;
                            }

                            cM = max(pWX->g(l,j-1) + hmm->probXM(l,k),
                                     small,
                                     pWM->g(l,j-1) + hmm->probMM(l,k)) + msr->fwdM(k);
                            if (cM > mM)
                            {
                                mM = cM;
                                sM = l*15+9+maxIndex;
                            }

                            l = hmm->transIndY(k,l+1);
                        }
                    }

                    if (seq1->fwdChildGapEndsNext( s1Beg + j ))   // flagged gap ends in seq1 child; Y-gap goes down so earlier
                    {

                        int l = hmm->transIndY(k,0);
                        while (l>=0)
                        {

                            cY = max(pWX->g(l,j) + hmm->probXY(l,k),
                                     small,
                                     pWM->g(l,j) + hmm->probMY(l,k)) + msr->indelY(k);
                            if (cY > mY)
                            {
                                mY = cY;
                                sY = l*15+9+maxIndex;
                            }

                            l = hmm->transIndY(k,l+1);
                        }
                    }


                    if (seq2->fwdGapStarts( s2bBeg + i ))   // flagged gap starts in seq2
                    {

                        cyY = pVY->g(k,j);
                        if (cyY > myY)
                        {
                            myY = cyY;
                            syY = k*15+1;
                        }

                        cyM = pVM->g(k,j);
                        if (cyM > myM)
                        {
                            myM = cyM;
                            syM = k*15+2;
                        }

                        if (seq1->fwdGapEndsNext( s1Beg + j ))   // .. and another closes in seq1
                        {

                            cyM = pXM->g(k,j);
                            if (cyM > myM)
                            {
                                myM = cyM;
                                syM = k*15+5;
                            }
                        }
                        if (seq1->fwdChildGapEndsNext( s1Beg + j ))   // .. and another closes in seq1 child
                        {

                            cyM = pWM->g(k,j);
                            if (cyM > myM)
                            {
                                myM = cyM;
                                syM = k*15+11;
                            }
                        }
                    }

                    if (seq2->fwdGapContinues( s2bBeg + i ))   // flagged gap continues in seq2
                    {

                        cyY = pYY->g(k,j);
                        if (cyY > myY)
                        {
                            myY = cyY;
                            syY = k*15+7;
                        }

                        cyM = pYM->g(k,j);
                        if (cyM > myM)
                        {
                            myM = cyM;
                            syM = k*15+8;
                        }
                    }

                    if (seq2->fwdGapEnds( s2bBeg + i ))   // flagged gap ends in seq2
                    {

                        int l = hmm->transIndY(k,0);
                        while (l>=0)
                        {
                            cY = max(small,
                                     pYY->g(l,j) + hmm->probYY(l,k),
                                     pYM->g(l,j) + hmm->probMY(l,k)) + msr->indelY(k);
                            if (cY > mY)
                            {
                                mY = cY;
                                sY = l*15+6+maxIndex;
                            }

                            cM = max(small,
                                     pYY->g(l,j-1) + hmm->probYM(l,k),
                                     pYM->g(l,j-1) + hmm->probMM(l,k)) + msr->fwdM(k);
                            if (cM > mM)
                            {
                                mM = cM;
                                sM = l*15+6+maxIndex;
                            }

                            l = hmm->transIndY(k,l+1);
                        }
                    }

                    if (seq2->fwdGapEndsNext( s2bBeg + i ))   // flagged gap ends in seq2; X-gap goes right so earlier
                    {

                        int l = hmm->transIndY(k,0);
                        while (l>=0)
                        {

                            cX = max(small,
                                     cfYY->g(l,j-1) + hmm->probYX(l,k),
                                     cfYM->g(l,j-1) + hmm->probMX(l,k)) + msr->indelX(k);
                            if (cX > mX)
                            {
                                mX = cX;
                                sX = l*15+6+maxIndex;
                            }

                            l = hmm->transIndY(k,l+1);
                        }
                    }


                    if (seq2->fwdChildGapStarts( s2bBeg + i ))   // flagged gap starts in seq2 child
                    {

                        czY = pVY->g(k,j);
                        if (czY > mzY)
                        {
                            mzY = czY;
                            szY = k*15+1;
                        }

                        czM = pVM->g(k,j);
                        if (czM > mzM)
                        {
                            mzM = czM;
                            szM = k*15+2;
                        }

                        if (seq1->fwdGapEndsNext( s1Beg + j ))   // .. and another closes in seq1
                        {

                            czM = pXM->g(k,j);
                            if (czM > mzM)
                            {
                                mzM = czM;
                                szM = k*15+5;
                            }
                        }
                        if (seq1->fwdChildGapEndsNext( s1Beg + j ))   // .. and another closes in seq1 child
                        {

                            czM = pWM->g(k,j);
                            if (czM > mzM)
                            {
                                mzM = czM;
                                szM = k*15+11;
                            }
                        }
                    }

                    if (seq2->fwdChildGapContinues( s2bBeg + i ))   // flagged gap continues in seq2
                    {

                        czY = pZY->g(k,j);
                        if (czY > mzY)
                        {
                            mzY = czY;
                            szY = k*15+13;
                        }

                        czM = pZM->g(k,j);
                        if (czM > mzM)
                        {
                            mzM = czM;
                            szM = k*15+14;
                        }
                    }

                    if (seq2->fwdChildGapEnds( s2bBeg + i ))   // flagged gap ends in seq2 child
                    {

                        int l = hmm->transIndY(k,0);
                        while (l>=0)
                        {
                            cY = max(small,
                                     pZY->g(l,j) + hmm->probYY(l,k),
                                     pZM->g(l,j) + hmm->probMY(l,k)) + msr->indelY(k);
                            if (cY > mY)
                            {
                                mY = cY;
                                sY = l*15+12+maxIndex;
                            }

                            cM = max(small,
                                     pZY->g(l,j-1) + hmm->probYM(l,k),
                                     pZM->g(l,j-1) + hmm->probMM(l,k)) + msr->fwdM(k);
                            if (cM > mM)
                            {
                                mM = cM;
                                sM = l*15+12+maxIndex;
                            }

                            l = hmm->transIndY(k,l+1);
                        }
                    }

                    if (seq2->fwdChildGapEndsNext( s2bBeg + i ))   // flagged gap ends in seq2 child; X-gap goes right so earlier
                    {

                        int l = hmm->transIndY(k,0);
                        while (l>=0)
                        {

                            cX = max(small,
                                     cfZY->g(l,j-1) + hmm->probYX(l,k),
                                     cfZM->g(l,j-1) + hmm->probMX(l,k)) + msr->indelX(k);
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

                        cX = max(cfVX->g(l,j-1) + hmm->probXX(l,k),
                                 cfVY->g(l,j-1) + hmm->probYX(l,k),
                                 cfVM->g(l,j-1) + hmm->probMX(l,k)) + msr->indelX(k);
                        if (cX > mX)
                        {
                            mX = cX;
                            sX = l*15+maxIndex;
                        }

                        cY = max(pVX->g(l,j) + hmm->probXY(l,k),
                                 pVY->g(l,j) + hmm->probYY(l,k),
                                 pVM->g(l,j) + hmm->probMY(l,k)) + msr->indelY(k);
                        if (cY > mY)
                        {
                            mY = cY;
                            sY = l*15+maxIndex;
                        }

                        cM = max(pVX->g(l,j-1) + hmm->probXM(l,k),
                                 pVY->g(l,j-1) + hmm->probYM(l,k),
                                 pVM->g(l,j-1) + hmm->probMM(l,k)) + msr->fwdM(k);
                        if (cM > mM)
                        {
                            mM = cM;
                            sM = l*15+maxIndex;
                        }

                        l = hmm->transIndY(k,l+1);
                    }

                    cfVX->s(mX,k,j);
                    cfVY->s(mY,k,j);
                    cfVM->s(mM,k,j);

                    ptVX->s(sX,k,j);
                    ptVY->s(sY,k,j);
                    ptVM->s(sM,k,j);

                    cfXX->s(mxX,k,j);
                    cfXM->s(mxM,k,j);

                    ptXX->s(sxX,k,j);
                    ptXM->s(sxM,k,j);

                    cfYY->s(myY,k,j);
                    cfYM->s(myM,k,j);

                    ptYY->s(syY,k,j);
                    ptYM->s(syM,k,j);

                    cfWX->s(mwX,k,j);
                    cfWM->s(mwM,k,j);

                    ptWX->s(swX,k,j);
                    ptWM->s(swM,k,j);

                    cfZY->s(mzY,k,j);
                    cfZM->s(mzM,k,j);

                    ptZY->s(szY,k,j);
                    ptZM->s(szM,k,j);
                }
            }
        } // FOR(j,mLen)


        // change the rows that are pointed
        tmpVX = pVX;
        tmpVY = pVY;
        tmpVM = pVM;
        tmpXX = pXX;
        tmpXM = pXM;
        tmpYY = pYY;
        tmpYM = pYM;
        tmpWX = pWX;
        tmpWM = pWM;
        tmpZY = pZY;
        tmpZM = pZM;

        pVX = cfVX;
        pVY = cfVY;
        pVM = cfVM;
        pXX = cfXX;
        pXM = cfXM;
        pYY = cfYY;
        pYM = cfYM;
        pWX = cfWX;
        pWM = cfWM;
        pZY = cfZY;
        pZM = cfZM;

        cfVX = tmpVX;
        cfVY = tmpVY;
        cfVM = tmpVM;
        cfXX = tmpXX;
        cfXM = tmpXM;
        cfYY = tmpYY;
        cfYM = tmpYM;
        cfWX = tmpWX;
        cfWM = tmpWM;
        cfZY = tmpZY;
        cfZM = tmpZM;

        if (NOISE>2)
        {
            printMatrix("fM",i,pVM);
            printMatrix("fX",i,pVX);
            printMatrix("fY",i,pVY);
        }
        if (NOISE>3)
        {
            printMatrix("fxM",i,pXM);
            printMatrix("fxX",i,pXX);
            printMatrix("fyM",i,pYM);
            printMatrix("fyY",i,pYY);
        }

    }


    // change the pointers back so "previous" can be recycled
    // and the mid-row calculation is correct
    cfVX = pVX;
    cfVY = pVY;
    cfVM = pVM;
    cfXX = pXX;
    cfXM = pXM;
    cfYY = pYY;
    cfYM = pYM;
    cfWX = pWX;
    cfWM = pWM;
    cfZY = pZY;
    cfZM = pZM;

    // Define pointers for current & previous row
    //
    cbVX = bVX1;
    cbVY = bVY1;
    cbVM = bVM1;
    cbXX = bXX1;
    cbXM = bXM1;
    cbYY = bYY1;
    cbYM = bYM1;
    cbWX = bWX1;
    cbWM = bWM1;
    cbZY = bZY1;
    cbZM = bZM1;

    pVX = bVX2;
    pVY = bVY2;
    pVM = bVM2;
    pXX = bXX2;
    pXM = bXM2;
    pYY = bYY2;
    pYM = bYM2;
    pWX = bWX2;
    pWM = bWM2;
    pZY = bZY2;
    pZM = bZM2;


    if (s2<e2)
    {


        // A loop through the second half of seq2
        //
        RFOR(i,s2eLen+1)
        {

            // A loop through seq1
            //
            RFOR(j,s1Len)
            {

                // Starting: set the corner values
                //
                if (i==s2eLen+1 && j==s1Len)
                {

                    FOR(k,nState)
                    {
                        // starting values
                        cbVX->s(bwdvX->g(k),k,j);
                        cbVY->s(bwdvY->g(k),k,j);
                        cbVM->s(bwdvM->g(k),k,j);

                        cbXX->s(bwdxX->g(k),k,j);
                        cbXM->s(bwdxM->g(k),k,j);
                        cbYY->s(bwdyY->g(k),k,j);
                        cbYM->s(bwdyM->g(k),k,j);

                        cbWX->s(bwdwX->g(k),k,j);
                        cbWM->s(bwdwM->g(k),k,j);
                        cbZY->s(bwdzY->g(k),k,j);
                        cbZM->s(bwdzM->g(k),k,j);

                    }
                    continue;
                }

                // Compute the substitution prices
                //
                msr->computeBwd( s1Beg+j, s2eBeg+i );


                if (i<s2eLen+1 && j==s1Len)   // y-gaps are possible
                {
                    FOR(k,nState)
                    {

                        mX=mY=mM=mxX=mxM=myY=myM=mwX=mwM=mzY=mzM=small;
                        cX=cY=cM=cxX=cxM=cyY=cyM=cwX=cwM=czY=czM=small;

                        if (seq2->bwdGapStarts( s2eBeg + i ))   // flagged gap starts in seq2
                        {

                            int l = hmm->transIndX(k,0);
                            while (l>=0)
                            {

                                cyY = hmm->probYY(k,l) + msr->indelY(l) + pVY->g(l,j);
                                if (cyY > myY)
                                {
                                    myY = cyY;
                                }

                                cyM = hmm->probMY(k,l) + msr->indelY(l) + pVY->g(l,j);
                                if (cyM > myM)
                                {
                                    myM = cyM;
                                }

                                l = hmm->transIndX(k,l+1);
                            }
                        }

                        if (seq2->bwdGapContinues( s2eBeg + i ))   // flagged gap continues in seq2
                        {

                            cyY = pYY->g(k,j);
                            if (cyY > myY)
                            {
                                myY = cyY;
                            }

                            cyM = pYM->g(k,j);
                            if (cyM > myM)
                            {
                                myM = cyM;
                            }
                        }

                        if (seq2->bwdGapEnds( s2eBeg + i ))   // flagged gap ends in seq2
                        {

                            cY =  pYY->g(k,j);
                            if (cY > mY)
                            {
                                mY = cY;
                            }

                            cM =  pYM->g(k,j);
                            if (cM > mM)
                            {
                                mM = cM;
                            }
                        }

                        if (seq2->bwdChildGapStarts( s2eBeg + i ))   // flagged gap starts in seq2 child
                        {

                            int l = hmm->transIndX(k,0);
                            while (l>=0)
                            {

                                czY = hmm->probYY(k,l) + msr->indelY(l) + pVY->g(l,j);
                                if (czY > mzY)
                                {
                                    mzY = czY;
                                }

                                czM = hmm->probMY(k,l) + msr->indelY(l) + pVY->g(l,j);
                                if (czM > mzM)
                                {
                                    mzM = czM;
                                }

                                l = hmm->transIndX(k,l+1);
                            }
                        }

                        if (seq2->bwdChildGapContinues( s2eBeg + i ))   // flagged gap continues in seq2 child
                        {

                            czY = pZY->g(k,j);
                            if (czY > mzY)
                            {
                                mzY = czY;
                            }

                            czM = pZM->g(k,j);
                            if (czM > mzM)
                            {
                                mzM = czM;
                            }
                        }


                        if (seq2->bwdChildGapEnds( s2eBeg + i ))   // flagged gap ends in seq2 child
                        {

                            cY =  pZY->g(k,j);
                            if (cY > mY)
                            {
                                mY = cY;
                            }

                            cM =  pZM->g(k,j);
                            if (cM > mM)
                            {
                                mM = cM;
                            }
                        }


                        if (seq1->bwdGapStarts( s1Beg + j ))   // flagged gap starts in seq1
                        {

                            int l = hmm->transIndX(k,0);
                            while (l>=0)
                            {

                                cxX = hmm->probXY(k,l) + msr->indelY(l) + pVY->g(l,j);
                                if (cxX > mxX)
                                {
                                    mxX = cxX;
                                }

                                cxM = hmm->probMY(k,l) + msr->indelY(l) + pVY->g(l,j);
                                if (cxM > mxM)
                                {
                                    mxM = cxM;
                                }

                                l = hmm->transIndX(k,l+1);
                            }
                        }

                        // flagged gap starts in seq1 and another closes in seq2
                        if (seq1->bwdGapStarts( s1Beg + j ) && seq2->bwdGapEnds( s2eBeg + i))
                        {

                            cxM = pYM->g(k,j);
                            if (cxM > mxM)
                            {
                                mxM = cxM;
                            }
                        }

                        if (seq1->bwdGapStarts( s1Beg + j ) && seq2->bwdChildGapEnds( s2eBeg + i))
                        {

                            cxM = pZM->g(k,j);
                            if (cxM > mxM)
                            {
                                mxM = cxM;
                            }

                        }


                        if (seq1->bwdChildGapStarts( s1Beg + j ))   // flagged gap starts in seq1 child
                        {

                            int l = hmm->transIndX(k,0);
                            while (l>=0)
                            {

                                cwX = hmm->probXY(k,l) + msr->indelY(l) + pVY->g(l,j);
                                if (cwX > mwX)
                                {
                                    mwX = cwX;
                                }

                                cwM = hmm->probMY(k,l) + msr->indelY(l) + pVY->g(l,j);
                                if (cwM > mwM)
                                {
                                    mwM = cwM;
                                }

                                l = hmm->transIndX(k,l+1);
                            }
                        }

                        // flagged gap starts in seq1 child and another closes in seq2
                        if (seq1->bwdChildGapStarts( s1Beg + j ) && seq2->bwdGapEnds( s2eBeg + i))
                        {

                            cwM = pYM->g(k,j);
                            if (cwM > mwM)
                            {
                                mwM = cwM;
                            }
                        }

                        if (seq1->bwdChildGapStarts( s1Beg + j ) && seq2->bwdChildGapEnds( s2eBeg + i))
                        {

                            cwM = pZM->g(k,j);
                            if (cwM > mwM)
                            {
                                mwM = cwM;
                            }
                        }

                        int l = hmm->transIndX(k,0);
                        while (l>=0)
                        {

                            cX = hmm->probXY(k,l) + msr->indelY(l) + pVY->g(l,j);
                            if (cX > mX)
                            {
                                mX = cX;
                            }
                            cY = hmm->probYY(k,l) + msr->indelY(l) + pVY->g(l,j);
                            if (cY > mY)
                            {
                                mY = cY;
                            }
                            cM = hmm->probMY(k,l) + msr->indelY(l) + pVY->g(l,j);
                            if (cM > mM)
                            {
                                mM = cM;
                            }

                            l = hmm->transIndX(k,l+1);
                        }

                        cbVX->s(mX,k,j);
                        cbVY->s(mY,k,j);
                        cbVM->s(mM,k,j);

                        cbXX->s(mxX,k,j);
                        cbXM->s(mxM,k,j);
                        cbYY->s(myY,k,j);
                        cbYM->s(myM,k,j);

                        cbWX->s(mwX,k,j);
                        cbWM->s(mwM,k,j);
                        cbZY->s(mzY,k,j);
                        cbZM->s(mzM,k,j);

                    }

                }
                else if (i==s2eLen+1 && j<s1Len)   // x-gaps are possible
                {

                    FOR(k,nState)
                    {

                        mX=mY=mM=mxX=mxM=myY=myM=mwX=mwM=mzY=mzM=small;
                        cX=cY=cM=cxX=cxM=cyY=cyM=cwX=cwM=czY=czM=small;


                        if (seq1->bwdGapStarts( s1Beg + j ))   // flagged gap starts in seq1
                        {

                            int l = hmm->transIndX(k,0);
                            while (l>=0)
                            {

                                cxX = hmm->probXX(k,l) + msr->indelX(l) + cbVX->g(l,j+1);
                                if (cxX > mxX)
                                {
                                    mxX = cxX;
                                }

                                cxM = hmm->probMX(k,l) + msr->indelX(l) + cbVX->g(l,j+1);
                                if (cxM > mxM)
                                {
                                    mxM = cxM;
                                }

                                l = hmm->transIndX(k,l+1);
                            }
                        }

                        if (seq1->bwdGapContinues( s1Beg + j ))   // flagged gap continues in seq1
                        {

                            cxX = cbXX->g(k,j+1);
                            if (cxX > mxX)
                            {
                                mxX = cxX;
                            }

                            cxM = cbXM->g(k,j+1);
                            if (cxM > mxM)
                            {
                                mxM = cxM;
                            }
                        }

                        if (seq1->bwdGapEnds( s1Beg + j ))     // flagged gap ends in seq1
                        {

                            cX = cbXX->g(k,j+1);
                            if (cX > mX)
                            {
                                mX = cX;
                            }

                            cM = cbXM->g(k,j+1);
                            if (cM > mM)
                            {
                                mM = cM;
                            }
                        }


                        if (seq1->bwdChildGapStarts( s1Beg + j ))   // flagged gap starts in seq1 child
                        {

                            int l = hmm->transIndX(k,0);
                            while (l>=0)
                            {

                                cwX = hmm->probXX(k,l) + msr->indelX(l) + cbVX->g(l,j+1);
                                if (cwX > mwX)
                                {
                                    mwX = cwX;
                                }

                                cwM = hmm->probMX(k,l) + msr->indelX(l) + cbVX->g(l,j+1);
                                if (cwM > mwM)
                                {
                                    mwM = cwM;
                                }

                                l = hmm->transIndX(k,l+1);
                            }
                        }

                        if (seq1->bwdChildGapContinues( s1Beg + j ))   // flagged gap continues in seq1 child
                        {

                            cwX = cbWX->g(k,j+1);
                            if (cwX > mwX)
                            {
                                mwX = cwX;
                            }

                            cwM = cbWM->g(k,j+1);
                            if (cwM > mwM)
                            {
                                mwM = cwM;
                            }
                        }

                        if (seq1->bwdChildGapEnds( s1Beg + j ))     // flagged gap ends in seq1 child
                        {

                            cX = cbWX->g(k,j+1);
                            if (cX > mX)
                            {
                                mX = cX;
                            }

                            cM = cbWM->g(k,j+1);
                            if (cM > mM)
                            {
                                mM = cM;
                            }
                        }


                        if (seq2->bwdGapStarts( s2eBeg + i ))   // flagged gap starts in seq2
                        {

                            int l = hmm->transIndX(k,0);
                            while (l>=0)
                            {
                                cyY = hmm->probYX(k,l) + msr->indelX(l) + cbVX->g(l,j+1);
                                if (cyY > myY)
                                {
                                    myY = cyY;
                                }

                                cyM = hmm->probMX(k,l) + msr->indelX(l) + cbVX->g(l,j+1);
                                if (cyM > myM)
                                {
                                    myM = cyM;
                                }

                                l = hmm->transIndX(k,l+1);
                            }
                        }

                        // flagged gap starts in seq2 and another closes in seq1
                        if (seq2->bwdGapStarts( s2eBeg + i ) && seq1->bwdGapEnds( s1Beg + j ))
                        {

                            cyM = cbXM->g(k,j+1);
                            if (cyM > myM)
                            {
                                myM = cyM;
                            }
                        }

                        if (seq2->bwdGapStarts( s2eBeg + i ) && seq1->bwdChildGapEnds( s1Beg + j ))
                        {

                            cyM = cbWM->g(k,j+1);
                            if (cyM > myM)
                            {
                                myM = cyM;
                            }
                        }


                        if (seq2->bwdChildGapStarts( s2eBeg + i ))   // flagged gap starts in seq2 child
                        {

                            int l = hmm->transIndX(k,0);
                            while (l>=0)
                            {
                                czY = hmm->probYX(k,l) + msr->indelX(l) + cbVX->g(l,j+1);
                                if (czY > mzY)
                                {
                                    mzY = czY;
                                }

                                czM = hmm->probMX(k,l) + msr->indelX(l) + cbVX->g(l,j+1);
                                if (czM > mzM)
                                {
                                    mzM = czM;
                                }

                                l = hmm->transIndX(k,l+1);
                            }
                        }

                        // flagged gap starts in seq2 child and another closes in seq1
                        if (seq2->bwdChildGapStarts( s2eBeg + i ) && seq1->bwdGapEnds( s1Beg + j ))
                        {

                            czM = cbXM->g(k,j+1);
                            if (czM > mzM)
                            {
                                mzM = czM;
                            }
                        }

                        if (seq2->bwdChildGapStarts( s2eBeg + i ) && seq1->bwdChildGapEnds( s1Beg + j ))
                        {

                            czM = cbWM->g(k,j+1);
                            if (czM > mzM)
                            {
                                mzM = czM;
                            }
                        }

                        int l = hmm->transIndX(k,0);
                        while (l>=0)
                        {

                            cX = hmm->probXX(k,l) + msr->indelX(l) + cbVX->g(l,j+1);
                            if (cX > mX)
                            {
                                mX = cX;
                            }
                            cY = hmm->probYX(k,l) + msr->indelX(l) + cbVX->g(l,j+1);
                            if (cY > mY)
                            {
                                mY = cY;
                            }
                            cM = hmm->probMX(k,l) + msr->indelX(l) + cbVX->g(l,j+1);
                            if (cM > mM)
                            {
                                mM = cM;
                            }

                            l = hmm->transIndX(k,l+1);
                        }

                        cbVX->s(mX,k,j);
                        cbVY->s(mY,k,j);
                        cbVM->s(mM,k,j);

                        cbXX->s(mxX,k,j);
                        cbXM->s(mxM,k,j);
                        cbYY->s(myY,k,j);
                        cbYM->s(myM,k,j);

                        cbWX->s(mwX,k,j);
                        cbWM->s(mwM,k,j);
                        cbZY->s(mzY,k,j);
                        cbZM->s(mzM,k,j);
                    }

                }
                else if (i<s2eLen+1 && j<s1Len)   // also matches are possible
                {

                    FOR(k,nState)
                    {

                        mX=mY=mM=mxX=mxM=myY=myM=mwX=mwM=mzY=mzM=small;
                        cX=cY=cM=cxX=cxM=cyY=cyM=cwX=cwM=czY=czM=small;

                        if (seq1->bwdGapStarts( s1Beg + j ))   // flagged gap starts in seq1
                        {

                            int l = hmm->transIndX(k,0);
                            while (l>=0)
                            {

                                cxX = max(hmm->probXX(k,l) + msr->indelX(l) + cbVX->g(l,j+1),
                                          hmm->probXY(k,l) + msr->indelY(l) + pVY->g(l,j),
                                          hmm->probXM(k,l) + msr->bwdM(l) + pVM->g(l,j+1));
                                if (cxX > mxX)
                                {
                                    mxX = cxX;
                                }

                                cxM = max(hmm->probMX(k,l) + msr->indelX(l) + cbVX->g(l,j+1),
                                          hmm->probMY(k,l) + msr->indelY(l) + pVY->g(l,j),
                                          hmm->probMM(k,l) + msr->bwdM(l) + pVM->g(l,j+1));
                                if (cxM > mxM)
                                {
                                    mxM = cxM;
                                }

                                l = hmm->transIndX(k,l+1);
                            }
                        }

                        if (seq1->bwdGapContinues( s1Beg + j ))   // flagged gap continues in seq1
                        {

                            cxX = cbXX->g(k,j+1);
                            if (cxX > mxX)
                            {
                                mxX = cxX;
                            }

                            cxM = cbXM->g(k,j+1);
                            if (cxM > mxM)
                            {
                                mxM = cxM;
                            }
                        }

                        if (seq1->bwdGapEnds( s1Beg + j ))     // flagged gap ends in seq1 or its child
                        {

                            cX = cbXX->g(k,j+1);
                            if (cX > mX)
                            {
                                mX = cX;
                            }

                            cM = cbXM->g(k,j+1);
                            if (cM > mM)
                            {
                                mM = cM;
                            }
                        }

                        // flagged gap starts in seq1 and another closes in seq2
                        if (seq1->bwdGapStarts( s1Beg + j ) && seq2->bwdGapEnds( s2eBeg + i ))
                        {

                            cxM = pYM->g(k,j);
                            if (cxM > mxM)
                            {
                                mxM = cxM;
                            }
                        }

                        if (seq1->bwdGapStarts( s1Beg + j ) && seq2->bwdChildGapEnds( s2eBeg + i ))
                        {

                            cxM = pZM->g(k,j);
                            if (cxM > mxM)
                            {
                                mxM = cxM;
                            }
                        }


                        if (seq1->bwdChildGapStarts( s1Beg + j ))   // flagged gap starts in seq1 child
                        {

                            int l = hmm->transIndX(k,0);
                            while (l>=0)
                            {

                                cwX = max(hmm->probXX(k,l) + msr->indelX(l) + cbVX->g(l,j+1),
                                          hmm->probXY(k,l) + msr->indelY(l) + pVY->g(l,j),
                                          hmm->probXM(k,l) + msr->bwdM(l) + pVM->g(l,j+1));
                                if (cwX > mwX)
                                {
                                    mwX = cwX;
                                }

                                cwM = max(hmm->probMX(k,l) + msr->indelX(l) + cbVX->g(l,j+1),
                                          hmm->probMY(k,l) + msr->indelY(l) + pVY->g(l,j),
                                          hmm->probMM(k,l) + msr->bwdM(l) + pVM->g(l,j+1));
                                if (cwM > mwM)
                                {
                                    mwM = cwM;
                                }

                                l = hmm->transIndX(k,l+1);
                            }
                        }

                        if (seq1->bwdChildGapContinues( s1Beg + j ))   // flagged gap continues in seq1 child
                        {

                            cwX = cbWX->g(k,j+1);
                            if (cwX > mwX)
                            {
                                mwX = cwX;
                            }

                            cwM = cbWM->g(k,j+1);
                            if (cwM > mwM)
                            {
                                mwM = cwM;
                            }
                        }

                        if (seq1->bwdChildGapEnds( s1Beg + j ))     // flagged gap ends in seq1 child
                        {

                            cX = cbWX->g(k,j+1);
                            if (cX > mX)
                            {
                                mX = cX;
                            }

                            cM = cbWM->g(k,j+1);
                            if (cM > mM)
                            {
                                mM = cM;
                            }
                        }

                        // flagged gap starts in seq1 child and another closes in seq2
                        if (seq1->bwdChildGapStarts( s1Beg + j ) && seq2->bwdGapEnds( s2eBeg + i ))
                        {

                            cwM = pYM->g(k,j);
                            if (cwM > mwM)
                            {
                                mwM = cwM;
                            }
                        }

                        if (seq1->bwdChildGapStarts( s1Beg + j ) && seq2->bwdChildGapEnds( s2eBeg + i ))
                        {

                            cwM = pZM->g(k,j);
                            if (cwM > mwM)
                            {
                                mwM = cwM;
                            }
                        }


                        if (seq2->bwdGapStarts( s2eBeg + i ))   // flagged gap starts in seq2
                        {

                            int l = hmm->transIndX(k,0);
                            while (l>=0)
                            {
                                cyY = max(hmm->probYX(k,l) + msr->indelX(l) + cbVX->g(l,j+1),
                                          hmm->probYY(k,l) + msr->indelY(l) + pVY->g(l,j),
                                          hmm->probYM(k,l) + msr->bwdM(l) + pVM->g(l,j+1));
                                if (cyY > myY)
                                {
                                    myY = cyY;
                                }

                                cyM = max(hmm->probMX(k,l) + msr->indelX(l) + cbVX->g(l,j+1),
                                          hmm->probMY(k,l) + msr->indelY(l) + pVY->g(l,j),
                                          hmm->probMM(k,l) + msr->bwdM(l) + pVM->g(l,j+1));
                                if (cyM > myM)
                                {
                                    myM = cyM;
                                }

                                l = hmm->transIndX(k,l+1);
                            }
                        }

                        if (seq2->bwdGapContinues( s2eBeg + i ))   // flagged gap continues in seq2
                        {

                            cyY = pYY->g(k,j);
                            if (cyY > myY)
                            {
                                myY = cyY;
                            }

                            cyM = pYM->g(k,j);
                            if (cyM > myM)
                            {
                                myM = cyM;
                            }
                        }

                        if (seq2->bwdGapEnds( s2eBeg + i ))   // flagged gap ends in seq2
                        {

                            cY = pYY->g(k,j);
                            if (cY > mY)
                            {
                                mY = cY;
                            }

                            cM = pYM->g(k,j);
                            if (cM > mM)
                            {
                                mM = cM;
                            }
                        }

                        // flagged gap starts in seq2 and another closes in seq1
                        if (seq2->bwdGapStarts( s2eBeg + i ) && seq1->bwdGapEnds( s1Beg + j ))
                        {

                            cyM = cbXM->g(k,j+1);
                            if (cyM > myM)
                            {
                                myM = cyM;
                            }
                        }

                        if (seq2->bwdGapStarts( s2eBeg + i ) && seq1->bwdChildGapEnds( s1Beg + j ))
                        {

                            cyM = cbWM->g(k,j+1);
                            if (cyM > myM)
                            {
                                myM = cyM;
                            }
                        }


                        if (seq2->bwdChildGapStarts( s2eBeg + i ))   // flagged gap starts in seq2 child
                        {

                            int l = hmm->transIndX(k,0);
                            while (l>=0)
                            {
                                czY = max(hmm->probYX(k,l) + msr->indelX(l) + cbVX->g(l,j+1),
                                          hmm->probYY(k,l) + msr->indelY(l) + pVY->g(l,j),
                                          hmm->probYM(k,l) + msr->bwdM(l) + pVM->g(l,j+1));
                                if (czY > mzY)
                                {
                                    mzY = czY;
                                }

                                czM = max(hmm->probMX(k,l) + msr->indelX(l) + cbVX->g(l,j+1),
                                          hmm->probMY(k,l) + msr->indelY(l) + pVY->g(l,j),
                                          hmm->probMM(k,l) + msr->bwdM(l) + pVM->g(l,j+1));
                                if (czM > mzM)
                                {
                                    mzM = czM;
                                }

                                l = hmm->transIndX(k,l+1);
                            }
                        }

                        if (seq2->bwdChildGapContinues( s2eBeg + i ))   // flagged gap continues in seq2 child
                        {

                            czY = pZY->g(k,j);
                            if (czY > mzY)
                            {
                                mzY = czY;
                            }

                            czM = pZM->g(k,j);
                            if (czM > mzM)
                            {
                                mzM = czM;
                            }
                        }

                        if (seq2->bwdChildGapEnds( s2eBeg + i ))   // flagged gap ends in seq2
                        {

                            cY = pZY->g(k,j);
                            if (cY > mY)
                            {
                                mY = cY;
                            }

                            cM = pZM->g(k,j);
                            if (cM > mM)
                            {
                                mM = cM;
                            }
                        }

                        // flagged gap starts in seq2 child and another closes in seq1
                        if (seq2->bwdChildGapStarts( s2eBeg + i ) && seq1->bwdGapEnds( s1Beg + j ))
                        {

                            czM = cbXM->g(k,j+1);
                            if (czM > mzM)
                            {
                                mzM = czM;
                            }
                        }

                        if (seq2->bwdChildGapStarts( s2eBeg + i ) && seq1->bwdChildGapEnds( s1Beg + j ))
                        {

                            czM = cbWM->g(k,j+1);
                            if (czM > mzM)
                            {
                                mzM = czM;
                            }
                        }

                        int l = hmm->transIndX(k,0);
                        while (l>=0)
                        {

                            cX = max(hmm->probXX(k,l) + msr->indelX(l) + cbVX->g(l,j+1),
                                     hmm->probXY(k,l) + msr->indelY(l) + pVY->g(l,j),
                                     hmm->probXM(k,l) + msr->bwdM(l) + pVM->g(l,j+1));
                            if (cX > mX)
                            {
                                mX = cX;
                            }

                            cY = max(hmm->probYX(k,l) + msr->indelX(l) + cbVX->g(l,j+1),
                                     hmm->probYY(k,l) + msr->indelY(l) + pVY->g(l,j),
                                     hmm->probYM(k,l) + msr->bwdM(l) + pVM->g(l,j+1));
                            if (cY > mY)
                            {
                                mY = cY;
                            }

                            cM = max(hmm->probMX(k,l) + msr->indelX(l) + cbVX->g(l,j+1),
                                     hmm->probMY(k,l) + msr->indelY(l) + pVY->g(l,j),
                                     hmm->probMM(k,l) + msr->bwdM(l) + pVM->g(l,j+1));
                            if (cM > mM)
                            {
                                mM = cM;
                            }

                            l = hmm->transIndX(k,l+1);
                        }

                        cbVX->s(mX,k,j);
                        cbVY->s(mY,k,j);
                        cbVM->s(mM,k,j);

                        cbXX->s(mxX,k,j);
                        cbXM->s(mxM,k,j);
                        cbYY->s(myY,k,j);
                        cbYM->s(myM,k,j);

                        cbWX->s(mwX,k,j);
                        cbWM->s(mwM,k,j);
                        cbZY->s(mzY,k,j);
                        cbZM->s(mzM,k,j);
                    }
                }
            }	/// RFOR(j,s1Len)

            // change the rows that are pointed
            tmpVX = pVX;
            tmpVY = pVY;
            tmpVM = pVM;
            tmpXX = pXX;
            tmpXM = pXM;
            tmpYY = pYY;
            tmpYM = pYM;
            tmpWX = pWX;
            tmpWM = pWM;
            tmpZY = pZY;
            tmpZM = pZM;

            pVX = cbVX;
            pVY = cbVY;
            pVM = cbVM;
            pXX = cbXX;
            pXM = cbXM;
            pYY = cbYY;
            pYM = cbYM;
            pWX = cbWX;
            pWM = cbWM;
            pZY = cbZY;
            pZM = cbZM;

            cbVX = tmpVX;
            cbVY = tmpVY;
            cbVM = tmpVM;
            cbXX = tmpXX;
            cbXM = tmpXM;
            cbYY = tmpYY;
            cbYM = tmpYM;
            cbWX = tmpWX;
            cbWM = tmpWM;
            cbZY = tmpZY;
            cbZM = tmpZM;

            if (NOISE>2)
            {
                printMatrix("bM",i,pVM);
                printMatrix("bX",i,pVX);
                printMatrix("bY",i,pVY);
            }
            if (NOISE>3)
            {
                printMatrix("bxM",i,pXM);
                printMatrix("bxX",i,pXX);
                printMatrix("byM",i,pYM);
                printMatrix("byY",i,pYY);
            }
        }

        // change the pointers back so the mid-row calculation is correct
        cbVX = pVX;
        cbVY = pVY;
        cbVM = pVM;
        cbXX = pXX;
        cbXM = pXM;
        cbYY = pYY;
        cbYM = pYM;
        cbWX = pWX;
        cbWM = pWM;
        cbZY = pZY;
        cbZM = pZM;

    }

    // Cases where only x-gaps possible
    //
    if (s2==e2)
    {

        // Starting: set the corner values
        //
        FOR(k,nState)
        {
            // starting values
            cbVM->s(bwdvM->g(k),k,s1Len);
            cbVX->s(bwdvX->g(k),k,s1Len);
            cbVY->s(bwdvY->g(k),k,s1Len);

            cbXM->s(bwdxM->g(k),k,s1Len);
            cbXX->s(bwdxX->g(k),k,s1Len);
            cbYM->s(bwdyM->g(k),k,s1Len);
            cbYY->s(bwdyY->g(k),k,s1Len);

            cbWM->s(bwdwM->g(k),k,s1Len);
            cbWX->s(bwdwX->g(k),k,s1Len);
            cbZM->s(bwdzM->g(k),k,s1Len);
            cbZY->s(bwdzY->g(k),k,s1Len);
        }

        RFOR(j,s1Len-1)
        {

            // Compute the substitution prices
            //
            msr->computeBwd( s1Beg+j, s2eBeg );

            FOR(k,nState)
            {

                // move into X-matrix
                //
                mX=mY=mM=mxX=mxM=myY=myM=mwX=mwM=mzY=mzM=small;
                cX=cY=cM=cxX=cxM=cyY=cyM=cwX=cwM=czY=czM=small;

                if (seq1->bwdGapStarts( s1Beg + j ))   // flagged gap starts in seq1
                {

                    int l = hmm->transIndX(k,0);
                    while (l>=0)
                    {

                        cxX = hmm->probXX(k,l) + msr->indelX(l) + cbVX->g(l,j+1);
                        if (cxX > mxX)
                        {
                            mxX = cxX;
                        }

                        cxM = hmm->probMX(k,l) + msr->indelX(l) + cbVX->g(l,j+1);
                        if (cxM > mxM)
                        {
                            mxM = cxM;
                        }

                        l = hmm->transIndX(k,l+1);
                    }
                }

                if (seq1->bwdGapContinues( s1Beg + j ))   // flagged gap continues in seq1
                {

                    cxX = cbXX->g(k,j+1);
                    if (cxX > mxX)
                    {
                        mxX = cxX;
                    }

                    cxM = cbXM->g(k,j+1);
                    if (cxM > mxM)
                    {
                        mxM = cxM;
                    }
                }

                if (seq1->bwdGapEnds( s1Beg + j ))     // flagged gap ends in seq1 or its child
                {

                    cX = cbXX->g(k,j+1);
                    if (cX > mX)
                    {
                        mX = cX;
                    }

                    cM = cbXM->g(k,j+1);
                    if (cM > mM)
                    {
                        mM = cM;
                    }
                }


                if (seq1->bwdChildGapStarts( s1Beg + j ))   // flagged gap starts in seq1 child
                {

                    int l = hmm->transIndX(k,0);
                    while (l>=0)
                    {

                        cwX = hmm->probXX(k,l) + msr->indelX(l) + cbVX->g(l,j+1);
                        if (cwX > mwX)
                        {
                            mwX = cwX;
                        }

                        cwM = hmm->probMX(k,l) + msr->indelX(l) + cbVX->g(l,j+1);
                        if (cwM > mwM)
                        {
                            mwM = cwM;
                        }

                        l = hmm->transIndX(k,l+1);
                    }
                }

                if (seq1->bwdChildGapContinues( s1Beg + j ))   // flagged gap continues in seq1
                {

                    cwX = cbWX->g(k,j+1);
                    if (cwX > mwX)
                    {
                        mwX = cwX;
                    }

                    cwM = cbWM->g(k,j+1);
                    if (cwM > mwM)
                    {
                        mwM = cwM;
                    }
                }

                if (seq1->bwdChildGapEnds( s1Beg + j ))     // flagged gap ends in seq1 or its child
                {

                    cX = cbWX->g(k,j+1);
                    if (cX > mX)
                    {
                        mX = cX;
                    }

                    cM = cbWM->g(k,j+1);
                    if (cM > mM)
                    {
                        mM = cM;
                    }
                }

                int l = hmm->transIndX(k,0);
                while (l>=0)
                {

                    cX = hmm->probXX(k,l) + msr->indelX(l) + cbVX->g(l,j+1);
                    if (cX > mX)
                    {
                        mX = cX;
                    }
                    cY = hmm->probYX(k,l) + msr->indelX(l) + cbVX->g(l,j+1);
                    if (cY > mY)
                    {
                        mY = cY;
                    }
                    cM = hmm->probMX(k,l) + msr->indelX(l) + cbVX->g(l,j+1);
                    if (cM > mM)
                    {
                        mM = cM;
                    }

                    l = hmm->transIndX(k,l+1);
                }

                cbVX->s(mX,k,j);
                cbVY->s(mY,k,j);
                cbVM->s(mM,k,j);

                cbXX->s(mxX,k,j);
                cbXM->s(mxM,k,j);
                cbYY->s(myY,k,j);
                cbYM->s(myM,k,j);

                cbWX->s(mwX,k,j);
                cbWM->s(mwM,k,j);
                cbZY->s(mzY,k,j);
                cbZM->s(mzM,k,j);
            }
        }

        if (NOISE>2)
        {
            printMatrix("BM",0,cbVM);
            printMatrix("BX",0,cbVX);
            printMatrix("BY",0,cbVY);
        }
        if (NOISE>3)
        {
            printMatrix("BxM",0,cbXM);
            printMatrix("BxX",0,cbXX);
            printMatrix("ByM",0,cbYM);
            printMatrix("ByY",0,cbYY);
        }
    }

    // Find k (i.e. the column through which the alignment path goes)
    //
    vector<Cell> maxCell;
    double maxScore = small;

    j=0;
//    if(s2==e2 && s1Len>1)

    if (s2==e2)
        j++;
    for (; j<mLen; j++)
    {

        double tmp;
        for (int k=0; k<nState; k++)
        {

            tmp = cfVY->g(k,j)+cbVY->g(k,j);
            if (tmp>maxScore)
            {

                maxScore = tmp;
                maxCell.clear();
                Cell c = {ptVY->g(k,j),k*15+1,j};
                maxCell.push_back(c);

            }
            else if (tmp==maxScore)
            {

                Cell c = {ptVY->g(k,j),k*15+1,j};
                maxCell.push_back(c);

            }

            tmp = cfYY->g(k,j)+cbYY->g(k,j);
            if (tmp>maxScore)
            {

                maxScore = tmp;
                maxCell.clear();
                Cell c = {ptYY->g(k,j),k*15+7,j};
                maxCell.push_back(c);

            }
            else if (tmp==maxScore)
            {

                Cell c = {ptYY->g(k,j),k*15+7,j};
                maxCell.push_back(c);

            }

            tmp = cfYM->g(k,j)+cbYM->g(k,j);
            if (tmp>maxScore)
            {

                maxScore = tmp;
                maxCell.clear();
                Cell c = {ptYM->g(k,j),k*15+8,j};
                maxCell.push_back(c);

            }
            else if (tmp==maxScore)
            {

                Cell c = {ptYM->g(k,j),k*15+8,j};
                maxCell.push_back(c);

            }

            tmp = cfZY->g(k,j)+cbZY->g(k,j);
            if (tmp>maxScore)
            {

                maxScore = tmp;
                maxCell.clear();
                Cell c = {ptZY->g(k,j),k*15+13,j};
                maxCell.push_back(c);

            }
            else if (tmp==maxScore)
            {

                Cell c = {ptZY->g(k,j),k*15+13,j};
                maxCell.push_back(c);

            }

            tmp = cfZM->g(k,j)+cbZM->g(k,j);
            if (tmp>maxScore)
            {

                maxScore = tmp;
                maxCell.clear();
                Cell c = {ptZM->g(k,j),k*15+14,j};
                maxCell.push_back(c);

            }
            else if (tmp==maxScore)
            {

                Cell c = {ptZM->g(k,j),k*15+14,j};
                maxCell.push_back(c);

            }

            if (s2<e2 || j>0)
            {

                tmp = cfVX->g(k,j)+cbVX->g(k,j);
                if (tmp>maxScore)
                {

                    maxScore = tmp;
                    maxCell.clear();
                    Cell c = {ptVX->g(k,j),k*15+0,j};
                    maxCell.push_back(c);

                }
                else if (tmp==maxScore)
                {

                    Cell c = {ptVX->g(k,j),k*15+0,j};
                    maxCell.push_back(c);

                }


                tmp = cfVM->g(k,j)+cbVM->g(k,j);
                if (tmp>maxScore)
                {

                    maxScore = tmp;
                    maxCell.clear();
                    Cell c = {ptVM->g(k,j),k*15+2,j};
                    maxCell.push_back(c);

                }
                else if (tmp==maxScore)
                {

                    Cell c = {ptVM->g(k,j),k*15+2,j};
                    maxCell.push_back(c);

                }

                tmp = cfXX->g(k,j)+cbXX->g(k,j);
                if (tmp>maxScore)
                {

                    maxScore = tmp;
                    maxCell.clear();
                    Cell c = {ptXX->g(k,j),k*15+3,j};
                    maxCell.push_back(c);

                }
                else if (tmp==maxScore)
                {

                    Cell c = {ptXX->g(k,j),k*15+3,j};
                    maxCell.push_back(c);

                }

                tmp = cfXM->g(k,j)+cbXM->g(k,j);
                if (tmp>maxScore)
                {

                    maxScore = tmp;
                    maxCell.clear();
                    Cell c = {ptXM->g(k,j),k*15+5,j};
                    maxCell.push_back(c);

                }
                else if (tmp==maxScore)
                {

                    Cell c = {ptXM->g(k,j),k*15+5,j};
                    maxCell.push_back(c);

                }

                tmp = cfWX->g(k,j)+cbWX->g(k,j);
                if (tmp>maxScore)
                {

                    maxScore = tmp;
                    maxCell.clear();
                    Cell c = {ptWX->g(k,j),k*15+9,j};
                    maxCell.push_back(c);

                }
                else if (tmp==maxScore)
                {

                    Cell c = {ptWX->g(k,j),k*15+9,j};
                    maxCell.push_back(c);

                }

                tmp = cfWM->g(k,j)+cbWM->g(k,j);
                if (tmp>maxScore)
                {

                    maxScore = tmp;
                    maxCell.clear();
                    Cell c = {ptWM->g(k,j),k*15+11,j};
                    maxCell.push_back(c);

                }
                else if (tmp==maxScore)
                {

                    Cell c = {ptWM->g(k,j),k*15+11,j};
                    maxCell.push_back(c);

                }

            }

        }
    }

//    if(count==2)
    maxFullScore = maxScore;

    int ms = maxCell.size();
    int rc = 0;
    if (ms>1)
    {
        rc = rndInt(ms);
        if (rc==ms)
        {
            cout<<"Random number error. Tell Tim (timm@ebi.ac.uk) that he was wrong."<<endl;
            exit(-1);
        }
    }
    Cell c = maxCell.at(rc);
    maxCell.clear();

    // Define the site
    //
    newsite->addNewSite();
    newsite->nullSite(false);
    newsite->isAnchor(false);

    msr->computeFwd( s1Beg+c.k , s2eBeg );

    int fState = c.prev/15;
    int fMatch = c.prev%15;

    int bState = c.curr/15;
    int bMatch = c.curr%15;

    newsite->currModelState(bState);
    newsite->currMatchState(bMatch);

    double forwardEnd = small;
    double backwardEnd = small;

    if (bMatch==0)
    {
        forwardEnd = cfVX->g(bState,c.k);
        backwardEnd = cbVX->g(bState,c.k);
    }
    else if (bMatch==1)
    {
        forwardEnd = cfVY->g(bState,c.k);
        backwardEnd = cbVY->g(bState,c.k);
    }
    else if (bMatch==2)
    {
        forwardEnd = cfVM->g(bState,c.k);
        backwardEnd = cbVM->g(bState,c.k);
    }
    else if (bMatch==3)
    {
        forwardEnd = cfXX->g(bState,c.k);
        backwardEnd = cbXX->g(bState,c.k);
    }
    else if (bMatch==5)
    {
        forwardEnd = cfXM->g(bState,c.k);
        backwardEnd = cbXM->g(bState,c.k);
    }
    else if (bMatch==7)
    {
        forwardEnd = cfYY->g(bState,c.k);
        backwardEnd = cbYY->g(bState,c.k);
    }
    else if (bMatch==8)
    {
        forwardEnd = cfYM->g(bState,c.k);
        backwardEnd = cbYM->g(bState,c.k);
    }
    else if (bMatch==9)
    {
        forwardEnd = cfWX->g(bState,c.k);
        backwardEnd = cbWX->g(bState,c.k);
    }
    else if (bMatch==11)
    {
        forwardEnd = cfWM->g(bState,c.k);
        backwardEnd = cbWM->g(bState,c.k);
    }
    else if (bMatch==13)
    {
        forwardEnd = cfZY->g(bState,c.k);
        backwardEnd = cbZY->g(bState,c.k);
    }
    else if (bMatch==14)
    {
        forwardEnd = cfZM->g(bState,c.k);
        backwardEnd = cbZM->g(bState,c.k);
    }
    else
    {
        cout<<"Hirschberg::error1 ("<<fMatch<<","<<bMatch<<")"<<endl;
    }

    // here, changes needed !!
    //
    if (bMatch==0 || bMatch==1 || bMatch==2)   // move penalty for unflagged sites
    {
        if (bMatch==0)
        {
            if (fMatch==0 || fMatch==3 || fMatch==9)
            {
                backwardEnd += hmm->probXX(fState,bState) + msr->indelX(bState);
            }
            else if (fMatch==1 || fMatch==7 || fMatch==13)
            {
                backwardEnd += hmm->probYX(fState,bState) + msr->indelX(bState);
            }
            else if (fMatch==2 || fMatch==5 || fMatch==8 || fMatch==11 || fMatch==14)
            {
                backwardEnd += hmm->probMX(fState,bState) + msr->indelX(bState);
            }
        }
        else if (bMatch==1)
        {
            if (fMatch==0 || fMatch==3 || fMatch==9)
            {
                backwardEnd += hmm->probXY(fState,bState) +msr->indelY(bState);
            }
            else if (fMatch==1 || fMatch==7 || fMatch==13)
            {
                backwardEnd += hmm->probYY(fState,bState) + msr->indelY(bState);
            }
            else if (fMatch==2 || fMatch==5 || fMatch==8 || fMatch==11 || fMatch==14)
            {
                backwardEnd += hmm->probMY(fState,bState) + msr->indelY(bState);
            }
        }
        else if (bMatch==2)
        {
            if (fMatch==0 || fMatch==3 || fMatch==9)
            {
                backwardEnd += hmm->probXM(fState,bState) + msr->fwdM(bState);
            }
            else if (fMatch==1 || fMatch==7 || fMatch==13)
            {
                backwardEnd += hmm->probYM(fState,bState) + msr->fwdM(bState);
            }
            else if (fMatch==2 || fMatch==5 || fMatch==8 || fMatch==11 || fMatch==14)
            {
                backwardEnd += hmm->probMM(fState,bState) + msr->fwdM(bState);
            }
        }
        else
        {
            cout<<"Hirschberg::error2 ("<<fMatch<<","<<bMatch<<")"<<endl;
        }
    }

    newsite->vitfM(bMatch);
    newsite->vitfS(bState);
    newsite->vitf(forwardEnd);

    newsite->vitbM(fMatch);
    newsite->vitbS(fState);
    newsite->vitb(backwardEnd);


    int K = s1Beg+c.k;

    if (newsite->currMatchState()==0)
    {
        newsite->cInd1(K);
        newsite->cInd2(-1);
        newsite->nInd1(K);
        newsite->nInd2(h);
        newsite->rInd1(K-1);
        newsite->rInd2(h);
        newsite->lInd1(K);
        newsite->lInd2(h); // char (starting!) on left hasn't changed
    }
    else if (newsite->currMatchState()==1)
    {
        newsite->cInd1(-1);
        newsite->cInd2(h);
        newsite->nInd1(K);
        newsite->nInd2(h);
        newsite->rInd1(K);
        newsite->rInd2(h-1);
        newsite->lInd1(K);
        newsite->lInd2(h);
    }
    else if (newsite->currMatchState()==2)
    {
        newsite->cInd1(K);
        newsite->cInd2(h);
        newsite->nInd1(K);
        newsite->nInd2(h);
        newsite->rInd1(K-1); // new char (one over!) on right
        newsite->rInd2(h-1);
        newsite->lInd1(K);   // new char (starting!) on left
        newsite->lInd2(h);
        countSites++;
    }
    else if (newsite->currMatchState()==3 || newsite->currMatchState()==5 || newsite->currMatchState()==9 || newsite->currMatchState()==11)
    {
        newsite->cInd1(K);
        newsite->cInd2(-1);
        newsite->nInd1(-1);
        newsite->nInd2(-1);
        newsite->rInd1(K-1);
        newsite->rInd2(h);
        newsite->lInd1(K);
        newsite->lInd2(h); // char (starting!) on left hasn't changed
        newsite->nullSite(true);
    }
    else if (newsite->currMatchState()==7 || newsite->currMatchState()==8 || newsite->currMatchState()==13 || newsite->currMatchState()==14)
    {
        newsite->cInd1(-1);
        newsite->cInd2(h);
        newsite->nInd1(-1);
        newsite->nInd2(-1);
        newsite->rInd1(K);
        newsite->rInd2(h-1);
        newsite->lInd1(K);
        newsite->lInd2(h);
        newsite->nullSite(true);
    }
    else
    {
        cout<<"Hirschberg: illegal matrix pointer "<<bMatch<<endl;
        exit(1);
    }

    countSites++; // counter to show the percentage aligned

    if (NOISE>1)
    {
        cout<<"Site: ("<<s1<<"-"<<e1<<" ; "<<s2<<"-"<<e2<<")("<<K<<" "<<h<<"); states";
        cout<<fState<<" "<<bState<<" ; "<<fMatch<<" "<<bMatch<<" : "<<maxScore<<endl;
        cout<<"vitf "<<newsite->vitf()<<" "<<newsite->vitfM()<<" "<<newsite->vitfS();
        cout<<"; vitb "<<newsite->vitb()<<" "<<newsite->vitbM()<<" "<<newsite->vitbS()<<endl;
        cout<<"cInd: "<<newsite->cInd1()<<" "<<newsite->cInd2()<<" rInd: "<<newsite->rInd1()<<" "<<newsite->rInd2();
        cout<<" lInd: "<<newsite->lInd1()<<" "<<newsite->lInd2()<<" ; "<<c.k<<"\n"<<endl;
    }
}




bool Hirschberg::rndBool()
{
    if(REPRODUCIBLE)
        srand(random_seed);

    double p = (double)rand()/(double)RAND_MAX;
    if (p>0.5)
        return true;
    else
        return false;
}

int Hirschberg::rndInt(int i)
{
    if(REPRODUCIBLE)
        srand(random_seed);

    return (int)(i*(rand()/(RAND_MAX+1.0)));
}


double Hirschberg::max(double a,double b)
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

double Hirschberg::max(double a,double b, double c)
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
            cout <<"Hirschberg::random number error: i="<<i<<endl;
            exit(1);
        }
    }
}

void Hirschberg::printMatrix(string n,int i,DbMatrix* m)
{
    cout<<n<<i<<": ";
    m->print();
//    m->print(i);
}

void Hirschberg::printMatrix(string n,int i,IntMatrix* m)
{
    cout<<n<<i<<": ";
    m->print();
}

