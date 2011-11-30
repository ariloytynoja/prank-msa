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
#include "config.h"
#include "pwhirschberg.h"
#include "exonerate_reads.h"
#include "pwsite.h"

using namespace std;

PwHirschberg::~PwHirschberg()
{
    cleanUp();
}

void PwHirschberg::cleanUp()
{

//     cout<<"cleanUp()"<<endl;
    // clean up
    delete fVM1;
    delete fVX1;
    delete fVY1;

    delete fVM2;
    delete fVX2;
    delete fVY2;

    delete bVM1;
    delete bVX1;
    delete bVY1;

    delete bVM2;
    delete bVX2;
    delete bVY2;

    delete ptVM;
    delete ptVX;
    delete ptVY;

    delete beg;
    delete end;
    delete pwsite;

}

PwHirschberg::PwHirschberg(int sl1)
{
    small = -100000000;
    sl1 = sl1+1;

    // Initialize matrices
    //
    fVX1 = new IntMatrix(sl1,"fVX1");  // matrices for fwd Viterbi scores
    fVY1 = new IntMatrix(sl1,"fVY1");
    fVM1 = new IntMatrix(sl1,"fVM1");

    fVX1->initialise(small);
    fVY1->initialise(small);
    fVM1->initialise(small);

    fVX2 = new IntMatrix(sl1,"fVX2");  // matrices for fwd Viterbi scores
    fVY2 = new IntMatrix(sl1,"fVY2");
    fVM2 = new IntMatrix(sl1,"fVM2");

    fVX2->initialise(small);
    fVY2->initialise(small);
    fVM2->initialise(small);

    bVM1 = new IntMatrix(sl1,"bVM1");  // matrices for bwd Viterbi scores
    bVX1 = new IntMatrix(sl1,"bVX1");
    bVY1 = new IntMatrix(sl1,"bVY1");

    bVX1->initialise(small);
    bVY1->initialise(small);
    bVM1->initialise(small);

    bVM2 = new IntMatrix(sl1,"bVM2");  // matrices for bwd Viterbi scores
    bVX2 = new IntMatrix(sl1,"bVX2");
    bVY2 = new IntMatrix(sl1,"bVY2");

    bVX2->initialise(small);
    bVY2->initialise(small);
    bVM2->initialise(small);

    ptVM = new IntMatrix(sl1,"ptVM");   // matrices for the backward pointers
    ptVX = new IntMatrix(sl1,"ptVX");
    ptVY = new IntMatrix(sl1,"ptVY");

    ptVX->initialise(-1);
    ptVY->initialise(-1);
    ptVM->initialise(-1);

    beg = new PwSite(0);
    end = new PwSite(1);
    pwsite = new PwSite();

}

void PwHirschberg::setModel( IntMatrix* scores,int delta, int epsilon )
{
    substScores = scores;
    deltaX1 = deltaX2 = deltaY1 = deltaY2 = delta;
    epsilonX = epsilonY = epsilon;

    if (scores->X() > 5)
    {
        sAlpha = 20;
        alpha = "ARNDCQEGHILKMFPSTWYV";
    }
    else
    {
        sAlpha = 4;
        alpha = "ACGT";
    }
}


int PwHirschberg::count = 2;
int PwHirschberg::depth = 0;

void PwHirschberg::setSequences(string* s1,string* s2)
{
    seq1 = s1;
    seq2 = s2;

    sl1 = s1->length();
    sl2 = s2->length();

    pwsite->resetCounter();

//    defineBegin();
//    defineEnd();

    count = 2;
}

void PwHirschberg::alignSeqs()
{

    totalSites = sl1+sl2;
    countSites = 0;

    defineBegin();

    if (EXONERATE)
    {
         this->getAnchors();
    }

    defineEnd();

    divideSeq();

}


void PwHirschberg::getAnchors()
{

    if (NOISE>0)
        cout<<"lengths: "<<sl1<<" and "<<sl2<<" ("<<anchSkipDist<<")"<<endl;

    if (sl1>anchSkipDist && sl2>anchSkipDist)
    {

        vector<hit> exonerate_hits;
        Exonerate_reads er;
        if (er.test_executable())
        {
            er.local_alignment(seq1,seq2,&exonerate_hits, true);
        }

        vector<pair<int,int> > anchor_pairs;

        for (int i=0; i<exonerate_hits.size(); i++)
        {
            hit h = exonerate_hits.at(i);
            if (NOISE>0)
                cout<<"e "<<h.query<<" "<<h.node<<" "<<h.score<<" "<<h.q_start<<" "<<h.q_end<<" "<<h.q_strand<<" "<<h.t_start<<" "<<h.t_end<<" "<<h.t_strand<<"\n";

            for (int j=0; j+h.q_start<h.q_end; j+=10)
            {
                if (j+h.q_start>5 && j+h.t_start>5)
                    anchor_pairs.push_back(make_pair(j+h.q_start,j+h.t_start));
            }
        }

        if (anchor_pairs.size()>0)
        {
            for (int i=0; i<anchor_pairs.size(); i++)
            {

                if (NOISE>0)
                    cout<<" ex anchor "<<anchor_pairs.at(i).first<<","<<anchor_pairs.at(i).second<<" *"<<endl;


                defineESite(anchor_pairs.at(i).first,anchor_pairs.at(i).second);

                if (NOISE>0)
                {
                    cout<<" beg: "<<beg->index()<<" "<<beg->lInd1()<<" "<<beg->lInd2()<<" | ";
                    cout<<" anc: "<<end->index()<<" "<<end->lInd1()<<" "<<end->lInd2()<<endl;
                }

                divideSeq();
            }
        }

        if (NOISE>1)
        {
            cout<<" beg: "<<beg->index()<<" "<<beg->lInd1()<<" "<<beg->lInd2()<<" | ";
            cout<<" end: "<<end->index()<<" "<<end->lInd1()<<" "<<end->lInd2()<<endl;
        }
    }

}

void PwHirschberg::defineBegin()
{
    beg->index(0);
    beg->currMatchState(-1);

    beg->cInd1(0);
    beg->cInd2(0);
    beg->rInd1(-1);
    beg->rInd2(-1); // before the start
    beg->lInd1(0);
    beg->lInd2(0);

    beg->vitfX(deltaX1);
    beg->vitfY(deltaX1);
    /*    beg->vitfX(0);
        beg->vitfY(0);*/
    beg->vitfM(0);

    beg->vitbX(small);
    beg->vitbY(small);
    beg->vitbM(small);

}

void PwHirschberg::defineESite(int l,int r)
{
    end->index(1);
//    end->isAnchor(true);
//    end->nullSite(true);

    end->cInd1(l);
    end->lInd1(l);
    end->cInd2(r);
    end->lInd2(r);
    end->rInd1(l-1);
    end->rInd2(r-1);

    end->vitfX(small);
    end->vitfY(small);
    end->vitfM(small);

    end->vitbX(deltaX1);
    end->vitbY(deltaX1);
    end->vitbM(0);
}

void PwHirschberg::defineEnd()
{
    end->index(1);
    end->cInd1(-1);
    end->cInd2(-1);
    end->rInd1(seq1->length());
    end->rInd2(seq2->length());
    end->lInd1(-1);
    end->lInd2(-1); // over the end

    end->vitfX(small);
    end->vitfY(small);
    end->vitfM(small);

    end->vitbX(deltaX1);
    end->vitbY(deltaX1);
    /*    end->vitbX(0);
        end->vitbY(0);*/
    end->vitbM(0);

}

void PwHirschberg::divideSeq()
{
//  	depth++;

    getMidSite(beg->lInd1(),end->rInd1(),beg->lInd2(),end->rInd2());
    pwsite->setNeighbours(beg,end);

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

// 	depth--;
}


void PwHirschberg::getMidSite(int s1,int e1,int s2,int e2)
{

    int h = (s2+e2)/2+1;   // midpoint
    if (s2==e2)            // exception for zero-length seq2
        h = (s2+e2)/2;

    int s1Beg = s1;        // seq1 start site
    int s1Len = e1-s1;     // seq1 length

    int s2bBeg = s2;       // seq2_begin start site
    int s2bLen = h-s2;     // seq2_begin length
    if (s2==e2)
        s2bLen=0;

    int s2eBeg = h;        // seq2_end start site
    int s2eLen = e2-(h+1); // seq2_end length

    mLen = s1Len+1;        // short cuts

    // Define pointers for current & previous row
    //
    cfVX = fVX1;
    cfVY = fVY1;
    cfVM = fVM1;

    pVX = fVX2;
    pVY = fVY2;
    pVM = fVM2;

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

                // set starting values
                cfVX->s(beg->vitfX(),j);
                cfVY->s(beg->vitfY(),j);
                cfVM->s(beg->vitfM(),j);

                ptVX->s(0,j);
                ptVY->s(1,j);
                ptVM->s(2,j);

                continue;
            }

            // Compute the substitution prices
            //
            this->computeFwd( s1Beg + j, s2bBeg + i );

            if (i==0 && j>0)   // only X-gaps are possible
            {

                sX=sY=sM=-1;
                cY=cM=small;


                cX = max(cfVX->g(j-1) + epsilonX,
                         cfVY->g(j-1) + deltaY2+deltaX1,
                         cfVM->g(j-1) + deltaX1);
                sX = maxIndex;


                cfVX->s(cX,j);
                cfVY->s(cY,j);
                cfVM->s(cM,j);

                ptVX->s(sX,j);
                ptVY->s(sY,j);
                ptVM->s(sM,j);

            }
            else if (i>0 && j==0)   // only Y-gaps are possible
            {

                sX=sY=sM=-1;
                cX=cM=small;

                cY = max(pVX->g(j) + deltaY1+deltaX2,
                         pVY->g(j) + epsilonY,
                         pVM->g(j) + deltaY1);
                sY = maxIndex;

                cfVX->s(cX,j);
                cfVY->s(cY,j);
                cfVM->s(cM,j);

                ptVX->s(sX,j);
                ptVY->s(sY,j);
                ptVM->s(sM,j);

            }
            else    // so far, the moves have been exceptional; from now on they are "normal"
            {


                cX = max(cfVX->g(j-1) + epsilonX,
                         cfVY->g(j-1) + deltaY2+deltaX1,
                         cfVM->g(j-1) + deltaX1);
                sX = maxIndex;

                cY = max(pVX->g(j) + deltaY1+deltaX2,
                         pVY->g(j) + epsilonY,
                         pVM->g(j) + deltaY1);
                sY = maxIndex;

                cM = max(pVX->g(j-1) + deltaX2,
                         pVY->g(j-1) + deltaY2,
                         pVM->g(j-1)) + matchScore;
                sM = maxIndex;

                cfVX->s(cX,j);
                cfVY->s(cY,j);
                cfVM->s(cM,j);

                ptVX->s(sX,j);
                ptVY->s(sY,j);
                ptVM->s(sM,j);

            }
        } // FOR(j,mLen)


        // change the rows that are pointed
        tmpVX = pVX;
        tmpVY = pVY;
        tmpVM = pVM;

        pVX = cfVX;
        pVY = cfVY;
        pVM = cfVM;

        cfVX = tmpVX;
        cfVY = tmpVY;
        cfVM = tmpVM;

        if (NOISE>2)
        {
            printMatrix("fM",mLen,pVM);
            printMatrix("fX",mLen,pVX);
            printMatrix("fY",mLen,pVY);
        }
    }

    // change the pointers back so "previous" can be recycled
    // and the mid-row calculation is correct
    cfVX = pVX;
    cfVY = pVY;
    cfVM = pVM;

    // Define pointers for current & previous row
    //
    cbVX = bVX1;
    cbVY = bVY1;
    cbVM = bVM1;

    pVX = bVX2;
    pVY = bVY2;
    pVM = bVM2;

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

                    // starting values
                    cbVX->s(end->vitbX(),j);
                    cbVY->s(end->vitbY(),j);
                    cbVM->s(end->vitbM(),j);

                    continue;
                }

                // Compute the substitution prices
                //
                this->computeBwd( s1Beg+j, s2eBeg+i );

                if (i<s2eLen+1 && j==s1Len)   // y-gaps are possible
                {

                    cX = deltaY1+deltaX2 + pVY->g(j);
                    cY = epsilonY + pVY->g(j);
                    cM = deltaY1 + pVY->g(j);

                    cbVX->s(cX,j);
                    cbVY->s(cY,j);
                    cbVM->s(cM,j);


                }
                else if (i==s2eLen+1 && j<s1Len)   // x-gaps are possible
                {

                    cX = epsilonX + cbVX->g(j+1);
                    cY = deltaY2+deltaX1 + cbVX->g(j+1);
                    cM = deltaX1 + cbVX->g(j+1);

                    cbVX->s(cX,j);
                    cbVY->s(cY,j);
                    cbVM->s(cM,j);

                }
                else if (i<s2eLen+1 && j<s1Len)   // also matches are possible
                {

                    cX = max(epsilonX + cbVX->g(j+1),
                             deltaY1+deltaX2 + pVY->g(j),
                             deltaX2 + matchScore + pVM->g(j+1));

                    cY = max(deltaY2+deltaX1 + cbVX->g(j+1),
                             epsilonY + pVY->g(j),
                             deltaY2 + matchScore + pVM->g(j+1));

                    cM = max(deltaX1 + cbVX->g(j+1),
                             deltaY1 + pVY->g(j),
                             matchScore + pVM->g(j+1));

                    cbVX->s(cX,j);
                    cbVY->s(cY,j);
                    cbVM->s(cM,j);

                }
            }	/// RFOR(j,s1Len)

            // change the rows that are pointed
            tmpVX = pVX;
            tmpVY = pVY;
            tmpVM = pVM;

            pVX = cbVX;
            pVY = cbVY;
            pVM = cbVM;

            cbVX = tmpVX;
            cbVY = tmpVY;
            cbVM = tmpVM;

            if (NOISE>2)
            {
                printMatrix("bM",mLen,pVM);
                printMatrix("bX",mLen,pVX);
                printMatrix("bY",mLen,pVY);
            }
        }

        // change the pointers back so the mid-row calculation is correct
        cbVX = pVX;
        cbVY = pVY;
        cbVM = pVM;
    }

    // Cases where only x-gaps possible
    //
    if (s2==e2)
    {

        // Starting: set the corner values
        //
        // starting values
        cbVX->s(end->vitbX(),s1Len);
        cbVY->s(end->vitbY(),s1Len);
        cbVM->s(end->vitbM(),s1Len);

        RFOR(j,s1Len-1)
        {

            // Compute the substitution prices
            //
            this->computeBwd( s1Beg+j, s2eBeg );

            // move into X-matrix
            //

            cX = epsilonX + cbVX->g(j+1);
            cY = deltaY2+deltaX1 + cbVX->g(j+1);
            cM = deltaX1 + cbVX->g(j+1);

            cbVX->s(cX,j);
            cbVY->s(cY,j);
            cbVM->s(cM,j);
        }

        if (NOISE>2)
        {
            printMatrix("BM",mLen,cbVM);
            printMatrix("BX",mLen,cbVX);
            printMatrix("BY",mLen,cbVY);
        }
    }

    // Find k (i.e. the column through which the alignment path goes)
    //
    vector<Cell> maxCell;
    int maxScore = small;

    j=0;

    if (s2==e2)
        j++;

    for (; j<mLen; j++)
    {

        int tmp;

        tmp = cfVY->g(j)+cbVY->g(j);
        if (tmp>maxScore)
        {

            maxScore = tmp;
            maxCell.clear();
            Cell c = {ptVY->g(j),1,j};
            maxCell.push_back(c);

        }
        else if (tmp==maxScore)
        {

            Cell c = {ptVY->g(j),1,j};
            maxCell.push_back(c);

        }

        if (s2<e2 || j>0)
        {

            tmp = cfVX->g(j)+cbVX->g(j);
            if (tmp>maxScore)
            {

                maxScore = tmp;
                maxCell.clear();
                Cell c = {ptVX->g(j),0,j};
                maxCell.push_back(c);

            }
            else if (tmp==maxScore)
            {

                Cell c = {ptVX->g(j),0,j};
                maxCell.push_back(c);

            }


            tmp = cfVM->g(j)+cbVM->g(j);
            if (tmp>maxScore)
            {

                maxScore = tmp;
                maxCell.clear();
                Cell c = {ptVM->g(j),2,j};
                maxCell.push_back(c);

            }
            else if (tmp==maxScore)
            {

                Cell c = {ptVM->g(j),2,j};
                maxCell.push_back(c);

            }

        }
    }

    if (maxScore == small)
    {
        cout<<endl<<"Pairwise alignment score below the minimum limit. Guidetree failed. Exiting."<<endl<<endl;
        exit(-1);
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
    pwsite->addNewSite();

    this->computeFwd( s1Beg+c.k , s2eBeg );

    int fMatch = c.prev;
    int bMatch = c.curr;

    pwsite->currMatchState(bMatch);

    int forwardEnd = small;
    int backwardEnd = small;

    if (bMatch==0)
    {
        forwardEnd = cfVX->g(c.k);
        backwardEnd = cbVX->g(c.k);
    }
    else if (bMatch==1)
    {
        forwardEnd = cfVY->g(c.k);
        backwardEnd = cbVY->g(c.k);
    }
    else if (bMatch==2)
    {
        forwardEnd = cfVM->g(c.k);
        backwardEnd = cbVM->g(c.k);
    }
    else
    {
        cout<<"PwHirschberg::error1 ("<<fMatch<<","<<bMatch<<")"<<endl;
    }

    if (bMatch==0 || bMatch==1 || bMatch==2)
    {
        if (bMatch==0)
        {
            if (fMatch==0)
            {
                backwardEnd += epsilonX;
            }
            else if (fMatch==1)
            {
                backwardEnd += deltaY2+deltaX1;
            }
            else if (fMatch==2)
            {
                backwardEnd += deltaX1;
            }
            else
            {
                cout<<"PwHirschberg::error2A ("<<fMatch<<","<<bMatch<<")"<<endl;
            }
        }
        else if (bMatch==1)
        {
            if (fMatch==0)
            {
                backwardEnd += deltaY2+deltaX1;
            }
            else if (fMatch==1)
            {
                backwardEnd += epsilonY;
            }
            else if (fMatch==2)
            {
                backwardEnd += deltaY1;
            }
            else
            {
                cout<<"PwHirschberg::error2B ("<<fMatch<<","<<bMatch<<")"<<endl;
            }
        }
        else if (bMatch==2)
        {
            if (fMatch==0)
            {
                backwardEnd += deltaX2 + matchScore;
            }
            else if (fMatch==1)
            {
                backwardEnd += deltaY2 + matchScore;
            }
            else if (fMatch==2)
            {
                backwardEnd += matchScore;
            }
            else
            {
                cout<<"PwHirschberg::error2C ("<<fMatch<<","<<bMatch<<")"<<endl;
            }
        }
        else
        {
            cout<<"PwHirschberg::error2 ("<<fMatch<<","<<bMatch<<")"<<endl;
        }
    }

    if (bMatch==0)
    {
        pwsite->vitfX(forwardEnd);
        pwsite->vitfY(small);
        pwsite->vitfM(small);
    }
    else if (bMatch==1)
    {
        pwsite->vitfX(small);
        pwsite->vitfY(forwardEnd);
        pwsite->vitfM(small);
    }
    else if (bMatch==2)
    {
        pwsite->vitfX(small);
        pwsite->vitfY(small);
        pwsite->vitfM(forwardEnd);
    }
    else
    {
        cout<<"PwHirschberg::error3 ("<<fMatch<<","<<bMatch<<")"<<endl;
    }

    if (fMatch==0)
    {
        pwsite->vitbX(backwardEnd);
        pwsite->vitbY(small);
        pwsite->vitbM(small);
    }
    else if (fMatch==1)
    {
        pwsite->vitbX(small);
        pwsite->vitbY(backwardEnd);
        pwsite->vitbM(small);
    }
    else if (fMatch==2)
    {
        pwsite->vitbX(small);
        pwsite->vitbY(small);
        pwsite->vitbM(backwardEnd);
    }
    else
    {
        cout<<"PwHirschberg::error4 ("<<fMatch<<","<<bMatch<<") ("<<s1<<","<<e1<<") ("<<s2<<","<<e2<<") ("<<s1Beg+c.k<<", "<<h<<")"<<" "<<count<<endl;
        cout<<fwdvM<<" "<<fwdvX<<" "<<fwdvY<<endl;
        cout<<bwdvM<<" "<<bwdvX<<" "<<bwdvY<<endl;
    }

    int K = s1Beg+c.k;

    if (pwsite->currMatchState()==0)
    {
        pwsite->cInd1(K);
        pwsite->cInd2(-1);
        pwsite->rInd1(K-1);
        pwsite->rInd2(h);
        pwsite->lInd1(K);
        pwsite->lInd2(h); // char (starting!) on left hasn't changed
    }
    else if (pwsite->currMatchState()==1)
    {
        pwsite->cInd1(-1);
        pwsite->cInd2(h);
        pwsite->rInd1(K);
        pwsite->rInd2(h-1);
        pwsite->lInd1(K);
        pwsite->lInd2(h);
    }
    else if (pwsite->currMatchState()==2)
    {
        pwsite->cInd1(K);
        pwsite->cInd2(h);
        pwsite->rInd1(K-1); // new char (one over!) on right
        pwsite->rInd2(h-1);
        pwsite->lInd1(K);   // new char (starting!) on left
        pwsite->lInd2(h);
        countSites++;
    }
    else
    {
        cout<<"PwHirschberg: illegal matrix pointer "<<bMatch<<endl;
        exit(1);
    }

    countSites++; // counter to show the percentage aligned

    if (NOISE>1)
    {
        cout<<"PwSite: ("<<s1<<"-"<<e1<<" ; "<<s2<<"-"<<e2<<")("<<K<<" "<<h<<"); states";
        cout<<" ; "<<fMatch<<" "<<bMatch<<" : "<<maxScore<<endl;
    }

}




bool PwHirschberg::rndBool()
{
    int p = (int)rand()/(int)RAND_MAX;
    if (p>0.5)
        return true;
    else
        return false;
}

int PwHirschberg::rndInt(int i)
{
    return (int)(i*(rand()/(RAND_MAX+1.0)));
}


int PwHirschberg::max(int a,int b)
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

int PwHirschberg::max(int a,int b, int c)
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
            cout <<"PwHirschberg::random number error: i="<<i<<endl;
            exit(1);
        }
    }
}

void PwHirschberg::printMatrix(string n,int i,IntMatrix* m)
{
    cout<<n<<i<<": ";
    m->print();
}

void PwHirschberg::computeFwd(int j,int i)
{
// cout<<"j="<<j<<",i="<<i<<endl;
    matchScore = 0;
    if (j==0 || i==0)
        return;

    int c1 = alpha.find(toupper(seq1->at(j-1)));
    if (c1<0)
        c1 = sAlpha;

    int c2 = alpha.find(toupper(seq2->at(i-1)));
    if (c2<0)
        c2 = sAlpha;

    matchScore = substScores->g(c1,c2);

    return;
}


void PwHirschberg::computeBwd(int j,int i)
{
    matchScore = 0;
    if (j==sl1 || i==sl2)
        return;

    int c1 = alpha.find(toupper(seq1->at(j)));
    if (c1<0)
        c1 = sAlpha;

    int c2 = alpha.find(toupper(seq2->at(i)));
    if (c2<0)
        c2 = sAlpha;

    matchScore = substScores->g(c1,c2);

    return;
}
