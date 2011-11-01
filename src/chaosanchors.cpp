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
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "config.h"
#include "chaosanchors.h"

using namespace std;

//
// This class needs serious re-thinking.....
// WORK HERE
//

ChaosAnchors::ChaosAnchors(string* s1,string* s2)
{
    left = s1;
    right= s2;

    subst = new IntMatrix(5,5,"chaos subst");

    alpha = "ACGTN";

    int tmp[] = {   91, -114,  -31, -123, -43,
                    -114,  100, -125,  -31, -43,
                    -31, -125,  100, -114, -43,
                    -123,  -31, -114,   91, -43,
                    -43,  -43,  -43,  -43, -43
                };

    FOR(i,5)
    {
        FOR(j,5)
        {
            subst->s( tmp[i*5+j], i, j );
        }
    }

    gOpen = -750;
    gExt = -25;

}


ChaosAnchors::~ChaosAnchors()
{
    delete subst;

    if (nanch>0)
    {
        delete anchps;
    }
}

IntMatrix* ChaosAnchors::getAnchors(int *na)
{

    ofstream seqFileL((outfile+".pwlt").c_str());
    seqFileL<<">left"<<endl<<*left<<endl;
    seqFileL.close();
    ofstream seqFileR((outfile+".pwrt").c_str());
    seqFileR<<">right"<<endl<<*right<<endl;
    seqFileR.close();

    vector<int> pairs;

    for (int co=25; co>=15; co-=5)
    {

        pairs.clear();
        char cutoff[4];
        sprintf(cutoff,"%i",co);
        string cmd = "rechaos.pl "+outfile+".pwlt "+outfile+".pwrt -chaos \"-co "+cutoff+"\" -out "+outfile+".anchs 2> /dev/null";
        if (NOISE>0)
            cout<<cmd<<endl;
        int tmp = system(cmd.c_str());

        if (NOISE>0)
            cout<<"cutoff: "<<co<<" ;reading anchors"<<endl;

        bool DISTOK = true;
        ifstream read((outfile+".anchs").c_str());
        string s;

        int j=0;
        while (getline(read,s))
        {
            int beg1 = atoi(s.substr(1,s.find(" ")).c_str());
            int end1 = atoi(s.substr(s.find(" ")+1,s.find(")=")).c_str());
            s=s.substr(s.find("=")+1);
            int beg2 = atoi(s.substr(1,s.find(" ")).c_str());
            int end2 = atoi(s.substr(s.find(" ")+1,s.find(") ")).c_str());

            int c1,c2;
            if (end1-beg1>2*minAnchDist && end2-beg2>2*minAnchDist)
            {
                int pos1,pos2;
                for (pos1=end1-minAnchDist,pos2=end2-minAnchDist; pos1>(beg1+minAnchDist) && pos2>(beg2+minAnchDist); pos1-=minAnchDist,pos2-=minAnchDist)
                {
                    alignRegions(&c1,&c2,pos1,pos1+minAnchDist,pos2,pos2+minAnchDist);
                    pairs.push_back(c1);
                    pairs.push_back(c2);
                    j++;
                    j++;
                    if (NOISE>0)
                        cout<<"fragment anchor*; "<<c1<<" "<<c2<<endl;
                }

                if (pos1>beg1 && pos2>beg2)
                {
                    alignRegions(&c1,&c2,beg1,pos1,beg2,pos2);
                    if (c1>minAnchDist && c2>minAnchDist)
                    {
                        pairs.push_back(c1);
                        pairs.push_back(c2);
                        j++;
                        j++;
                        if (NOISE>0)
                            cout<<"last anchors*; "<<c1<<" "<<c2<<endl;
                    }
                }
            }
            else
            {
                alignRegions(&c1,&c2,beg1,end1,beg2,end2);
                if (j==0 && (int)left->length()-c1>minAnchDist && (int)right->length()-c2>minAnchDist)
                {
                    pairs.push_back(c1);
                    pairs.push_back(c2);
                    j++;
                    j++;
                    if (NOISE>0)
                        cout<<"first anchor; "<<c1<<" "<<c2<<endl;
                }
                else if (j>0 && pairs.at(j-2)-c1>minAnchDist && pairs.at(j-1)-c2>minAnchDist)
                {
                    pairs.push_back(c1);
                    pairs.push_back(c2);
                    j++;
                    j++;
                    if (NOISE>0)
                        cout<<"other anchors; "<<c1<<" "<<c2<<endl;
                }
                else
                {
                    if (NOISE>0)
                        cout<<"skipped anchor; "<<c1<<" "<<c2<<endl;
                }
            }

            if (j>2 && pairs.at(j-4)-pairs.at(j-2)>maxAnchDist && pairs.at(j-3)-pairs.at(j-1)>maxAnchDist)
            {
                if (NOISE>0)
                    cout<<"anchors too distant: "<<pairs.at(j-2)<<"; "<<pairs.at(j-1)<<endl;
                DISTOK=false;
            }
        }
        if (j>0 && pairs.at(j-2)>maxAnchDist && pairs.at(j-1)>maxAnchDist)
        {
            if (NOISE>0)
                cout<<"first anchors too distant: "<<pairs.at(j-2)<<"; "<<pairs.at(j-1)<<endl;
            DISTOK=false;
        }

        if (j>0 && (int)left->length()-pairs.at(0)>maxAnchDist && (int)right->length()-pairs.at(1)>maxAnchDist)
        {
            if (NOISE>0)
                cout<<"last anchors too distant: "<<pairs.at(1)<<"; "<<pairs.at(0)<<endl;
            DISTOK=false;
        }

        *na=j/2;

        // only one round of anchoring
//	if(DISTOK)
        break;

    }

    if (fopen((outfile+".anchs").c_str(),"r")>0)
        rename((outfile+".anchs").c_str(),(outfile+".anchs_old").c_str());
    if (fopen((outfile+".pwlt").c_str(),"r")>0)
        rename((outfile+".pwlt").c_str(),(outfile+".pwlt_old").c_str());
    if (fopen((outfile+".pwrt").c_str(),"r")>0)
        rename((outfile+".pwrt").c_str(),(outfile+".pwrt_old").c_str());

    if (*na>0)
    {
        anchps = new IntMatrix(2,*na,"chaos anchorpairs");

        j=0;
        FOR(i,*na)
        {
            anchps->s( pairs.at(j), 0, i );
            j++;
            anchps->s( pairs.at(j), 1, i );
            j++;
            if (NOISE>0)
                cout<<"anchor "<<i<<" "<<anchps->g(0,i)<<" "<<anchps->g(1,i)<<endl;
        }
        pairs.clear();
    }

    nanch = *na;

    return anchps;
}

void ChaosAnchors::alignRegions(int *c1,int *c2,int beg1,int end1,int beg2,int end2)
{
    //small = (int)-HUGE_VAL;
    small = -1000000;

    int len1 = end1-beg1+1;
    int len2 = end2-beg2+1;

    IntMatrix* m1X = new IntMatrix(len1);
    IntMatrix* m1Y = new IntMatrix(len1);
    IntMatrix* m1M = new IntMatrix(len1);
    IntMatrix* m2X = new IntMatrix(len1);
    IntMatrix* m2Y = new IntMatrix(len1);
    IntMatrix* m2M = new IntMatrix(len1);
    IntMatrix* m3X = new IntMatrix(len1);
    IntMatrix* m3Y = new IntMatrix(len1);
    IntMatrix* m3M = new IntMatrix(len1);

    IntMatrix* tX;
    IntMatrix* tY;
    IntMatrix* tM;

    IntMatrix* pX = m1X;
    IntMatrix* pY = m1Y;
    IntMatrix* pM = m1M;
    IntMatrix* cfX = m2X;
    IntMatrix* cfY = m2Y;
    IntMatrix* cfM = m2M;
//    IntMatrix* cbX = m3X;
//    IntMatrix* cbY = m3Y;
//    IntMatrix* cbM = m3M;

    int match;
    int maxScore = small;
    int maxi =0;
    int maxj =0;

    int t;

    FOR(j,len2)
    {

        FOR(i,len1)
        {

            if (i==0 && j==0)
            {
                cfX->s(0,i);
                cfY->s(0,i);
                cfM->s(0,i);
                continue;
            }

            if (i==0)
            {
                cfY->s(0,i);
                cfX->s(small,i);
                cfM->s(small,i);
                continue;
            }

            if (j==0)
            {
                cfX->s(0,i);
                cfY->s(small,i);
                cfM->s(small,i);
                continue;
            }

            match=matchScore(beg1+i-1,beg2+j-1);

            t = max(cfX->g(i-1)+gExt,
                    cfY->g(i-1)+gOpen,
                    cfM->g(i-1)+gOpen);
            cfX->s(t,i);
            t = max(pX->g(i)+gOpen,
                    pY->g(i)+gExt,
                    pM->g(i)+gOpen);
            cfY->s(t,i);
            t = max(pX->g(i-1),
                    pY->g(i-1),
                    pM->g(i-1))+match;
            cfM->s(t,i);

            if (i==len1-1)
            {
                t = max(pX->g(i)+gOpen,
                        pY->g(i),
                        pM->g(i)+gOpen);
                cfY->s(t,i);
            }

            if (j==len2-1)
            {
                t = max(cfX->g(i-1),
                        cfY->g(i-1)+gOpen,
                        cfM->g(i-1)+gOpen);
                cfX->s(t,i);
            }

            if (cfM->g(i)>maxScore)
            {
                maxScore=cfM->g(i);
                maxi=i;
                maxj=j;
            }
        }

        tX = pX;
        tY = pY;
        tM = pM;

        pX=cfX;
        pY=cfY;
        pM=cfM;

        cfX=tX;
        cfY=tY;
        cfM=tM;

    }


    *c1 = beg1+maxi-1;
    *c2 = beg2+maxj-1;

    if (NOISE>0)
        cout<<"pair "<<*c1<<" "<<*c2<<"; "<<left->at(*c1)<<" "<<right->at(*c2)<<" ("<<maxScore<<"); "<<beg1<<" "<<end1<<" "<<len1<<"; "<<beg2<<" "<<end2<<" "<<len2<<endl;

    delete m1X;
    delete m1Y;
    delete m1M;
    delete m2X;
    delete m2Y;
    delete m2M;
    delete m3X;
    delete m3Y;
    delete m3M;

}

void ChaosAnchors::reverseAlignRegions(int *c1,int *c2,int beg1,int end1,int beg2,int end2)
{
    small = (int)-HUGE_VAL;

    int len1 = end1-beg1+1;
    int len2 = end2-beg2+1;

    IntMatrix* m1X = new IntMatrix(len1);
    IntMatrix* m1Y = new IntMatrix(len1);
    IntMatrix* m1M = new IntMatrix(len1);
    IntMatrix* m2X = new IntMatrix(len1);
    IntMatrix* m2Y = new IntMatrix(len1);
    IntMatrix* m2M = new IntMatrix(len1);
    IntMatrix* m3X = new IntMatrix(len1);
    IntMatrix* m3Y = new IntMatrix(len1);
    IntMatrix* m3M = new IntMatrix(len1);

    IntMatrix* tX;
    IntMatrix* tY;
    IntMatrix* tM;

    IntMatrix* pX = m1X;
    IntMatrix* pY = m1Y;
    IntMatrix* pM = m1M;
    IntMatrix* cfX = m2X;
    IntMatrix* cfY = m2Y;
    IntMatrix* cfM = m2M;
//    IntMatrix* cbX = m3X;
//    IntMatrix* cbY = m3Y;
//    IntMatrix* cbM = m3M;

    int match;
    int maxScore = small;
    int maxi =0;
    int maxj =0;

    int t;

    RFOR(j,len2-1)
    {

        RFOR(i,len1-1)
        {

            if (i==len1-1 && j==len2-1)
            {
                cfX->s(0,i);
                cfY->s(0,i);
                cfM->s(0,i);
                continue;
            }

            if (i==len1-1)
            {
                cfY->s(0,i);
                cfX->s(small,i);
                cfM->s(small,i);
                continue;
            }

            if (j==len2-1)
            {
                cfX->s(0,i);
                cfY->s(small,i);
                cfM->s(small,i);
                continue;
            }

            match=matchScore(beg1+i,beg2+j);

            t = max(cfX->g(i+1)+gExt,
                    cfY->g(i+1)+gOpen,
                    cfM->g(i+1)+gOpen);
            cfX->s(t,i);
            t = max(pX->g(i)+gOpen,
                    pY->g(i)+gExt,
                    pM->g(i)+gOpen);
            cfY->s(t,i);
            t = max(pX->g(i+1),
                    pY->g(i+1),
                    pM->g(i+1))+match;
            cfM->s(t,i);

            if (i==0)
            {
                t = max(pX->g(i)+gOpen,
                        pY->g(i),
                        pM->g(i)+gOpen);
                cfY->s(t,i);
            }

            if (j==0)
            {
                t = max(cfX->g(i+1),
                        cfY->g(i+1)+gOpen,
                        cfM->g(i+1)+gOpen);
                cfX->s(t,i);
            }

            if (cfM->g(i)>maxScore)
            {
                maxScore=cfM->g(i);
                maxi=i;
                maxj=j;
            }
        }

        tX = pX;
        tY = pY;
        tM = pM;

        pX=cfX;
        pY=cfY;
        pM=cfM;

        cfX=tX;
        cfY=tY;
        cfM=tM;

    }


    *c1 = beg1+maxi-1;
    *c2 = beg2+maxj-1;

    cout<<"rev.pair "<<*c1<<" "<<*c2<<"; "<<left->at(*c1)<<" "<<right->at(*c2)<<" ("<<maxScore<<")"<<endl;

    delete m1X;
    delete m1Y;
    delete m1M;
    delete m2X;
    delete m2Y;
    delete m2M;
    delete m3X;
    delete m3Y;
    delete m3M;

}

int ChaosAnchors::matchScore(int site1,int site2)
{
    if (site1<0 || site1 > (int)left->length())
        cout<<"ChaosAnchors: site1 "<<site1<<" ("<<left->length()<<")"<<endl;
    int i1=alpha.find(left->at(site1));
    if (i1<0 || i1>3)
    {
        i1=4;
    }

    if (site2<0 || site2> (int)right->length())
        cout<<"ChaosAnchors: site2 "<<site2<<" ("<<right->length()<<")"<<endl;
    int i2=alpha.find(right->at(site2));
    if (i2<0 || i2>3)
    {
        i2=4;
    }

    return subst->g(i1,i2);
}


int ChaosAnchors::max(int a,int b)
{
    if (a<=0 && b<=0)
    {
        return 0;
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

int ChaosAnchors::max(int a,int b, int c)
{
    if (a<=0 && b<=0 && c<=0)
    {
        return 0;
    }
    else if (a>b && a>c)
    {
        return a;
    }
    else if (a<b && b>c)
    {
        return b;
    }
    else if (a<c && b<c)
    {
        return c;
    }
    else if (a>b && a==c)
    {
        if (rndBool())
        {
            return a;
        }
        else
        {
            return c;
        }
    }
    else if (a>c && a==b)
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
    else if (a<b && b==c)
    {
        if (rndBool())
        {
            return b;
        }
        else
        {
            return c;
        }
    }
    else
    {
        int i = rndInt(3);
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
            cout <<"ChaosAnchors::random number error: i="<<i<<endl;
            exit(1);
        }
    }
}

bool ChaosAnchors::rndBool()
{
    double p = (double)rand()/(double)RAND_MAX;
    if (p>0.5)
        return true;
    else
        return false;
}

int ChaosAnchors::rndInt(int i)
{
    return (int)(rand()/(double)RAND_MAX*(float)i);
}
