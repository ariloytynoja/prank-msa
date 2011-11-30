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

#include <cstdio>
#include "guidetree.h"
#include "pwhirschberg.h"
#include "pwsite.h"
#include "config.h"
#include "translatesequences.h"
#include <cstdlib>

using namespace std;

GuideTree::GuideTree(vector<string>* sequences,vector<string>* names,IntMatrix* substScores)
{
    bool isDna = (substScores->X()<=5);
    int ns = sequences->size();

    string full_alphabet = "ARNDCQEGHILKMFPSTWYVX";
    if(isDna)
        full_alphabet = "ACGTURYMKSWHBVDN";

    if (NOISE>=0)
        cout<<"Generating approximate guidetree."<<endl;

    vector<string>::iterator si = sequences->begin();
    int longest = 0;
    int slongest = 0;

    vector<string> local_seqs;

    for (; si!=sequences->end(); si++)
    {
        string seq = *si;

        string::iterator ci = seq.begin();
        for (;ci != seq.end();ci++)
        {
            char c = *ci;
            switch (c)
            {
            case '-':
                seq.erase(ci);
                ci--;
                break;
            default:
                // Remove characters not in full alphabet
                if(full_alphabet.find(c) == string::npos) {
                    seq.erase(ci);
                    ci--;
                }
            }
        }

        local_seqs.push_back(seq);

        if ((int)seq.length()>longest)
        {
            slongest = longest;
            longest = seq.length();
        }
        else if ((int)seq.length()>slongest)
        {
            slongest = seq.length();
        }
    }

    vector<string> *seqs = &local_seqs;

    if (isDna && TRANSLATE)
    {
        TranslateSequences *trseq = new TranslateSequences();
        if (!trseq->translateProtein(names,seqs))
        {
            cout<<"Translation failed. Exiting."<<endl<<endl;
            delete trseq;
            exit(-1);
        }

        delete trseq;
    }

    int delta = int(log( 1-exp(-1.0*pwGapRate*pwDist) )*1000);
    int epsilon = int(log( pwGapExt )*1000);

//     cout<<delta<<" "<<epsilon<<endl;
//     substScores->print();

    if (NOISE>1)
    {
        cout<<"Pairwise gap scoring penalties"<<endl;
        cout<<"open: "<<delta<<", extension: "<<epsilon<<endl;
    }

    PwSite *pws = new PwSite();
    pws->setMatrices(longest,slongest);

    PwHirschberg* pwh = new PwHirschberg(longest);
    pwh->setModel( substScores,delta,epsilon );

    FlMatrix* distance = new FlMatrix(ns,ns,"pw distances");
    distance->initialise(0);

    si = seqs->begin();
    vector<string>::iterator se = seqs->end();
    se--;

    int total = ns*(ns-1)/2;
    int done = 1;

    int i = 0;
    for (; si!=se; si++)
    {
        vector<string>::iterator si2 = si;
        si2++;

        int j = i+1;
        for (; si2!=seqs->end(); si2++)
        {

// 	    cout<<endl<<"unaligned"<<endl;
// 	    cout<<">1"<<endl<<*si<<endl<<">2"<<endl<<*si2<<endl;

            if (SCREEN)
            {
                unsigned int m;
                FOR(m,message.length())
                {
                    cout<<'\b';
                }

                char prop[20];
                sprintf(prop,"%i to %i",i,j);
                message = "aligning ";
                message += prop;
                sprintf(prop,"%i",done*100/total);
                message += " (";
                message += prop;
                message += "%)                    ";

                cout<<message;
                cout.flush();
            }

            pwh->setSequences(&(*si),&(*si2));
            pwh->alignSeqs();

            pws->index(0);
            pws->next();
            int l1 = si->length();
            int l2 = si2->length();

            string a = "";
            string b = "";
            int s = 0;
            int m = 0;
            int s1,s2;
            char c1,c2;

            while (pws->index()!=1)
            {

                c1 = '-';
                c2 = '-';

                s1 = pws->cInd1()-1;
                s2 = pws->cInd2()-1;
                if ( s1>=0 && s1<l1)
                    c1 = si->at(s1);

                if ( s2>=0 && s2<l2)
                    c2 = si2->at(s2);

                if (c1!='-' && c2!='-')
                {
                    s++;
                    if (c1==c2)
                        m++;
                }

// 				cout<<"output:"<<pws->index()<<" "<<c1<<" "<<c2<<endl;

                a+=c1;
                b+=c2;

                pws->next();
            }

            if (NOISE>1)
                cout<<endl<<">1"<<endl<<a<<endl<<">2"<<endl<<b<<endl;

            float p = 1-(float)m/(float)s;
//             float p = 1-(float)m/(float)min(l1,l2);

//               cout<<m<<" "<<s<<" "<<p<<" ";//<<min(l1,l2)<<" "<<1-(float)m/(float)min(l1,l2)<<" ";
            if (CORRECTP)
            {
                if (isDna)
                {
                    if (p>0.7)
                        p=0.9;
                    else
                        p = -0.75*log(1-4/3*p);
                }
                else
                {
                    if (p>0.85)
                        p=2.26;
                    else
                        p = -1*log(1-p-0.2*p*p);
                }
            }

            if (CORRECTP)
            {
                if (isDna && p!=p)
                    p=0.9;
                else if (p!=p)
                    p=2.26;
            }
            else
            {
                if (isDna && p!=p)
                    p=1;
                else if (p!=p)
                    p=1;
            }

//    	    cout<<p<<endl;
            distance->s(p,i,j);
            distance->s(p,j,i);

            j++;
            done++;
        }
        i++;
    }

    delete pwh;
    pws->deleteMatrices();
    delete pws;

    if (SCREEN)
    {
        unsigned int m;
        FOR(m,message.length())
        {
            cout<<'\b';
        }
    }

    if (NOISE>=1)
    {
        cout<<"Pairwise distances"<<endl;
        distance->print();
    }

    this->makeTree(distance,names);
    delete distance;
}

GuideTree::GuideTree(vector<string>* seqs,vector<string>* names,bool isDna)
{

    int ns = seqs->size();
    FlMatrix* distance = new FlMatrix(ns,ns,"pw distances");
    distance->initialise(0);

    if (NOISE>=0)
        cout<<"Generating improved guidetree."<<endl;

    vector<string>::iterator si = seqs->begin();
    vector<string>::iterator se = seqs->end();
    se--;

    int i = 0;
    for (; si!=se; si++)
    {
        vector<string>::iterator si2 = si;
        si2++;

        int j = i+1;
        for (; si2!=seqs->end(); si2++)
        {

            int s = 0;
            int m = 0;
//            int s1,s2;
            char c1,c2;

            for (unsigned int k=0; k<si->length(); k++)
            {

                c1 = si->at(k);
                c2 = si2->at(k);

                if (c1!='-' && c2!='-')
                {
                    s++;
                    if (c1==c2)
                        m++;
                }
            }

            float p = 1-(float)m/(float)s;

            if (CORRECTP)
            {
                if (isDna)
                {
                    if (p>0.7)
                        p=0.9;
                    else
                        p = -0.75*log(1-4/3*p);
                }
                else
                {
                    if (p>0.85)
                        p=2.26;
                    else
                        p = -1*log(1-p-0.2*p*p);
                }
            }

            if (CORRECTP)
            {
                if (isDna && p!=p)
                    p=0.9;
                else if (p!=p)
                    p=2.26;
            }
            else
            {
                if (isDna && p!=p)
                    p=1;
                else if (p!=p)
                    p=1;
            }

            distance->s(p,i,j);
            distance->s(p,j,i);

            j++;
        }
        i++;
    }


    if (NOISE>=1)
    {
        cout<<"Pairwise distances"<<endl;
        distance->print();
    }

    this->makeTree(distance,names);

    delete distance;
}


void GuideTree::makeTree(FlMatrix* distance,vector<string>* nms)
{
    int no = distance->X();

    string* names = new string[no];
    string* newNames = new string[no];

    FlMatrix* newDistance = new FlMatrix(no,no,"pw new distances");
    newDistance->initialise(0);
    FlMatrix* rDist = new FlMatrix(no,"pw rDist"); // sum of distances d_i,j
    rDist->initialise(0);

    vector<string>::iterator ir = nms->begin();
    int i = 0;
    for (; ir!=nms->end(); ir++)
    {
        names[i++] = *ir;
    }

    while (no>2)
    {
        joinNeighbors(distance,names,newDistance,newNames,rDist,&no);
    }

    if (names[0].at(names[0].length()-1) == ')')
    {
        char dist[10];
        sprintf(dist,"%.5f",abs(distance->g(0,0)+distance->g(0,1)) );
        tree = names[0].substr(0,names[0].length()-1)+','+names[1]+':'+dist+");";
    }
    else
    {
        char dist[10];
        sprintf(dist,"%.5f",abs(distance->g(0,0)+distance->g(0,1))/2 );
        tree = '('+names[0]+':'+dist+','+names[1]+':'+dist+");";
    }

    delete[] names;
    delete[] newNames;
    delete newDistance;
    delete rDist;
}

void GuideTree::joinNeighbors(FlMatrix* distance, string* names,FlMatrix* newDistance, string* newNames,FlMatrix* rDist,int* no)
{

    int otu1=0, otu2=0;
    float minM = HUGE_VAL;


    int i,j;

    rDist->initialise(0);
    FOR(i,*no)
    {
        FOR(j,*no)
        {
            rDist->a( distance->g(i,j), i);
        }
    }
    FOR(i,*no)
    {
        FOR(j,*no)
        {
            if (j!=i)
            {
                float mDist = distance->g(i,j)-( rDist->g(i)+rDist->g(j) )/( (*no)-2);
                if (mDist<minM)
                {
                    minM = mDist;
                    otu1 = min(i,j);
                    otu2 = max(i,j);
                }
            }
        }
    }


    double brl1 = distance->g(otu1,otu2)/2+(rDist->g(otu1)-rDist->g(otu2))/(2*(*no-2));
    double brl2 = distance->g(otu1,otu2)-brl1;

    int ci=0;
    int cj=0;
    float v;

    FOR(i,*no)
    {
        FOR(j,*no)
        {
            if (i==j)
            {
                continue;
            }

            if (j<otu2)
            {
                cj=j;
            }
            else if (j==otu2)
            {
                continue;
            }
            else
            {
                cj=j-1;
            }

            if (i<otu2)
            {
                ci=i;
            }
            else if (i==otu2)
            {
                continue;
            }
            else
            {
                ci=i-1;
            }

            if (i==otu1)
            {
                v = ( distance->g(otu1,j)+distance->g(otu2,j)-distance->g(otu1,otu2) )/2;
                newDistance->s(v,ci,cj);
            }
            else if (j==otu1)
            {
                v = ( distance->g(otu1,i)+distance->g(otu2,i)-distance->g(otu1,otu2) )/2;
                newDistance->s(v,ci,cj);
            }
            else
            {
                v = distance->g(i,j);
                newDistance->s(v,ci,cj);
            }
        }
    }

    string s;
    FOR(i,*no)
    {
        if (i==otu1)
        {
            char l1[10];
            sprintf(l1,"%.5f",abs(brl1));
            char l2[10];
            sprintf(l2,"%.5f",abs(brl2));
            newNames[i] = '('+names[otu1]+':'+l1+','+names[otu2]+':'+l2+')';
//	    cout<<*no<<" "<<'('+names[otu1]+':'+l1+','+names[otu2]+':'+l2+')'<<endl;
            continue;
        }
        else if (i<otu2)
        {
            ci = i;
        }
        else if (i==otu2)
        {
            continue;
        }
        else
        {
            ci=i-1;
        }
        newNames[ci] = names[i];
    }

    (*no)--;
    FOR(i,*no)
    {
        names[i] = newNames[i];
        FOR(j,*no)
        {
            distance->s(newDistance->g(i,j),i,j);
        }
    }
}


GuideTree::~GuideTree()
{
}


