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

#ifndef RFOR
#define RFOR(i,n) for(i=n; i>=0; i--)
#endif
#ifndef FOR
#define FOR(i,n) for(i=0; i<n; i++)
#endif

#ifndef HMMODEL_H
#define HMMODEL_H

#include <iostream>
#include <string>
#include "dbmatrix.h"
#include "flmatrix.h"
#include "intmatrix.h"
#include "ancestralnode.h"

extern bool LOGVALUES;

class HMModel
{
    int as,fas;           // size of alphabet / full alphabet
    std::string alphabet, fullAlphabet; // alphabet as string
    int sn;           // number of structure states

    DbMatrix* cPi;     // character background frequencies
    DbMatrix* nPi;     // character background frequencies for null model
    DbMatrix* logcPi;  // character background frequencies
    DbMatrix* lognPi;  // character background frequencies for null model
    DbMatrix* wRoot;   // substitution matrix eigen values
    DbMatrix* wU;      // substitution matrix eigen vector 1
    DbMatrix* wV;      // substitution matrix eigen vector 2

    DbMatrix* cQ;

    DbMatrix* sbf;     // structure background frequencies
    DbMatrix* lsbf;     // structure background frequencies
    DbMatrix* stp;     // structure transition probabilities
    DbMatrix* sir;     // state indel rates
    DbMatrix* gep;     // gap extension probabilities
    DbMatrix* mep;     // match extension probabilities
    IntMatrix* codon;      // position in a codon state
    IntMatrix* drawPt;     // draw pattern for kav
    IntMatrix* drawCl;     // draw color for kav
    IntMatrix* drawOf;     // draw offset for kav
    std::string* stNames; // state names
    bool* stShow;      // state show/not

    DbMatrix* cPl;     // substitution probabilities left
    DbMatrix* cPr;     // substitution probabilities right
    DbMatrix* logcPl;  // substitution probabilities left
    DbMatrix* logcPr;  // substitution probabilities right

    DbMatrix* trp;     // state transition probabilities
    DbMatrix* pba;     // probabilities to begin alignment
    DbMatrix* pea;     // probabilities to end alignment
    IntMatrix* tiX;    // non-zero structure transition index
    IntMatrix* tiY;    // non-zero structure transition index

    int end;         // index for file reader

    AncestralNode *node;     // current tree node
    int i,j,k,l,m;
public:
    HMModel();
    ~HMModel();

    void alignmentModel(AncestralNode* tn);
    void readModel(const char* filename);
    void proteinModel();
    void codonModel();
    void dnaModel(float* pi,bool isRna);
    void buildModel();
    void pairwiseModel(IntMatrix* scores,float dist);

    // general functions
    int getNStates()
    {
        return sn;
    }
    int getASize()
    {
        return as;
    }
    int getFullASize()
    {
        return fas;
    }
    std::string getAlphabet()
    {
        return alphabet;
    }
    std::string getFullAlphabet()
    {
        return fullAlphabet;
    }

    int getCodon(int i)
    {
        return codon->g(i);
    }

    // function for sequence alignment
    double charBgFreq(int k,int j)
    {
        return cPi->g(k,j);
    }
    double nullBgFreq(int j)
    {
        return nPi->g(j);
    }
    double charSbProbL(int k,int i,int j)
    {
        return cPl->g(k,i,j);
    }
    double charSbProbR(int k,int i,int j)
    {
        return cPr->g(k,i,j);
    }

    double logCharBgFreq(int k,int j)
    {
        return logcPi->g(k,j);
    }
    double logNullBgFreq(int j)
    {
        return lognPi->g(j);
    }
    double logCharSbProbL(int k,int i,int j)
    {
        return logcPl->g(k,i,j);
    }
    double logCharSbProbR(int k,int i,int j)
    {
        return logcPr->g(k,i,j);
    }

    double structBgFreq(int k)
    {
        return lsbf->g(k);
    }

    double probWX(int k)
    {
        return pba->g(k,0);
    }
    double probWY(int k)
    {
        return pba->g(k,1);
    }
    double probWM(int k)
    {
        return pba->g(k,2);
    }
    double probXW(int k)
    {
        return pea->g(k,0);
    }
    double probYW(int k)
    {
        return pea->g(k,1);
    }
    double probMW(int k)
    {
        return pea->g(k,2);
    }

    double probXX(int k,int l)
    {
        return trp->g(0,k,0,l);
    }
    double probXY(int k,int l)
    {
        return trp->g(0,k,1,l);
    }
    double probXM(int k,int l)
    {
        return trp->g(0,k,2,l);
    }

    double probYX(int k,int l)
    {
        return trp->g(1,k,0,l);
    }
    double probYY(int k,int l)
    {
        return trp->g(1,k,1,l);
    }
    double probYM(int k,int l)
    {
        return trp->g(1,k,2,l);
    }

    double probMX(int k,int l)
    {
        return trp->g(2,k,0,l);
    }
    double probMY(int k,int l)
    {
        return trp->g(2,k,1,l);
    }
    double probMM(int k,int l)
    {
        return trp->g(2,k,2,l);
    }

    int transIndX(int k,int i)
    {
        return tiX->g(k,i);
    }
    int transIndY(int k,int i)
    {
        return tiY->g(i,k);
    }

    std::string getDrawPt(int i);
    std::string getDrawCl(int i);
    int getDrawOf(int i)
    {
        return drawOf->g(i);
    }
    std::string getStName(int i)
    {
        return stNames[i];
    }
    bool getStShow(int i)
    {
        return stShow[i];
    }

    // functions for reading model file
    std::string nextNotComment(std::ifstream* in);
    std::string getString(std::string row,std::string chars);
    int nextInt(std::string row);
    double nextDouble(std::string row);

};

#endif
