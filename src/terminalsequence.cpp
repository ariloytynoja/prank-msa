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
#include <string>
#include <iostream>
#include "terminalsequence.h"
#include "config.h"

using namespace std;

TerminalSequence::~TerminalSequence()
{
    delete seqvec;

}

// Define a terminal sequence (matrix) from a non-gapped or gapped character string
// non-gapped = plain alignment
// gapped = re-alignment or posterior probability computation
//
TerminalSequence::TerminalSequence(string* s)
        : Sequence()
{
    terminal = true;

    string alpha = hmm->getAlphabet();
    sAlpha = alpha.length();

    charseq = "";
    int ci;
    string fullAlpha = hmm->getFullAlphabet();
    int sFullAlpha = fullAlpha.length();

//    cout<<"sAlpha "<<sAlpha<<", sFullAlpha "<<sFullAlpha<<endl;

    map<string,int> codons;

    if (PREALIGNED || PARTLYALIGNED || UPDATE)
        gappedseq = *s;

    if (CODON)
    {

        bool gaps_removed = false;

        if (! (PREALIGNED || PARTLYALIGNED || UPDATE) )
            this->removeGaps(s);

        if (s->size()%3!=0)
        {
            this->removeGaps(s);
            gaps_removed = true;

            if (s->size()%3!=0)
            {
//                cout<<s->size()<<" "<<*s<<endl;
                cout<<"codon sequence length is not multiple of three!"<<endl;
                exit(0);
            }
        }

        for (int i=0; i<sFullAlpha; i+=3)
        {
            codons.insert(make_pair(alpha.substr(i,3),i/3));
        }
        sAlpha /= 3;

        codons.insert(make_pair("---",sAlpha));


        bool stop_masked = false;

        string S;
        for (int i=0; i<(int)s->length(); i++)
        {
            S += toupper(s->at(i));
        }

        for (int i=0; i<(int)S.length(); i+=3)
        {
            ci=-1;
            if(codons.find(S.substr(i,3)) != codons.end())
                ci = codons.find(S.substr(i,3))->second;

            if (ci>=0 && ci<sAlpha)
                charseq += S.substr(i,3);
            else if (S.substr(i,3)=="---" || S.substr(i,3)=="...")
                ;
            else
            {
                // An unknown codon -- which includes every STOP codon, since
                // the codon alphabet is 61-state and stops are not in it.
                //
                // Mid-sequence these were already masked as NNN, preserving
                // the column so a caller can restore the real bases.  A
                // TERMINAL one was DROPPED instead, silently shortening the
                // sequence by three: a 3822 nt CDS ending in TAA came back as
                // 3819 ending ...CATTACACA, and the caller had no column left
                // to restore into.  Downstream that reads as "the stop codon
                // was lost" rather than "the model cannot score it".
                //
                // Mask it like any other unknown codon.  The alignment is
                // unaffected -- NNN is what the model already scores for the
                // mid-sequence case -- but the length is now preserved, so the
                // terminal stop stays addressable.  It is still reported at
                // NOISE>0 below, as a masking rather than a removal.
                charseq += "NNN";
                if(!(i+3<(int)S.length() || PREALIGNED))
                    stop_masked = true;
            }
        }

        seqLength = realLength = charseq.size()/3;

        if(NOISE>0 && stop_masked)
            cout<<"Note: terminal stop codon was masked as 'NNN'\n";

        if(NOISE>0 && gaps_removed)
            cout<<"Note: gaps were removed\n";

//        cout<<charseq<<endl;
    }

    else   // Protein or DNA
    {

        if (sAlpha==20)
        {
            for (int i=0; i<(int)s->length(); i++)
            {
                ci=-1;
                if(fullAlpha.find(toupper(s->at(i))) != string::npos)
                    ci = fullAlpha.find(toupper(s->at(i)));

                if (ci>=0 && ci<sFullAlpha)
                    charseq += s->at(i);
                else
                {
                    if (s->at(i)!='-' && s->at(i)!='.')
                        charseq += 'X';
                }
            }
        }

        else
        {
            for (int i=0; i<(int)s->length(); i++)
            {
                ci=-1;
                if(fullAlpha.find(toupper(s->at(i))) != string::npos)
                    ci = fullAlpha.find(toupper(s->at(i)));
                if (ci>=0 && ci<sFullAlpha)
                    charseq += s->at(i);
                else
                {
                    if (s->at(i)!='-' && s->at(i)!='.')
                        charseq += 'N';
                }
            }
        }
        seqLength = realLength = charseq.size();
    }


    // Store the sequence as a probability matrix; note gapped vs non-gapped
    //
    seqvec = new IntMatrix(seqLength,"seqvec");
    seqvec->initialise(-1);

// Note: "NNN" defined as 62nd codon

    if (CODON)
    {
        FOR(i,seqLength)
        {
            ci = codons.find(charseq.substr(i*3,3))->second;
            if (ci>=0 && ci<sAlpha)
            {
                seqvec->s(ci,i);
            }
            else
            {
                seqvec->s(sAlpha,i);
            }
        }
    }

    else if (sAlpha==20)
    {
        FOR(i,seqLength)
        {
            ci = fullAlpha.find(toupper(charseq.at(i)));
            if (ci>=0 && ci<sFullAlpha)
            {
                seqvec->s(ci,i);
            }
            else
            {
                seqvec->s(sAlpha+1,i);
            }
        }
    }

    else
    {
        FOR(i,seqLength)
        {
            ci = fullAlpha.find(toupper(charseq.at(i)));
            if (ci>=0 && ci<sFullAlpha)
            {
                seqvec->s(ci,i);
            }
            else
            {
                seqvec->s(sAlpha+1,i);
            }
        }
    }
}

