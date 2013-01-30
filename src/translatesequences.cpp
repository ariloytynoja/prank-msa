/***************************************************************************
 *   Copyright (C) 2008 by Ari Loytynoja   *
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
#include "translatesequences.h"
#include "config.h"

using namespace std;

TranslateSequences::TranslateSequences()
{

    string codon[66] = {"TTT", "TTC", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG",
                        "ATT", "ATC", "ATA", "ATG", "GTT", "GTC", "GTA", "GTG",
                        "TCT", "TCC", "TCA", "TCG", "CCT", "CCC", "CCA", "CCG",
                        "ACT", "ACC", "ACA", "ACG", "GCT", "GCC", "GCA", "GCG",
                        "TAT", "TAC", "TAA", "TAG", "CAT", "CAC", "CAA", "CAG",
                        "AAT", "AAC", "AAA", "AAG", "GAT", "GAC", "GAA", "GAG",
                        "TGT", "TGC", "TGA", "TGG", "CGT", "CGC", "CGA", "CGG",
                        "AGT", "AGC", "AGA", "AGG", "GGT", "GGC", "GGA", "GGG",
                        "NNN", "---"
                       };
    string unaa[66] = {"F", "F", "L", "L", "L", "L", "L", "L", "I", "I", "I", "M", "V", "V", "V", "V",
                       "S", "S", "S", "S", "P", "P", "P", "P", "T", "T", "T", "T", "A", "A", "A", "A",
                       "Y", "Y", "X", "X", "H", "H", "Q", "Q", "N", "N", "K", "K", "D", "D", "E", "E",
                       "C", "C", "X", "W", "R", "R", "R", "R", "S", "S", "R", "R", "G", "G", "G", "G",
                       "X", "-"
                      };
    string mtaa[66] = {"F", "F", "L", "L", "L", "L", "L", "L", "I", "I", "M", "M", "V", "V", "V", "V",
                       "S", "S", "S", "S", "P", "P", "P", "P", "T", "T", "T", "T", "A", "A", "A", "A",
                       "Y", "Y", "X", "X", "H", "H", "Q", "Q", "N", "N", "K", "K", "D", "D", "E", "E",
                       "C", "C", "W", "W", "R", "R", "R", "R", "S", "S", "X", "X", "G", "G", "G", "G",
                       "X", "-"
                      };

    if (MTTABLE)
    {
        for (int i=0; i<66; i++)
        {
            codonToAa.insert(make_pair(codon[i],mtaa[i]));
            aaToCodon.insert(make_pair(mtaa[i],codon[i]));
        }
    }
    else
    {
        for (int i=0; i<66; i++)
        {
            codonToAa.insert(make_pair(codon[i],unaa[i]));
            aaToCodon.insert(make_pair(unaa[i],codon[i]));
        }
    }
}


TranslateSequences::~TranslateSequences()
{
}


bool TranslateSequences::translateProtein(vector<string> *names,vector<string> *sequences,map<string,string> *dnaSequences)
{

    vector<string>::iterator nit = names->begin();
    vector<string>::iterator sit = sequences->begin();
//    dnaSeqs.clear();
    dnaSequences->clear();
    bool replaced = false;
    string full_alphabet = "ACGTN";

    bool inFrame = true;

    for (; sit!=sequences->end(); sit++)
    {
        for (unsigned int j=0; j<sit->length(); j+=3)
        {
            string codon = sit->substr(j,3);
            if (codonToAa.find(codon)==codonToAa.end())
            {
                inFrame = false;
                if(UPDATE || PREALIGNED)
                    cout<<"Input alignment not in frame. Gaps removed and realignment needed.\n";
                UPDATE = false;
                PREALIGNED = false;
                break;
            }
        }
        if(not inFrame)
            break;
    }

    sit = sequences->begin();
    for (; sit!=sequences->end(); sit++)
    {
        if(not inFrame)
        {
            string seq = *sit;
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
            *sit = seq;
        }

        for (unsigned int j=0; j<sit->length(); j+=3)
        {
            string codon = sit->substr(j,3);
            if (codonToAa.find(codon)==codonToAa.end())
            {
                sit->replace(j,3,"NNN");
                replaced = true;
            }
        }
    }
    if (replaced)
    {
        cout<<"Warning: Unknown codons replaced with 'NNN'."<<endl;
    }

    sit = sequences->begin();

    for (; sit!=sequences->end(); sit++,nit++)
    {
//        dnaSeqs.insert(make_pair(*nit,*sit));
        string seq = *sit;
        for (string::iterator ci = seq.begin();ci != seq.end();ci++)
        {
            if(*ci == '-')
            {
                seq.erase(ci);
                ci--;
            }
        }
        dnaSequences->insert(make_pair(*nit,seq));
        string tmp ="";
        for (unsigned int j=0; j<sit->length(); j+=3)
        {
            string codon = sit->substr(j,3);
            tmp+=codonToAa.find(codon)->second;
        }
        if (NOISE>1)
            cout<<tmp<<endl;

        *sit=tmp;
    }
    return true;
}

bool TranslateSequences::translateDNA(std::vector<std::string> *names,std::vector<std::string> *protein,std::vector<std::string> *dna,map<string,string> *dnaSequences)
{

    vector<string>::iterator nit = names->begin();
    vector<string>::iterator pit = protein->begin();

    for (; pit!=protein->end(); pit++,nit++)
    {

//        string dnaSeq = dnaSeqs.find(*nit)->second;
        string dnaSeq = dnaSequences->find(*nit)->second;

        string nuc = "";

        for (unsigned int j=0,i=0; j<pit->length(); j++)
        {
            string aa = pit->substr(j,1);
            if (aa=="-")
            {
                nuc+="---";
            }
            else
            {
                string codon = dnaSeq.substr(i,3);
                i+=3;
                if (aa !=codonToAa.find(codon)->second)
                {
                    cout<<"Mismatch in backtranslation: ("<<*nit<<";"<<j<<") "<<codon<<" "<<codonToAa.find(codon)->second<<" != "<<aa<<"."<<endl;
                    return false;
                }
                nuc+=codon;
//         cout<<nuc<<endl;
            }
        }
        dna->push_back(nuc);
    }

    return true;

}
