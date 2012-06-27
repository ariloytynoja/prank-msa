/***************************************************************************
 *   Copyright (C) 2005 by Ari Loytynoja                                   *
 *   ari@ebi.ac.uk                                                         *
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
#include <iomanip>
#include <iostream>
#include <fstream>
#include "writefile.h"

extern bool DOPOST;
extern bool CODON;

using namespace std;

WriteFile::WriteFile()
{
    error = "";
    chars_by_line = 60;
}


WriteFile::~WriteFile()
{
}


void WriteFile::writeSeqs(const char* outputfile, std::vector<std::string>* names, std::vector<std::string> *seqs, int outform, bool isDna,AncestralNode *root, bool translate)
{
    if (outform == 17)
        this->writeNexus(outputfile,names,seqs,isDna,root,translate);
    else
        this->writeSeqs(outputfile,names,seqs,outform);
}

void WriteFile::writeSeqs(const char* outputfile, std::vector<std::string>* names, std::vector<std::string> *seqs, int outform )
{
    if (outform == 8)
        this->writeFasta(outputfile,names,seqs);
    else if (outform == 12)
        this->writeInterleaved(outputfile,names,seqs);
    else if (outform == 11)
        this->writeSequential(outputfile,names,seqs,true);
    else if (outform == 17)
        this->writeSimpleNexus(outputfile,names,seqs);
    else if (outform == 18)
        this->writeSequential(outputfile,names,seqs,false);
    else if (outform == 19)
        this->writeLongSequential(outputfile,names,seqs);
    else
        this->writeFasta(outputfile,names,seqs);

}


void WriteFile::writeFasta(const char* outputfile, vector<string> * names, vector<string> * seqs)
{
    ofstream output( outputfile );

    // Checking the existence of specified file, and possibility to open it in write mode
    if (! output)
    {
        cout<<"Failed to open output file "<<outputfile<<". Exiting.\n\n";
        exit(-1);
    }

    string seq, temp = "";  // Initialization

    vector<string>::iterator si = seqs->begin();
    vector<string>::iterator ni = names->begin();

    // Main loop : for all sequences in vector container
    for (; ni != names->end(); ni++,si++)
    {
        output << ">" << *ni;
        output << endl;

        // Sequence cutting to specified characters number per line
        seq = *si;
//        cout<<*ni<<"\n"<<seq<<"\n\n";
        while (seq != "")
        {
            if ((int)seq.size() > chars_by_line)
            {
                temp = string(seq.begin(), seq.begin() + chars_by_line);
                output << temp  << endl;
                seq.erase(seq.begin(), seq.begin() + chars_by_line);
            }
            else
            {
                output << seq << endl;
                seq = "";
            }
        }
    }

    output.close();
}

void WriteFile::writeInterleaved(const char* outputfile,std::vector<std::string> *names,std::vector<std::string> *seqs)
{
    ofstream output( outputfile );

    // Checking the existence of specified file, and possibility to open it in write mode
    if (! output)
    {
        cout<<"Failed to open output file "<<outputfile<<". Exiting.\n\n";
        exit(-1);
    }

    string seq, temp = "";  // Initialization

    int length = seqs->begin()->length();


    output<<seqs->size()<<" "<<length<<endl;


    for (int offset = 0; offset<length; offset+=chars_by_line)
    {
        vector<string>::iterator si = seqs->begin();
        vector<string>::iterator ni = names->begin();

        for (; ni!=names->end(); ni++,si++)
        {
            string tmp = ni->substr(0,10)+"          ";
            if (offset > 0)
            {
                tmp = "           ";
            }
            output << tmp.substr(0,10)<<" ";

            output<<si->substr(offset,chars_by_line)<<endl;
        }
    }
}

void WriteFile::writeSequential(const char* outputfile,std::vector<std::string> *names,std::vector<std::string> *seqs,bool truncate)
{
    ofstream output( outputfile );

    // Checking the existence of specified file, and possibility to open it in write mode
    if (! output)
    {
        cout<<"Failed to open output file "<<outputfile<<". Exiting.\n\n";
        exit(-1);
    }

    string seq, temp = "";  // Initialization

    vector<string>::iterator si = seqs->begin();
    vector<string>::iterator ni = names->begin();

    output<<seqs->size()<<" "<<si->length()<<endl;

    // Main loop : for all sequences in vector container
    for (; ni != names->end(); ni++,si++)
    {
        if (truncate)
            output << (*ni+"          ").substr(0,10)<<" "<< endl;
        else
            output << *ni<< endl;

        // Sequence cutting to specified characters number per line
        seq = *si;
        while (seq != "")
        {
            if ((int)seq.size() > chars_by_line)
            {
                temp = string(seq.begin(), seq.begin() + chars_by_line);
                output << temp  << endl;
                seq.erase(seq.begin(), seq.begin() + chars_by_line);
            }
            else
            {
                output << seq << endl;
                seq = "";
            }
        }
    }

    output.close();
}

void WriteFile::writeLongSequential(const char* outputfile,std::vector<std::string> *names,std::vector<std::string> *seqs)
{
    ofstream output( outputfile );

    // Checking the existence of specified file, and possibility to open it in write mode
    if (! output)
    {
        cout<<"Failed to open output file "<<outputfile<<". Exiting.\n\n";
        exit(-1);
    }

    string seq, temp = "";  // Initialization

    vector<string>::iterator si = seqs->begin();
    vector<string>::iterator ni = names->begin();

    output<<seqs->size()<<" "<<si->length()<<endl;

    // Main loop : for all sequences in vector container
    for (; ni != names->end(); ni++,si++)
    {
        output << *ni<< endl<<*si<<endl;
    }

    output.close();
}

void WriteFile::writeNexus(const char* outputfile,std::vector<std::string> *names,std::vector<std::string> *seqs, bool isDna, AncestralNode *root, bool translate)
{

    ofstream output( outputfile );

    // Checking the existence of specified file, and possibility to open it in write mode
    if (! output)
    {
        cout<<"Failed to open output file "<<outputfile<<". Exiting.\n\n";
        exit(-1);
    }

    string datatype = "protein";
    if (isDna)
        datatype = "dna";

    int length = seqs->begin()->length();

    output<<"#NEXUS\nbegin data;\ndimensions ntax="<<seqs->size()<<" nchar="<<length<<";\nformat datatype="<<datatype<<" interleave=yes gap=-;\nmatrix\n"<<endl;


    for (int offset = 0; offset<length; offset+=chars_by_line)
    {
        output<<endl;

        vector<string>::iterator si = seqs->begin();
        vector<string>::iterator ni = names->begin();

        for (; ni!=names->end(); ni++,si++)
        {
            string tmp = ni->substr(0,20)+"'                    ";
            output << "'"<<tmp.substr(0,21)<<"     ";

            output<<si->substr(offset,chars_by_line)<<endl;
        }
    }
    output<<";\nend;\nbegin trees;\n translate\n";
    for (int i=0; i<(int)names->size(); i++)
    {
        output<<"  "<<i+1<<" '"<<names->at(i).substr(0,20)<<"'";
        if (i<(int)names->size()-1)
            output<<",\n";
        else
            output<<"\n";
    }

    string tree = "";
    int count = 1;
    root->getNexusTree(&tree,&count);

    output<<" ;\n tree * PRANK =\n"<<tree<<";\nend;\n";

    if (DOPOST)
    {
        output << "begin assumptions;\n wtset PrankMinimum (VECTOR) = \n ";

        int l = root->getSequence()->length();

        root->setSiteLength(l);
        for (int i=0; i<l; i++)
        {
            root->setSiteIndex(i,i);
        }

        for (int offset = 0; offset<l; offset++)
        {
            double p = 1.0;
            root->getLowestAlignmentPostProbAt(&p,offset);
            output<<fixed<<setprecision(2);
            output<<" "<<p;
            if (translate || CODON)
                output<<" "<<p<<" "<<p;
        }

        output<<";\nend;\n";
    }

    output.close();
}


void WriteFile::writeSimpleNexus(const char* outputfile, vector<string> * names, vector<string> * seqs)
{
    ofstream output( outputfile );

    // Checking the existence of specified file, and possibility to open it in write mode
    if (! output)
    {
        cout<<"Failed to open output file "<<outputfile<<". Exiting.\n\n";
        exit(-1);
    }

    string datatype = "protein";

    if (this->dnaSeqs(seqs))
        datatype = "dna";

    int length = seqs->begin()->length();

    output<<"#NEXUS\nbegin data;\ndimensions ntax="<<seqs->size()<<" nchar="<<length<<";\nformat datatype="<<datatype<<" interleave=yes gap=-;\nmatrix\n"<<endl;


    for (int offset = 0; offset<length; offset+=chars_by_line)
    {
        output<<endl;

        vector<string>::iterator si = seqs->begin();
        vector<string>::iterator ni = names->begin();

        for (; ni!=names->end(); ni++,si++)
        {
            string tmp = ni->substr(0,20)+"'                    ";
            output << "'"<<tmp.substr(0,21)<<"     ";

            output<<si->substr(offset,chars_by_line)<<endl;
        }
    }
    output<<";\nend;";

    output.close();
}

bool WriteFile::dnaSeqs(vector<string> * seqs)
{

    string nucs = "ACGTUN";

    int match=0;
    int total1=0;
    int total2=0;
    int pos;
    vector<string>::iterator si = seqs->begin();
    for (; si!=seqs->end(); si++)
    {
        total1 += (*si).length();
        for (unsigned int i=0; i<(*si).length(); i++)
        {
            pos= nucs.find((*si).at(i));
            if (pos>=0 && pos<=(int)nucs.length())
                match++;
//            if((*si).at(i) != '-')
            if ((*si).at(i) != '-' && (*si).at(i) != '?')
                total2++;
        }
    }
    return (float)match/(float)total2 > 0.95;
}
