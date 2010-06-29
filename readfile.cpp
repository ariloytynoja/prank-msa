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

#include <string>
#include <map>
#include <sstream>
#include <set>
#include <iostream>
#include <fstream>
#include "readfile.h"
#include <cstdlib>
#include <cctype>
#include <algorithm>

using namespace std;

ReadFile::ReadFile()
{
}

ReadFile::~ReadFile()
{
    names.clear();
    seqs.clear();
}

bool ReadFile::dnaSeqs() {

    string nucs = "ACGTUN";

    int match=0;
    int total1=0;
    int total2=0;
    int pos;
    vector<string>::iterator si = seqs.begin();
    for (;si!=seqs.end();si++){
        total1 += (*si).length();
        for (unsigned int i=0;i<(*si).length();i++){
            pos= nucs.find((*si).at(i));
            if (pos>=0 && pos<=(int)nucs.length())
                match++;
            if((*si).at(i) != '-')
                total2++;
        }
    }
    return (float)match/(float)total2 > 0.95;
}

void ReadFile::countDnaFreqs(float* freqs)
{
    isRna=false;
    string nucs = "ACGTU";

    int pos;
    int nt = 0;
    int nu = 0;
    vector<string>::iterator si = seqs.begin();
    for (;si!=seqs.end();si++){
        for (unsigned int i=0;i<(*si).length();i++){
            pos= nucs.find((*si).at(i));
            if (pos<0 || pos>4)
                continue;
            if (pos==3)
                nt++;
            if (pos==4)
                nu++;
            if (pos==4)
                pos--;
            freqs[pos]++;
        }
    }
    if (nu>nt)
        isRna=true;
}

int ReadFile::readFile(const char* inputfile)
{

    ifstream input(inputfile, ios::in);


    if (!input) {
        cout<<"Failed to open sequence file "<<inputfile<<". Exiting.\n\n";
        exit(-1);
    }


    char c = input.get();
    while(c==' ' || c=='\n')
    {
        c = input.get();
    }

    if(c=='>')
    {
        input.unget();
        this->readFasta(input);
    }
    else if(c=='#')
    {
        input.unget();
        this->readNexus(input);
    }
    else if(isdigit( c ))
    {
        input.unget();
        this->readPhylip(input);
    }

    else
    {
        cout<<"Input file format unrecognized. Only FASTA format supported. Exiting.\n\n";
        exit(-1);
    }

    set<string> nameset;
    for(int i=0;i<(int)seqs.size();i++)
    {
        string name = names.at(i);
        while(nameset.find(name) != nameset.end())
        {
            cout<<"Sequence name "<<name<<" is defined more than once! Adding suffix '.1'.\n";
            names.at(i) += ".1";
            name += ".1";
        }
        nameset.insert(name);
    }


    if (names.size() == seqs.size())
        return names.size();
    return 0;

}


void ReadFile::readFasta(istream & input)
{

    string temp, name, sequence = "";  // Initialization

    while(!input.eof())
    {
        getline(input, temp, '\n');  // Copy current line in temporary string


        // If first character is >
        if(temp[0] == '>')
        {
            temp = this->remove_last_whitespaces(temp);

            // If a name and a sequence were found
            if((name != "") && (sequence != ""))
            {
                names.push_back(name);
                sequence = this->remove_whitespaces(sequence);
                seqs.push_back(sequence);
                name = "";
                sequence = "";
            }
            name = temp;
            name.erase(name.begin());  // Character > deletion
        }
        else
        {
            sequence += temp;  // Sequence isolation
        }
    }

    // Addition of the last sequence in file
    if((name != "") && (sequence != ""))
    {
        names.push_back(name);
        sequence = this->remove_whitespaces(sequence);
        seqs.push_back(sequence);
    }

}

void ReadFile::readNexus(std::istream & input)
{
    string temp, name, sequence = "";  // Initialization
    getline(input, temp, '\n');  // Copy current line in temporary string

    transform( temp.begin(), temp.end(), temp.begin(), (int(*)(int))toupper );
    temp = this->remove_whitespaces(temp);

    if(temp != "#NEXUS")
    {
        cout<<"Input file starts with '#' but not with '#NEXUS'. Reading the file failed. Exiting.\n";
        exit(-1);
    }

    int ntax = -1;
    int length = -1;

    while(!input.eof())
    {
        getline(input, temp, '\n');  // Copy current line in temporary string

        string::size_type loc = temp.find("ntax");
        if(loc!=string::npos)
        {
            string str = temp.substr(loc+5,temp.find_first_of(" ;",loc+5)-(loc+5));
            stringstream ss(str);
            ss>>ntax;
        }

        loc = temp.find("nchar");
        if(loc!=string::npos)
        {
            string str = temp.substr(loc+6,temp.find_first_of(" ;",loc+6)-(loc+6));
            stringstream ss(str);
            ss>>length;
        }

        if(temp.find("matrix")!=string::npos)
            break;
    }

    if(ntax<1 || length<1)
    {
        cout<<"Failed to read the dimensions of the Nexus alignment. Exiting.\n";
        exit(-1);
    }

    stringstream rows;
    map<string,string> data;

    while(!input.eof())
    {
        getline(input, temp, '\n');  // Copy current line in temporary string
        if(temp.find("end;")!=string::npos)
            break;

        rows.clear();
        rows.str(temp);

        string name,seq = "";
        rows>>name>>seq;

        name = this->remove_last_whitespaces(name);
        seq = this->remove_whitespaces(seq);

        if(name.length()>0 && seq.length()>0)
        {
            map<string,string>::iterator mi = data.find(name);
            if(mi==data.end())
            {
                data.insert(make_pair(name,seq));
                names.push_back(name);
//                cout<<"new "<<name<<endl;
            }
            else
            {
                mi->second += seq;
//                cout<<"old "<<name<<endl;
            }
        }
    }

    for(int i=0;i<ntax;i++)
    {
        map<string,string>::iterator mi = data.find(names.at(i));
        seqs.push_back(mi->second);

        if(mi->second.length()!=length)
        {
            cout<<"Reading may have failed: sequences are not equally long!\n";
        }
    }

}

void ReadFile::readPhylip(std::istream & input)
{

    int nseq = -1;
    int length = -1;

    string temp, name, sequence = "";  // Initialization
    getline(input, temp, '\n');  // Copy current line in temporary string

    stringstream nums(temp);
    nums>>nseq>>length;

    if(nseq<1 || length<1)
    {
        cout<<"Input file starts with a digit but not with two positive digits. Reading the file failed. Exiting.\n";
        exit(-1);
    }

    getline(input, temp, '\n');

    stringstream rows(temp);
    rows>>name>>sequence;

    name = this->remove_last_whitespaces(name);
    sequence = this->remove_whitespaces(sequence);

    if((int) name.length()>0 && (int) name.length()<=10 && ( ( (int) sequence.length()>=50 && (int) sequence.length()<=60 ) || (int) sequence.length()==length )  )
        readInterleaved(temp,input,nseq,length);
    else
        readSequential(temp,input,nseq,length);
}

void ReadFile::readInterleaved(string temp,istream & input,int nseq, int length)
{
    stringstream rows;

    for(int i=0;i<nseq;i++)
    {
        rows.clear();
        rows.str(temp);
        string name,seq = "";
        rows>>name>>seq;

        name = this->remove_last_whitespaces(name);
        seq = this->remove_whitespaces(seq);

        names.push_back(name);
        seqs.push_back(seq);

        getline(input, temp, '\n');
    }

    int i=0;
    do
    {
        temp = this->remove_whitespaces(temp);
        seqs.at(i++) += temp;
        if(i==nseq)
            i=0;
    }
    while(getline(input, temp, '\n'));

    for(i=0;i<nseq;i++)
    {
        if((int) seqs.at(i).length()!=length)
        {
            cout<<"Reading may have failed: interleaved sequences are not equally long!\n";
        }
    }

}

void ReadFile::readSequential(string temp,istream & input,int nseq, int length)
{
    for(int i=0;i<nseq;i++)
    {
        string name = temp;
        name = this->remove_last_whitespaces(name);
        names.push_back(name);

        string seq = "";
        while((int) seq.length()<length)
        {
            getline(input, temp, '\n');
            temp = this->remove_whitespaces(temp);
            seq += temp;
        }
        seqs.push_back(seq);

        getline(input, temp, '\n');
    }

    for(int i=0;i<nseq;i++)
    {
        if((int) seqs.at(i).length()!=length)
        {
            cout<<"Reading may have failed: sequential sequences are not equally long!\n";
            cout<<seqs.at(i)<<" "<<seqs.at(i).length()<<endl;
        }
    }
}

string ReadFile::remove_last_whitespaces(const string & s)
{
  // Copy sequence
  string st (s);

  while(st.size() > 0 && this->is_whitespace_character(st[st.size() - 1]))
  {
    st.erase(st.end() - 1);
  }

  // Send result
  return st;
}

string ReadFile::remove_whitespaces(const string & s)
{
  string st="";

  for (unsigned int i = 0; i < s.size(); i++)
  {
    if(!this->is_whitespace_character(s[i]))
    {
      st+=s[i];
    }
  }
  return st;
}

bool ReadFile::is_whitespace_character(char c)
{
    return (c == ' ')
        || (c == '\t')
        || (c == '\n')
        || (c == '\r')
        || (c == '\f');
}
