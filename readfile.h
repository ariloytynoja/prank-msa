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
#ifndef READFILE_H
#define READFILE_H

#include <string>
#include <iostream>
#include <vector>

using namespace std;

class ReadFile{
    vector<string> names;
    vector<string> seqs;
    
public:
    ReadFile();
    ~ReadFile();
 
    int readFile(const char* filename);
    void readFasta(istream & input);
    void readNexus(istream & input);
    void readPhylip(istream & input);

    void readInterleaved(string temp,istream & input,int nseq, int length);
    void readSequential(string temp,istream & input,int nseq, int length);

    string remove_last_whitespaces(const string & s);
    string remove_whitespaces(const string & s);
    bool is_whitespace_character(char c);

    vector<string> getNames(){ return names; }
    vector<string> getSeqs(){ return seqs; }

    bool dnaSeqs();
    void countDnaFreqs(float* freqs);
    bool isRna;
};

#endif
