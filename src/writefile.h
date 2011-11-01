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
#ifndef WRITEFILE_H
#define WRITEFILE_H

#include <vector>
#include <string>
#include "ancestralnode.h"

class WriteFile
{
    std::string error;
    int chars_by_line;
public:
    WriteFile();
    ~WriteFile();
    void writeSeqs(const char* outputfile, std::vector<std::string>* names, std::vector<std::string> *seqs, int outform, bool isDna, AncestralNode *root, bool translate);
    void writeSeqs(const char* outputfile, std::vector<std::string> *names, std::vector<std::string> *seqs,int outform);
    void writeFasta(const char* outputfile, std::vector<std::string> * names, std::vector<std::string> * seqs);
    void writeInterleaved(const char* outputfile,std::vector<std::string> *names,std::vector<std::string> *seqs);
    void writeSequential(const char* outputfile,std::vector<std::string> *names,std::vector<std::string> *seqs,bool truncate);
    void writeNexus(const char* outputfile,std::vector<std::string> *names,std::vector<std::string> *seqs, bool isDna, AncestralNode *root, bool translate);
    void writeSimpleNexus(const char* outputfile, std::vector<std::string> * names, std::vector<std::string> * seqs);

    bool dnaSeqs(std::vector<std::string> * seqs);

    bool hasError()
    {
        return error != "";
    }
    std::string getError()
    {
        return error;
    }

};

#endif
