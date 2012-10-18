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
#ifndef TERMINALNODE_H
#define TERMINALNODE_H

#include <treenode.h>
#include <terminalsequence.h>

class TerminalNode : public TreeNode
{
    TerminalSequence* seq;
public:
    ~TerminalNode();
    TerminalNode(std::string s,float l);

    TerminalSequence* getSequence()
    {
        return seq;
    }

    int getTerminalNodeNumber();
    int getInternalNodeNumber();

    void getNames(std::vector<std::string>* nms);
    void getTerminalNames(std::vector<std::string>* nms);
    void getInternalNames(std::vector<std::string>* nms);

    void setCharString(std::vector<std::string>* sns,std::vector<std::string>* sqs,int* count);
    void setCharString(std::vector<std::string>* sns,std::vector<std::string>* sqs);
    void getCharStrings(std::vector<std::string>* sqs);

    void alignSequences();

    void getCleanNewick(std::string* tree);
    void getLowestAlignmentPostProbAt(double*,int);
    void outputXml(std::ofstream* out,std::map<std::string,std::string> *anc_seqs,bool triple);

    void writeNewick(std::string* tree,int* sInd);
    void getNewick(std::string* tree);
    void getLabelledNewickBrl(std::string* tree)
    {
        this->getNewickBrl(tree);
    };
    void getNewickBrl(std::string* tree);
    void getNexusTree(std::string* tree, int *count);
    void writeAncCharacters(int *parSite,int iteration);

    void getMLAncestralSeqs(std::vector<std::string>* sqs,std::vector<std::string>* nms);
    void setSiteLength(int ) {}
    void setSiteIndex(int ,int ) {}

    void getAllCharactersAt(std::vector<std::string>* col,int i,bool parentIns,bool parentPermIns) { this->getCharactersAt(col,i,parentPermIns); }
    void getAncCharactersAt(std::vector<std::string>* col,int i,bool parentIns);
    void getCharactersAt(std::vector<std::string>* col,int i,bool parentPermIns=false);

};

#endif
