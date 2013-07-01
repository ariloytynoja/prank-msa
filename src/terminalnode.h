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

    void concatenateTerminalNames(std::string *s)
    {
        s->append(nodeName);
    }

    void getNames(std::vector<std::string>* nms);
    void getTerminalNames(std::vector<std::string>* nms);
    void getInternalNames(std::vector<std::string>* nms);

    void setCharString(std::vector<std::string>* sns,std::vector<std::string>* sqs,int* count);
    void setCharString(std::vector<std::string>* sns,std::vector<std::string>* sqs);
    void getCharStrings(std::vector<std::string>* sqs);

    void getAllSubtrees(std::map<std::string,float> *subtrees) {}
    void getAllSubtreesWithNodename(std::map<std::string,std::string> *subtrees) {}
    void getSubtreeBelow(std::string *subtree) { *subtree = nodeName; }
    void markRealignSubtrees(std::map<std::string,float> *subtrees) {}

    bool anyChildNodeRealigned() { return false; }

    void alignSequences();

    void getCleanNewick(std::string* tree);
    void getLowestAlignmentPostProbAt(double*,int);
    void outputXml(std::ofstream* out,std::map<std::string,std::string> *anc_seqs,bool triple);

    void writeNewick(std::string* tree,int* sInd);
    void writeLabelledNewick(std::string* tree,int* sInd) {writeNewick(tree,sInd);}
    void getNewick(std::string* tree);
    void getLabelledNewickBrl(std::string* tree)
    {
        this->getNewickBrl(tree);
    }
    void getLabelledNewick(std::string* tree)
    {
        this->getNewick(tree);
    }

    void getNHXBrl(std::string* tree,int *nodeNumber)
    {
        this->getNewickBrl(tree);
    }

    void getNewickBrl(std::string* tree);
    void getNexusTree(std::string* tree, int *count);

    void getMLAncestralSeqs(std::vector<std::string>* sqs,std::vector<std::string>* nms);
    void setSiteLength(int ) {}
    void setSiteIndex(int ,int ) {}

    void getAllCharactersAt(std::vector<std::string>* col,int i,bool parentIns,bool parentPermIns) { this->getCharactersAt(col,i,parentPermIns); }
    void getAncCharactersAt(std::vector<std::string>* col,int i,bool parentIns);
    void getCharactersAt(std::vector<std::string>* col,int i,bool parentPermIns=false);
    void getIndelEvents(std::vector<indelEvent> *indels){}
    void getSubstEvents(std::vector<substEvent> *substs){}

    void setAncSequenceStrings(std::vector<std::string>*){}
    void getAncSequenceStrings(std::vector<std::string>*){}

    void setAncSequenceGaps(std::vector<std::string>*){}
    void setAncSequenceStrings(std::map<std::string,std::string>*){}

    void setAlignedSequenceStrings(std::vector<std::string>* aseqs)
    {
        alignedseqstr = aseqs->at(0);
        aseqs->erase(aseqs->begin());
    }

    void getAlignedSequenceStrings(std::vector<std::string>* aseqs)
    {
        aseqs->push_back(alignedseqstr);
    }

    void fixTerminalNodenames()
    {
        if(nodeName.find('_') != std::string::npos)
        {
            nodeName = nodeName.substr(nodeName.find('_')+1);
        }
    }

    bool updateInsertionSite(int i);

};

#endif
