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

#ifndef ANCESTRALNODE_H
#define ANCESTRALNODE_H

#include <treenode.h>
#include <ancestralsequence.h>

class AncestralNode : public TreeNode
{
    AncestralSequence* seq;
public:
    AncestralNode(std::string s);

    ~AncestralNode();

    std::string left_nhx_tag;
    std::string right_nhx_tag;

    AncestralSequence* getSequence()
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

    void getAllSubtrees(std::map<std::string,float> *subtrees);
//    void getAllSubtrees(std::set<std::string> *subtrees);
    void getSubtreeBelow(std::string *subtree);
    void markRealignSubtrees(std::map<std::string,float> *subtrees);
//    void markRealignSubtrees(std::set<std::string> *subtrees);

    bool anyChildNodeRealigned()
    {
        if(getLChild()->anyChildNodeRealigned())
            return true;
        if(getRChild()->anyChildNodeRealigned())
            return true;
        return realignNode;
    }

    void getThisAlignmentPostProbAt(double* p,int i);
    void getLowestAlignmentPostProbAt(double* p,int i);

    void alignSequences( );
    void alignThisNode( );
    void readAlignment();
    void readThisNode();
    void printDebugNodes();

    void partlyAlignSequences();
    void updateAlignedSequences();

    void getCleanNewick(std::string* tree);
    void outputXml(std::ofstream* out,std::map<std::string,std::string> *anc_seqs,bool triple);

    void writeNewick(std::string* tree,int* sInd);
    void getNewick(std::string* tree);
    void getLabelledNewickBrl(std::string* tree);
    void getNewickBrl(std::string* tree);
    void getNexusTree(std::string* tree, int *count);
    void writeAncCharacters(int *parSite,int iteration);

    void getMLAncestralSeqs(std::vector<std::string>* sqs,std::vector<std::string>* nms);
    void setSiteLength(int l);
    void setSiteIndex(int site,int index);

    void getAllCharactersAt(std::vector<std::string>* col,int i,bool parentIns,bool parentPermIns);
    void getAncCharactersAt(std::vector<std::string>* col,int i,bool parentIns,bool parentPermIns);
    void getCharactersAt(std::vector<std::string>* col,int i,bool parentPermIns=false);

    void setPermanentInsertion(int i);
    void printChildAlignment(TreeNode *node,std::string filename);
};

#endif
