/***************************************************************************
 *   Copyright (C) 2005 by Ari Loytynoja                                   *
 *   ari@sink                                                              *
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
#ifndef TREENODE_H
#define TREENODE_H

/**
 * A node in a hierarchical tree.
 */

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include "sequence.h"
#include "site.h"
#include "flmatrix.h"

extern float minBrL;

class TreeNode {
protected:
    std::string nodeName; // name
    float branchLength;   // length of the branch below
    bool terminal;        // is/not terminal node
    bool root;            // is/not root node

    std::string ln,rn,n3; // names of left, right and third (unrooted) branch
    float ld,rd,d3;       // lengths of those branches
    TreeNode *lChild;     // left child
    TreeNode *rChild;
    bool lInternal;       // is/not left internal
    bool rInternal;

    static bool rooted;   // tree is/not rooted

    float tot;
    float left;
    float right;

    std::string groupName;

    std::string charString;  // unaligned sequence
    Sequence* seq;

    static int totalNodes;
    static int alignedNodes;

    int siteLength;
    int* siteIndex;         // index for site output

public:
    virtual ~TreeNode();

    virtual void alignSequences(int ) {}
    virtual void readAlignment() {}
    virtual void partlyAlignSequences() {}

    virtual int getTerminalNodeNumber() = 0;
    virtual int getInternalNodeNumber() = 0;

    virtual std::string getGroupName() { return groupName; }

    void setTotalNodes() { totalNodes = this->getTerminalNodeNumber(); alignedNodes = 1; }

    void setBranchLength(float l) { if(l<minBrL){ branchLength = minBrL;} else { branchLength = l;} } //set the length of the  branch below
    float getBranchLength() { return branchLength; }


    std::string getNodeName(){ return nodeName; }    // get node name
    void setNodeName(std::string s){ nodeName = s; } // set node name

    bool isRoot(){ return root; }
    void setRoot(){ root = true; }

    bool isRooted(){ return rooted; }
    bool isLInternal(){ return lInternal; }
    bool isRInternal(){ return rInternal; }

    bool isTerminal(){ return terminal; }

    float getLeftBrL(){ return ld; }
    float getRightBrL(){ return rd; }

    int gappedLength();

    virtual void getNames(std::vector<std::string>* nms) = 0;
    virtual void getTerminalNames(std::vector<std::string>* nms) = 0;
    virtual void getInternalNames(std::vector<std::string>* nms) = 0;

    virtual void setCharString(std::vector<std::string>* sns,std::vector<std::string>* sqs,int* count) = 0;
    virtual void setCharString(std::vector<std::string>* sns,std::vector<std::string>* sqs) = 0;
    virtual void getCharStrings(std::vector<std::string>* sqs) = 0;

    virtual void setAnnotation(std::map<std::string,FlMatrix*>* annotation) = 0;

    virtual Sequence* getSequence() = 0;

    void setLChild(TreeNode* tn){ lChild = tn; } // set child nodes
    void setRChild(TreeNode* tn){ rChild = tn; }

    std::string getLName(){ return ln; }      // get child node names
    std::string getRName(){ return rn; }

    virtual void getCharactersAt(std::vector<std::string>* ,int ){}
    virtual void getAncCharactersAt(std::vector<std::string>* ,int ,bool ){}

    virtual void setSiteLength(int ){}
    virtual void setSiteIndex(int ,int ){};

    virtual void getLowestAlignmentPostProbAt(double*,int) = 0;

    virtual void outputXml(std::ofstream* out,bool triple) = 0;

    virtual void writeNewick(std::string* ,int* ){}
    virtual void getNewick(std::string* tree) = 0;
    virtual void getLabelledNewickBrl(std::string* tree) = 0;
    virtual void getNewickBrl(std::string* tree) = 0;
    virtual void getNexusTree(std::string* tree, int *count) = 0;
    virtual void writeAncCharacters(int *parSite,int iteration) = 0;

    void getCleanNewick(std::string* tree);
    virtual void getMLAncestralSeqs(std::vector<std::string>* ,std::vector<std::string>* ){}

    virtual void setPermanentInsertion(int ){}

};

#endif
