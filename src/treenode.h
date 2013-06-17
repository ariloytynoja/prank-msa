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
#include <set>
#include <string>
#include <iostream>
#include "sequence.h"
#include "site.h"
#include "flmatrix.h"

extern float minBrL;
extern int rnd_seed;

struct substEvent
{
    std::string branch;
    int realPos;
    int alignedPos;
    int pChar;
    int dChar;
};

struct indelEvent
{
    std::string branch;
    int realStart;
    int realEnd;
    int alignedStart;
    int alignedEnd;
    int length;
    bool isTerminal;
    bool isInsertion;
};


class TreeNode
{
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

    std::string alignedseqstr; // aligned sequence
    std::vector<int> alignedstates;
    int alignedStartSite;
    int alignedEndSite;

    static int totalNodes;
    static int alignedNodes;

    int siteLength;
    int* siteIndex;         // index for site output

public:
    virtual ~TreeNode();

    std::string nhx_tag;

    virtual void alignSequences() {}
    virtual bool readAlignment() { return true; }
    virtual bool partlyAlignSequences() { return true; }
    virtual void updateAlignedSequences() {}

    virtual int getTerminalNodeNumber() = 0;
    virtual int getInternalNodeNumber() = 0;

    bool realignNode;
    bool LRealign;
    bool RRealign;

    std::string getAlignedSeqStr() { return alignedseqstr; }
    std::vector<int> *getAlignedStates() { return &alignedstates; }

    virtual std::string getGroupName()
    {
        return groupName;
    }

    void setTotalNodes()
    {
        totalNodes = this->getTerminalNodeNumber();
        alignedNodes = 1;
    }

    void setBranchLength(float l)
    {
        if (l<minBrL)
        {
            branchLength = minBrL;    //set the length of the  branch below
        }
        else
        {
            branchLength = l;
        }
    }
    float getBranchLength()
    {
        return branchLength;
    }


    std::string getNodeName()
    {
        return nodeName;    // get node name
    }
    void setNodeName(std::string s)
    {
        nodeName = s;    // set node name
    }

    bool isRoot()
    {
        return root;
    }
    void setRoot()
    {
        root = true;
    }

    bool isRooted()
    {
        return rooted;
    }
    bool isLInternal()
    {
        return lInternal;
    }
    bool isRInternal()
    {
        return rInternal;
    }

    bool isTerminal()
    {
        return terminal;
    }

    float getLeftBrL()
    {
        return ld;
    }
    float getRightBrL()
    {
        return rd;
    }

    int gappedLength();

    virtual void concatenateTerminalNames(std::string *s) = 0;

    virtual void getNames(std::vector<std::string>* nms) = 0;
    virtual void getTerminalNames(std::vector<std::string>* nms) = 0;
    virtual void getInternalNames(std::vector<std::string>* nms) = 0;

    virtual void setCharString(std::vector<std::string>* sns,std::vector<std::string>* sqs,int* count) = 0;
    virtual void setCharString(std::vector<std::string>* sns,std::vector<std::string>* sqs) = 0;
    virtual void getCharStrings(std::vector<std::string>* sqs) = 0;

    virtual void getAllSubtrees(std::map<std::string,float> *subtrees) = 0;
    virtual void getAllSubtreesWithNodename(std::map<std::string,std::string> *subtrees) = 0;
    virtual void getSubtreeBelow(std::string *subtree) = 0;
    virtual void markRealignSubtrees(std::map<std::string,float> *subtrees) = 0;

    virtual bool anyChildNodeRealigned() = 0;

    virtual Sequence* getSequence() = 0;

    void setLChild(TreeNode* tn)
    {
        lChild = tn;    // set child nodes
    }
    void setRChild(TreeNode* tn)
    {
        rChild = tn;
    }

    TreeNode* getLChild()
    {
        return lChild;    // get child nodes
    }
    TreeNode* getRChild()
    {
        return rChild;    // get child nodes
    }

    std::string getLName()
    {
        return ln;    // get child node names
    }
    std::string getRName()
    {
        return rn;
    }


    int hash(const char *str)
    {
        unsigned hash = rnd_seed;
        while (*str)
        {
            hash = hash * 101  +  *str++;
        }
        return hash;
    }

    virtual void getCharactersAt(std::vector<std::string>* ,int,bool t=false ) {}
    virtual void getAncCharactersAt(std::vector<std::string>* ,int ,bool,bool ) {}
    virtual void getAllCharactersAt(std::vector<std::string>* ,int ,bool, bool ) {}

    virtual void setSiteLength(int ) {}
    virtual void setSiteIndex(int ,int ) {}

    virtual void getLowestAlignmentPostProbAt(double*,int) = 0;

    virtual void outputXml(std::ofstream* out,std::map<std::string,std::string> *anc_seqs,bool triple) = 0;

    virtual void writeNewick(std::string* ,int* ) {}
    virtual void writeLabelledNewick(std::string* tree,int* sInd) {}
    virtual void getNewick(std::string* tree) = 0;
    virtual void getLabelledNewickBrl(std::string* tree) = 0;
    virtual void getLabelledNewick(std::string* tree) = 0;
    virtual void getNewickBrl(std::string* tree) = 0;
    virtual void getNexusTree(std::string* tree, int *count) = 0;
    virtual void getNHXBrl(std::string* tree,int *nodeNumber) = 0;

    void getCleanNewick(std::string* tree);
    virtual void getMLAncestralSeqs(std::vector<std::string>* ,std::vector<std::string>* ) {}

    virtual void setPermanentInsertion(int ) {}
    virtual void setAncSequenceStrings(std::vector<std::string>*){}
    virtual void setAncSequenceStrings(std::map<std::string,std::string>*){}
    virtual void getAncSequenceStrings(std::vector<std::string>*){}
    virtual void setAlignedSequenceStrings(std::vector<std::string>*){}
    virtual void getAlignedSequenceStrings(std::vector<std::string>*){}

    virtual void setAncSequenceGaps(std::vector<std::string>*){}

    void getAllSequenceStrings(std::vector<std::string>* aseqs)
    {
        if(!isTerminal())
            lChild->getAllSequenceStrings(aseqs);
        aseqs->push_back(alignedseqstr);
        if(!isTerminal())
            rChild->getAllSequenceStrings(aseqs);
    }
    std::string getThisSequenceString()
    {
        return alignedseqstr;
    }

    virtual void fixTerminalNodenames() = 0;

    virtual void getIndelEvents(std::vector<indelEvent> *indels) = 0;
    virtual void getSubstEvents(std::vector<substEvent> *substs) = 0;

    void getColumnParsimonyScore(int position,int *stateChanges)
    {
        this->getColumnParsimonyScoreAt(position,stateChanges,this->alignedstates.at(position));
    }

    void getColumnParsimonyScoreAt(int position,int *stateChanges,int parentState)
    {
        int thisState = alignedstates.at(position);

        if(thisState != parentState)
        {
//            std::cout<<nodeName<<" "<<parentState<<" "<<thisState<<std::endl;

            if(thisState == -1 || parentState == -1)
            {


            }
            else
            {
                (*stateChanges)++;
            }
        }

        if(!isTerminal())
        {
            lChild->getColumnParsimonyScoreAt(position,stateChanges,thisState);
            rChild->getColumnParsimonyScoreAt(position,stateChanges,thisState);
        }
    }

    void setAlignedStates(std::map<std::string,int> *alphabet,int wordsize)
    {
        if(!terminal)
        {
            lChild->setAlignedStates(alphabet,wordsize);
            rChild->setAlignedStates(alphabet,wordsize);
        }
        alignedstates.clear();
        for(int i=0;i<alignedseqstr.length();i+=wordsize)
        {
            std::string c = alignedseqstr.substr(i,wordsize);
            std::map<std::string,int>::iterator it = alphabet->find(c);
            if(it!=alphabet->end())
                alignedstates.push_back(it->second);
            else
                alignedstates.push_back(-3);
        }
    }

    virtual void deleteAncestralSeqs() {}

};

#endif
