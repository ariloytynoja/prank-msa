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
#include <cstdio>
#include <cstdlib>
#include "terminalnode.h"
#include "config.h"

using namespace std;

extern float fixedBranchLength;
extern float branchScalingFactor;
extern bool MAXBRANCH;
extern bool FIXEDBRANCH;

TerminalNode::~TerminalNode()
{
    if (siteLength > 0)
        delete []siteIndex;

    delete seq;

}


TerminalNode::TerminalNode(string s,float l)
        : TreeNode()
{

    l *= branchScalingFactor;

    if (MAXBRANCH)
    {
        if (l>fixedBranchLength)
            l=fixedBranchLength;
    }

    if (FIXEDBRANCH)
    {
        l=fixedBranchLength;
    }

    if (l<minBrL)
    {
        cout<<"Branch length <"<<minBrL<<". Set to "<<minBrL<<"."<<endl;
        l=minBrL;
    }

    root = false;
    terminal = true;
    seq = 0;
    siteLength = 0;

    groupName = "null";
    nodeName = s;

    branchLength = l;
    ld = rd = 0.0;
}


void TerminalNode::alignSequences(int )
{
    return;
}

int TerminalNode::getTerminalNodeNumber()
{
    return 1;
}

int TerminalNode::getInternalNodeNumber()
{
    return 0;
}

void TerminalNode::getNames(vector<string>* nms)
{
    nms->push_back(nodeName);
}

void TerminalNode::getTerminalNames(vector<string>* nms)
{
    nms->push_back(nodeName);
}

void TerminalNode::getInternalNames(vector<string>* )
{
    return;
}

// ClustaW tree - no names
void TerminalNode::setCharString(vector<string>* sns,vector<string>* sqs)
{
    int index = atoi(nodeName.c_str());
    this->setNodeName(sns->at(index));
//     charString = sqs->at(index);
    seq = new TerminalSequence(&sqs->at(index));
    charString = *seq->getMLsequence();
    if (NOISE>1)
        cout<<nodeName<<endl<<charString<<endl;
}

// user-defined - number of tree nodes and sequences may not match
void TerminalNode::setCharString(vector<string>* sns,vector<string>* sqs,int* count)
{
    vector<string>::iterator ni = sns->begin();
    vector<string>::iterator si = sqs->begin();

    for (; ni!=sns->end(); si++,ni++)
    {
        string seqname = (*ni);
        string tmpGroup = "null";
        if (PARTLYALIGNED)
        {
            size_t pos = seqname.find("_group_");
            if (pos != string::npos)
            {
                tmpGroup = seqname.substr((int)pos);
                seqname = seqname.substr(0,(int)pos);
            }
        }

        if (seqname==nodeName)
        {
            if (NOISE>1)
                cout<<(*ni)<<"\n"<<(*si)<<endl;

            seq = new TerminalSequence(&(*si));
            charString = *seq->getMLsequence();
            groupName = tmpGroup;
            (*count)++;
        }
    }
}


void TerminalNode::getCharStrings(vector<string>* sqs)
{
    sqs->push_back(charString);
}

void TerminalNode::setAnnotation(map<string,FlMatrix*>* annotation)
{
    map<string,FlMatrix*>::iterator smi;
    smi = annotation->find(nodeName);
    if (smi!=annotation->end())
    {
        seq->setAnnotation(smi->second);
    }

}

void TerminalNode::getLowestAlignmentPostProbAt(double* ,int )
{
    return;
}

void TerminalNode::outputXml(std::ofstream* ,bool )
{
    return;
}


void TerminalNode::writeNewick(std::string* tree,int* sInd)
{
    char str[25];
    sprintf(str,"seq%i:%.5f",*sInd,branchLength);
    (*sInd)++;
    *tree += str;

    return;
}

void TerminalNode::getNewickBrl(string* tree)
{
    *tree += nodeName;
    char str[10];
    sprintf(str,":%.5f",branchLength);
    *tree += str;

    return;
}

void TerminalNode::getNexusTree(std::string* tree, int *count)
{
    *tree += itos(*count);
    char str[10];
    sprintf(str,":%.5f",branchLength);
    *tree += str;

    (*count)++;

    return;
}

void TerminalNode::getNewick(string* tree)
{
    *tree += nodeName;

    return;
}

void TerminalNode::getMLAncestralSeqs(vector<string>* ,vector<string>* )
{
    return;
}

void TerminalNode::writeAncCharacters(int *,int )
{
    return;
}

void TerminalNode::getAncCharactersAt(vector<string>* ,int ,bool )
{
    return;
}

void TerminalNode::getCharactersAt(vector<string>* col,int i)
{
    if (i<0)
    {
        if (CODON)
        {
            col->push_back("---");
        }
        else
        {
            col->push_back("-");
        }
    }
    else if (i<getSequence()->length())
    {
        if (CODON)
        {
            col->push_back(charString.substr(i*3,3));
        }
        else
        {
            col->push_back(charString.substr(i,1));
        }
    }
    else
    {
        cout<<nodeName<<": index out of scope ("<<i<<")"<<endl;
    }
}
