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
#include <cstdlib>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <ctime>
#include "config.h"
#include "writefile.h"
#include "treenode.h"
#include "hirschberg.h"
#include "site.h"
#include "fullprobability.h"
#include "postprobability.h"
#include "characterprobability.h"
#include "ancestralnode.h"
#include "terminalnode.h"
#include "readalignment.h"

using namespace std;

extern float defaultBranchLength;
extern float fixedBranchLength;
extern float branchScalingFactor;
extern bool MAXBRANCH;
extern bool FIXEDBRANCH;
extern bool FOREVER_FOR_PA;

string tmpNodeName;

AncestralNode::~AncestralNode()
{

    delete lChild;
    delete rChild;

    if (siteLength > 0)
        delete []siteIndex;

    delete seq;
}

AncestralNode::AncestralNode(string s)
        : TreeNode()
{
    root = false;
    terminal = false;
    seq = 0;
    siteLength = 0;
    branchLength = 0;

    rooted = true;
    ln = s.substr(0,s.find(","));
    rn = s.substr(s.find(",")+1);
    ld = defaultBranchLength;
    rd = defaultBranchLength;

    if(ln.find(":") != string::npos)
    {
        ld = atof(ln.substr(ln.find(":")+1).c_str());
        ln = ln.substr(0,ln.find(":"));
    }
    if(rn.find(":") != string::npos)
    {
        rd = atof(rn.substr(rn.find(":")+1).c_str());
        rn = rn.substr(0,rn.find(":"));
    }
    if (ln.at(0)=='#' && ln.at(ln.length()-1)=='#')
    {
        lInternal = true;
    }
    else
    {
        lInternal = false;
        lChild = new TerminalNode(ln,ld);
    }
    if (rn.at(0)=='#' && rn.at(rn.length()-1)=='#')
    {
        rInternal = true;
    }
    else
    {
        rInternal = false;
        rChild = new TerminalNode(rn,rd);
    }

    ld *= branchScalingFactor;

    if (MAXBRANCH)
    {
        if (ld>fixedBranchLength)
            ld=fixedBranchLength;
    }

    if (FIXEDBRANCH)
    {
        ld=fixedBranchLength;
    }

    rd *= branchScalingFactor;

    if (MAXBRANCH)
    {
        if (rd>fixedBranchLength)
            rd=fixedBranchLength;
    }

    if (FIXEDBRANCH)
    {
        rd=fixedBranchLength;
    }

    if (ld<minBrL)
        ld=minBrL;
    if (rd<minBrL)
        rd=minBrL;

    groupName = "null";
}

void AncestralNode::partlyAlignSequences()
{

    lChild->partlyAlignSequences();
    rChild->partlyAlignSequences();

    if (lChild->getGroupName() !="null" && lChild->getGroupName() == rChild->getGroupName())
    {
        FOREVER = false;
        this->readThisNode();
        groupName = lChild->getGroupName();
        FOREVER = FOREVER_FOR_PA;
    }
    else
    {
        this->alignThisNode();
    }
}

/*
 * Recursions for an alignment from scratch:
 *  - HirschbergAlignment is either exact or guided
 * Recursions for full probability
 *  - FullBand is either exact or within a band and needs max l1*k space
 */
void AncestralNode::alignSequences()
{
//    cout<<endl<<nodeName+": children "+lChild->getNodeName()+" and "+rChild->getNodeName()<<endl;

    lChild->alignSequences();
    rChild->alignSequences();
    this->alignThisNode();
}

void AncestralNode::alignThisNode()
{

    char prop[20];
    sprintf(prop,"(%i/%i)",alignedNodes,totalNodes-1);

    currentNode = nodeName+prop;
    if (NOISE>0)
        cout<<endl<<nodeName+prop+": aligning "+lChild->getNodeName()+" and "+rChild->getNodeName()<<endl;

    hmm->alignmentModel(this);

    PhyloMatchScore *pms = new PhyloMatchScore(lChild->getSequence(),rChild->getSequence());
    int time1 = time(0);

    Hirschberg* hp = new Hirschberg();
    hp->alignSeqs(lChild->getSequence(),rChild->getSequence(),pms);
    if (NOISE>0)
        cout<<"Hirschberg: "<<hp->getMaxScore()<<"; time "<<(time(0)-time1)<<"s"<<endl;

    delete hp;
    delete pms;

    if (!lChild->isTerminal())
    {
        AncestralSequence *a1 = static_cast<AncestralSequence*>(lChild->getSequence());
        a1->setRealIndex(true);
    }
    if (!rChild->isTerminal())
    {
        AncestralSequence *a2 = static_cast<AncestralSequence*>(rChild->getSequence());
        a2->setRealIndex(false);
    }
    time1 = time(0);

    if (DOPOST)
    {

        PhyloMatchScore *pms = new PhyloMatchScore(lChild->getSequence(),rChild->getSequence());

        if (NOISE>=0 && SCREEN)
        {
            for (unsigned int i=0; i<message.length(); i++)
            {
                cout<<'\b';
            }

            message = currentNode+": computing full probability               ";

            cout<<message;
            cout.flush();
        }
        else if (NOISE>0)
        {
            cout<<currentNode+": computing full probability"<<endl;
        }

        FullProbability* fp = new FullProbability(lChild->getSequence(),rChild->getSequence(),pms);

        if (FULLBAND)
            fp->alignBand();
        else
            fp->alignSeqs();

        if (NOISE>0)
            cout <<"FullProbability: "<< fp->getMaxFwdScore()<<" "<<fp->getMaxBwdScore()<<" "<<fp->getMaxFwdScore()-fp->getMaxBwdScore()<<"; time "<<(time(0)-time1)<<"s"<<endl;

        time1 = time(0);

        PostProbability* pp = new PostProbability(lChild->getSequence(),rChild->getSequence(),fp->getMaxFwdScore(),pms);

        if (NOISE>0)
            cout<<"PostProbability: time "<<(time(0)-time1)<<"s"<<endl;

        delete fp;
        delete pp;

        time1 = time(0);

        delete pms;
    }

    CharacterProbability *cp = new CharacterProbability(lChild->getSequence(),rChild->getSequence());

    if (NOISE>0)
        cout<<"CharacterProbability: "<< cp->getFwdScore()<<" "<<cp->getBwdScore()<<"; time "<<(time(0)-time1)<<"s"<<endl;

    delete cp;


    seq = new AncestralSequence();
    seq->setChildGaps(lChild->getSequence(),rChild->getSequence());


    if (DOTS)
    {
        int l = getSequence()->length();
        for (int i=0; i<l; i++)
        {
            if (seq->isPermInsertion(i))
            {
                if (seq->getLIndex(i)<0)
                    rChild->setPermanentInsertion(seq->getRIndex(i));
                if (seq->getRIndex(i)<0)
                    lChild->setPermanentInsertion(seq->getLIndex(i));
            }
        }
    }

    // debugging: print each intermediate MA
    if (PRINTNODES)
    {

        int n = getTerminalNodeNumber();
        int l = getSequence()->length();
        int nState = hmm->getNStates();

        vector<string> nms;
        this->getTerminalNames(&nms);

        vector<string> sqs;
        for (int i=0; i<n; i++)
        {
            string s = "";
            sqs.push_back(s);
        }

        vector<string>::iterator si = sqs.begin();


        vector<string> col;

        char* alignment;
        if (CODON)
        {
            alignment = new char[n*l*3];
        }
        else
        {
            alignment = new char[n*l];
        }

        for (int i=0; i<l; i++)
        {
            col.clear();
            this->getCharactersAt(&col,i);
            vector<string>::iterator cb = col.begin();
            vector<string>::iterator ce = col.end();

            si = sqs.begin();
            int j=0;
            for (; cb!=ce; cb++,si++,j++)
            {

                *si+=*cb;

                if (CODON)
                {
                    alignment[j*l*3+i*3] = cb->at(0);
                    alignment[j*l*3+i*3+1] = cb->at(1);
                    alignment[j*l*3+i*3+2] = cb->at(2);
                }
                else
                {
                    alignment[j*l+i] = cb->at(0);
                }
            }

        }

        if (CODON)
            l*=3;


        WriteFile* wfa = new WriteFile();
        wfa->writeSeqs((outfile+"_"+nodeName).c_str(),&nms,&sqs,8);
        delete wfa;

        l = getSequence()->length();

        if (WRITEXML)
        {
            ofstream seqout((outfile+"_"+nodeName+".xml").c_str());

            si = nms.begin();

            // header
            seqout<<"<ms_alignment>"<<endl;
            // tree
            string* treeStr = new string();
            int sInd = 1;
            this->writeNewick(treeStr,&sInd);
            seqout<<"<newick>"<<*treeStr<<"</newick>"<<endl;
            delete treeStr;

            // nodes
            seqout<<"<nodes>"<<endl;
            // terminal nodes

            int ll = l;
            if (CODON)
                ll*=3;

            for (int j=0; j<n; j++)
            {
                seqout<<"<leaf id=\"seq"<<j+1<<"\" name=\""<<(*si++)<<"\">"<<endl;
                seqout<<"  <sequence>"<<endl<<"    ";
                for (int i=0; i<ll; i++)
                {
                    seqout<<alignment[j*ll+i];
                }
                seqout<<endl;
                seqout<<"  </sequence>"<<endl<<"</leaf>"<<endl;
            }

            sqs.clear();
            nms.clear();

            // internal nodes
            this->setSiteLength(l);
            for (int i=0; i<l; i++)
            {
                this->setSiteIndex(i,i);
            }

            map<string,string> anc_seqs;
            this->outputXml(&seqout,&anc_seqs,false);
            seqout<<"</nodes>"<<endl<<"<model>"<<endl;

            // model
            for (int k=0; k<nState; k++)
            {
                seqout<<"  <probability id=\""<<k+1<<"\" name=\""<<hmm->getStName(k)<<"\" ";
                seqout<<"color=\""<<hmm->getDrawCl(k)<<"\" style=\""<<hmm->getDrawPt(k)<<"\" ";
                seqout<<"offset=\""<<hmm->getDrawOf(k)<<"\" show=\"yes\"/>"<<endl;
            }
            seqout<<"  <probability id=\""<<nState+1<<"\" name=\"postprob\" color=\"gray\" ";
            if (DOPOST)
                seqout<<"style=\"bar\" show=\"yes\"/>"<<endl;
            else
                seqout<<"style=\"bar\" show=\"no\"/>"<<endl;
            seqout<<"</model>"<<endl<<"</ms_alignment>"<<endl;
        }

        delete []alignment;

    }

    alignedNodes++;

    lChild->getSequence()->cleanSpace();
    rChild->getSequence()->cleanSpace();

}

void AncestralNode::readAlignment()
{
    lChild->readAlignment();
    rChild->readAlignment();
    this->readThisNode();
}

void AncestralNode::readThisNode()
{
//    cout<<"AncestralNode::readThisNode() "<<nodeName<<endl;
    vector<int> path;
    Sequence *seq1 = lChild->getSequence();
    Sequence *seq2 = rChild->getSequence();

    if (NOISE>1)
        cout<<*(seq1->getGappedSeq())<<endl<<*(seq2->getGappedSeq())<<endl;

    string* ancSeq = new string();
    int i;
    if (!CODON)
    {
        FOR(i,seq1->gappedLength())
        {
            bool c1 = seq1->prealignedGapAt(i);
            bool c2 = seq2->prealignedGapAt(i);
            if (NOISE>1)
            {
                cout<<i<<"/"<<seq1->gappedLength()<<" ";
                cout<<c1<<" "<<c2<<endl;
            }
            if (c1 && c2)
            {
                ancSeq->append("-");
            }
            else if (!c1 && c2)
            {
                path.push_back(0);
                ancSeq->append("A");
            }
            else if (c1 && !c2)
            {
                path.push_back(1);
                ancSeq->append("A");
            }
            else if (!c1 && !c2)
            {
                path.push_back(2);
                ancSeq->append("A");
            }
        }
    }
    else
    {
        for (i=0; i<seq1->gappedLength(); i+=3)
        {
            bool c1a = seq1->prealignedGapAt(i);
            bool c1b = seq1->prealignedGapAt(i+1);
            bool c1c = seq1->prealignedGapAt(i+2);
            bool c2a = seq2->prealignedGapAt(i);
            bool c2b = seq2->prealignedGapAt(i+1);
            bool c2c = seq2->prealignedGapAt(i+2);

            if ( ( (c1a && c1b && c1c) || (!c1a && !c1b && !c1c) ) &&
                    ( (c2a && c2b && c2c) || (!c2a && !c2b && !c2c) ) )
            {
                //ok
            }
            else
            {
                cout<<"ReadAlignment: Error reading the alignment. Gaps not following codon structure. Exiting.\n\n";
                exit(-1);
            }

            if (NOISE>1)
            {
                cout<<i<<"/"<<seq1->gappedLength()<<" ";
                cout<<c1a<<c1b<<c1c<<" "<<c2a<<c2b<<c2c<<endl;
            }
            if (c1a && c2a)
            {
                ancSeq->append("---");
            }
            else if (!c1a && c2a)
            {
                path.push_back(0);
                ancSeq->append("AAA");
            }
            else if (c1a && !c2a)
            {
                path.push_back(1);
                ancSeq->append("AAA");
            }
            else if (!c1a && !c2a)
            {
                path.push_back(2);
                ancSeq->append("AAA");
            }
        }
    }

    ////////

    char prop[20];
    sprintf(prop,"(%i/%i)",alignedNodes,totalNodes-1);

    currentNode = nodeName+prop;
    if (NOISE>0)
        cout<<endl<<nodeName+prop+": reading "+lChild->getNodeName()+" and "+rChild->getNodeName()<<endl;


    hmm->alignmentModel(this);

    PhyloMatchScore *pms = new PhyloMatchScore(lChild->getSequence(),rChild->getSequence());
    int time1 = time(0);


    ReadAlignment* ra = new ReadAlignment();
    ra->readSeqs(lChild->getSequence(),rChild->getSequence(),pms,this,&path);
    if (NOISE>0)
        cout<<"ReadAlignment: "<<ra->getMaxScore()<<"; time "<<(time(0)-time1)<<"s"<<endl;

    delete ra;
    delete pms;


    if (!lChild->isTerminal())
    {
        AncestralSequence *a1 = static_cast<AncestralSequence*>(lChild->getSequence());
        a1->setRealIndex(true);
    }
    if (!rChild->isTerminal())
    {
        AncestralSequence *a2 = static_cast<AncestralSequence*>(rChild->getSequence());
        a2->setRealIndex(false);
    }
    time1 = time(0);

    if (DOPOST)
    {

        PhyloMatchScore *pms = new PhyloMatchScore(lChild->getSequence(),rChild->getSequence());

        if (NOISE>=0 && SCREEN)
        {
            for (unsigned int i=0; i<message.length(); i++)
            {
                cout<<'\b';
            }

            message = currentNode+": computing full probability               ";

            cout<<message;
            cout.flush();
        }
        else if (NOISE>0)
        {
            cout<<currentNode+": computing full probability"<<endl;
        }

        FullProbability* fp = new FullProbability(lChild->getSequence(),rChild->getSequence(),pms);

        if (FULLBAND)
            fp->alignBand();
        else
            fp->alignSeqs();

        if (NOISE>0)
            cout <<"FullProbability: "<< fp->getMaxFwdScore()<<" "<<fp->getMaxBwdScore()<<" "<<fp->getMaxFwdScore()-fp->getMaxBwdScore()<<"; time "<<(time(0)-time1)<<"s"<<endl;

        time1 = time(0);

        PostProbability* pp = new PostProbability(lChild->getSequence(),rChild->getSequence(),fp->getMaxFwdScore(),pms);

        if (NOISE>0)
            cout<<"PostProbability: time "<<(time(0)-time1)<<"s"<<endl;

        delete fp;
        delete pp;

        time1 = time(0);

        delete pms;
    }

    CharacterProbability *cp = new CharacterProbability(lChild->getSequence(),rChild->getSequence());

    if (NOISE>0)
        cout<<"CharacterProbability: "<< cp->getFwdScore()<<" "<<cp->getBwdScore()<<"; time "<<(time(0)-time1)<<"s"<<endl;

    delete cp;

    seq = new AncestralSequence();
    seq->setChildGaps(lChild->getSequence(),rChild->getSequence());
    seq->setGappedSeq(ancSeq);


    alignedNodes++;

    lChild->getSequence()->cleanSpace();
    rChild->getSequence()->cleanSpace();

    delete ancSeq;
}




void AncestralNode::setPermanentInsertion(int i)
{
    if (i<0 || (seq->getLIndex(i)>=0 && seq->getRIndex(i)>=0))
        return;

    this->getSequence()->setPermInsertion(i);

    if (isLInternal())
    {
        lChild->setPermanentInsertion(seq->getLIndex(i));
    }
    if (isRInternal())
    {
        rChild->setPermanentInsertion(seq->getRIndex(i));
    }
}
int AncestralNode::getTerminalNodeNumber()
{
    int n = 0;
    n += lChild->getTerminalNodeNumber();
    n += rChild->getTerminalNodeNumber();
    return n;
}

int AncestralNode::getInternalNodeNumber()
{
    int n = 1;
    n += lChild->getInternalNodeNumber();
    n += rChild->getInternalNodeNumber();
    return n;
}

void AncestralNode::getNames(vector<string>* nms)
{
    lChild->getNames(nms);
    nms->push_back(nodeName);
    rChild->getNames(nms);
}

void AncestralNode::getTerminalNames(vector<string>* nms)
{
    lChild->getTerminalNames(nms);
    rChild->getTerminalNames(nms);
}

void AncestralNode::getInternalNames(vector<string>* nms)
{
    lChild->getInternalNames(nms);
    rChild->getInternalNames(nms);
    nms->push_back(nodeName);
}

// ClustaW tree - no names
void AncestralNode::setCharString(vector<string>* sns,vector<string>* sqs)
{
    lChild->setCharString(sns,sqs);
    rChild->setCharString(sns,sqs);
}

// user-defined - number of tree nodes and sequences may not match
void AncestralNode::setCharString(vector<string>* sns,vector<string>* sqs,int* count)
{
    lChild->setCharString(sns,sqs,count);
    rChild->setCharString(sns,sqs,count);
}


void AncestralNode::getCharStrings(vector<string>* sqs)
{
    lChild->getCharStrings(sqs);
    rChild->getCharStrings(sqs);
}


void AncestralNode::getThisAlignmentPostProbAt(double* p,int i)
{
    if (i>=0)
    {
        if (seq->postProbAt(i)>=0)
            (*p) = seq->postProbAt(i);
    }
}

void AncestralNode::getLowestAlignmentPostProbAt(double* p,int i)
{
    if (seq->getLIndex(i)>=0)
    {
        lChild->getLowestAlignmentPostProbAt(p,seq->getLIndex(i));
    }

    if (seq->getRIndex(i)>=0)
    {
        rChild->getLowestAlignmentPostProbAt(p,seq->getRIndex(i));
    }

    double tp = 1.0;

    this->getThisAlignmentPostProbAt(&tp,i);

    if (tp<(*p))
        (*p)=tp;

}

void AncestralNode::outputXml(std::ofstream* out,map<string,string> *anc_seqs,bool triple)
{
    lChild->outputXml(out,anc_seqs,triple);
    rChild->outputXml(out,anc_seqs,triple);

    (*out)<<"<node id=\""<<nodeName<<"\">"<<endl;

    map<string,string>::iterator it = anc_seqs->find(nodeName);

    if(it != anc_seqs->end())
    {
        (*out)<<"  <sequence>"<<endl<<"    "<<it->second<<endl<<"  </sequence>"<<endl;
    }

    int nState = hmm->getNStates();

    if(nState>1)
    {
        for (int k=0; k<hmm->getNStates(); k++)
        {
            (*out)<<"  <probability id=\""<<k+1<<"\">"<<endl<<"    ";

            for (int m=0; m<siteLength; m++)
            {

                if (m>0)
                {
                    (*out)<<",";
                }

                int i = siteIndex[m];

                if (CODON || triple)
                {
                    if (i<0 || (SKIPINS && getSequence()->isInsertion(i)) )
                    {
                        (*out)<<"-1,-1,-1";
                    }
                    else
                    {
                        if (seq->stateProbAt(k,i)>=0)
                        {
                            int t = (int)(seq->stateProbAt(k,i)*100+0.5);
                            (*out)<<t<<","<<t<<","<<t;
                        }
                        else
                        {
                            (*out)<<"0,0,0";
                        }
                    }
                }
                else
                {
                    if (i<0 || (SKIPINS && getSequence()->isInsertion(i)) )
                    {
                        (*out)<<"-1";
                    }
                    else
                    {
                        if (seq->stateProbAt(k,i)>=0)
                        {
                            (*out)<<(int)(seq->stateProbAt(k,i)*100+0.5);
                        }
                        else
                        {
                            (*out)<<"0";
                        }
                    }
                }
            }
            (*out)<<endl<<"  </probability>"<<endl;
        }
    }

    if (DOPOST)
    {
        (*out)<<"  <probability id=\""<<nState+1<<"\">"<<endl<<"    ";

        double a;
        double* ap = &a;

        for (int m=0; m<siteLength; m++)
        {
            if (m>0)
            {
                (*out)<<",";
            }
            int i = siteIndex[m];
            if (CODON || triple)
            {
                if (i<0 || getSequence()->isInsertion(i))
                {
                    (*out)<<"-1,-1,-1";
                }
                else
                {
                    *ap = 0.0;
                    getThisAlignmentPostProbAt(ap,i);
                    *ap=(int)((*ap)*100+0.5);
                    (*out)<<*ap<<","<<*ap<<","<<*ap;
                }
            }
            else
            {
                if (i<0 || getSequence()->isInsertion(i))
                {
                    (*out)<<"-1";
                }
                else
                {
                    *ap = 0.0;
                    getThisAlignmentPostProbAt(ap,i);
                    *ap=(int)((*ap)*100+0.5);
                    (*out)<<*ap;
                }
            }
        }

        (*out)<<endl<<"  </probability>"<<endl;
    }
    (*out)<<"</node>"<<endl;
}

void AncestralNode::writeNewick(std::string* tree,int* sInd)
{
    *tree += '(';
    lChild->writeNewick(tree,sInd);
    *tree += ',';
    rChild->writeNewick(tree,sInd);
    *tree += + ')';
    *tree += nodeName;
    char str[10];
    sprintf(str,":%.5f",branchLength);
    *tree += str;

    return;
}

void AncestralNode::getNewickBrl(string* tree)
{
    *tree += '(';
    lChild->getNewickBrl(tree);
    *tree += ',';
    rChild->getNewickBrl(tree);
    *tree += + ')';
    char str[10];
    sprintf(str,":%.5f",branchLength);
    *tree += str;

    return;
}

void AncestralNode::getLabelledNewickBrl(string* tree)
{
    *tree += '(';
    lChild->getLabelledNewickBrl(tree);
    *tree += ',';
    rChild->getLabelledNewickBrl(tree);
    *tree += + ')';
    *tree += nodeName;
    char str[10];
    sprintf(str,":%.5f",branchLength);
    *tree += str;

    return;
}

void AncestralNode::getCleanNewick(string* tree)
{
    *tree += "(";
    this->lChild->getNewickBrl(tree);
    *tree += ",";
    this->rChild->getNewickBrl(tree);
    *tree += + ");";
    return;
}

void AncestralNode::getNewick(string* tree)
{
    *tree += '(';
    lChild->getNewick(tree);
    *tree += ',';
    rChild->getNewick(tree);
    *tree += + ')';
    *tree += nodeName;

    return;
}

void AncestralNode::getNexusTree(std::string* tree, int *count)
{
    *tree += '(';
    lChild->getNexusTree(tree,count);
    *tree += ',';
    rChild->getNexusTree(tree,count);
    *tree += + ')';
    char str[10];
    sprintf(str,":%.5f",branchLength);
    *tree += str;

    return;
}

void AncestralNode::getMLAncestralSeqs(vector<string>* sqs,vector<string>* nms)
{
    lChild->getMLAncestralSeqs(sqs,nms);
    rChild->getMLAncestralSeqs(sqs,nms);
    sqs->push_back(*(this->getSequence()->getMLsequence()));
    nms->push_back(nodeName);
}

void AncestralNode::writeAncCharacters(int *parSite,int iteration)
{

    int *insSite = new int[siteLength];
    int h;
    FOR(h,siteLength)
    {
        insSite[h]=0;
        int i = siteIndex[h];
        if (i>=0)
        {
            insSite[h] = this->getSequence()->isInsertion(i);
        }
    }

    lChild->writeAncCharacters(insSite,iteration);
    rChild->writeAncCharacters(insSite,iteration);

    delete []insSite;

    FILE *ancPro = fopen((outfile+"."+itos(iteration)+".ancprof").c_str(),"a");

    fprintf(ancPro,"## %s ##\n",nodeName.c_str());
    string alpha = hmm->getAlphabet();
    int sAlpha = alpha.length();

    int nState = hmm->getNStates();
    int maxState = -1;
    float maxProb = -HUGE_VAL;
    int j,k;


    FOR(h,siteLength)
    {
        int i = siteIndex[h];

        if (i<0)
        {
            fprintf(ancPro,"- %i (-1) ;",h+1);
            FOR(k,nState)
            {
                FOR(j,sAlpha-1)
                {
                    fprintf(ancPro," 0,");
                }
                fprintf(ancPro," 0;");
            }
        }
        else
        {

            if (this->getSequence()->isInsertion(i))
            {
                fprintf(ancPro,"+ %i",h+1);
            }
            else if (this->getSequence()->isGap(i) && parSite[h])
            {
                fprintf(ancPro,"* %i",h+1);
            }
            else
            {
                fprintf(ancPro,"  %i",h+1);
            }

            FOR(k,nState)
            {
                if (this->getSequence()->stateProbAt(k,i)>maxProb)
                {
                    maxProb = this->getSequence()->stateProbAt(k,i);
                    maxState = k;
                }
            }

            fprintf(ancPro," (%i) ;",maxState+1);


            FOR(k,nState)
            {
                double sum = 0;

                FOR(j,sAlpha)
                {
                    double t = this->getSequence()->mlCharProbAt(j,i,k);
                    if (LOGVALUES)
                        sum += exp(t);
                    else
                        sum += t;
                }

                FOR(j,sAlpha-1)
                {
                    double t = this->getSequence()->mlCharProbAt(j,i,k);
                    if (LOGVALUES)
                        fprintf(ancPro," %.3f,",exp(t)/sum);
                    else
                        fprintf(ancPro," %.3f,",t/sum);
                }
                double t = this->getSequence()->mlCharProbAt(sAlpha-1,i,k);

                if (LOGVALUES)
                    fprintf(ancPro," %.3f;",exp(t)/sum);
                else
                    fprintf(ancPro," %.3f;",t/sum);
            }
        }
        fprintf(ancPro,"\n");
    }
    fprintf(ancPro,"\n");
    fclose(ancPro);
}

void AncestralNode::getAncCharactersAt(vector<string>* col,int i,bool parentIns)
{
    if (i<0)
    {
        for (int j=0; j<getInternalNodeNumber(); j++)
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

    }
    else
    {

        lChild->getAncCharactersAt(col,seq->getLIndex(i),this->getSequence()->isInsertion(i));
        rChild->getAncCharactersAt(col,seq->getRIndex(i),this->getSequence()->isInsertion(i));

        if (this->getSequence()->isInsertion(i) || ( parentIns && this->getSequence()->isGap(i) ) )   /*e090626*/
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
        else
        {
            string alpha = hmm->getAlphabet();
            int sAlpha = alpha.length();
            if (CODON)
                sAlpha /= 3;

            int nState = hmm->getNStates();
            int maxState = -1;
            float maxProb = -HUGE_VAL;
            int j,k;

            FOR(k,nState)
            {
                if (this->getSequence()->stateProbAt(k,i)>maxProb)
                {
                    maxProb = this->getSequence()->stateProbAt(k,i);
                    maxState = k;
                }
            }

            if (LOGVALUES)
            {
                float ms = -HUGE_VAL;

                int mi = -1;
                FOR(j,sAlpha)
                {
                    if (this->getSequence()->mlCharProbAt(j,i,maxState)>= ms)
                    {
                        ms = this->getSequence()->mlCharProbAt(j,i,maxState);
                        mi = j;
                    }
                }

                if (mi>=0)
                {
                    if (CODON)
                    {
                        col->push_back(alpha.substr(mi*3,3));
                    }
                    else
                    {
                        col->push_back(string(1,alpha.at(mi)));
                    }
                }
                else
                {
                    cout<<"impossible index: site "<<i<<", "<<mi<<endl;
                }

            }
            else
            {
                float ms = 0;

                int mi = -1;
                FOR(j,sAlpha)
                {
                    if (this->getSequence()->mlCharProbAt(j,i,maxState) >= ms)
                    {
                        ms = this->getSequence()->mlCharProbAt(j,i,maxState);
                        mi = j;
                    }
                }
                if (mi>=0)
                {
                    if (CODON)
                    {
                        col->push_back(alpha.substr(mi*3,3));
                    }
                    else
                    {
                        col->push_back(string(1,alpha.at(mi)));
                    }
                }
                else
                {
                    cout<<"impossible index: site "<<i<<", "<<mi<<endl;
                }
            }
        }
    }
}

void AncestralNode::getAllCharactersAt(vector<string>* col,int i,bool parentIns)
{
    if (i<0)
    {
        for (int j=0; j<getInternalNodeNumber()+getTerminalNodeNumber(); j++)
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

    }
    else
    {

        lChild->getAllCharactersAt(col,seq->getLIndex(i),this->getSequence()->isInsertion(i));

        if (this->getSequence()->isInsertion(i) || ( parentIns && this->getSequence()->isGap(i) ) )   /*e090626*/
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
        else
        {
            string alpha = hmm->getAlphabet();
            int sAlpha = alpha.length();
            if(CODON)
                sAlpha /= 3;

            int nState = hmm->getNStates();
            int maxState = -1;
            float maxProb = -HUGE_VAL;
            int j,k;

            FOR(k,nState)
            {
                if (this->getSequence()->stateProbAt(k,i)>maxProb)
                {
                    maxProb = this->getSequence()->stateProbAt(k,i);
                    maxState = k;
                }
            }

            if (LOGVALUES)
            {
                float ms = -HUGE_VAL;

                int mi = -1;
                FOR(j,sAlpha)
                {
                    if (this->getSequence()->mlCharProbAt(j,i,maxState)>= ms)
                    {
                        ms = this->getSequence()->mlCharProbAt(j,i,maxState);
                        mi = j;
                    }
                }

                if (mi>=0)
                {
                    if (CODON)
                    {
                        col->push_back(alpha.substr(mi*3,3));
                    }
                    else
                    {
                        col->push_back(string(1,alpha.at(mi)));
                    }
                }
                else
                {
                    cout<<"impossible index: site "<<i<<", "<<mi<<endl;
                }

            }
            else
            {
                float ms = 0;

                int mi = -1;
                FOR(j,sAlpha)
                {
                    if (this->getSequence()->mlCharProbAt(j,i,maxState) >= ms)
                    {
                        ms = this->getSequence()->mlCharProbAt(j,i,maxState);
                        mi = j;
                    }
                }
                if (mi>=0)
                {
                    if (CODON)
                    {
                        col->push_back(alpha.substr(mi*3,3));
                    }
                    else
                    {
                        col->push_back(string(1,alpha.at(mi)));
                    }
                }
                else
                {
                    cout<<"impossible index: site "<<i<<", "<<mi<<endl;
                }
            }
        }
        rChild->getAllCharactersAt(col,seq->getRIndex(i),this->getSequence()->isInsertion(i));
    }
}

void AncestralNode::getCharactersAt(vector<string>* col,int i)
{
    if (seq->getLIndex(i)<0)
    {
        if (DOTS && this->getSequence()->isPermInsertion(i))
        {
            for (int j=0; j<lChild->getTerminalNodeNumber(); j++)
            {
                if (CODON)
                {
                    col->push_back("...");
                }
                else
                {
                    col->push_back(".");
                }
            }
        }
        else
        {
            for (int j=0; j<lChild->getTerminalNodeNumber(); j++)
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
        }
    }
    else
    {
        lChild->getCharactersAt(col,seq->getLIndex(i));
    }
    if (seq->getRIndex(i)<0)
    {
        if (DOTS && this->getSequence()->isPermInsertion(i))
        {
            for (int j=0; j<rChild->getTerminalNodeNumber(); j++)
            {
                if (CODON)
                {
                    col->push_back("...");
                }
                else
                {
                    col->push_back(".");
                }
            }
        }
        else
        {
            for (int j=0; j<rChild->getTerminalNodeNumber(); j++)
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
        }
    }
    else
    {
        rChild->getCharactersAt(col,seq->getRIndex(i));
    }
}



void AncestralNode::setSiteLength(int l)
{
    lChild->setSiteLength(l);
    rChild->setSiteLength(l);

    if (siteLength>0)
        delete []siteIndex;

    siteIndex = new int[l];
    siteLength = l;
}

void AncestralNode::setSiteIndex(int site,int index)
{
    siteIndex[site] = index;

    if (index>=0)
    {
        lChild->setSiteIndex(site,getSequence()->getLIndex(index));
    }
    else
    {
        lChild->setSiteIndex(site,-1);
    }
    if (index>=0)
    {
        rChild->setSiteIndex(site,getSequence()->getRIndex(index));
    }
    else
    {
        rChild->setSiteIndex(site,-1);
    }
}


void AncestralNode::printChildAlignment(TreeNode *node,string filename)
{

    int n = node->getTerminalNodeNumber();
    int l = node->getSequence()->length();

    vector<string> nms;
    node->getTerminalNames(&nms);

    vector<string> sqs;
    for (int i=0; i<n; i++)
    {
        string s = "";
        sqs.push_back(s);
    }

    vector<string>::iterator si = sqs.begin();
    vector<string> col;

    char* alignment;
    if (CODON)
    {
        alignment = new char[n*l*3];
    }
    else
    {
        alignment = new char[n*l];
    }

    for (int i=0; i<l; i++)
    {
        col.clear();
        node->getCharactersAt(&col,i);
        vector<string>::iterator cb = col.begin();
        vector<string>::iterator ce = col.end();

        si = sqs.begin();
        int j=0;
        for (; cb!=ce; cb++,si++,j++)
        {

            *si+=*cb;

            if (CODON)
            {
                alignment[j*l*3+i*3] = cb->at(0);
                alignment[j*l*3+i*3+1] = cb->at(1);
                alignment[j*l*3+i*3+2] = cb->at(2);
            }
            else
            {
                alignment[j*l+i] = cb->at(0);
            }
        }

    }

    if (CODON)
        l*=3;


    WriteFile* wfa = new WriteFile();
    wfa->writeSeqs((filename).c_str(),&nms,&sqs,8);
    delete wfa;

    delete []alignment;

}
