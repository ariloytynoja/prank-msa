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
#include <map>
#include <fstream>
#include <sstream>
#include <list>
#include "readnewick.h"
#include "readfile.h"
#include "writefile.h"
#include "guidetree.h"
#include "progressivealignment.h"
#include "hirschberg.h"
#include "readalignment.h"
#include "node.h"
#include "exonerate_reads.h"
#include "mafft_alignment.h"

using namespace std;

ProgressiveAlignment::~ProgressiveAlignment()
{

}

ProgressiveAlignment::ProgressiveAlignment(string treefile,string seqfile,string dnafile)
{

    // Write general info unless silenced
    //
    if (NOISE>=0)
        this->showInfo();

    Exonerate_reads er;
    if (!CONVERT && EXONERATE && !er.test_executable())
    {
        cout<<"The executable for Exonerate not found. Fast alignment anchoring is not used.\n";
        EXONERATE = false;
    }

    // Backtranslate predefined protein alignment to DNA and exit.
    //
    if (BACKTRANSLATE)
        this->backTranslate();



    // Get the sequence data
    //
    vector<string> names;
    vector<string> sequences;
    bool isDna;

    this->getSequenceData(&names,&sequences,&isDna);


    // Convert predefined alignment to required format and exit.
    //
    if (CONVERT)
        this->convertSequencesOnly(&names,&sequences,&isDna);


    this->cleanupSeqNames(&names);


    if (isDna && TRANSLATE)
        this->translateSequences(&names,&sequences,&isDna);


    this->makeSettings(isDna);


    // Get the guidetree -- or generate one
    //
    string tree;
    this->getGuideTree(&names,&sequences,&tree,isDna);


    vector<string>::iterator si = sequences.begin();
    int longest = 0; int slongest = 0;
    for (; si!=sequences.end(); si++)
    {
        if ((int)si->length()>longest)
        {
            slongest = longest;
            longest = (int)si->length();
        }
    }

    Site *sites = new Site();
    sites->setASize(hmm->getASize());
    sites->setNState(hmm->getNStates());
    sites->setMatrices(longest,slongest);


    map<string,TreeNode*> nodes;

    ReadNewick rn;
    rn.buildTree(tree,&nodes);

    AncestralNode* root = static_cast<AncestralNode*>(nodes[rn.getRoot()]);

    int nsqs = 0;
    root->setCharString(&names,&sequences,&nsqs);

    if (nsqs!=root->getTerminalNodeNumber())
    {
        cout<<"Names in sequence file "<<seqfile<<" and guidetree "<<treefile<<" do not match!"<<endl;
        exit(-1);
    }

    if (PREALIGNED)
    {
        if (TREEFROMALIGNMENT || PRINTTREE)
            this->printNewickTree(root,outfile+".0.dnd");

        ReadAlignment ra;
        ra.initialiseMatrices(longest+2);

        if (NOISE>=0)
            cout<<"Reading multiple alignment."<<endl;

        root->setTotalNodes();
        root->readAlignment();

        printAlignment(root,&names,&sequences,0,isDna);
        ra.cleanUp();
    }

    else if (PARTLYALIGNED)
    {
        ReadAlignment ra;
        ra.initialiseMatrices(longest+2);

        Hirschberg hir;
        hir.initialiseMatrices((int)(((float)longest+2)*initialMatrixSize));

        if (NOISE>=0)
            cout<<"Finishing partially aligned alignment."<<endl;

        root->setTotalNodes();
        root->partlyAlignSequences();

        printAlignment(root,&names,&sequences,0,isDna);

        ra.cleanUp();
        hir.cleanUp();
    }

    else
    {
        Hirschberg hir;
        if (EXONERATE)
            hir.initialiseMatrices(initialAnchDist);
        else
            hir.initialiseMatrices((int)(((float)longest+2)*initialMatrixSize));

        root->setTotalNodes();

        if (PRINTTREE)
            this->printNewickTree(root,outfile+".1.dnd");

        bool saveDOPOST = DOPOST;
        if (TWICE)
            DOPOST = false;

        if (NOISE>=0)
            cout<<"Generating multiple alignment."<<endl;

        root->alignSequences();

        if (TWICE)
        {

            if (NOISE>=0)
                cout<<endl;

            printAlignment(root,&names,&sequences,1,isDna);

            this->getNewSequences(root,&names,&sequences);

            GuideTree gt;
            gt.computeTree(&sequences,&names,isDna);
            tree = gt.getTree();

            if (NOISE>0)
                cout<<tree<<endl;

            delete root;
            nodes.clear();

            rn.buildTree(tree,&nodes);


            root = static_cast<AncestralNode*>(nodes[rn.getRoot()]);
            nsqs = 0;

            this->removeGaps(&sequences);
            root->setCharString(&names,&sequences,&nsqs);

            root->setTotalNodes();

            if (PRINTTREE)
                this->printNewickTree(root,outfile+".2.dnd");

            DOPOST = saveDOPOST;

            if (NOISE>=0)
                cout<<"Generating improved multiple alignment."<<endl;

            root->alignSequences();
        }

        hir.cleanUp();
        printAlignment(root,&names,&sequences,2,isDna);
    }

    sites->deleteMatrices();
    delete sites;

    nodes.clear();

    delete root;

}

void ProgressiveAlignment::printAlignment(AncestralNode *root,vector<string> *nms,vector<string> *seqs,int iteration, bool isDna)
{

    int l = root->getSequence()->length();

    nms->clear();
    root->getTerminalNames(nms);

    vector<string>::iterator si = seqs->begin();
    for (; si!=seqs->end(); si++)
    {
        si->clear();
    }

    vector<string> col;

    for (int i=0; i<l; i++)
    {

        col.clear();
        root->getCharactersAt(&col,i);
        vector<string>::iterator cb = col.begin();
        vector<string>::iterator ce = col.end();

        si = seqs->begin();
        for (; cb!=ce; cb++,si++)
        {
            *si+=*cb;
        }
    }

    if (CODON)
        l*=3;

    if (!TRANSLATE)
    {

        WriteFile wfa;
        string file = outfile+"."+itos(iteration)+formatExtension(format);
        wfa.writeSeqs(file.c_str(),nms,seqs,format,isDna,root,false);

        if (WRITEXML)
            printXml(root,iteration,false);

    }
    else
    {

        WriteFile wfa;
        TranslateSequences trseq;
        string file = outfile+".pep."+itos(iteration)+formatExtension(format);
        wfa.writeSeqs(file.c_str(),nms,seqs,format,false,root,false);

        printXml(root,iteration,false);

        vector<string> dnaSeqs;
        if (!trseq.translateDNA(nms,seqs,&dnaSeqs))
        {
            cout<<"Backtranslation failed. Exiting."<<endl;
            exit(-1);
        }


        file = outfile+".nuc."+itos(iteration)+formatExtension(format);
        wfa.writeSeqs(file.c_str(),nms,&dnaSeqs,format,true,root,true);

        if (WRITEXML)
            printXml(root,iteration,true);
    }

     if (WRITEANC || WRITEANCSEQ)
        printAncestral(root,nms,seqs,iteration);

}


void ProgressiveAlignment::printXml(AncestralNode *root,int iteration,bool translate)
{

    int n = root->getTerminalNodeNumber();
    int l = root->getSequence()->length();
    int nState = hmm->getNStates();

    char* alignment;
    if (CODON || translate)
    {
        alignment = new char[n*l*3];
    }
    else
    {
        alignment = new char[n*l];
    }

    this->getAlignmentMatrix(root,alignment,translate);

    string type = "";

    if (TRANSLATE && !translate)
        type="pep.";
    else if (TRANSLATE && translate)
        type="nuc.";

    ofstream seqout((outfile+"."+type+itos(iteration)+".xml").c_str());

    vector<string> nms;
    root->getTerminalNames(&nms);
    vector<string>::iterator si = nms.begin();

    // header
    seqout<<"<ms_alignment>"<<endl;
    // tree
    string* treeStr = new string();
    int sInd = 1;
    root->writeNewick(treeStr,&sInd);
    seqout<<"<newick>"<<*treeStr<<"</newick>"<<endl;
    delete treeStr;

    // nodes
    seqout<<"<nodes>"<<endl;
    // terminal nodes
    for (int j=0; j<n; j++)
    {
        seqout<<"<leaf id=\"seq"<<j+1<<"\" name=\""<<(*si++)<<"\">"<<endl;
        seqout<<"  <sequence>"<<endl<<"    ";
        if (CODON || translate)
        {
            for (int i=0; i<l*3; i++)
            {
                seqout<<alignment[j*l*3+i];
            }
        }
        else
        {
            for (int i=0; i<l; i++)
            {
                seqout<<alignment[j*l+i];
            }
        }
        seqout<<endl;
        seqout<<"  </sequence>"<<endl<<"</leaf>"<<endl;
    }
    nms.clear();

    // internal nodes
    map<string,string> anc_seqs;
    this->getAncestralAlignmentSeqs(root,&anc_seqs);


    root->setSiteLength(l);
    for (int i=0; i<l; i++)
    {
        root->setSiteIndex(i,i);
    }
    root->outputXml(&seqout,&anc_seqs,translate);
    seqout<<"</nodes>"<<endl;

    if(nState>1 || DOPOST)
        seqout<<"<model>"<<endl;

    // model
    if(nState>1)
    {
        for (int k=0; k<nState; k++)
        {
            seqout<<"  <probability id=\""<<k+1<<"\" name=\""<<hmm->getStName(k)<<"\" ";
            seqout<<"color=\""<<hmm->getDrawCl(k)<<"\" style=\""<<hmm->getDrawPt(k)<<"\" ";
            seqout<<"offset=\""<<hmm->getDrawOf(k)<<"\"";
            if (nState>1)
            {
                seqout<<" show=\"yes\"/>"<<endl;
            }
            else
            {
                seqout<<" show=\"no\"/>"<<endl;
            }
        }
    }
    if (DOPOST)
        seqout<<"  <probability id=\""<<nState+1<<"\" name=\"postprob\" color=\"gray\" style=\"bar\" show=\"yes\"/>"<<endl;

    if(nState>1 || DOPOST)
        seqout<<"</model>"<<endl;

    seqout<<"</ms_alignment>"<<endl;

    delete []alignment;

}

void ProgressiveAlignment::getAlignmentMatrix(AncestralNode *root,char* alignment,bool translate)
{
    int n = root->getTerminalNodeNumber();
    int l = root->getSequence()->length();

    if (!translate)
    {

        vector<string> col;
        for (int i=0; i<l; i++)
        {
            col.clear();
            root->getCharactersAt(&col,i);
            vector<string>::iterator cb = col.begin();
            vector<string>::iterator ce = col.end();

            int j=0;
            for (; cb!=ce; cb++)
            {
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
                j++;
            }
        }
    }
    else
    {

        char *tmp = new char[n*l];

        vector<string> col;
        for (int i=0; i<l; i++)
        {
            col.clear();
            root->getCharactersAt(&col,i);
            vector<string>::iterator cb = col.begin();
            vector<string>::iterator ce = col.end();

            int j=0;
            for (; cb!=ce; cb++,j++)
            {
                tmp[j*l+i] = cb->at(0);
            }
        }

        vector<string> names;
        root->getTerminalNames(&names);
        vector<string>::iterator si = names.begin();
        vector<string> prot;

        for (int j=0; j<n; j++)
        {
            string seq = "";
            for (int i=0; i<l; i++)
            {
                seq+=tmp[j*l+i];
            }
            prot.push_back(seq);
        }

        si = prot.begin();

        vector<string> dna;

        TranslateSequences trseq;

        if (!trseq.translateDNA(&names,&prot,&dna))
        {
            cout<<"Backtranslation failed. Exiting."<<endl;
            exit(-1);
        }

        si = dna.begin();
        for (int j=0; j<n; j++)
        {
            for (int i=0; i<l*3; i++)
            {
                alignment[j*l*3+i] = si->at(i);
            }
            si++;
        }
    }

}

void ProgressiveAlignment::printAncestral(AncestralNode *root,vector<string> *nms,vector<string> *sqs,int iteration)
{

    string tree = "";
    if (WRITEANCSEQ)
    {
        root->getLabelledNewickBrl(&tree);
        tree += ';';


        ofstream ancTre((outfile+"."+itos(iteration)+".anc.dnd").c_str());
        ancTre<<tree<<endl;
        ancTre.close();

        ofstream ancSeq((outfile+"."+itos(iteration)+".anc.fas").c_str());

        int n = root->getInternalNodeNumber()+root->getTerminalNodeNumber();
        int l = root->getSequence()->length();

        char* alignment;
        if (CODON)
            alignment = new char[n*l*3];
        else
            alignment = new char[n*l];


        this->getFullAlignmentMatrix(root,alignment);

        vector<string> anms;
        root->getNames(&anms);

        vector<string>::iterator ni = anms.begin();
        int j=0;
        for (; ni!=anms.end(); j++,ni++)
        {
            ancSeq<<">"<<*ni<<endl;
            for (int i=0; i<l; i++)
            {
                ancSeq<<alignment[j*l+i];
            }
            ancSeq<<endl;
        }

        delete []alignment;
        ancSeq.close();

    }

    if (WRITEANC)
    {

        root->getNewick(&tree);


        vector<string> anms;
        root->getInternalNames(&anms);

        ofstream ancSeq((outfile+"."+itos(iteration)+".ancseq").c_str());

        ancSeq<<"# "<<tree<<endl;

        vector<string>::iterator si = sqs->begin();
        vector<string>::iterator ni = nms->begin();
        for (; ni!=nms->end(); si++,ni++)
        {
            ancSeq<<">"<<*ni<<endl;
            ancSeq<<*si<<endl;
        }


        int n = root->getInternalNodeNumber();
        int l = root->getSequence()->length();

        char* alignment;
        if (CODON)
            alignment = new char[n*l*3];
        else
            alignment = new char[n*l];


        this->getAncestralAlignmentMatrix(root,alignment);


        ni = anms.begin();
        int j=0;
        for (; ni!=anms.end(); j++,ni++)
        {
            ancSeq<<">"<<*ni<<endl;
            for (int i=0; i<l; i++)
            {
                ancSeq<<alignment[j*l+i];
            }
            ancSeq<<endl;
        }

        delete []alignment;
        ancSeq.close();


        FILE *ancPro = fopen((outfile+"."+itos(iteration)+".ancprof").c_str(),"w");
        fclose(ancPro);

        int *insSite = new int[l];
        int i;
        FOR(i,l)
        {
            insSite[i]=0;
        }
        root->writeAncCharacters(insSite,iteration);

        delete []insSite;

    }
    return;
}

void ProgressiveAlignment::getAncestralAlignmentMatrix(AncestralNode *root,char* alignment)
{
    vector<string> col;
    int n = root->getInternalNodeNumber();
    int l = root->getSequence()->length();

    int i;
    FOR(i,l)
    {

        col.clear();
        root->getAncCharactersAt(&col,i,0);
        vector<string>::iterator cb = col.begin();
        vector<string>::iterator ce = col.end();

        int j=0;
        for (; cb!=ce; cb++)
        {
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
            j++;
        }
    }
}

void ProgressiveAlignment::getFullAlignmentMatrix(AncestralNode *root,char* alignment)
{
    vector<string> col;
    int n = root->getInternalNodeNumber()+root->getTerminalNodeNumber();
    int l = root->getSequence()->length();

    int i;
    FOR(i,l)
    {

        col.clear();
        root->getAllCharactersAt(&col,i,0);
        vector<string>::iterator cb = col.begin();
        vector<string>::iterator ce = col.end();

        int j=0;
        for (; cb!=ce; cb++)
        {
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
            j++;
        }
    }
}

void ProgressiveAlignment::getAncestralAlignmentSeqs(AncestralNode *root,map<string,string> *anc_seqs)
{
    int l = root->getSequence()->length();
    int n = root->getInternalNodeNumber();

    char* anc_alignment;
    if (CODON)
        anc_alignment = new char[n*l*3];
    else
        anc_alignment = new char[n*l];

    this->getAncestralAlignmentMatrix(root,anc_alignment);

    vector<string> anms;
    root->getInternalNames(&anms);

    vector<string>::iterator ni = anms.begin();
    int j=0;
    for (; ni!=anms.end(); j++,ni++)
    {
        stringstream ss;
        for (int i=0; i<l; i++)
        {
            ss<<anc_alignment[j*l+i];
        }
        anc_seqs->insert(pair<string,string>(*ni,ss.str()));
    }

    delete []anc_alignment;

}

string ProgressiveAlignment::formatExtension(int format)
{
    if (format==1)
    {
        return ".igs";
    }
    else if (format==2)
    {
        return ".gen";
    }
    else if (format==3)
    {
        return ".nbr";
    }
    else if (format==4)
    {
        return ".emb";
    }
    else if (format==6)
    {
        return ".dst";
    }
    else if (format==7)
    {
        return ".fch";
    }
    else if (format==8)
    {
        return ".fas";
    }
    else if (format==11 || format==12 || format==18 || format==19)
    {
        return ".phy";
    }
    else if (format==14)
    {
        return ".pir";
    }
    else if (format==15)
    {
        return ".msf";
    }
    else if (format==17)
    {
        return ".nex";
    }
    else
    {
        return "";
    }
}
