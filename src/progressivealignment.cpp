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
#include <set>
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
#include "bppancestors.h"

using namespace std;

ProgressiveAlignment::~ProgressiveAlignment(){}

ProgressiveAlignment::ProgressiveAlignment(string treefile,string seqfile,string dnafile)
{

    // Write general info unless silenced
    //
    if (NOISE>=0)
        this->showInfo();

    Exonerate_reads er;
    if (!CONVERT && EXONERATE && !er.test_executable())
    {
//        cout<<"The executable for Exonerate not found. Fast alignment anchoring is not used.\n";
        EXONERATE = false;
    }

    // Backtranslate predefined protein alignment to DNA and exit.
    //
    if (BACKTRANSLATE)
    {
        this->backTranslate();
        exit(0);
    }

    // Convert predefined alignment to required format and exit.
    //
    if (CONVERT)
    {
        this->convertSequencesOnly();
        exit(0);
    }


    // Get the sequence data
    //
    vector<string> names;
    vector<string> sequences;
    bool isDna;

    this->getSequenceData(&names,&sequences,&isDna);

    this->cleanupSeqNames(&names);


    if (isDna && TRANSLATE)
    {
        this->translateSequences(&names,&sequences);
        isDna = false;
    }


//    map<string,string> org_stuff;
//    for(int i=0;i<names.size();i++)
//        org_stuff.insert(org_stuff.begin(),pair<string,string>(names.at(i),sequences.at(i)));
//    cout<<"org "<<org_stuff.size()<<endl;


//    cout<<"check 1\n";
//    this->checkStuff(&org_stuff,&names,&sequences);
//    cout<<"check 1\n";


    // Make setting and set the alignment model
    //
    this->makeSettings(isDna);
    this->setHMModel(&sequences,isDna);

     // Find the lengths and reserve space
    //
    int longest = 0; int slongest = 0;
    sites = new Site();
    this->findLongestSeq(&sequences,&longest,&slongest,sites);

     // Get the guidetree -- or generate one
    //
    string tree;
    this->getGuideTree(&names,&sequences,&tree,isDna);

    if(TREEONLY)
    {
        cout<<"\n\nWriting\n";
        cout<<" - estimated tree to '"<<outfile<<".dnd'\n\n";

        string name = outfile+".dnd";
        ofstream seqout(name.c_str());
        seqout<<tree<<endl;
        seqout.close();

        sites->deleteMatrices();
        delete sites;

        exit(0);
    }
    // Build the tree structure and get its root
    //
    map<string,TreeNode*> nodes;

    ReadNewick rn;
    rn.buildTree(tree,&nodes);

    AncestralNode* root = static_cast<AncestralNode*>(nodes[rn.getRoot()]);

//    string tmpstr;
//    root->getNewickBrl(&tmpstr);
//    cout<<tmpstr<<endl;

    // If an old tree is provided, mark the shared sub-trees
    //
    this->checkOldTree(root,&sequences);


    // Now set the sequences ...
    //
    int nsqs = 0;
    root->setCharString(&names,&sequences,&nsqs);

    // and check that the sequence names match
    //
    this->checkMatchingNames(root,&names,nsqs);


    /////////////////////////////////
    // Different alignment options //
    /////////////////////////////////


    // Prealigned data: compute ancestral sequences or convert to xml
    //
    if (PREALIGNED)
    {
        this->readAlignment(root,&names,&sequences,isDna,longest);

        int nSubst; int nIns; int nDel; int nInsDel; bool noSuffix=true;
        int bestScore = this->computeParsimonyScore(root,isDna,-1,&nSubst,&nIns,&nDel,&nInsDel,noSuffix);

        cout<<"\nAlignment score: "<<bestScore;
        if(PRINTSCOREONLY)
        {
            if(nInsDel>0)
                cout<<" [ "<<nSubst<<" subst., "<<nIns<<" ins., "<<nDel<<" del., "<<nInsDel<<" indel. ]"<<endl<<endl;
            else
                cout<<" [ "<<nSubst<<" subst., "<<nIns<<" ins., "<<nDel<<" del. ]"<<endl<<endl;
             exit(0);
        }
        else
            cout<<endl;
    }

    // Partly aligned data: just do the unaligned nodes
    //
    else if (PARTLYALIGNED)
    {
        this->partlyAlign(root,&names,&sequences,isDna,longest);

        int nSubst; int nIns; int nDel; int nInsDel; bool noSuffix=true;
        int bestScore = this->computeParsimonyScore(root,isDna,-1,&nSubst,&nIns,&nDel,&nInsDel,noSuffix);
        cout<<"\nAlignment score: "<<bestScore<<endl;
    }

    // Alignment based on an old tree: re-align nodes that have changed
    //
    else if (UPDATE)
    {
        this->updateAlignment(root,&names,&sequences,isDna,longest);

        int nSubst; int nIns; int nDel; int nInsDel; bool noSuffix=true;
        int bestScore = this->computeParsimonyScore(root,isDna,-1,&nSubst,&nIns,&nDel,&nInsDel,noSuffix);
        cout<<"\nAlignment score: "<<bestScore<<endl;
    }

    // Regular alignment
    //
    else
    {
        int thisIteration = 1;

        Hirschberg hir;
        if (EXONERATE)
            hir.initialiseMatrices(initialAnchDist);
        else
            hir.initialiseMatrices((int)(((float)longest+2)*initialMatrixSize));

        root->setTotalNodes();


        if (NOISE>=0)
            cout<<"\nGenerating multiple alignment: iteration 1."<<endl;

        root->alignSequences();

        this->updateIndelSites(root);

        if(iterations>1 && WRITEITER)
        {
            cout<<"\n\nWriting\n";
            if (PRINTTREE)
                this->printNewickTree(root,outfile+".1.dnd",true);

            this->printAlignment(root,&names,&sequences,outfile+".1",isDna);
        }

        // Write best so far..
        //
        if (PRINTTREE)
            this->printNewickTree(root,outfile+".best.dnd",false);

        this->printAlignment(root,&names,&sequences,outfile+".best",isDna,false);

        int bestScore = this->computeParsimonyScore(root,isDna);
        cout<<"\nAlignment score: "<<bestScore<<endl;

        thisIteration++;

        while (thisIteration<=iterations)
        {

            map<string,float> subtreesOld;
            if (UPDATESECOND)
                root->getAllSubtrees(&subtreesOld);

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

            if (UPDATESECOND)
            {
                UPDATE = true;
                root->setCharString(&names,&sequences,&nsqs);
                UPDATE = false;
            }
            else
            {
                this->removeGaps(&sequences);
                root->setCharString(&names,&sequences,&nsqs);
            }
            root->setTotalNodes();

            if (NOISE>=0)
                cout<<"\nGenerating multiple alignment: iteration "<<thisIteration<<"."<<endl;

            if (UPDATESECOND)
            {
                longest = sequences.at(0).length();

                ReadAlignment ra;
                ra.initialiseMatrices(longest+2);

                // mark changed subtrees
                root->markRealignSubtrees(&subtreesOld);
                root->updateAlignedSequences();

                ra.cleanUp();
            }
            else
            {
                root->alignSequences();
            }

            this->updateIndelSites(root);


            if(WRITEITER)
            {
                string fname = outfile+"."+itos(thisIteration);

                cout<<"\n\nWriting\n";
                if (PRINTTREE)
                    this->printNewickTree(root,fname+".dnd",true);

                this->printAlignment(root,&names,&sequences,fname,isDna);
            }

            int thisScore = this->computeParsimonyScore(root,isDna,bestScore);
            cout<<"\nAlignment score: "<<thisScore<<endl;

            if(thisScore<bestScore)
            {
                bestScore = thisScore;

                if (PRINTTREE)
                    this->printNewickTree(root,outfile+".best.dnd",false);

                this->printAlignment(root,&names,&sequences,outfile+".best",isDna,false);
            }

            thisIteration++;
        }

        hir.cleanUp();

        //************************************************************************//

        string filename = outfile+".best";

        cout<<"\n\nWriting\n";
        if (PRINTTREE)
            cout<<" - alignment guide tree to '"<<filename<<".dnd'\n";

        if(TRANSLATE)
        {
            if (WRITEXML)
                cout<<" - alignment to '"<<filename<<".[nuc|pep]"<<this->formatExtension(format)<<"' and '"<<filename<<".[nuc|pep]"<<".xml'.\n";
            else
                cout<<" - alignment to '"<<filename<<".[nuc|pep]"<<this->formatExtension(format)<<"'.\n";
        }
        else
        {
            if (WRITEXML)
                cout<<" - alignment to '"<<filename<<this->formatExtension(format)<<"' and '"<<filename<<".xml'\n";
            else
                cout<<" - alignment to '"<<filename<<this->formatExtension(format)<<"'\n";
        }

        if (WRITEANCSEQ)
            cout<<" - ancestors to '"<<filename<<".anc"<<this->formatExtension(format)<<"' and '"<<filename<<".anc.dnd'\n";

        if(LISTEVENTS)
            cout<<" - inferred events to file '"<<filename<<".events'\n";

        //************************************************************************//

    }

    sites->deleteMatrices();
    delete sites;

    nodes.clear();
    delete root;

}


void ProgressiveAlignment::updateIndelSites(AncestralNode *root)
{

    for(int i=0;i<root->getSequence()->length();i++)
            root->updateInsertionSite(i,not root->getSequence()->isInsertion(i));

}


void ProgressiveAlignment::printAlignment(AncestralNode *root,vector<string> *nms,vector<string> *seqs,string filename, bool isDna,bool verbose)
{
    if(verbose)
    {
        if(TRANSLATE)
        {
            if (WRITEXML)
                cout<<" - alignment to '"<<filename<<".[nuc|pep]"<<this->formatExtension(format)<<"' and '"<<filename<<".xml'.\n";
            else
                cout<<" - alignment to '"<<filename<<".[nuc|pep]"<<this->formatExtension(format)<<"'.\n";
        }
        else
        {
            if (WRITEXML)
                cout<<" - alignment to '"<<filename<<this->formatExtension(format)<<"' and '"<<filename<<".xml'.\n";
            else
                cout<<" - alignment to '"<<filename<<this->formatExtension(format)<<"'.\n";
        }
    }
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
        string file = filename+formatExtension(format);
        wfa.writeSeqs(file.c_str(),nms,seqs,format,isDna,root,false);

        if (WRITEXML)
            this->printXml(root,filename,false);
    }
    else
    {
        WriteFile wfa;
        TranslateSequences trseq;
        string file = filename+".pep"+formatExtension(format);
        wfa.writeSeqs(file.c_str(),nms,seqs,format,false,root,false);

        if (WRITEXML)
            this->printXml(root,filename,false);

        vector<string> dSeqs;
        if (!trseq.translateDNA(nms,seqs,&dSeqs,&dnaSeqs))
        {
            cout<<"Backtranslation failed. Exiting."<<endl;
            exit(-1);
        }


        bool tmpCODON = CODON;
        bool tmpisDna = isDna;
        bool tmpPREALIGNED = PREALIGNED;

        AncestralNode* tmpRoot = root;

        if(MTTABLE)
        {
            cout<<"Warning: DNA ancestor reconstruction based on the standard codon model.\n";
        }

        // Make setting and set the alignment model
        //
        CODON = true;
        isDna = true;
        PREALIGNED = true;
        this->makeSettings(isDna);
        this->setHMModel(&dSeqs,isDna);

         // Find the lengths and reserve space
        //
        int longest = 0; int slongest = 0;
        this->findLongestSeq(&dSeqs,&longest,&slongest,sites);

         // Get the guidetree -- or generate one
        //
        string tree = "";
        root->getCleanNewick(&tree);

        // Build the tree structure and get its root
        //
        map<string,TreeNode*> nodes;

        ReadNewick rn;
        rn.buildTree(tree,&nodes);

        AncestralNode* codonRoot = static_cast<AncestralNode*>(nodes[rn.getRoot()]);

        // Now set the sequences ...
        //
        int nsqs = 0;
        codonRoot->setCharString(nms,&dSeqs,&nsqs);

        //
        if(!this->sequencesAligned(&dSeqs))
        {
            cout<<"Sequences don't seem to be aligned. Exiting.\n\n";
            exit(0);
        }

        ReadAlignment ra;
        ra.initialiseMatrices(longest+2);

        codonRoot->setTotalNodes();

        bool success = codonRoot->readAlignment();

        if(not success)
        {
            codonRoot->deleteAncestralSeqs();

            cout<<"\nReading the alignment failed. Trying without option '+F'.\n";
            FOREVER = false;
            ra.cleanUp();
            ra.initialiseMatrices(longest+2);

            codonRoot->setTotalNodes();
            success = codonRoot->readAlignment();
            if(not success)
            {
                cout<<"Reading the alignment failed. Terminating.\n";
                exit(-1);
            }
        }

        this->updateIndelSites(codonRoot);

        nms->clear();
        codonRoot->getTerminalNames(nms);

        vector<string>::iterator si = dSeqs.begin();
        for (; si!=dSeqs.end(); si++)
        {
            si->clear();
        }

        vector<string> col;

        for (int i=0; i<l; i++)
        {

            col.clear();
            codonRoot->getCharactersAt(&col,i);
            vector<string>::iterator cb = col.begin();
            vector<string>::iterator ce = col.end();

            si = dSeqs.begin();
            for (; cb!=ce; cb++,si++)
            {
                *si+=*cb;
            }
        }

        l*=3;

        tmpRoot = codonRoot;
        ra.cleanUp();


        file = filename+".nuc"+formatExtension(format);
        wfa.writeSeqs(file.c_str(),nms,&dSeqs,format,true,tmpRoot,true);

        if (WRITEXML)
            this->printXml(tmpRoot,filename,true);

        CODON = tmpCODON;
        isDna = tmpisDna;
        PREALIGNED = tmpPREALIGNED;

        this->makeSettings(isDna);
        this->setHMModel(seqs,isDna);

        delete tmpRoot;


    }

    if (WRITEANCSEQ)
        this->printAncestral(root,filename,isDna,verbose);

}


void ProgressiveAlignment::printXml(AncestralNode *root,string filename,bool translate)
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

    if (TRANSLATE && !translate)
        filename+=".pep";
    else if (TRANSLATE && translate)
        filename+=".nuc";

    ofstream seqout((filename+".xml").c_str());

    vector<string> nms;
    root->getTerminalNames(&nms);
    vector<string>::iterator si = nms.begin();

    // header
    seqout<<"<ms_alignment>"<<endl;
    // tree
    string* treeStr = new string();
    int sInd = 1;
    root->writeLabelledNewick(treeStr,&sInd);
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


void ProgressiveAlignment::reconstructAncestors(AncestralNode *root,bool isDna)
{

    BppAncestors bppa;
    if(bppa.testExecutable())
    {
        if(NOISE>0)
            cout<<"Using BppAncestor to infer ancestral sequences\n";

        map<string,string> aseqs;
        string atree;
        bppa.inferAncestors(root,&aseqs,&atree,isDna);

        ReadNewick rn;
        /*
        // This is not needed with the fixed bppancestor
        //
        map<string,TreeNode*> anodes;
        rn.buildTree(atree,&anodes);

        AncestralNode* aroot = static_cast<AncestralNode*>(anodes[rn.getRoot()]);
        aroot->fixTerminalNodenames();


        map<string,string> subtrees;
        root->getAllSubtreesWithNodename(&subtrees);

        map<string,string> asubtrees;
        aroot->getAllSubtreesWithNodename(&asubtrees);

        map<string,string> asequences;

        map<string,string>::iterator oldit = subtrees.begin();
        for(;oldit!=subtrees.end();oldit++)
        {
            string name = oldit->second;
            if(name == root->getNodeName())
                continue;

            string aname = asubtrees.find(oldit->first)->second;

            aname.erase(aname.begin());
            aname.erase(aname.length()-1);
            asequences.insert(asequences.begin(),pair<string,string>(name,aseqs.find(aname)->second));
        }

        root->setAncSequenceStrings(&asequences);
        */

        root->setAncSequenceStrings(&aseqs);

        // BppAncestor works on unrooted trees; root needs to be done separately
        //
        vector<string> twoseqs;
        twoseqs.push_back(root->getLChild()->getThisSequenceString());
        twoseqs.push_back(root->getRChild()->getThisSequenceString());
        vector<string> twonms;
        twonms.push_back("left");
        twonms.push_back("right");

        stringstream tree;
        tree << "(left:" << root->getLeftBrL()<<",right:"<<root->getRightBrL()<<");";

        map<string,TreeNode*> twonodes;
        rn.buildTree(tree.str(),&twonodes);

        bool tmpPREALIGNED = PREALIGNED;
        PREALIGNED = true;

        bool tmpFOREVER = FOREVER;
        FOREVER = false;

        bool tmpSCREEN = SCREEN;
        SCREEN = false;

        int tmpNOISE = NOISE;
        NOISE = -1;

//        cout<<twoseqs.at(0)<<endl<<twoseqs.at(1)<<endl;
        AncestralNode* tworoot = static_cast<AncestralNode*>(twonodes[rn.getRoot()]);
        int nsqs = 0;
        tworoot->setCharString(&twonms,&twoseqs,&nsqs);

        ReadAlignment ra;
        ra.initialiseMatrices(twoseqs.at(0).length()+2);

        tworoot->setTotalNodes();
        tworoot->readAlignment();

        FOREVER = tmpFOREVER;
        PREALIGNED = tmpPREALIGNED;
        NOISE = tmpNOISE;
        SCREEN = tmpSCREEN;

        string rootstr;
        int len = twoseqs.at(0).length();
        if(CODON)
            len /= 3;
        for(int i=0;i<len;i++)
            rootstr += tworoot->getThisAncCharactersAt(i);

        root->setThisAncSequenceString(rootstr);


        vector<string> aseqs2;
        this->getAncestralAlignmentMatrix(root,&aseqs2);
        root->setAncSequenceGaps(&aseqs2);

        if(NOISE>1)
            cout<<"BppAncestor done\n";

    }
    else
    {
        vector<string> aseqs;
        this->getAncestralAlignmentMatrix(root,&aseqs);

        root->getLChild()->setAncSequenceStrings(&aseqs);
        root->setThisAncSequenceString(&aseqs);
        root->getRChild()->setAncSequenceStrings(&aseqs);

    }
}

void ProgressiveAlignment::setAlignedSequences(AncestralNode *root)
{
    vector<string> aseqs;
    this->getAlignmentMatrix(root,&aseqs,false);
    root->setAlignedSequenceStrings(&aseqs);
}

int ProgressiveAlignment::computeParsimonyScore(AncestralNode *root,bool isDna,int bestScore,int *nSubst,int *nIns,int *nDel,int *nInsDel,bool noSuffix)
{

    this->setAlignedSequences(root);
    this->reconstructAncestors(root,isDna);

    string alpha = hmm->getFullAlphabet();
    int sAlpha = alpha.length();

    map<string,int> alphabet;

    int wordsize = 1;
    if(CODON)
    {
        for(int i=0;i<sAlpha/3;i++)
            alphabet.insert(alphabet.begin(),pair<string,int>(alpha.substr(i*3,3),i));

        alphabet.insert(alphabet.begin(),pair<string,int>("---",-1));
        alphabet.insert(alphabet.begin(),pair<string,int>("...",-1));
        wordsize = 3;
    }
    else
    {
        for(int i=0;i<sAlpha;i++)
            alphabet.insert(alphabet.begin(),pair<string,int>(alpha.substr(i,1),i));

        alphabet.insert(alphabet.begin(),pair<string,int>("-",-1));
        alphabet.insert(alphabet.begin(),pair<string,int>(".",-1));
    }

    root->setAlignedStates(&alphabet,wordsize);

    vector<indelEvent> indels;
    root->getIndelEvents(&indels);

    int substScore = 0;
    for(int i=0;i<root->getSequence()->length();i++)
    {
        int thisScore = 0;
        root->getColumnParsimonyScore(i,&thisScore);
        substScore += thisScore;
    }

    int score = substScore;

    int insCount = 0;
    int delCount = 0;
    int insdelCount = 0;
    int idLength = 0;

    int idscore_1 = 6;
    int idscore_2 = 8;
    int idscore_3 = 9;
    int idscore_4 = 10;

    if(INDELSCORE != "")
    {
        char c;
        stringstream ids(INDELSCORE);
        ids>>idscore_1>>c>>idscore_2>>c>>idscore_3>>c>>idscore_4;
    }

    if(NOISE>0)
        cout<<"Using indel scores:\n  "<<idscore_1<<" for length 1\n  "<<idscore_2<<" for length 2\n  "<<idscore_3<<" for length 3\n  "<<idscore_4<<" for longer than 3\n";

    vector<indelEvent>::iterator ite = indels.begin();
    for(;ite!=indels.end();ite++)
    {

        if(ite->branch == root->getLChild()->getNodeName() || ite->branch == root->getRChild()->getNodeName() )
        {
//            cout<<ite->branch<<endl;
            insdelCount++;
        }
        else
        {
            if(ite->isInsertion)
                insCount++;
            else
                delCount++;
        }

        if(NOTGAP && ite->isTerminal)
        {
            score += idscore_1;
//            cout<<ite->branch<<" "<<ite->alignedStart<<" "<<ite->alignedEnd<<"\n";
        }
        else
        {
            if(ite->length == 1)
                score += idscore_1;
            else if(ite->length == 2)
                score += idscore_2;
            else if(ite->length == 3)
                score += idscore_3;
            else
                score += idscore_4;
        }
        idLength += ite->length;
    }

    if(nSubst)
    {
        *nSubst = substScore;
        *nIns = insCount;
        *nDel = delCount;
        *nInsDel = insdelCount;
    }

    if(NOISE>0)
        cout<<"\nInferred events: "<<substScore<<" subst, "<<insCount<<" ins, "<<delCount<<" dels, "<<insdelCount<<" insdels\n";

    if(LISTEVENTS && ( bestScore<0 || score < bestScore) )
    {

//        cout<<"events "<<score<<" "<<bestScore<<endl;

        vector<substEvent> substs;
        root->getSubstEvents(&substs);

        vector<string> ralphabet;

        if(CODON)
        {
            for(int i=0;i<sAlpha/3;i++)
                ralphabet.push_back(alpha.substr(i*3,3));
        }
        else
        {
            for(int i=0;i<sAlpha;i++)
                ralphabet.push_back(alpha.substr(i,1));
        }

        vector<string> names;
        root->getTerminalNames(&names);
        root->getInternalNames(&names);

        stringstream alloutput;

        for(int n=0;n<names.size();n++)
        {
            stringstream output;
            bool printOut = false;

            string tName = names.at(n);
//            cout<<tName<<endl;
            output<<"\nbranch "<<tName<<endl;

            for(int i=0;i<substs.size();i++)
            {
                if(substs.at(i).branch==tName)
                {
                    output<<" "<<substs.at(i).alignedPos+1<<" "//<<"("<<substs.at(i).realPos<<") "
                        <<ralphabet.at(substs.at(i).pChar)<<" -> "<<ralphabet.at(substs.at(i).dChar)<<"\n";

                    printOut = true;
                }
            }
            for(int i=0;i<indels.size();i++)
            {
                if(indels.at(i).branch==tName)
                {
                    output<<" "<<indels.at(i).alignedStart+1<<".."<<indels.at(i).alignedEnd+1;//<<" ("
//                        <<indels.at(i).realStart<<".."<<indels.at(i).realEnd<<") ";
                    if(indels.at(i).isInsertion)
                        output<<" insertion";
                    else
                        output<<" deletion";
                    if(indels.at(i).isTerminal)
                        output<<" (terminal)";
                    output<<"\n";

                    printOut = true;
                }
            }
            if(printOut)
            {
                alloutput<<output.str();
            }
        }

        if(alloutput.str().length()>0)
        {
            string outname = outfile+".best.events";
            if(noSuffix)
                outname = outfile+".events";

//            if (NOISE>=0)
//                cout<<" - inferred events to file '"<<outname<<"'.\n";

            ofstream seqout(outname.c_str());
            string tree;
            root->getLabelledNewick(&tree);

            seqout<<"\nAlignment topology with node labels:\n\n"<<tree<<endl<<endl;

            seqout<<"Inferred evolutionary events per branch:\n";
            seqout<<alloutput.str()<<"\n";
            seqout.close();
        }
    }

    return score;

}

void ProgressiveAlignment::printAncestral(AncestralNode *root,string filename, bool isDna, bool verbose)
{
    this->setAlignedSequences(root);
    this->reconstructAncestors(root,isDna);

    if(verbose)
        cout<<" - ancestors to '"<<filename<<".anc"<<this->formatExtension(format)<<"' and '"<<filename<<".anc.dnd'.\n";

    string tree = "";
    root->getLabelledNewickBrl(&tree);
    tree += ";";

    ofstream ancTre((filename+".anc.dnd").c_str());
    ancTre<<tree<<endl;
    ancTre.close();

    ofstream ancSeq((filename+".anc.fas").c_str());

    vector<string> anms;
    root->getNames(&anms);

    vector<string> aseqs;
    root->getAllSequenceStrings(&aseqs);


    vector<string>::iterator ni = anms.begin();
    vector<string>::iterator si = aseqs.begin();

    for (; ni!=anms.end(); si++,ni++)
        ancSeq<<">"<<*ni<<endl<<*si<<endl;

   return;
}

void ProgressiveAlignment::getAlignmentMatrix(AncestralNode *root,char* alignment,bool translate)
{
    int n = root->getTerminalNodeNumber();
    int l = root->getSequence()->length();

    if (1) //!translate)
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
        if (!trseq.translateDNA(&names,&prot,&dna,&dnaSeqs))
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

void ProgressiveAlignment::getAlignmentMatrix(AncestralNode *root,vector<string> *aseqs,bool translate)
{
    int n = root->getTerminalNodeNumber();
    int l = root->getSequence()->length();


    for (int i=0; i<n; i++)
        aseqs->push_back(string(""));

    if (!translate)
    {
        vector<string> col;
        for (int i=0; i<l; i++)
        {
            col.clear();
            root->getCharactersAt(&col,i);
            vector<string>::iterator cb = col.begin();
            vector<string>::iterator ce = col.end();
            vector<string>::iterator si = aseqs->begin();

            for (; cb!=ce; cb++,si++)
            {
                if (CODON)
                {
                    *si += cb->at(0);
                    *si += cb->at(1);
                    *si += cb->at(2);
                }
                else
                {
                    *si += cb->at(0);
                }
            }
        }
    }
    else
    {

        vector<string> prot;
        for (int i=0; i<n; i++)
            prot.push_back(string(""));

        vector<string> col;
        for (int i=0; i<l; i++)
        {
            col.clear();
            root->getCharactersAt(&col,i);
            vector<string>::iterator cb = col.begin();
            vector<string>::iterator ce = col.end();
            vector<string>::iterator si = prot.begin();

            for (; cb!=ce; cb++,si++)
            {
                *si += cb->at(0);
            }
        }

        vector<string> names;
        root->getTerminalNames(&names);

        TranslateSequences trseq;
        if (!trseq.translateDNA(&names,&prot,aseqs,&dnaSeqs))
        {
            cout<<"Backtranslation failed. Exiting."<<endl;
            exit(-1);
        }
    }

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
        root->getAncCharactersAt(&col,i,0,root->getSequence()->isPermInsertion(i));
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

void ProgressiveAlignment::getAncestralAlignmentMatrix(AncestralNode *root,vector<string> *aseqs)
{
    vector<string> col;
    int n = root->getInternalNodeNumber();
    int l = root->getSequence()->length();

    int i;
    FOR(i,n)
            aseqs->push_back(string(""));

    FOR(i,l)
    {

        col.clear();
        root->getAncCharactersAt(&col,i,0,root->getSequence()->isPermInsertion(i));

        vector<string>::iterator cb = col.begin();
        vector<string>::iterator ce = col.end();
        vector<string>::iterator si = aseqs->begin();

        for (; cb!=ce; cb++,si++)
        {
            if (CODON)
            {
                *si += cb->at(0);
                *si += cb->at(1);
                *si += cb->at(2);
            }
            else
            {
                *si += cb->at(0);
            }
        }
    }
}

void ProgressiveAlignment::getFullAlignmentMatrix(AncestralNode *root,char* alignment)
{
    vector<string> col;
    int l = root->getSequence()->length();
    int sl = l;
    if(CODON)
        sl *= 3;
    int i;
    FOR(i,l)
    {

        col.clear();
        root->getAllCharactersAt(&col,i,0,root->getSequence()->isPermInsertion(i));
        vector<string>::iterator cb = col.begin();
        vector<string>::iterator ce = col.end();

        int j=0;
        for (; cb!=ce; cb++)
        {
            if (CODON)
            {
                alignment[j*sl+i*3] = cb->at(0);
                alignment[j*sl+i*3+1] = cb->at(1);
                alignment[j*sl+i*3+2] = cb->at(2);
            }
            else
            {
                alignment[j*sl+i] = cb->at(0);
            }
            j++;
        }
    }
}


void ProgressiveAlignment::getFullAlignmentMatrix(AncestralNode *root,vector<string> *aseqs)
{
    vector<string> col;
    int l = root->getSequence()->length();
    int sl = l;
    if(CODON)
        sl *= 3;

    int n = root->getInternalNodeNumber()+root->getTerminalNodeNumber();

    int i;
    FOR(i,n)
        aseqs->push_back(string(""));

    FOR(i,l)
    {

        col.clear();
        root->getAllCharactersAt(&col,i,0,root->getSequence()->isPermInsertion(i));
        vector<string>::iterator cb = col.begin();
        vector<string>::iterator ce = col.end();
        vector<string>::iterator si = aseqs->begin();

        for (; cb!=ce; cb++)
        {
            if (CODON)
            {
                *si += cb->at(0);
                *si += cb->at(1);
                *si += cb->at(2);
            }
            else
            {
                *si += cb->at(0);
            }
        }
    }
}

void ProgressiveAlignment::getAncestralAlignmentSeqs(AncestralNode *root,map<string,string> *anc_seqs)
{
    int l = root->getSequence()->length();
    int n = root->getInternalNodeNumber();

    char* anc_alignment;
    if (CODON)
        l*=3;

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
