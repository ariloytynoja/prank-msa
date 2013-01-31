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
#ifndef PROGRESSIVEALIGNMENT_H
#define PROGRESSIVEALIGNMENT_H

/**
 * Wrapper for sequence loading and multiple alignment.
 */

#include <string>
#include <map>
#include <list>
#include "config.h"
#include "ancestralnode.h"
#include "treenode.h"
#include "readfile.h"
#include "writefile.h"
#include "readnewick.h"
#include "guidetree.h"
#include "translatesequences.h"
#include "node.h"
#include "mafft_alignment.h"
#include "exonerate_reads.h"
#include "bppancestors.h"
#include "readalignment.h"
#include "hirschberg.h"

class ProgressiveAlignment
{
public:
    ProgressiveAlignment(std::string treefile,std::string seqfile,std::string dnafile);
    ~ProgressiveAlignment();

private:
    void getAlignmentMatrix(AncestralNode *root,char* alignment,bool translate);
    void getAlignmentMatrix(AncestralNode *root,vector<string> *aseqs,bool translate);

    void getAncestralAlignmentMatrix(AncestralNode *root,char* alignment);
    void getAncestralAlignmentMatrix(AncestralNode *root,vector<string> *aseqs);
    void getFullAlignmentMatrix(AncestralNode *root,char* alignment);
    void getFullAlignmentMatrix(AncestralNode *root,vector<string> *aseqs);
    void getAncestralAlignmentSeqs(AncestralNode *root,map<string,string> *anc_seqs);

    void readAlignment(AncestralNode *root,vector<string> *names,vector<string> *sequences,bool isDna,int longest,bool verbose=true,bool writeOutput=true)
    {
        if(!this->sequencesAligned(sequences))
        {
            cout<<"Sequences don't seem to be aligned. Exiting.\n\n";
            exit(0);
        }

        ReadAlignment ra;
        ra.initialiseMatrices(longest+2);

        if (NOISE>=0)
            cout<<"Reading multiple alignment."<<endl;

        root->setTotalNodes();
        root->readAlignment();

        if(verbose)
            cout<<"\n\nWriting\n";

        if(writeOutput)
        {
            if (PRINTTREE)
                this->printNewickTree(root,outfile+".dnd",verbose);

            this->printAlignment(root,names,sequences,outfile,isDna,verbose);
        }

        ra.cleanUp();
    }

    void partlyAlign(AncestralNode *root,vector<string> *names,vector<string> *sequences,bool isDna,int longest)
    {
        ReadAlignment ra;
        ra.initialiseMatrices(longest+2);

        Hirschberg hir;
        hir.initialiseMatrices((int)(((float)longest+2)*initialMatrixSize));

        if (NOISE>=0)
            cout<<"Finishing partially aligned alignment."<<endl;

        root->setTotalNodes();
        root->partlyAlignSequences();

        cout<<"\n\nWriting\n";
        if (PRINTTREE)
            this->printNewickTree(root,outfile+".dnd");

        printAlignment(root,names,sequences,0,isDna);

        ra.cleanUp();
        hir.cleanUp();
    }

    void updateAlignment(AncestralNode *root,vector<string> *names,vector<string> *sequences,bool isDna,int longest)
    {
        ReadAlignment ra;
        ra.initialiseMatrices(longest+2);

        Hirschberg hir;
        hir.initialiseMatrices((int)(((float)longest+2)*initialMatrixSize));

        if (NOISE>=0)
            cout<<"Updating partially aligned alignment."<<endl;

        root->setTotalNodes();
        root->updateAlignedSequences();

        cout<<"\n\nWriting\n";
        if (PRINTTREE)
            this->printNewickTree(root,outfile+".dnd");

        printAlignment(root,names,sequences,0,isDna);

        ra.cleanUp();
        hir.cleanUp();
    }

    void reconstructAncestors(AncestralNode *root,bool isDna);
    void setAlignedSequences(AncestralNode *root);
    int computeParsimonyScore(AncestralNode *root,bool isDna,int bestScore=-1);

    void printAlignment(AncestralNode *root,std::vector<std::string> *nms,std::vector<std::string> *sqs,string filename, bool isDna, bool verbose=true);
    void printAncestral(AncestralNode *root,string filename,bool isDna, bool verbose=true);
    void printXml(AncestralNode *root,string filename,bool translate);
    std::string formatExtension(int format);

    std::map<std::string,std::string> dnaSeqs;

    /********************************************/

    void showInfo()
    {
        cout<<endl<<"-----------------\n"<<" PRANK v."<<version<<":\n"<<"-----------------\n\n"<<"Input for the analysis\n";

        if (CONVERT)
        {
            cout<<" - converting '"<<seqfile<<"' to '"<<outfile<<this->formatExtension(format)<<"'"<<endl<<endl;
        }
        else
        {
            if (BACKTRANSLATE)
            {
                cout<<" - creating a DNA alignment '"<<outfile<<this->formatExtension(format)
                    <<"' based on a protein alignment '"<<seqfile<<"' and DNA sequences in '"<<dnafile<<"'"<<endl<<endl;
            }
            else
            {
                if(MERGE)
                {
                    cout<<" - merging alignments in '"<<seqfile1<<"' and '"<<seqfile2<<"'";
                    if (treefile1!="" & treefile2!="")
                        cout<<"\n - using alignment guide trees '"<<treefile1<<"'' and '"<<treefile2<<"'\n";
                    else
                        cout<<"\n - using inferred alignment guide trees\n";
                }
                else
                {
                    if(PREALIGNED)
                        cout<<" - reading alignment '"<<seqfile<<"'";
                    else
                        cout<<" - aligning sequences in '"<<seqfile<<"'";
                    if (treefile!="" & hmmname!="")
                        cout<<"\n - using alignment guide tree '"<<treefile<<"'' and model '"<<hmmname<<"'\n";
                    else if (treefile!="" && oldtreefile!="")
                        cout<<"\n - using alignment guide trees '"<<treefile<<"' (new) and '"<<oldtreefile<<"' (old)\n";
                    else if (treefile!="")
                        cout<<"\n - using alignment guide tree '"<<treefile<<"'\n";
                    else if (hmmname!="")
                        cout<<"\n - using alignment model '"<<hmmname<<"'\n";
                    else
                        cout<<"\n - using inferred alignment guide tree\n";
                }
                if(!FOREVER && !PREALIGNED)
                    cout<<" - option '+F' is not used; it can be enabled with '+F'"<<endl;
                else if(PREALIGNED)
                    cout<<" - option '+F' is not compatible with '-keep'"<<endl;

                Mafft_alignment ma;
                bool mafftOK = ma.test_executable();

                Exonerate_reads er;
                bool exonerateOK = er.test_executable();

                BppAncestors ba;
                bool bppaOK = ba.testExecutable();

                if(mafftOK || exonerateOK || bppaOK)
                {
                    cout<<" - external tools available:\n";
                    if(mafftOK)
                        cout<<"    MAFFT for initial alignment\n";
                    if(exonerateOK)
                        cout<<"    Exonerate for alignment anchoring\n";
                    if(bppaOK)
                        cout<<"    BppAncestor for ancestral state reconstruction\n";
                }
                cout<<endl;
            }
        }
    }

    /********************************************/

    void checkOldTree(AncestralNode *root,vector<string> *sequences)
    {
        // If an old tree is provided, mark the shared sub-trees
        //
        if(oldtreefile!="")
        {
            if(!this->sequencesAligned(sequences))
            {
                cout<<"Sequences don't seem to be aligned. Realignment needed.\n";
                UPDATE = false;

                return;
            }

            ReadNewick rn;
            string oldtree = rn.readFile(oldtreefile.c_str());
            map<string,TreeNode*> oldnodes;

            rn.buildTree(oldtree,&oldnodes);
            AncestralNode* oldroot = static_cast<AncestralNode*>(oldnodes[rn.getRoot()]);

            map<string,float> subtreesOld;
            oldroot->getAllSubtrees(&subtreesOld);

            root->markRealignSubtrees(&subtreesOld);
            UPDATE = true;
        }
    }

    /********************************************/

    void checkMatchingNames(AncestralNode *root,vector<string> *names,int nsqs)
    {

        if (nsqs!=root->getTerminalNodeNumber())
        {
            cout<<"Names in sequence file "<<seqfile<<" and guide tree "<<treefile<<" do not match!"<<endl;
            exit(-1);
        }

        // .. and make sure that the sequence data and the tree match!
        //
        if(nsqs != names->size() && !PRUNEDATA)
        {
            cout<<"Of the "<<names->size()<<" sequences, only "<<nsqs<<" match the tree leaves.\n"
                  "The data can be pruned to match the tree using the flag '-prunedata'. Now exiting.\n\n";

            vector<string> leaf_nms;
            root->getTerminalNames(&leaf_nms);

            set<string> all_nms;
            for(vector<string>::iterator it=leaf_nms.begin();it!=leaf_nms.end();it++)
                all_nms.insert(*it);

            if(names->size()-nsqs>10)
            {
                cout<<"First ten unmatched sequences are:\n";
                int count =0;
                for(int i=0;i<names->size();i++)
                {
                    if(all_nms.find(names->at(i))==all_nms.end())
                    {
                        cout<<"  "<<names->at(i)<<endl;
                        count++;
                    }
                    if(count>=10)
                        break;
                }
                cout<<"\n\n";
            }
            else
            {
                cout<<"The unmatched sequences are:\n";
                for(int i=0;i<names->size();i++)
                {
                    if(all_nms.find(names->at(i))==all_nms.end())
                    {
                        cout<<"  "<<names->at(i)<<endl;
                    }
                }
                cout<<"\n\n";
            }


            exit(-1);
        }
    }

    /********************************************/

    void backTranslate()
    {
        ReadFile rfa;
        TranslateSequences trseq;

        int ns = rfa.readFile(dnafile.c_str());
        if (ns>0)
        {
            vector<string> dnaNames = rfa.getNames();
            vector<string> dnaSequences = rfa.getSeqs();
            bool isDna = rfa.dnaSeqs();
            if (!isDna)
            {
                cout<<"Not DNA in "<<dnafile<<". Exiting."<<endl;
                exit(-1);
            }
            if (!trseq.translateProtein(&dnaNames,&dnaSequences,&dnaSeqs))
            {
                cout<<"Translation of DNA sequences to protein failed. Exiting."<<endl<<endl;
                exit(-1);
            }
        }


        ns = rfa.readFile(seqfile.c_str());
        if (ns>0)
        {
            vector<string> protNames = rfa.getNames();
            vector<string> protSequences = rfa.getSeqs();
            bool isDna = rfa.dnaSeqs();
            if (isDna)
            {
                cout<<"Not protein in "<<seqfile<<". Exiting."<<endl;
                exit(-1);
            }

            vector<string> codons;
            if (!trseq.translateDNA(&protNames,&protSequences,&codons,&dnaSeqs))
            {
                cout<<"Backtranslation of protein sequences to DNA failed. Exiting."<<endl;
                exit(-1);
            }

            WriteFile wfa;
            string file = outfile+formatExtension(format);
            wfa.writeSeqs(file.c_str(),&protNames,&codons,format);

            exit(0);

        }
    }

    /********************************************/

    void getSequenceData(vector<string> *names,vector<string> *sequences,bool *isDna)
    {

        ReadFile rfa;

        if(MERGE)
        {
            int ns1 = rfa.readFile(seqfile1.c_str());

            vector<string> names1;
            vector<string> sequences1;
            bool isDna1;

            if (ns1>=1)
            {
                names1 = rfa.getNames();
                sequences1 = rfa.getSeqs();
                isDna1 = rfa.dnaSeqs();
            }

            vector<string> names2;
            vector<string> sequences2;
            bool isDna2;

            int ns2 = rfa.readFile(seqfile2.c_str());
            if (ns2>=1)
            {
                names2 = rfa.getNames();
                sequences2 = rfa.getSeqs();
                isDna2 = rfa.dnaSeqs();
            }

            if(isDna1 != isDna2 || ns1<1 || ns2<1)
            {
                cout<<"Problem reading sequence files to be merged. Please check '"<<seqfile1<<"' and '"<<seqfile2<<"'\n. Exiting.\n\n";
                exit(0);
            }

            vector<string>::iterator si = sequences1.begin();
            vector<string>::iterator ni = names1.begin();
            for(;si!=sequences1.end();si++,ni++)
            {
                sequences->push_back(*si);
                names->push_back(*ni+"_group_1");
            }

            si = sequences2.begin();
            ni = names2.begin();
            for(;si!=sequences2.end();si++,ni++)
            {
                sequences->push_back(*si);
                names->push_back(*ni+"_group_2");
            }
            *isDna = isDna1;
        }
        else
        {
            int ns = rfa.readFile(seqfile.c_str());
            if (ns>=2)
            {
                *names = rfa.getNames();
                *sequences = rfa.getSeqs();
                *isDna = rfa.dnaSeqs();
            }
            else if (ns==1)
            {
                cout<<"Only one sequence in "<<seqfile<<"!"<<endl;
                *names = rfa.getNames();
                *sequences = rfa.getSeqs();

                WriteFile wfa;
                string file = outfile+formatExtension(format);
                wfa.writeSeqs(file.c_str(),names,sequences,format);

                exit(0);
            }
            else
            {
                cout<<"No sequences found in "<<seqfile<<"!"<<endl;
                exit(-1);
            }
        }


        if (DNA && !isDna)
        {
            cout<<"Warning autodetection suggests protein but DNA model forced.\n";
            *isDna = true;
        }
        if (CODON && !isDna)
        {
            cout<<"Warning autodetection suggests protein but codon model forced.\n";
            *isDna = true;
        }
        if (PROTEIN && isDna)
        {
            cout<<"Warning autodetection suggests DNA but protein model forced.\n";
            *isDna = false;
        }

        if(isDna)
        {
            DNA = true;
            PROTEIN = false;
        }
        else
        {
            DNA = false;
            PROTEIN = true;
        }

        list<string> sorted_names;
        sorted_names.assign(names->begin(),names->end());
        sorted_names.sort();

        string pname = "";
        bool unique = true;

        if (SHORTNAMES)
        {
            for (vector<string>::iterator lit = names->begin(); lit!=names->end(); lit++)
            {

                string n = *lit;
                n = n.substr(0,n.find_first_of(' '));
                *lit = n;
            }
        }

        for (list<string>::iterator lit = sorted_names.begin(); lit!=sorted_names.end(); lit++)
        {

            if (pname==*lit)
            {
                cout<<"Sequence name "<<pname<<" appears more than once!"<<endl;
                unique = false;
            }
            pname=*lit;
        }
        if (!unique)
        {
            cout<<"Sequence names not unique. Exiting."<<endl<<endl;
            exit(-1);
        }


        if (PREALIGNED)
        {
            vector<string>::iterator si = sequences->begin();
            for (; si!=sequences->end(); si++)
            {
                string ts = *si;
                string us = "-";
                for (unsigned int i=0; i<ts.length(); i++)
                {
                    if (ts.at(i)=='.')
                    {
                        ts.replace(i,1,us);
                    }
                }
                *si = ts;
            }
        }

    }

    /********************************************/

    void convertSequencesOnly()
    {

        // Get the sequence data
        //
        vector<string> names;
        vector<string> sequences;
        bool isDna;

        this->getSequenceData(&names,&sequences,&isDna);

        WriteFile wfa;

        if (!TRANSLATE)
        {
            string file = outfile+formatExtension(format);
            wfa.writeSeqs(file.c_str(),&names,&sequences,format);
        }

        else if (isDna && TRANSLATE)
        {
            TranslateSequences trseq;

            if (!trseq.translateProtein(&names,&sequences,&dnaSeqs))
            {
                cout<<"Translation to protein failed. Exiting."<<endl<<endl;

                exit(-1);
            }

            string file = outfile+formatExtension(format);
            wfa.writeSeqs(file.c_str(),&names,&sequences,format);
        }

        else
        {
            cout<<"Inconsistent options: conversion between formats failed. Exiting.\n\n";
        }
    }

    /********************************************/

    void findLongestSeq(vector<string> *sequences,int *longest,int *slongest,Site *sites)
    {
        // Find the lengths and reserve space
        //
        vector<string>::iterator si = sequences->begin();

        for (; si!=sequences->end(); si++)
        {
            if ((int)si->length()>=*longest)
            {
                *slongest = *longest;
                *longest = (int)si->length();
            }
        }

        sites->setASize(hmm->getASize());
        sites->setNState(hmm->getNStates());
        sites->setMatrices(*longest,*slongest);
    }

    /********************************************/

    void cleanupSeqNames(vector<string> *names)
    {
        bool warning = false;
        vector<string>::iterator tit = names->begin();

        for (; tit!=names->end(); tit++)
        {
            string ts = *tit;
            string us = "_";
            for (unsigned int i=0; i<ts.length(); i++)
            {
                if (ts.at(i)==' ' || ts.at(i)==':' || ts.at(i)=='\t' || ts.at(i)=='\n'
                        || ts.at(i)==')' || ts.at(i)=='(' || ts.at(i)==',')
                {
                    ts.replace(i,1,us);
                    warning = true;
                }
            }
//            cout<<*tit<<" "<<ts<<endl;
            *tit = ts;
        }

        if (warning)
        {
            cout<<"Warning: sequence names changed."<<endl<<endl;
        }
    }

    /********************************************/

    void translateSequences(vector<string> *names,vector<string> *sequences)
    {
        TranslateSequences trseq;

        if (!trseq.translateProtein(names,sequences,&dnaSeqs))
        {
            cout<<"Failed to translate DNA sequences to protein. Exiting."<<endl<<endl;
            exit(-1);
        }
    }

    /********************************************/

    void setHMModel(vector<string> *sequences,bool isDna,bool isRna=false)
    {
        // Check that sequences match the model
        if (HASHMM)
        {
            if (isDna && hmm->getAlphabet().length()==20)
            {
                cout<<"Sequences in "<<seqfile<<" seem to be DNA but the alphabet in "<<hmmname<<" is for protein!"<<endl;
                exit(-1);
            }
            else if (!isDna && hmm->getAlphabet().length()==4)
            {
                cout<<"Sequences in "<<seqfile<<" seem to be protein but the alphabet in "<<hmmname<<" is for DNA!"<<endl;
                exit(-1);
            }

            // Or make a new model
        }
        else
        {
            // protein model
            if (isDna)
            {
                // codon model
                if (CODON)
                {
                    hmm = new HMModel;
                    hmm->codonModel();

                }
                // dna model
                else
                {
                    hmm = new HMModel;

                    float freqs[4];
                    bool isRna = false;
                    this->getDNABaseFreqs(sequences,freqs,&isRna);

                    hmm->dnaModel(freqs,isRna);

                }
            }
            else
            {
                hmm = new HMModel;
                hmm->proteinModel();
            }
        }

    }

    void getDNABaseFreqs(vector<string> *sequences,float *freqs,bool isRna=false)
    {
        ReadFile rfa;

        // get nucleotide frequencies - either empirical (NOTE! all seqs used!) or user-defined
        //
        freqs[0]=freqs[1]=freqs[2]=freqs[3]=0;

        // user-defined nucleotide frequencies
        if (dnaFreqs!="")
        {
            int i=0; int j=0;
            while (dnaFreqs.find(",",i)<=dnaFreqs.length())
            {
                freqs[j++]=atof(dnaFreqs.substr(i,dnaFreqs.find(",",i)-i).c_str());
                i=dnaFreqs.find(",",i)+1;
            }
            freqs[3]=atof(dnaFreqs.substr(i).c_str());
            float total=freqs[0]+freqs[1]+freqs[2]+freqs[3];
            freqs[0]/=total;
            freqs[1]/=total;
            freqs[2]/=total;
            freqs[3]/=total;

        }
        // empirical nucleotide frequencies
        else
        {
            rfa.countDnaFreqs(freqs,sequences);
            float total=freqs[0]+freqs[1]+freqs[2]+freqs[3];
            freqs[0]/=total;
            freqs[1]/=total;
            freqs[2]/=total;
            freqs[3]/=total;
        }
    }

    void getPairwiseSubstScores(vector<string> *sequences,IntMatrix *substScores,bool isDna)
    {
        if (isDna)
        {
            HMModel *tmp_hmm = new HMModel;

            float freqs[4];
            bool isRna = false;
            this->getDNABaseFreqs(sequences,freqs,&isRna);

            tmp_hmm->dnaModel(freqs,isRna);
            tmp_hmm->pairwiseModel(substScores,pwDist);

            delete tmp_hmm;
        }
        else
        {
            HMModel *tmp_hmm = new HMModel;

            tmp_hmm->proteinModel();
            tmp_hmm->pairwiseModel(substScores,pwDist);

            delete tmp_hmm;
        }
    }

    void getGuideTree(vector<string> *names,vector<string> *sequences,string *tree,bool isDna)
    {

        ReadNewick rn;

        if (treefile=="" && !MERGE || (MERGE && (treefile1=="" || treefile2=="") ) )
        {
            int time1 = time(0);

            GuideTree gt;

            if(MERGE)
            {
                string tree1;
                string tree2;

                vector<string> sequences1;
                vector<string> sequences2;
                vector<string> names1;
                vector<string> names2;

                vector<string>::iterator si=sequences->begin();
                vector<string>::iterator ni=names->begin();

                for(;si!=sequences->end();si++,ni++)
                {
                    size_t pos = ni->find("group_1");
                    if(pos != string::npos)
                    {
                        sequences1.push_back(*si);
                        names1.push_back(ni->substr(0,pos-1));
                        continue;
                    }

                    pos = ni->find("group_2");
                    if(pos != string::npos)
                    {
                        sequences2.push_back(*si);
                        names2.push_back(ni->substr(0,pos-1));
                    }
                }

                if(this->sequencesAligned(&sequences1) && this->sequencesAligned(&sequences2))
                {

                    gt.computeTree(&sequences1,&names1,isDna);
                    tree1 = gt.getTree();
                    if(sequences1.size()==1)
                    {
                        if(tree1.at(0)=='(')
                            tree1.erase(0,1);

                        if(tree1.at(tree1.size()-2)==')')
                            tree1.erase(tree1.size()-2,1);
                    }
                    else if(sequences1.size()>2)
                    {
                        Node* n1 = new Node(tree1);
                        tree1=n1->rootedTree();
                        delete n1;
                    }

                    gt.computeTree(&sequences2,&names2,isDna);
                    tree2 = gt.getTree();
                    if(sequences2.size()==1)
                    {
                        if(tree2.at(0)=='(')
                            tree2.erase(0,1);

                        if(tree2.at(tree2.size()-2)==')')
                            tree2.erase(tree2.size()-2,1);
                    }
                    else if(sequences2.size()>2)
                    {
                        Node* n2 = new Node(tree2);
                        tree2=n2->rootedTree();
                        delete n2;
                    }
                }
                else
                {
                    cout<<"Sequences don't seem to be aligned. Exiting.\n\n";
                    exit(0);
                }


                float dist1 = -1; float dist2 = -1;

                string merge_tree = "";

                if (TREESTRING)
                    merge_tree = treefile;
                else if(treefile!="")
                    merge_tree = rn.readFile(treefile.c_str());

                if(merge_tree != "")
                {
                    char p,c1,c2,co1,cm,c3,c4,co2;
                    string end;

                    stringstream merge_ss(merge_tree);
                    merge_ss >> p >> c1 >> c2  >> co1 >> dist1 >> cm >> c3  >> c4  >> co2 >> dist2 >> end;

                    if(c2=='2' && c4=='1')
                    {
                        float tmp = dist1;
                        dist1 = dist2;
                        dist2 = tmp;
                    }
                }

                float dist = defaultBranchLength;
                if(mergeBranchLength>0)
                    dist = mergeBranchLength;

                dist /= 2;

                if(dist1<0 || dist2<0)
                {
                    dist1 = dist;
                    dist2 = dist;
                }

                if(tree1.at(tree1.size()-1)==';')
                    tree1.erase(tree1.size()-1);

                if(tree2.at(tree2.size()-1)==';')
                    tree2.erase(tree2.size()-1);

                if (NOISE>0)
                    cout<<"GuideTree; time "<<(time(0)-time1)<<"s"<<endl;

                stringstream ss;
                ss<<"("<<tree1<<":"<<dist1<<","<<tree2<<":"<<dist2<<");";
                *tree = ss.str();
            }


            else
            {
                Mafft_alignment ma;

                if (TREEFROMALIGNMENT)
                {
                    if(!this->sequencesAligned(sequences))
                    {
                        cout<<"Sequences don't seem to be aligned. Exiting.\n\n";
                        exit(0);
                    }
                    gt.computeTree(sequences,names,isDna);
                }
                else if(MAFFTALIGNMENT && ma.test_executable())
                {

                    vector<string> tmp_names;
                    vector<string> tmp_seqs;
                    vector<string>::iterator ni = names->begin();
                    vector<string>::iterator si = sequences->begin();

                    for(;si!=sequences->end();si++,ni++)
                    {
                        tmp_names.push_back(*ni);
                        tmp_seqs.push_back(*si);
                    }

                    this->removeGaps(&tmp_seqs);

                    bool tmp_isDna = isDna;

                    if(isDna && CODON)
                    {
                        this->translateSequences(&tmp_names,&tmp_seqs);
                        tmp_isDna = false;
                    }


                    ma.align_sequences(&tmp_names,&tmp_seqs);


                    if(!this->sequencesAligned(&tmp_seqs))
                    {
                        cout<<"Sequences don't seem to be aligned. Exiting.\n\n";
                        exit(0);
                    }

                    gt.computeTree(&tmp_seqs,&tmp_names,tmp_isDna);

                    if(isDna && CODON)
                    {
                        TranslateSequences trseq;
                        vector<string> dSeqs;
                        if (!trseq.translateDNA(&tmp_names,&tmp_seqs,&dSeqs,&dnaSeqs))
                        {
                            cout<<"Backtranslation failed. Exiting."<<endl;
                            exit(-1);
                        }
                        tmp_seqs.clear();
                        for(int i=0;i<dSeqs.size();i++)
                            tmp_seqs.push_back(dSeqs.at(i));

                        tmp_isDna = true;
                    }

                    if(true)
                    {
                        // Get the tree
                        //
                        string tmp_tree = gt.getTree();

                        // Build the tree    structure and get its root
                        //
                        map<string,TreeNode*> tmp_nodes;

                        ReadNewick rn;
                        rn.buildTree(tmp_tree,&tmp_nodes);

                        AncestralNode* tmp_root = static_cast<AncestralNode*>(tmp_nodes[rn.getRoot()]);

                        // Now set the sequences ...
                        //
                        bool tmpPREALIGNED = PREALIGNED;
                        PREALIGNED = true;
                        int nsqs = 0;
                        tmp_root->setCharString(&tmp_names,&tmp_seqs,&nsqs);
                        PREALIGNED = tmpPREALIGNED;

                        // and check that the sequence names match
                        //
                        this->checkMatchingNames(tmp_root,&tmp_names,nsqs);

                        this->readAlignment(tmp_root,&tmp_names,&tmp_seqs,tmp_isDna,tmp_seqs.at(0).length());

                        int bestScore = this->computeParsimonyScore(tmp_root,tmp_isDna);
                        cout<<"\nInitial alignment score: "<<bestScore;
                        if(CODON)
                            cout<<" (based on protein alignment)"<<endl<<endl;
                        else
                            cout<<endl<<endl;

                        delete tmp_root;
                    }
                }

                // Do not use MAFFT; make own pairwise alignments
                //
                else
                {
                    IntMatrix *substScores;

                    if (isDna)
                        substScores = new IntMatrix(5,5,"pwSubstScores");
                    else
                        substScores = new IntMatrix(21,21,"pwSubstScores");

                    this->getPairwiseSubstScores(sequences,substScores,isDna);

                    gt.computeTree(sequences,names,substScores);

                    delete substScores;
                }

                if (NOISE>0)
                    cout<<"GuideTree; time "<<(time(0)-time1)<<"s"<<endl;

                *tree = gt.getTree();
            }


            if (NOISE>0)
                cout<<"tree "<<*tree<<endl;

        }
        else
        {
            if (TREESTRING && !MERGE)
                *tree = treefile;

            else if(MERGE)
            {

                float dist1 = -1; float dist2 = -1;
                string merge_tree = "";

                if (TREESTRING)
                    merge_tree = treefile;
                else if(treefile!="")
                    merge_tree = rn.readFile(treefile.c_str());

                if(merge_tree != "")
                {
                    char p,c1,c2,co1,cm,c3,c4,co2;
                    string end;

                    stringstream merge_ss(merge_tree);
                    merge_ss >> p >> c1 >> c2  >> co1 >> dist1 >> cm >> c3  >> c4  >> co2 >> dist2 >> end;

                    if(c2=='2' && c4=='1')
                    {
                        float tmp = dist1;
                        dist1 = dist2;
                        dist2 = tmp;
                    }
                }

                float dist = defaultBranchLength;
                if(mergeBranchLength>0)
                    dist = mergeBranchLength;

                dist /= 2;

                if(dist1<0 || dist2<0)
                {
                    dist1 = dist;
                    dist2 = dist;
                }

                string tree1 = rn.readFile(treefile1.c_str());
                string tree2 = rn.readFile(treefile2.c_str());

                if(tree1.at(tree1.size()-1)==';')
                    tree1.erase(tree1.size()-1);

                if(tree2.at(tree2.size()-1)==';')
                    tree2.erase(tree2.size()-1);

                stringstream ss;
                ss<<"("<<tree1<<":"<<dist1<<","<<tree2<<":"<<dist2<<");";
                *tree = ss.str();

//                cout<<tree1<<"\n"<<tree2<<"\n"<<*tree<<endl;
            }

            else
                *tree = rn.readFile(treefile.c_str());

            if (*tree=="")
            {
                cout<<"No tree found in "<<treefile<<"! Exiting."<<endl;
                exit(-1);
            }
        }

//        this->setHMModel(sequences,isDna);


        Node* n = new Node(*tree);
        n->mark_sequences(names);
        int treeleaves = 0;
        int treematches = 0;
        n->countMatchingLeaves(&treeleaves,&treematches);

        if (treeleaves!=treematches)
        {
            if (PRUNETREE)
            {
                n->prune_tree();
                *tree=n->print_tree();
            }
            else if (PARTLYALIGNED)
            {
//                cout<<"Warning: The guide tree names do not match those in the sequence data.\n\n";
            }
            else
            {
                cout<<"The guide tree has "<<treeleaves<<" leaves but only "<<treematches<<" match with the names in the sequence data.\n"
                    "The tree can be pruned to match the data using the flag '-prunetree'. Now exiting.\n\n";

                vector<string> unmatching;
                n->collectUnmatchingLeaves(&unmatching);

                if(unmatching.size()>10)
                {
                    cout<<"First ten unmatched tree leaves are:\n";
                    for(int i=0;i<10;i++)
                        cout<<"  "<<unmatching.at(i)<<endl;

                    cout<<"\n\n";
                }
                else
                {
                    cout<<"The unmatched tree leaves are:\n";
                    for(int i=0;i<unmatching.size();i++)
                        cout<<"  "<<unmatching.at(i)<<endl;

                    cout<<"\n\n";
                }

                exit(-1);
            }
        }
        delete n;
    }

    /********************************************/

    void printNewickTree(AncestralNode* root,string name,bool verbose=false)
    {
        if(verbose)
            cout<<" - alignment guide tree to '"<<name<<"'.\n";

        string tmpTree = "";
        root->getCleanNewick(&tmpTree);

        ofstream seqout(name.c_str());
        seqout<<tmpTree<<endl;
        seqout.close();
    }

    /********************************************/

    void getNewSequences(AncestralNode* root,vector<string> *names,vector<string> *sequences)
    {

        names->clear();
        root->getTerminalNames(names);

        sequences->erase(sequences->begin(),sequences->end());

        vector<string>::iterator ni = names->begin();
        for (; ni!=names->end(); ni++)
        {
            sequences->push_back(string());
        }

        vector<string> col;

        for (int i=0; i<root->getSequence()->length(); i++)
        {
            col.clear();
            root->getCharactersAt(&col,i,false);
            vector<string>::iterator cb = col.begin();
            vector<string>::iterator ce = col.end();

            vector<string>::iterator si = sequences->begin();
            for (; cb!=ce; cb++,si++)
            {
                *si+=*cb;
            }
        }
    }

    /********************************************/

    void removeGaps(vector<string> *sequences)
    {
        vector<string>::iterator si = sequences->begin();
        for (; si!=sequences->end(); si++)
        {
            string s = "";
            for (int i=0; i<(int)si->length(); i++)
            {
                char c = si->at(i);
                if (c!='-')
                {
                    s+=c;
                }
            }
            *si = s;
        }
    }

    /********************************************/

    bool sequencesAligned(vector<string> *seqs)
    {
        int length = -1;
        vector<string>::iterator si = seqs->begin();
        for(;si!=seqs->end();si++)
        {
            if(length<0)
                length = si->length();
            else if(length != si->length())
                return false;
        }
        return true;
    }

    /********************************************/

    void makeSettings(bool isDna)
    {
        if (isDna)
        {
            if (gapRate<0)
                gapRate = dnaGapRate;

            if (gapExt<0)
                gapExt = dnaGapExt;

            if (pwDist<0)
                pwDist = pwDnaDist;

            if (pwGapRate<0)
                pwGapRate = dnaGapRate;

            if (pwGapExt<0)
                pwGapExt = dnaGapExt;

            if (branchScalingFactor<0)
                branchScalingFactor = dnaBranchScalingFactor;

        }
        else
        {
            if (gapRate<0)
                gapRate = protGapRate;

            if (gapExt<0)
                gapExt = protGapExt;

            if (pwDist<0)
                pwDist = pwProtDist;

            if (pwGapRate<0)
                pwGapRate = pwProtGapRate;

            if (pwGapExt<0)
                pwGapExt = pwProtGapExt;

            if (branchScalingFactor<0)
                branchScalingFactor = protBranchScalingFactor;

        }

    }
    void checkStuff(map<string,string> *org_stuff,vector<string> *names,vector<string> *seqs)
    {
        int c=0;
        for(int i=0;i<names->size();i++)
        {
            string name = names->at(i);
            string tseq = seqs->at(i);

            string oseq = org_stuff->find(name)->second;

            for(string::iterator it = tseq.begin();it!=tseq.end();)
            {
                if(*it=='-')
                    tseq.erase(it);
                else
                    it++;
            }
            for(string::iterator it = oseq.begin();it!=oseq.end();)
            {
                if(*it=='-')
                    oseq.erase(it);
                else
                    it++;
            }


            if(oseq!=tseq)
            {
                cout<<"difference: "<<name<<"\n"<<tseq<<"\n"<<oseq<<"\n";
            }
        }
    }

};


#endif
