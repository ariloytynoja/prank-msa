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
#include "readannotation.h"
#include "hirschberg.h"
#include "readalignment.h"
#include "node.h"

using namespace std;

ProgressiveAlignment::~ProgressiveAlignment()
{

}

ProgressiveAlignment::ProgressiveAlignment(string treefile,string seqfile,string dnafile,int aMethod)
{

    if (NOISE>=0) {
      if (CONVERT) {
        cout<<endl<<"PRANK: converting '"<<seqfile<<"' to '"<<outfile<<this->formatExtension(format)<<"'."<<endl<<endl;
      }
      else
      {
          if (BACKTRANSLATE) {
              cout<<endl<<"PRANK: creating a DNA alignment '"<<outfile<<this->formatExtension(format)
                      <<"' based on a protein alignment '"<<seqfile<<"' and DNA sequences in '"<<dnafile<<"'."<<endl<<endl;
          } else {
            cout<<endl<<"PRANK: aligning sequences in '"<<seqfile<<"', writing results to '"<<outfile<<".*' ("<<this->formatExtension(format)<<"/.xml/.dnd).\n";
            if (treefile!="" & hmmname!="")
              cout<<"Using alignment guidetree '"<<treefile<<"'' and model '"<<hmmname<<"'.\n";
            else if (treefile!="")
              cout<<"Using alignment guidetree '"<<treefile<<"'.\n";
            else if (hmmname!="")
              cout<<"Using alignment model '"<<hmmname<<"'.\n";
            cout<<endl;
          }
      }
    }

    vector<string> names;
    vector<string> sequences;
    bool isDna;


    if(BACKTRANSLATE){
      ReadFile* rfa = new ReadFile();
      int ns = rfa->readFile(dnafile.c_str());
      if (ns>0) {
        vector<string> dnaNames = rfa->getNames();
        vector<string> dnaSequences = rfa->getSeqs();
        isDna = rfa->dnaSeqs();
        if(!isDna){
          cout<<"Not DNA in "<<dnafile<<". Exiting."<<endl;
          exit(-1);
        }
        trseq = new TranslateSequences();
        if(!trseq->translateProtein(&dnaNames,&dnaSequences)){
          cout<<"Translation failed. Exiting."<<endl<<endl;
          delete trseq;
        }
      }
      delete rfa;

      rfa = new ReadFile();
      ns = rfa->readFile(seqfile.c_str());
      if (ns>0) {
        vector<string> protNames = rfa->getNames();
        vector<string> protSequences = rfa->getSeqs();
        isDna = rfa->dnaSeqs();
        if(isDna){
          cout<<"Not protein in "<<dnafile<<". Exiting."<<endl;
          exit(-1);
        }

        vector<string> codons;
        if(!trseq->translateDNA(&protNames,&protSequences,&codons)){
          cout<<"Backtranslation failed. Exiting."<<endl;
          exit(-1);
        }

        ofstream seqout((outfile+".fas").c_str());
        vector<string>::iterator ir;
        vector<string>::iterator nr = protNames.begin();
        for(ir=codons.begin();ir!=codons.end();ir++,nr++){
          seqout<<">"<<*nr<<endl;
          string seq = *ir;
          for(unsigned int i=0;i<seq.length();i+=60){
            seqout<<seq.substr(i,60)<<endl;
          }
        }
        exit(0);

      }
      exit(0);
    }

    // Get the sequence data
    //
    ReadFile* rfa = new ReadFile();
    int ns = rfa->readFile(seqfile.c_str());
    if (ns>=2) {

        names = rfa->getNames();
        sequences = rfa->getSeqs();
        isDna = rfa->dnaSeqs();
        if(DNA && !isDna)
        {
            cout<<"Warning autodetection suggests protein but DNA model forced.\n";
            isDna = true;
        }
        if(CODON && !isDna)
        {
            cout<<"Warning autodetection suggests protein but codon model forced.\n";
            isDna = true;
        }
        if(PROTEIN && isDna)
        {
            cout<<"Warning autodetection suggests DNA but protein model forced.\n";
            isDna = false;
        }

        list<string> sorted_names;
        sorted_names.assign(names.begin(),names.end());
        sorted_names.sort();

        string pname = "";
        bool unique = true;

        if(SHORTNAMES) {
          for(vector<string>::iterator lit = names.begin();lit!=names.end();lit++){

            string n = *lit;
            n = n.substr(0,n.find_first_of(' '));
            *lit = n;
          }
        }

        for(list<string>::iterator lit = sorted_names.begin();lit!=sorted_names.end();lit++){

          if(pname==*lit){
            cout<<"Sequence name "<<pname<<" appears more than once!"<<endl;
            unique = false;
          }
          pname=*lit;
        }
        if(!unique){
          cout<<"Sequence names not unique. Exiting."<<endl<<endl;
          exit(-1);
        }
    } else if (ns==1) {

        cout<<"Only one sequence in "<<seqfile<<"!"<<endl;
        names = rfa->getNames();
        sequences = rfa->getSeqs();

        WriteFile* wfa = new WriteFile();
        string file = outfile+formatExtension(format);
        wfa->writeSeqs(file.c_str(),&names,&sequences,format);
        delete wfa;
        exit(0);

    } else {
        cout<<"No sequences found in "<<seqfile<<"!"<<endl;
        exit(-1);
    }


    if (CONVERT && !TRANSLATE) {
        WriteFile* wfa = new WriteFile();
        string file = outfile+formatExtension(format);
        wfa->writeSeqs(file.c_str(),&names,&sequences,format);
        delete wfa;
        exit(0);
    }

    bool warning = false;
    vector<string>::iterator tit = names.begin();
    for(;tit!=names.end();tit++){
      string ts = *tit;
      string us = "_";
      for(unsigned int i=0;i<ts.length();i++) {
        if(ts.at(i)==' ' || ts.at(i)==':' || ts.at(i)=='\t' || ts.at(i)=='\n'
           || ts.at(i)==')' || ts.at(i)=='(' || ts.at(i)==','){
          ts.replace(i,1,us);
          warning = true;
        }
      }
      *tit = ts;
    }

    if(warning){
        cout<<"Warning: sequence names changed."<<endl<<endl;
    }


    if(isDna && TRANSLATE && CONVERT){
      trseq = new TranslateSequences();
      if(trseq->translateProtein(&names,&sequences)){
        isDna=false;
      } else {
        cout<<"Exiting."<<endl<<endl;
        delete trseq;
        exit(-1);
      }

      WriteFile* wfa = new WriteFile();
      string file = outfile+formatExtension(format);
      wfa->writeSeqs(file.c_str(),&names,&sequences,format);
      delete wfa;
      exit(0);

    }


    if(isDna && TRANSLATE){
      trseq = new TranslateSequences();
      if(trseq->translateProtein(&names,&sequences)){
        isDna=false;
      } else {
        cout<<"Exiting."<<endl<<endl;
        delete trseq;
        exit(-1);
      }
    }


    if (CONVERT) {
      unsigned int outputlength = sequences.at(0).length();
      vector<string>::iterator ir;
      for(ir=sequences.begin();ir!=sequences.end();ir++){
        if(ir->length()>outputlength)
          outputlength = ir->length();
      }
      WriteFile* wfa = new WriteFile();
      string file = outfile+formatExtension(format);
      wfa->writeSeqs(file.c_str(),&names,&sequences,format);
      delete wfa;
      exit(0);
    }


    this->makeSettings(isDna);


    // Get the guidetree -- or generate one
    //
    ReadNewick *rn = new ReadNewick();
    string tree;

    if (treefile==""){

        IntMatrix *substScores;

        if (isDna) {

            float* freqs = new float[4];
            freqs[0]=freqs[1]=freqs[2]=freqs[3]=0;

            rfa->countDnaFreqs(freqs);
            float total=freqs[0]+freqs[1]+freqs[2]+freqs[3];
            freqs[0]/=total;freqs[1]/=total;freqs[2]/=total;freqs[3]/=total;

            hmm = new HMModel;
            hmm->dnaModel(freqs,rfa->isRna);

            substScores = new IntMatrix(5,5,"pwSubstScores");
            hmm->pairwiseModel(substScores,pwDist);

            delete[] freqs;
            delete hmm;

        } else {

            hmm = new HMModel;
            hmm->proteinModel();

            substScores = new IntMatrix(21,21,"pwSubstScores");
            hmm->pairwiseModel(substScores,pwDist);

            delete hmm;
        }

        int time1 = time(0);

        GuideTree* gt;
        if(TREEFROMALIGNMENT)
            gt = new GuideTree(&sequences,&names,isDna);
        else
            gt = new GuideTree(&sequences,&names,substScores);

            if (NOISE>0)
                    cout<<"GuideTree; time "<<(time(0)-time1)<<"s"<<endl;

        tree = gt->getTree();
        delete gt;
        delete substScores;

        if (NOISE>0)
            cout<<tree<<endl;

    } else {
        if(PWGENOMIC){
          if (ns!=2) {
            cout<<"Pairwise alignment only valid for two sequences!"<<endl;
            exit(-1);
          }
          stringstream treest;
          treest << "("+names.at(0)<<":"<<pwgendist/2.0<<","<<names.at(1)<<":"<<pwgendist/2.0<<");";
          tree = treest.str();
          cout<<tree<<endl;
        } else {
          if(TREESTRING){
           tree = treefile;
          } else {
            tree = rn->readFile(treefile.c_str());
          }
          if (tree==""){
              cout<<"No tree found in "<<treefile<<"!"<<endl;
              exit(-1);
          }
        }
    }

    // Check that sequences match the model
    if (HASHMM){

        if (isDna && hmm->getAlphabet().length()==20) {
            cout<<"Sequences in "<<seqfile<<" seem to be DNA but the alphabet in "<<hmmname<<" is for protein!"<<endl;
            exit(-1);
        } else if (!isDna && hmm->getAlphabet().length()==4) {
            cout<<"Sequences in "<<seqfile<<" seem to be protein but the alphabet in "<<hmmname<<" is for DNA!"<<endl;
            exit(-1);
        }

        // Or make a new model
    } else {

        // protein model
        if (!isDna){
            hmm = new HMModel;
            hmm->proteinModel();

        } else {

            // codon model
            if (CODON) {
                hmm = new HMModel;
                hmm->codonModel();

                // dna model
            } else {

                // get nucleotide frequencies - either empirical (NOTE! all seqs used!) or user-defined
                float* freqs = new float[4];
                freqs[0]=freqs[1]=freqs[2]=freqs[3]=0;

                // user-defined nucleotide frequencies
                if (dnaFreqs!=""){
                    rfa->countDnaFreqs(freqs);
                    int i=0;
                    int j=0;
                    while (dnaFreqs.find(",",i)<=dnaFreqs.length()){
                        freqs[j++]=atof(dnaFreqs.substr(i,dnaFreqs.find(",",i)-i).c_str());
                        i=dnaFreqs.find(",",i)+1;
                    }
                    freqs[3]=atof(dnaFreqs.substr(i).c_str());
                    float total=freqs[0]+freqs[1]+freqs[2]+freqs[3];
                    freqs[0]/=total;freqs[1]/=total;freqs[2]/=total;freqs[3]/=total;

                    // empirical nucleotide frequencies
                } else {
                    rfa->countDnaFreqs(freqs);
                    float total=freqs[0]+freqs[1]+freqs[2]+freqs[3];
                    freqs[0]/=total;freqs[1]/=total;freqs[2]/=total;freqs[3]/=total;
                }

                //cout<<"freqs "<<freqs[0]<<" "<<freqs[1]<<" "<<freqs[2]<<" "<<freqs[3]<<endl;
                // make dna model
                hmm = new HMModel;
                hmm->dnaModel(freqs,rfa->isRna);
                delete[] freqs;
            }
        }
    }

    delete rfa;

    vector<string>::iterator si = sequences.begin();
    unsigned int longest = 0;
    unsigned int slongest = 0;
    for (;si!=sequences.end();si++){
            if (si->length()>longest) {
                    slongest = longest;
                    longest = si->length();
            } else if(si->length()>slongest) {
                    slongest = si->length();
            }
    }

    Site *sites = new Site();
    sites->setASize(hmm->getASize());
    sites->setNState(hmm->getNStates());
    sites->setMatrices(longest,slongest);

    Node* n = new Node(tree);
    n->mark_sequences(&names);
    int treeleaves = 0;
    int treematches = 0;
    n->countMatchingLeaves(&treeleaves,&treematches);

    if(treeleaves!=treematches)
    {
        if(PRUNETREE)
        {
            n->prune_tree();
            tree=n->print_tree();
        } else {
            cout<<"The guide tree has "<<treeleaves<<" leaves but only "<<treematches<<" match with the names in the sequence data.\n"
                    "The tree can be pruned to match the data using the flag '-prunetree'. Now exiting.\n\n";
            exit(-1);
        }
    }

    map<string,TreeNode*> nodes;
    rn->buildTree(tree,&nodes);

    AncestralNode* root = static_cast<AncestralNode*>(nodes[rn->getRoot()]);

    if(PREALIGNED) {
      si = sequences.begin();
      for (;si!=sequences.end();si++){
        string ts = *si;
        string us = "-";
        for(unsigned int i=0;i<ts.length();i++) {
          if(ts.at(i)=='.'){
            ts.replace(i,1,us);
          }
        }
        *si = ts;
      }
    }

    int nsqs = 0;
    root->setCharString(&names,&sequences,&nsqs);


    if(PREALIGNED) {

      if(TREEFROMALIGNMENT) {
        string tmpTree = "";
        root->getCleanNewick(&tmpTree);

        ofstream seqout((outfile+".0.dnd").c_str());
        seqout<<tmpTree<<endl;
        seqout.close();
      }

      ReadAlignment ra;
      ra.initialiseMatrices(longest+2);

      if (NOISE>=0)
        cout<<"Reading multiple alignment."<<endl;

      root->setTotalNodes();

      root->readAlignment();

      sites->deleteMatrices();
      delete sites;

      nodes.clear();
      delete rn;

      printAlignment(root,&names,&sequences,0,isDna);

      delete root;
      ra.cleanUp();

       return;

    } else if(PARTLYALIGNED) {
      ReadAlignment ra;
      ra.initialiseMatrices(longest+2);

      Hirschberg hir;
      hir.initialiseMatrices((int)(((float)longest+2)*initialMatrixSize));

      if (NOISE>=0)
        cout<<"Finishing partially aligned alignment."<<endl;

      root->setTotalNodes();

      root->partlyAlignSequences();

      sites->deleteMatrices();
      delete sites;

      nodes.clear();
      delete rn;

      printAlignment(root,&names,&sequences,0,isDna);

      delete root;
      ra.cleanUp();
      hir.cleanUp();
       return;

    } else {

      Hirschberg hir;
      if (ANCHORS) {
          hir.initialiseMatrices(initialAnchDist);
      }
      else {
          hir.initialiseMatrices((int)(((float)longest+2)*initialMatrixSize));
      }


      map<string,FlMatrix*> annotation;
      if (PRIORS){
          ReadAnnotation* ra = new ReadAnnotation(annofile,&annotation);
          root->setAnnotation(&annotation);
          delete ra;
      }

      if (nsqs!=root->getTerminalNodeNumber()){
          cout<<"Names in sequence file "<<seqfile<<" and guidetree "<<treefile<<" do not match!"<<endl;
          exit(-1);
      }

      root->setTotalNodes();

      string tmpTree = "";
      root->getCleanNewick(&tmpTree);

      if(PRINTTREE) {
        ofstream seqout((outfile+".1.dnd").c_str());
        seqout<<tmpTree<<endl;
        seqout.close();
      }

      if (NOISE>=0)
          cout<<"Generating multiple alignment."<<endl;

      bool saveDOPOST = DOPOST;
      if (TWICE) {
        DOPOST = false;
      }

      // HirschbergAlignment
      if (aMethod==0) {
          root->alignSequences(0);
      } else {
          cout<<"wrong method"<<endl;
          exit(-1);
      }

      if (PRIORS) {
          map<string,FlMatrix*>::iterator it;
          for ( it=annotation.begin();it!=annotation.end(); it++ ) {
              delete it->second;
          }
      }

      int iteration = 1;

      if (TWICE) {

          if (NOISE>=0)
              cout<<endl;

          printAlignment(root,&names,&sequences,iteration,isDna);
          iteration++;

          int l = root->getSequence()->length();

          names.clear();
          root->getTerminalNames(&names);

          vector<string>::iterator si = sequences.begin();
          for (;si!=sequences.end();si++) {
              si->clear();
          }

          vector<string> col;

          for (int i=0;i<l;i++) {
              col.clear();
              root->getCharactersAt(&col,i);
              vector<string>::iterator cb = col.begin();
              vector<string>::iterator ce = col.end();

              si = sequences.begin();
              for (;cb!=ce;cb++,si++) {
                  *si+=*cb;
              }
          }

          GuideTree* gt = new GuideTree(&sequences,&names,isDna);
          tree = gt->getTree();
          if (NOISE>0)
              cout<<tree<<endl;
          delete gt;

          delete root;
          nodes.clear();

          rn->buildTree(tree,&nodes);

          if(CODON)
            l*=3;

          si = sequences.begin();
          for (;si!=sequences.end();si++) {
              string s = "";
              for (int i=0;i<(int)si->length();i++) {
                  char c = si->at(i);
                  if (c!='-') {
                      s+=c;
                  }
              }
              *si = s;
          }

          root = static_cast<AncestralNode*>(nodes[rn->getRoot()]);
          nsqs = 0;
          root->setCharString(&names,&sequences,&nsqs);

          root->setTotalNodes();

          string tmpTree = "";
          root->getCleanNewick(&tmpTree);
          if(PRINTTREE) {
            ofstream seqout((outfile+".2.dnd").c_str());
            seqout<<tmpTree<<endl;
            seqout.close();
          }
          DOPOST = saveDOPOST;

          if (NOISE>=0)
              cout<<"Generating improved multiple alignment."<<endl;
          root->alignSequences(0);

      }

      hir.cleanUp();

      printAlignment(root,&names,&sequences,iteration,isDna);
    }

    sites->deleteMatrices();
    delete sites;

    nodes.clear();
    delete rn;

    delete root;

    if(TRANSLATE){
      delete trseq;
    }
}

void ProgressiveAlignment::printAlignment(AncestralNode *root,vector<string> *nms,vector<string> *seqs,int iteration, bool isDna)
{
    // done - write output
    //
//    int n = root->getTerminalNodeNumber();
    int l = root->getSequence()->length();

    nms->clear();
    root->getTerminalNames(nms);

    vector<string>::iterator si = seqs->begin();
    for (;si!=seqs->end();si++) {
        si->clear();
    }

    vector<string> col;

    for (int i=0;i<l;i++) {

        col.clear();
        root->getCharactersAt(&col,i);
        vector<string>::iterator cb = col.begin();
        vector<string>::iterator ce = col.end();

        si = seqs->begin();
        for (;cb!=ce;cb++,si++) {
            *si+=*cb;
        }
    }

    if (CODON)
        l*=3;

    if(!TRANSLATE) {

      WriteFile* wfa = new WriteFile();
      string file = outfile+"."+itos(iteration)+formatExtension(format);
      wfa->writeSeqs(file.c_str(),nms,seqs,format,isDna,root,false);
      delete wfa;

      printXml(root,iteration,false);

    } else {

      WriteFile* wfa = new WriteFile();
      string file = outfile+".pep."+itos(iteration)+formatExtension(format);
      wfa->writeSeqs(file.c_str(),nms,seqs,format,false,root,false);
      delete wfa;

      printXml(root,iteration,false);

      vector<string> dnaSeqs;
      trseq->translateDNA(nms,seqs,&dnaSeqs);

      wfa = new WriteFile();
      file = outfile+".nuc."+itos(iteration)+formatExtension(format);
      wfa->writeSeqs(file.c_str(),nms,&dnaSeqs,format,true,root,true);
      delete wfa;

      printXml(root,iteration,true);
    }

    printAncestral(root,nms,seqs,iteration);

}


void ProgressiveAlignment::printXml(AncestralNode *root,int iteration,bool translate) {

    if (!WRITEXML)
        return;

    int n = root->getTerminalNodeNumber();
    int l = root->getSequence()->length();
    int nState = hmm->getNStates();

    char* alignment;
    if (CODON || translate){
        alignment = new char[n*l*3];
    } else {
        alignment = new char[n*l];
    }

    if(!translate) {

      vector<string> col;
      for (int i=0;i<l;i++) {
          col.clear();
          root->getCharactersAt(&col,i);
          vector<string>::iterator cb = col.begin();
          vector<string>::iterator ce = col.end();

          int j=0;
          for (;cb!=ce;cb++) {
              if (CODON){
                  alignment[j*l*3+i*3] = cb->at(0);
                  alignment[j*l*3+i*3+1] = cb->at(1);
                  alignment[j*l*3+i*3+2] = cb->at(2);
              } else {
                  alignment[j*l+i] = cb->at(0);
              }
              j++;
          }
      }
    } else {

      char *tmp = new char[n*l];

      vector<string> col;
      for (int i=0;i<l;i++) {
        col.clear();
        root->getCharactersAt(&col,i);
        vector<string>::iterator cb = col.begin();
        vector<string>::iterator ce = col.end();

        int j=0;
        for (;cb!=ce;cb++,j++) {
            tmp[j*l+i] = cb->at(0);
        }
      }

      vector<string> names;
      root->getTerminalNames(&names);
      vector<string>::iterator si = names.begin();
      vector<string> prot;

      for (int j=0;j<n;j++) {
        string seq = "";
        for (int i=0;i<l;i++) {
          seq+=tmp[j*l+i];
        }
        prot.push_back(seq);
      }

      si = prot.begin();

      vector<string> dna;

      if(!trseq->translateDNA(&names,&prot,&dna)){
        cout<<"Backtranslation failed. Exiting."<<endl;
        exit(-1);
      }

      si = dna.begin();
      for (int j=0;j<n;j++) {
        for (int i=0;i<l*3;i++) {
          alignment[j*l*3+i] = si->at(i);
        }
        si++;
      }
    }


    string type = "";

    if(TRANSLATE && !translate)
      type="pep.";
    else if(TRANSLATE && translate)
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
    for (int j=0;j<n;j++) {
        seqout<<"<leaf id=\"seq"<<j+1<<"\" name=\""<<(*si++)<<"\">"<<endl;
        seqout<<"  <sequence>"<<endl<<"    ";
        if (CODON || translate){
          for (int i=0;i<l*3;i++) {
            seqout<<alignment[j*l*3+i];
          }
        } else {
          for (int i=0;i<l;i++) {
              seqout<<alignment[j*l+i];
          }
        }
        seqout<<endl;
        seqout<<"  </sequence>"<<endl<<"</leaf>"<<endl;
    }
    nms.clear();

    // internal nodes
    root->setSiteLength(l);
    for (int i=0;i<l;i++) {
        root->setSiteIndex(i,i);
    }
    root->outputXml(&seqout,translate);
    seqout<<"</nodes>"<<endl<<"<model>"<<endl;

    // model
    for (int k=0;k<nState;k++) {
        seqout<<"  <probability id=\""<<k+1<<"\" name=\""<<hmm->getStName(k)<<"\" ";
        seqout<<"color=\""<<hmm->getDrawCl(k)<<"\" style=\""<<hmm->getDrawPt(k)<<"\" ";
        seqout<<"offset=\""<<hmm->getDrawOf(k)<<"\"";
        if(nState>1){
          seqout<<" show=\"yes\"/>"<<endl;
        } else {
          seqout<<" show=\"no\"/>"<<endl;
        }
    }
    if (DOPOST)
        seqout<<"  <probability id=\""<<nState+1<<"\" name=\"postprob\" color=\"gray\" style=\"bar\" show=\"yes\"/>"<<endl;

    seqout<<"</model>"<<endl<<"</ms_alignment>"<<endl;

    delete []alignment;

}

void ProgressiveAlignment::printAncestral(AncestralNode *root,vector<string> *nms,vector<string> *sqs,int iteration)
{
    if (!WRITEANC)
        return;

    string tree = "";
    root->getNewick(&tree);

    vector<string> anms;
    root->getInternalNames(&anms);

    ofstream ancSeq((outfile+"."+itos(iteration)+".ancseq").c_str());

    ancSeq<<"# "<<tree<<endl;

    vector<string>::iterator si = sqs->begin();
    vector<string>::iterator ni = nms->begin();
    for (;ni!=nms->end();si++,ni++){
        ancSeq<<">"<<*ni<<endl;
        ancSeq<<*si<<endl;
    }


    int n = root->getInternalNodeNumber();
    int l = root->getSequence()->length();

    char* alignment;
    if (CODON){
        alignment = new char[n*l*3];
    } else {
        alignment = new char[n*l];
    }

    vector<string> col;

    int i;
    FOR(i,l) {

        col.clear();
        root->getAncCharactersAt(&col,i,0);
        vector<string>::iterator cb = col.begin();
        vector<string>::iterator ce = col.end();

        int j=0;
        for (;cb!=ce;cb++) {
            if (CODON){
                alignment[j*l*3+i*3] = cb->at(0);
                alignment[j*l*3+i*3+1] = cb->at(1);
                alignment[j*l*3+i*3+2] = cb->at(2);
            } else {
                alignment[j*l+i] = cb->at(0);
            }
            j++;
        }
    }

    ni = anms.begin();
    int j=0;
    for (;ni!=anms.end();j++,ni++){
        ancSeq<<">"<<*ni<<endl;
        for (int i=0;i<l;i++) {
            ancSeq<<alignment[j*l+i];
        }
        ancSeq<<endl;
    }

    delete []alignment;
    ancSeq.close();


    FILE *ancPro = fopen((outfile+"."+itos(iteration)+".ancprof").c_str(),"w");
    fclose(ancPro);

    int *insSite = new int[l];
    FOR(i,l){
        insSite[i]=0;
    }
    root->writeAncCharacters(insSite,iteration);

    delete []insSite;

    return;
}

string ProgressiveAlignment::formatExtension(int format)
{
  if(format==1){
    return ".igs";
  } else if(format==2){
    return ".gen";
  } else if(format==3){
    return ".nbr";
  } else if(format==4){
    return ".emb";
  } else if(format==6){
    return ".dst";
  } else if(format==7){
    return ".fch";
  } else if(format==8){
    return ".fas";
  } else if(format==11 || format==12 || format==18){
    return ".phy";
  } else if(format==14){
    return ".pir";
  } else if(format==15){
    return ".msf";
  } else if(format==17){
    return ".nex";
  } else {
    return "";
  }
}
