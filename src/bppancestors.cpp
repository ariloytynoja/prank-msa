/***************************************************************************
 *   Copyright (C) 2013 by Ari Loytynoja   *
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

#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>
#include <exception>
#include <unistd.h>
#include "bppancestors.h"
#include "readfile.h"
#include "readnewick.h"

#if defined (__APPLE__)
#include <mach-o/dyld.h>
#endif

using namespace std;

BppAncestors::BppAncestors()
{
}

bool BppAncestors::testExecutable()
{

    #if defined (__CYGWIN__)
    char path[200] = "";
    int length = readlink("/proc/self/exe",path,200-1);

    string epath = string(path).substr(0,length);
    if (epath.find("/")!=std::string::npos)
        epath = epath.substr(0,epath.rfind("/")+1);
    bppdistpath = epath;
    epath = epath+"bppancestor >/dev/null 2>/dev/null";
    int status = system(epath.c_str());

    return WEXITSTATUS(status) == 0;

    # else

    char path[200] = "";
    string epath;

    #if defined (__APPLE__)
    uint32_t size = sizeof(path);
    _NSGetExecutablePath(path, &size);
    epath = string(path);
    if (epath.find("/")!=std::string::npos)
        epath = epath.substr(0,epath.rfind("/")+1);
    //epath = "DYLD_LIBRARY_PATH="+epath+" "+epath;

    #else
    int length = readlink("/proc/self/exe",path,200-1);
    epath = string(path).substr(0,length);
    if (epath.find("/")!=std::string::npos)
        epath = epath.substr(0,epath.rfind("/")+1);

    #endif

    bppdistpath = epath;
    epath = epath+"bppancestor >/dev/null 2>/dev/null";
    int status = system(epath.c_str());

    if(WEXITSTATUS(status) == 0)
        return true;

    bppdistpath = "";
    status = system("bppancestor >/dev/null 2>/dev/null");

    return WEXITSTATUS(status) == 0;

    #endif
}

bool BppAncestors::inferAncestors(AncestralNode *root,map<string,string> *aseqs,string *atree,bool isDna)
{
    stringstream f_name;
    stringstream t_name;
    stringstream o_name;

    int r = rand();
    while(true)
    {

        f_name.str("");
        t_name.str("");
        o_name.str("");

        f_name <<tmp_dir<<"/f"<<r<<".fas";
        ifstream f_file(f_name.str().c_str());

        t_name <<tmp_dir<<"/t"<<r<<".tre";
        ifstream t_file(t_name.str().c_str());

        o_name <<tmp_dir<<"/o"<<r<<".fas";
        ifstream o_file(t_name.str().c_str());


        if(!f_file && !t_file && !o_file )
        {
            ofstream f_tmp;
            f_tmp.open(f_name.str().c_str(), (ios::out) );
            ofstream t_tmp;
            t_tmp.open(t_name.str().c_str(), (ios::out) );
            ofstream o_tmp;
            o_tmp.open(o_name.str().c_str(), (ios::out) );

            break;
        }
        r = rand();
    }



    ////////////

    int l = root->getSequence()->length();

    vector<string> names;
    vector<string> sequences;

    root->getTerminalNames(&names);
    for (int i=0; i<names.size(); i++)
    {
        sequences.push_back(string(""));
    }

    vector<string> col;

    bool tmpDOTS = DOTS;
    DOTS = false;
    for (int i=0; i<l; i++)
    {
        col.clear();
        root->getCharactersAt(&col,i);
        vector<string>::iterator cb = col.begin();
        vector<string>::iterator ce = col.end();
        vector<string>::iterator si = sequences.begin();

        for (; cb!=ce; cb++,si++)
        {
            *si+=*cb;
        }
    }
    DOTS = tmpDOTS;

    vector<string>::iterator si = sequences.begin();
    vector<string>::iterator ni = names.begin();

    ofstream f_output;
    f_output.open( f_name.str().c_str(), (ios::out) );
    int count = 0;
    map<string,string> tmp_names;
    for(;si!=sequences.end();si++,ni++)
    {
        string name1 = *ni+":";
        stringstream name2;
        name2<<"seq"<<count++;
        tmp_names.insert(make_pair(name1,name2.str()+":"));
        f_output<<">"<<name2.str()<<"\n"<<*si<<"\n";
    }
    f_output.close();


    string tree = "";

    int nodeNum = root->getTerminalNodeNumber();
    root->getNHXBrl(&tree,&nodeNum);

    for(map<string,string>::iterator it = tmp_names.begin();it != tmp_names.end(); it++)
    {
        size_t pos = 0;
        if((pos = tree.find("("+it->first)) != std::string::npos)
            tree.replace(pos+1, it->first.length(), it->second);
        if((pos = tree.find(","+it->first)) != std::string::npos)
            tree.replace(pos+1, it->first.length(), it->second);

    }

    stringstream tag;
    tag << root->getNodeName();
    char b,e; int num;
    tag >> b >> num >> e;
    tag.clear();
    tag.str("");
    tag << num;
    tree += "[&&NHX:ND="+tag.str()+"];";

    ofstream t_output;
    t_output.open( t_name.str().c_str(), (ios::out) );
    t_output<<tree<<endl;
    t_output.close();


    stringstream command;
    command << bppdistpath<<"bppancestor input.sequence.file="<<f_name.str()<<" input.sequence.format=Fasta input.sequence.sites_to_use=all input.tree.file="<<t_name.str()<<
            " input.tree.format=NHX input.sequence.max_gap_allowed=100% initFreqs=observed output.sequence.file="<<o_name.str()<<" output.sequence.format=Phylip";
    if(!isDna)
        command << " alphabet=Protein model=WAG01";
    else
    {
        if(CODON)
            command << " alphabet=Codon\\(letter=DNA,type=Standard\\) model=YN98\\(kappa=2,omega=0.5\\)";
        else
            command << " alphabet=DNA model=HKY85";
    }

    if(NOISE>0)
        cout<<"cmd: "<<command.str()<<endl;


    FILE *fpipe;
    char line[256];

    if ( !(fpipe = (FILE*)popen(command.str().c_str(),"r")) )
    {
        cout<<"Problems with bppancestor pipe.\nExiting.\n";
        exit(1);
    }
    while ( fgets( line, sizeof line, fpipe))
    {
        if(NOISE>1)
            cout<<"BppAncestor: "+string(line);
    }
    pclose(fpipe);

    ReadFile rf;
    int rv = rf.readBppPhylip(o_name.str().c_str());

    if(rv>0)
    {
        vector<string> s = rf.getSeqs();
        vector<string> n = rf.getNames();

        for(int i=0;i<n.size();i++)
            aseqs->insert(aseqs->begin(),pair<string,string>("#"+n.at(i)+"#",s.at(i)));
    }

    //delete files 
    if ( remove( f_name.str().c_str() ) != 0)
      perror("Error deleting temporary file in BppAncestors::inferAncestors");
    if ( remove( t_name.str().c_str()) != 0 )
      perror("Error deleting temporary file in BppAncestors::inferAncestors");
    if ( remove( o_name.str().c_str()) != 0 )
      perror("Error deleting temporary file in BppAncestors::inferAncestors");
    
    return (rv>0);
}
