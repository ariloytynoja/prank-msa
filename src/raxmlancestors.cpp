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
#include <algorithm>
#include <unistd.h>
#include "raxmlancestors.h"
#include "readfile.h"
#include "readnewick.h"

#if defined (__APPLE__)
#include <mach-o/dyld.h>
#endif

using namespace std;

RaxmlAncestors::RaxmlAncestors()
{
}

bool RaxmlAncestors::testExecutable()
{

    #if defined (__CYGWIN__)
    char path[200] = "";
    int length = readlink("/proc/self/exe",path,200-1);

    string epath = string(path).substr(0,length);
    if (epath.find("/")!=std::string::npos)
        epath = epath.substr(0,epath.rfind("/")+1);
    bppdistpath = epath;
    epath = epath+"raxml -h >/dev/null 2>/dev/null";
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

    raxmlpath = epath;
    epath = epath+"raxml -h >/dev/null 2>/dev/null";
    int status = system(epath.c_str());


    if(WEXITSTATUS(status) == 0)
        return true;

    raxmlpath = "";
    status = system("raxml -h >/dev/null 2>/dev/null");


    return WEXITSTATUS(status) == 0;

    #endif
}

bool RaxmlAncestors::inferAncestors(AncestralNode *root,map<string,string> *aseqs,string *atree,bool isDna)
{

    stringstream f_name;
    stringstream t_name;
    stringstream i_name;
    stringstream o_name;
    stringstream m_name;
    stringstream p_name;

    int r = rand();
    while(true)
    {

        f_name.str("");
        t_name.str("");
        i_name.str("");
        o_name.str("");
        m_name.str("");
        p_name.str("");

        f_name <<tmp_dir<<"/f"<<r<<".phy";
        ifstream f_file(f_name.str().c_str());

        t_name <<tmp_dir<<"/t"<<r<<".tre";
        ifstream t_file(t_name.str().c_str());

        i_name <<tmp_dir<<"/RAxML_info."<<r;
        ifstream i_file(i_name.str().c_str());

        o_name <<tmp_dir<<"/RAxML_marginalAncestralStates."<<r;
        ifstream o_file(o_name.str().c_str());

        m_name <<tmp_dir<<"/RAxML_nodeLabelledRootedTree."<<r;
        ifstream m_file(m_name.str().c_str());

        p_name <<tmp_dir<<"/RAxML_marginalAncestralProbabilities."<<r;
        ifstream p_file(p_name.str().c_str());

        if(!f_file && !t_file && !i_file && !o_file && !m_file && !p_file)
        {
            ofstream f_tmp;
            f_tmp.open(f_name.str().c_str(), (ios::out) );
            ofstream t_tmp;
            t_tmp.open(t_name.str().c_str(), (ios::out) );
            ofstream o_tmp;
            o_tmp.open(o_name.str().c_str(), (ios::out) );
            ofstream m_tmp;
            m_tmp.open(m_name.str().c_str(), (ios::out) );
            ofstream p_tmp;
            p_tmp.open(p_name.str().c_str(), (ios::out) );

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

    vector<int> masked_columns;

    for(int i=0;i<sequences.at(0).length();i++) {
        int ic=0;
        for(si = sequences.begin();si!=sequences.end();si++)
        {
            char c = si->at(i);
            if(c != 'N' && c != '-' && c != 'X' )
            {
                ic++;
                continue;
            }
        }
        if(ic==0) {
            for(si = sequences.begin();si!=sequences.end();si++)
            {
                char c = si->at(i);
                if(c == 'N')
                    si->at(i)='A';
                if(c == 'X')
                    si->at(i)='G';

                masked_columns.push_back(i);
            }
        }
    }

    ofstream f_output;    
    f_output.open( f_name.str().c_str(), (ios::out) );
    f_output<<sequences.size()<<" "<<sequences.at(0).length()<<"\n";
    int count = 0;
    map<string,string> tmp_names;
    for(si = sequences.begin();si!=sequences.end();si++,ni++)
    {
        string name1 = *ni+":";
        stringstream name2;
        name2<<"seq"<<count++;
        tmp_names.insert(make_pair(name1,name2.str()+":"));
        f_output<<name2.str()<<"\n"<<*si<<"\n";
    }
    f_output.close();


    string tree = "";

    root->getLabelledNewickBrl(&tree);
    tree.erase(std::remove(tree.begin(), tree.end(), '#'), tree.end());
    for(map<string,string>::iterator it = tmp_names.begin();it != tmp_names.end(); it++)
    {
        size_t pos = 0;
        if((pos = tree.find("("+it->first)) != std::string::npos)
        {
            tree.replace(pos+1, it->first.length(), it->second);
        }
        if((pos = tree.find(","+it->first)) != std::string::npos)
        {
            tree.replace(pos+1, it->first.length(), it->second);
        }
    }

    stringstream tag;
    tag << root->getNodeName();
    char b,e; int num;
    tag >> b >> num >> e;
    tag.clear();
    tag.str("");
    tag << num;
    tree += tag.str()+";";

    ofstream t_output;
    t_output.open( t_name.str().c_str(), (ios::out) );
    t_output<<tree<<endl;
    t_output.close();


    stringstream command;
    command << raxmlpath<<"raxml -s "<<f_name.str()<<" -t "<<t_name.str()<<" -w "+string(tmp_dir)+" -f A -T 2 -n "<<r;
    if(!isDna)
        command << " -m PROTCATJTT";
    else
        command << " -m GTRCAT";

    if(NOISE>0)
        cout<<"cmd: "<<command.str()<<endl;


    FILE *fpipe;
    char line[256];

    if ( !(fpipe = (FILE*)popen(command.str().c_str(),"r")) )
    {
        cout<<"Problems with raxml pipe.\nExiting.\n";
        exit(1);
    }
    while ( fgets( line, sizeof line, fpipe))
    {
        if(NOISE>1)
            cout<<"RAxML: "+string(line);
    }
    pclose(fpipe);

    o_name.str("");
    o_name <<tmp_dir<<"/RAxML_marginalAncestralStates."<<r;
    ifstream o_file(o_name.str().c_str());


    ReadNewick rn;
    string rtree = rn.readFile(m_name.str().c_str());

    if (rtree.length()==0)
    {

        //delete files
        if( remove( f_name.str().c_str() ) != 0 )
            perror("Error deleting temporary file in RaxmlAncestors::inferAncestors");
        if( remove( t_name.str().c_str() ) != 0 )
            perror("Error deleting temporary file in RaxmlAncestors::inferAncestors");
        if( remove( i_name.str().c_str() ) != 0 )
            perror("Error deleting temporary file in RaxmlAncestors::inferAncestors");
        if( remove( o_name.str().c_str() ) != 0 )
            perror("Error deleting temporary file in RaxmlAncestors::inferAncestors");
        if( remove( m_name.str().c_str() ) != 0 )
            perror("Error deleting temporary file in RaxmlAncestors::inferAncestors");
        if( remove( p_name.str().c_str() ) != 0 )
            perror("Error deleting temporary file in RaxmlAncestors::inferAncestors");

         stringstream r_name;
         r_name <<tmp_dir<<"/f"<<r<<".phy.reduced";
         ifstream r_file(r_name.str().c_str());
         if(r_file)
             if( remove( r_name.str().c_str() ) != 0 )
                 perror("Error deleting temporary file in RaxmlAncestors::inferAncestors");


        return false;
    }

    map<string,TreeNode*> rtnodes;
    rn.buildTree(rtree,&rtnodes);

    AncestralNode* rtroot = static_cast<AncestralNode*>(rtnodes[rn.getRoot()]);

    map<string,string> rtsubtrees;
    rtroot->getAllFullSubtreesWithNodename(&rtsubtrees,false);
    delete rtroot;


    map<string,TreeNode*> innodes;
    rn.buildTree(tree,&innodes);

    AncestralNode* inroot = static_cast<AncestralNode*>(innodes[rn.getRoot()]);

    map<string,string> subtrees;
    inroot->getAllFullSubtreesWithNodename(&subtrees,true);
    delete inroot;


    string temp, name, sequence = "";

    while (!o_file.eof())
    {
        getline(o_file, temp, '\n');  // Copy current line in temporary string
        stringstream tmpstr(temp);
        name="";sequence="";

        tmpstr >> name >> sequence;

        if(sequence.length()>0) {
            vector<int>::iterator mi;
            for(mi = masked_columns.begin();mi!=masked_columns.end();mi++) {
                for(si = sequences.begin();si!=sequences.end();si++)
                {
                    char c = sequence.at(*mi);
                    if(c == 'A')
                        sequence.at(*mi)='N';
                    if(c == 'G')
                        sequence.at(*mi)='X';
                }
            }
        }

        map<string,string>::iterator it = rtsubtrees.find("#"+name+"#");
        if(it!=rtsubtrees.end())
        {
            string subtree = it->second;
            it = subtrees.find(subtree);
            if(it!=subtrees.end())
            {
                std::replace(sequence.begin(),sequence.end(),'?','-');
                aseqs->insert(aseqs->begin(),pair<string,string>(it->second,sequence));
            }
        }
    }

    //delete files
    if( remove( f_name.str().c_str() ) != 0 )
        perror("Error deleting temporary file in RaxmlAncestors::inferAncestors");
    if( remove( t_name.str().c_str() ) != 0 )
        perror("Error deleting temporary file in RaxmlAncestors::inferAncestors");
    if( remove( i_name.str().c_str() ) != 0 )
        perror("Error deleting temporary file in RaxmlAncestors::inferAncestors");
    if( remove( o_name.str().c_str() ) != 0 )
        perror("Error deleting temporary file in RaxmlAncestors::inferAncestors");
    if( remove( m_name.str().c_str() ) != 0 )
        perror("Error deleting temporary file in RaxmlAncestors::inferAncestors");
    if( remove( p_name.str().c_str() ) != 0 )
        perror("Error deleting temporary file in RaxmlAncestors::inferAncestors");

    stringstream r_name;
    r_name <<tmp_dir<<"/f"<<r<<".phy.reduced";
    ifstream r_file(r_name.str().c_str());
    if(r_file)
        if( remove( r_name.str().c_str() ) != 0 )
            perror("Error deleting temporary file in RaxmlAncestors::inferAncestors");


    return aseqs->size()>0;
}
