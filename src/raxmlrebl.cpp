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
#include "raxmlrebl.h"
#include "readfile.h"
#include "readnewick.h"
#include "node.h"

#if defined (__APPLE__)
#include <mach-o/dyld.h>
#endif

using namespace std;

RaxmlRebl::RaxmlRebl()
{
}

bool RaxmlRebl::testExecutable()
{

    #if defined (__CYGWIN__)
    char path[200];
    int length = readlink("/proc/self/exe",path,200-1);

    string epath = string(path).substr(0,length);
    if (epath.find("/")!=std::string::npos)
        epath = epath.substr(0,epath.rfind("/")+1);
    bppdistpath = epath;
    epath = epath+"raxml -h >/dev/null 2>/dev/null";
    int status = system(epath.c_str());

    return WEXITSTATUS(status) == 0;

    # else

    char path[200];
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

bool RaxmlRebl::inferBranchLengths(AncestralNode *root,vector<string> *names,vector<string> *sequences,bool isDna)
{

    stringstream f_name;
    stringstream t_name;
    stringstream i_name;
    stringstream o_name;
    stringstream m_name;

    int r = rand();
    while(true)
    {

        f_name.str("");
        t_name.str("");
        i_name.str("");
        o_name.str("");
        m_name.str("");

        f_name <<tmp_dir<<"/f"<<r<<".phy";
        ifstream f_file(f_name.str().c_str());

        t_name <<tmp_dir<<"/t"<<r<<".tre";
        ifstream t_file(t_name.str().c_str());

        i_name <<tmp_dir<<"/RAxML_info."<<r;
        ifstream i_file(i_name.str().c_str());

        o_name <<tmp_dir<<"/RAxML_result."<<r;
        ifstream o_file(o_name.str().c_str());

        if(!f_file && !t_file && !i_file && !o_file)
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


    vector<string>::iterator si = sequences->begin();
    vector<string>::iterator ni = names->begin();

    ofstream f_output;    
    f_output.open( f_name.str().c_str(), (ios::out) );
    f_output<<sequences->size()<<" "<<sequences->at(0).length()<<"\n";
    int count = 0;
    map<string,string> tmp_names;
    map<string,string> rev_names;
    for(;si!=sequences->end();si++,ni++)
    {
        string name1 = *ni+":";
        stringstream name2;
        name2<<"seq"<<count++;
        tmp_names.insert(make_pair(name1,name2.str()+":"));
        rev_names.insert(make_pair(name2.str(),*ni));
        f_output<<name2.str()<<"\n"<<*si<<"\n";
//        cout<<name2.str()<<"\n"<<*si<<"\n";
    }
    f_output.close();


    string tree = "";
    root->getLabelledNewickBrl(&tree);
    tree.erase(std::remove(tree.begin(), tree.end(), '#'), tree.end());
    tree += ";";
    for(map<string,string>::iterator it = tmp_names.begin();it != tmp_names.end(); it++)
    {
        size_t pos = 0;
        if((pos = tree.find("("+it->first)) != std::string::npos)
            tree.replace(pos+1, it->first.length(), it->second);
        if((pos = tree.find(","+it->first)) != std::string::npos)
            tree.replace(pos+1, it->first.length(), it->second);
    }


    ofstream t_output;
    t_output.open( t_name.str().c_str(), (ios::out) );
    t_output<<tree<<endl;
    t_output.close();


    stringstream command;
    command << raxmlpath<<"raxml -s "<<f_name.str()<<" -t "<<t_name.str()<<" -w "+string(tmp_dir)+" -f e -T 2 -n "<<r;
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

    ReadNewick rn;
    string rtree = rn.readFile(o_name.str().c_str());


    Node* n1 = new Node(tree,false);
    Node* n2 = new Node(rtree,false);

    int h0,h1;
    string s = n2->markNodes();
    s = n1->markNodes();
    n1->getRootNamehashes(&h0,&h1);

    bool found = n2->rootAtMarkedNode(h0,h1);
    if(found)
    {
        string ntree = n2->rootedTree();

        Node* n3 = new Node(ntree,false);
        s = n3->markNodes();

        ReadNewick rn;
        map<string,TreeNode*> tmpnodes;
        rn.buildTree(tree,&tmpnodes);
        AncestralNode* tmproot = static_cast<AncestralNode*>(tmpnodes[rn.getRoot()]);

        map<int,string> ids;
        tmproot->markNodes();
        tmproot->getNamehashNames(&ids);

        map<int,float> map1;
        map<int,float> map3;

        n1->getNamehashesDistance(&map1);
        n3->getNamehashesDistance(&map3);

        map<string,float> brls;
        int nbrls = 0;

        for(map<int,float>::iterator it1=map1.begin();it1!=map1.end();it1++)
        {
            map<int,string>::iterator it=ids.find(it1->first);
            if(it != ids.end())
            {
                map<int,float>::iterator it3=map3.find(it1->first);
                if(it3 != map3.end())
                {
                    map<string,string>::iterator nit = rev_names.find(it->second);
                    if(nit != rev_names.end() )
                    {
                        brls.insert(make_pair(nit->second,it3->second));
                        nbrls++;
//                        cout<<nit->second<<" "<<it1->first<<" "<<it1->second<<" "<<it3->second<<endl;
                    }
                    else if(it->second.at(0)=='#')
                    {
                        brls.insert(make_pair(it->second,it3->second));
                        nbrls++;
//                        cout<<it->second<<" "<<it1->first<<" "<<it1->second<<" "<<it3->second<<endl;
                    }
                }
            }

        }

        delete tmproot;
        delete n3;

        if(nbrls == n1->getNodeNumber())
        {
            int nubrl = 0;
            root->updateBranchLengths(&brls,&nubrl);
            if(nubrl != nbrls)
                cout<<nubrl<<" Updating branch lengths failed!\n\n";
        }
    }

    delete n1;
    delete n2;

    //delete files
    if( remove( f_name.str().c_str() ) != 0 )
        perror("Error deleting temporary file in RaxmlRebl::inferAncestors");
    if( remove( t_name.str().c_str() ) != 0 )
        perror("Error deleting temporary file in RaxmlRebl::inferAncestors");
    if( remove( i_name.str().c_str() ) != 0 )
        perror("Error deleting temporary file in RaxmlRebl::inferAncestors");
    if( remove( o_name.str().c_str() ) != 0 )
        perror("Error deleting temporary file in RaxmlRebl::inferAncestors");
    if( remove( m_name.str().c_str() ) != 0 )
        perror("Error deleting temporary file in RaxmlRebl::inferAncestors");


    return true;
}
