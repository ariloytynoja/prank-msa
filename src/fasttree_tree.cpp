/***************************************************************************
 *   Copyright (C) 2010-2014 by Ari Loytynoja                              *
 *   ari.loytynoja@gmail.com                                               *
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

#include "fasttree_tree.h"
#include <netdb.h>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

#if defined (__APPLE__)
#include <mach-o/dyld.h>
#endif


FastTree_tree::FastTree_tree()
{
}

bool FastTree_tree::test_executable()
{
    #if defined (__CYGWIN__)
    char path[200] = "";
    int length = readlink("/proc/self/exe",path,200-1);
    path[length] = '\0';

    string epath = string(path).substr(0,length);
    if (epath.find("/")!=std::string::npos)
        epath = epath.substr(0,epath.rfind("/")+1);
    progpath = epath;
    epath = epath+"fasttree >/dev/null 2>/dev/null";
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

    #else
    int length = readlink("/proc/self/exe",path,200-1);
    path[length] = '\0';
    epath = string(path).substr(0,length);
    if (epath.find("/")!=std::string::npos)
        epath = epath.substr(0,epath.rfind("/")+1);

    #endif

    char hostname[1024];
    hostname[1023] = '\0';
    gethostname(hostname, 1023);

    progpath = epath;
    epath = epath+"fasttree >/dev/null 2>/dev/null";
    int status = system(epath.c_str());
    if(WEXITSTATUS(status) == 0)
        return true;

    if(WEXITSTATUS(status) == 1 && strcmp(hostname, "wasabi2")==0)
        return true;

    progpath = "";
    status = system("fasttree >/dev/null 2>/dev/null");

    if(WEXITSTATUS(status) == 1 && strcmp(hostname, "wasabi2")==0)
        return true;

    return WEXITSTATUS(status) == 0;

    #endif
}

string FastTree_tree::infer_phylogeny(std::vector<string> *names,std::vector<string> *sequences,bool is_protein)
{

    stringstream f_name;

    int r = rand();
    while(true)
    {

        f_name <<tmp_dir<<"/d"<<r<<".fas";
        ifstream f_file(f_name.str().c_str());

        if(!f_file)
        {
            break;
        }
        r = rand();
    }

    ofstream f_output;
    f_output.open( f_name.str().c_str(), (ios::out) );

    vector<string>::iterator si = sequences->begin();
    vector<string>::iterator ni = names->begin();
    for(;si!=sequences->end();si++,ni++)
    {
        f_output<<">"<<*ni<<"\n"<<*si<<"\n";
    }
    f_output.close();

    stringstream command;
    command << progpath<<"fasttree -quiet -nopr -nosupport ";
    if(is_protein)
        command << f_name.str() << " 2>/dev/null";
    else
        command << "-nt "<<f_name.str() << " 2>/dev/null";

    if(NOISE>0)
        cout<<"cmd: "<<command.str()<<endl;

    FILE *fpipe;
    if ( !(fpipe = (FILE*)popen(command.str().c_str(),"r")) )
    {
        cout<<"Problems with fasttree pipe.\nExiting.\n";
        exit(1);
    }

    char line[256];
     string tree = "";
    while ( fgets( line, sizeof line, fpipe))
    {
        if(NOISE>1)
            cout<<"FastTree: "<<line;
        tree += line;
    }
    pclose(fpipe);

    if(NOISE>1)
        cout<<"FastTree_tree: "<<tree<<"\n";


    //delete files
    if( remove( f_name.str().c_str() ) != 0 )
      perror("Error deleting temporary file in FastTree_tree::infer_phylogeny");

    return tree;
}

