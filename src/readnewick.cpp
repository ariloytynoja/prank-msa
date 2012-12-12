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
#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include "readnewick.h"
#include "ancestralnode.h"
#include "node.h"

using namespace std;

extern int NOISE;

ReadNewick::ReadNewick()
{
}

ReadNewick::~ReadNewick()
{
}

string ReadNewick::readFile(const char* filename)
{
    string t;
    s = "";
    ifstream in(filename);
    while (getline(in,t))
    {
        s += t;
    }

    return s;
}

void ReadNewick::buildTree(string s,map<string,TreeNode*>* nodes)
{
    string::iterator b = s.begin();
    string::iterator e = s.end();

    for (unsigned int i = 0; i < s.size(); i++)
    {
      if(s[i] == ' ')
      {
        s.erase(s.begin() + i);
        i--;
      }
    }


    int open=0;
    int end=0;
    int comma = 0;

    for (; b!=e; b++)
    {
        if ((*b)=='(')
            open++;
        if ((*b)==')')
            end++;
        if ((*b)==',')
            comma++;
    }

    if (open!=end)
    {
        cout<<"brackets do not match: "<<open<<" opening and "<<end<<" closing!"<<endl;
        exit(1);
    }

    // unrooted
    if (comma==open+1)
    {
        if (NOISE>0)
            cout<<"Unrooted tree, using midpoint rooting."<<endl;

        Node* n = new Node(s);
        s=n->rootedTree();

        b = s.begin();
        e = s.end();

        open=0;
        end=0;
        comma = 0;

        for (; b!=e; b++)
        {
            if ((*b)=='(')
                open++;
            if ((*b)==')')
                end++;
            if ((*b)==',')
                comma++;
        }

        if (NOISE>0)
            cout<<s<<endl;
        delete n;
    }
    else if (comma==open && comma==end)
    {

    }
    else
    {
        cout<<"Problem with the guidetee: brackets ("<<open<<","<<end<<") and commas ("<<comma<<") don't match)"<<endl<<endl;
        exit(-1);
    }

    int count = 1;
    string r;


    string n = "";

    do
    {
        r = "";
        b = s.begin();
        e = s.end();
        bool hasText = false;
        bool isOpen = false;

//        cout<<s<<endl;

        while (b!=e)
        {
            if ((*b)==' ' || (*b)=='\t' || (*b)=='\n')
            {
                b++;
                continue;
            }

            if ((*b)=='(')
            {
                isOpen = true;
                hasText = false;
                r += n+(*b);
                n = "";
                b++;
            }
            else if ((*b)==')')
            {
                if (hasText && isOpen)
                {
                    char tc[15];
                    sprintf(tc,"#%i#",count++);

                    AncestralNode *tn = new AncestralNode(n);
                    tn->setNodeName(tc);
                    if (tn->isLInternal())
                    {
                        tn->setLChild(nodes->find(tn->getLName())->second);
                        nodes->find(tn->getLName())->second->setBranchLength(tn->getLeftBrL());
                    }
                    if (tn->isRInternal())
                    {
                        tn->setRChild(nodes->find(tn->getRName())->second);
                        nodes->find(tn->getRName())->second->setBranchLength(tn->getRightBrL());
                    }

                    if(tn->LRealign){
                        tn->getLChild()->realignNode = true;
                        tn->realignNode = true;
                    }
                    if(tn->RRealign){
                        tn->getRChild()->realignNode = true;
                        tn->realignNode = true;
                    }

                    nodes->insert(make_pair(tc,tn));
                    r = r.substr(0,r.length()-1);
                    r += tc;
                    root = tc;
                    n = "";
                }
                else
                {
                    r += n+")";
                    n = "";
                }
                isOpen = false;
                hasText = false;
                b++;
            }
            else
            {
                hasText = true;
                n += (*b);
                b++;
            }
        }

        s = r;

        s+=n;
    }
    while (s.find(",")>0 && s.find(",")<s.length());

    if(n.find(":XN=realign") != string::npos)
    {
        n = n.substr(0,n.find(":XN=realign"))+n.substr(n.find(":XN=realign")+string(":XN=realign").length());
        nodes->find(root)->second->realignNode = true;
    }

    if(n.find("[&&NHX]") != string::npos)
    {
        n = n.substr(0,n.find("[&&NHX]"))+n.substr(n.find("[&&NHX]")+string("[&&NHX]").length());
    }

    nodes->find(root)->second->nhx_tag = n;
}

