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

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include "node.h"

using namespace std;

std::string mpTree;
float halfLength;
float maxSpan;
extern float defaultBranchLength;

Node::~Node()
{

//    if (!isLast){
//        if (child0!=0){
//            delete child0;
//        }
//        if (child1!=0){
//            delete child1;
//        }
//    }

    if (node_has_left_child)
    {
        delete child0;
        node_has_left_child = false;
    }
    if (node_has_right_child)
    {
        delete child1;
        node_has_right_child = false;
    }

//    cout<<"delete "<<this->get_name()<<endl;
}

int Node::count = 1;

Node::Node(string t)
{
    mpTree = "";
    maxLength = 0.0;
    maxSpan = 0.0;
    isLast = false;
    isFirst = true;

    this->count = 1;
    stringstream ss;
    ss<<Node::count;
    this->set_name(ss.str());
    Node::count++;

    node_has_sequence = false;
    node_has_left_child = true;
    node_has_right_child = true;

    has_missing_branch_lengths = false;

    tree = t;

    for (unsigned int i = 0; i < tree.size(); i++)
    {
      if(tree[i] == ' ')
      {
        tree.erase(tree.begin() + i);
        i--;
      }
    }

    while(tree.find("[&&NHX:") != string::npos)
    {
        int start = tree.find("[&&NHX:");
        int stop = tree.find("]",start);
        tree.erase(start,stop-start+1);
    }

    divideTree(tree,subTrees,subDistances);

    subDistances[0] = abs(subDistances[0]);
    subDistances[1] = abs(subDistances[1]);

    float tot = subDistances[0]+subDistances[1];

    char num[10];
    sprintf(num,"%.5f",tot);

    revTrees[0] = subTrees[1]+":"+num;
    revTrees[1] = subTrees[0]+":"+num;

    child0 = new Node(subTrees[0],this,0);
    child1 = new Node(subTrees[1],this,1);


    float currPair = subDistances[0]+child0->maxLength+subDistances[1]+child1->maxLength;
    if (currPair > maxSpan)
    {
        maxSpan = currPair;
    }


    findMiddlePoint();

}

string Node::rootedTree()
{
    if(has_missing_branch_lengths)
    {
        cout<<"The guide tree with missing branch lengths should be rooted. Exiting."<<endl<<endl;
        exit(0);
    }
    return mpTree;
}

void Node::checkBifurcation(string *t)
{

    string::iterator b = t->begin();
    string::iterator e = t->end();

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

    if (comma!=open || comma!=end)
    {
            cout<<"Correcting (arbitrarily) for multifurcating nodes."<<endl;

            if(open==end+1 && open == comma)
            {
                t->append(":0.0)");
            }
    }
}

Node::Node(string t,Node* p,int branch)
{
//    cout<<t<<endl;

    tree = t;
    parent = p;
    maxLength = 0;
    isLast = true;
    isFirst = false;

    node_has_sequence = false;
    node_has_left_child = false;
    node_has_right_child = false;
    this->set_name(tree);

    if (tree.find(",",0)>0 && tree.find(",",0)<tree.length())
    {
        isLast = false;

        node_has_left_child = true;
        node_has_right_child = true;

        stringstream ss;
        ss<<Node::count;
        this->set_name(ss.str());
        Node::count++;

        divideTree(tree,subTrees,subDistances);

        subDistances[0] = abs(subDistances[0]);
        subDistances[1] = abs(subDistances[1]);

        char num0[10];
        sprintf(num0,"%.5f",subDistances[0]);
        char num1[10];
        sprintf(num1,"%.5f",subDistances[1]);

        revTrees[0] = "("+parent->revTrees[branch]+","+subTrees[1]+":"+num1+"):"+num0;
        revTrees[1] = "("+parent->revTrees[branch]+","+subTrees[0]+":"+num0+"):"+num1;

        child0 = new Node(subTrees[0],this,0);
        child1 = new Node(subTrees[1],this,1);


        float currPair = subDistances[0]+child0->maxLength+subDistances[1]+child1->maxLength;
        if (currPair > maxSpan)
        {
            maxSpan = currPair;
        }

        if (subDistances[0]+child0->maxLength > subDistances[1]+child1->maxLength)
        {
            maxLength = subDistances[0]+child0->maxLength;
        }
        else
        {
            maxLength = subDistances[1]+child1->maxLength;
        }
    }
}

void Node::findMiddlePoint()
{
    halfLength = maxSpan/2;

    if (halfLength >= child0->maxLength && halfLength <= child0->maxLength+subDistances[0]+subDistances[1])
    {

        float b0 = halfLength-child0->maxLength;
        float b1 = subDistances[0]+subDistances[1]-b0;

        char num0[10];
        sprintf(num0,"%.5f",b0);
        char num1[10];
        sprintf(num1,"%.5f",b1);

        mpTree = "("+child0->tree+":"+num0+","+child1->tree+":"+num1+");";

        return;
    }

    child0->findMiddle(1);
    child1->findMiddle(0);

}

void Node::findMiddle(int branch)
{
    if (!isLast)
    {
        if (branch==0)
        {

            if (halfLength >= child0->maxLength && halfLength <= child0->maxLength+subDistances[0])
            {

                float b0 = halfLength-child0->maxLength;
                float b1 = subDistances[0]-b0;

                char num0[10];
                sprintf(num0,"%.5f",b0);
                char num1[10];
                sprintf(num1,"%.5f",b1);

                char num[10];
                sprintf(num,"%.5f",subDistances[1]);

                mpTree = "("+child0->tree+":"+num0+",("+parent->revTrees[1]+","+subTrees[1]+":"+num+"):"+num1+");";

                return;
            }

            if (halfLength >= child1->maxLength && halfLength <= child1->maxLength+subDistances[1])
            {

                float b0 = halfLength-child1->maxLength;
                float b1 = subDistances[1]-b0;

                char num0[10];
                sprintf(num0,"%.5f",b0);
                char num1[10];
                sprintf(num1,"%.5f",b1);

                char num[10];
                sprintf(num,"%.5f",subDistances[0]);

                mpTree = "("+child1->tree+":"+num0+",("+parent->revTrees[1]+","+subTrees[0]+":"+num+"):"+num1+");";

                return;
            }
            child0->findMiddle(1);
            child1->findMiddle(0);

        }
        else
        {

            if (halfLength >= child0->maxLength && halfLength <= child0->maxLength+subDistances[0])
            {

                float b0 = halfLength-child0->maxLength;
                float b1 = subDistances[0]-b0;

                char num0[10];
                sprintf(num0,"%.5f",b0);
                char num1[10];
                sprintf(num1,"%.5f",b1);

                char num[10];
                sprintf(num,"%.5f",subDistances[1]);

                mpTree = "("+child0->tree+":"+num0+",("+parent->revTrees[0]+","+subTrees[1]+":"+num+"):"+num1+");";

                return;
            }


            if (halfLength >= child1->maxLength && halfLength <= child1->maxLength+subDistances[1])
            {

                float b0 = halfLength-child1->maxLength;
                float b1 = subDistances[1]-b0;

                char num0[10];
                sprintf(num0,"%.5f",b0);
                char num1[10];
                sprintf(num1,"%.5f",b1);

                char num[10];
                sprintf(num,"%.5f",subDistances[0]);

                mpTree = "("+child1->tree+":"+num0+",("+parent->revTrees[0]+","+subTrees[0]+":"+num+"):"+num1+");";

                return;
            }
            child0->findMiddle(1);
            child1->findMiddle(0);
        }
    }
}

void Node::divideTree(string tree,string* trees,float* distance)
{

    checkBifurcation(&tree);

    trees[0] = "";

    if ((tree.substr(tree.length()-1)).compare(";")==0)
    {
        tree = tree.substr(0,tree.find_last_of(")")+1); // remove last ';' and anything after the last bracket
    }
    tree = tree.substr(1,tree.length()-2);     // remove first & last '('

    if (tree.at(0)!='(')   // only one taxon before midpoint comma
    {

        string tmp = tree.substr(0,tree.find(",",0));

        trees[0] = tmp;
        distance[0] = defaultBranchLength;
        if(tmp.find(":")!=string::npos)
        {
            trees[0] = tmp.substr(0,tmp.find(":",0));
            distance[0] = atof((tmp.substr(tmp.find(":",0)+1).c_str()));
        }
        else
        {
            has_missing_branch_lengths = true;
        }

        tree = tree.substr(tree.find(",",0)+1);

        bool trifurc = false;
        int open = 0;
        for (unsigned int j = 0; j<tree.length(); j++)
        {
            if (tree.at(j)=='(')
            {
                open++;
            }
            else if (tree.at(j)==')')
            {
                open--;
            }
            if (j>0 && open==0 && tree.substr(j).find(",",0)<=tree.length())
            {
                trifurc = true;
            }
        }

        // correction for trifurcating root
        if (trifurc)
        {
            isUnrooted = true;
            trees[1] = "("+tree+")";
            distance[0] = distance[0]/2;
            distance[1] = distance[0];
        }
        else
        {
            trees[1] = tree;
            distance[1] = defaultBranchLength;
            if(tree.find(":")!=string::npos)
            {
                trees[1] = tree.substr(0,tree.find_last_of(":"));
                tmp = tree.substr(tree.find_last_of(":")+1);

                distance[1] = atof(tmp.c_str());
            }
            else
            {
                has_missing_branch_lengths = true;
            }
        }

    }
    else
    {

        int open = 0;

        for (unsigned int i=0; i<tree.length(); i++)
        {

            // count parentheses that are "open"
            if (tree.at(i)=='(')
            {
                open++;
            }
            else if (tree.at(i)==')')
            {
                open--;
            }
            trees[0].append(tree.substr(i,1));

            if (open<=0)
            {
                distance[0] = defaultBranchLength;
                if(tree.at(i+1) == ':') {

                    string tmp = tree.substr(i+2,tree.find(",",i+2));
                    distance[0] = atof(tmp.c_str());
                }
                else
                {
                    has_missing_branch_lengths = true;
                }

                tree = tree.substr(tree.find(",",i)+1);

                bool trifurc = false;
                open = 0;
                for (unsigned int j = 0; j<tree.length(); j++)
                {
                    if (tree.at(j)=='(')
                    {
                        open++;
                    }
                    else if (tree.at(j)==')')
                    {
                        open--;
                    }

                    if (open==0 && tree.find(",",j)<=tree.length())
                    {
                        trifurc = true;
                    }
                }

                // correction for trifurcating root
                if (trifurc)
                {

                    isUnrooted = true;
                    trees[1] = "("+tree+")";
                    distance[0] = distance[0]/2;
                    distance[1] = distance[0];

                }
                else
                {

                    trees[1] = tree;
                    distance[1] = defaultBranchLength;
                    if(tree.find(":")!=string::npos)
                    {

                        trees[1] = tree.substr(0,tree.find_last_of(":"));

                        string tmp = tree.substr(tree.find_last_of(":")+1);
                        distance[1] = atof(tmp.c_str());
                    }
                    else
                    {
                        has_missing_branch_lengths = true;
                    }
                }
                break;
            }
        }
    }
}



void Node::mark_sequences(vector<string> *names)
{
    if (is_leaf())
    {
        for (int i=0; i<(int)names->size(); i++)
        {
            if (names->at(i)==tree)
            {
                this->has_sequence(true);
                break;
            }
        }
    }
    else
    {
        child0->mark_sequences(names);
        child1->mark_sequences(names);
    }
}

void Node::prune_down()
{
//    cout<<"prune down in "<<this->get_name()<<endl;

    if (this->is_leaf())
        return;

    child0->set_distance_to_parent(subDistances[0]);
    child0->prune_down();
    child1->set_distance_to_parent(subDistances[1]);
    child1->prune_down();

    if (!child0->has_sequence())
    {
        this->delete_left_child();
        this->has_left_child(false);
    }

    if (!child1->has_sequence())
    {
        this->delete_right_child();
        this->has_right_child(false);
    }

    if (this->has_left_child() && !child0->is_leaf())
    {
        if (!child0->has_left_child() && child0->has_right_child())
        {
            Node *new_child = child0->child1;
            new_child->set_distance_to_parent (child0->get_distance_to_parent()+
                                               child0->child1->get_distance_to_parent());

            child0->has_right_child(false);
            this->delete_left_child();
            this->add_left_child(new_child);
        }
        else if (child0->has_left_child() && !child0->has_right_child())
        {
            Node *new_child = child0->child0;
            new_child->set_distance_to_parent (child0->get_distance_to_parent()+
                                               child0->child0->get_distance_to_parent());

            child0->has_left_child(false);
            this->delete_left_child();
            this->add_left_child(new_child);
        }
    }

    if (this->has_right_child() && !child1->is_leaf())
    {
        if (!child1->has_left_child() && child1->has_right_child())
        {
            Node *new_child = child1->child1;
            new_child->set_distance_to_parent (child1->get_distance_to_parent()+
                                               child1->child1->get_distance_to_parent() );

            child1->has_right_child(false);
            this->delete_right_child();
            this->add_right_child(new_child);
        }
        else if (child1->has_left_child() && !child1->has_right_child())
        {
            Node *new_child = child1->child0;
            new_child->set_distance_to_parent (child1->get_distance_to_parent()+
                                               child1->child0->get_distance_to_parent());
            child1->has_left_child(false);
            this->delete_right_child();
            this->add_right_child(new_child);
        }
    }

    if (this->has_left_child() && child0->has_sequence())
        this->has_sequence(true);

    if (this->has_right_child() && child1->has_sequence())
        this->has_sequence(true);

//    cout<<"prune down out "<<this->get_name()<<endl;
}

void Node::prune_up()
{
//    cout<<"prune up in "<<this->get_name()<<endl;

    if (!this->is_leaf() && !this->has_left_child() && this->has_right_child())
    {
        Node* tmp_child = child1;

        child0 = tmp_child->child0;
        child1 = tmp_child->child1;

        tmp_child->has_left_child(false);
        tmp_child->has_right_child(false);

        this->has_left_child(true);
        this->has_right_child(true);

        delete tmp_child;
    }

    if (!this->is_leaf() && this->has_left_child() && !this->has_right_child())
    {
        Node* tmp_child = child0;

        child0 = tmp_child->child0;
        child1 = tmp_child->child1;

        tmp_child->has_left_child(false);
        tmp_child->has_right_child(false);

        this->has_left_child(true);
        this->has_right_child(true);

        delete tmp_child;
    }
//    cout<<"prune up out "<<this->get_name()<<endl;
}
