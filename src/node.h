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
#ifndef NODE_H
#define NODE_H

/*
 * Clumsy way of re-rooting
 */

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include <stdlib.h>

using namespace std;
class Node
{
    std::string tree;
    std::string subTrees[2];
    std::string revTrees[2];

    std::string reverseTree;

    Node* parent;
    Node* child0;
    Node* child1;

    float subDistances[2];
    float maxLength;

    bool isLast;
    bool isFirst;
    bool isUnrooted;

    Node(std::string t,Node* p,int branch);

    void checkBifurcation(std::string *t);
    void findMiddlePoint();
    void findMiddle(int branch);
    void divideTree(std::string tree,std::string* trees,float* distances);

    bool has_missing_branch_lengths;

    bool node_has_sequence;
    bool node_has_left_child;
    bool node_has_right_child;

    bool is_leaf()
    {
        return isLast;
    }
    void is_leaf(bool i)
    {
        isLast = i;
    }

    bool has_left_child()
    {
        return node_has_left_child;
    }
    void has_left_child(bool h)
    {
        node_has_left_child = h;
    }

    bool has_right_child()
    {
        return node_has_right_child;
    }
    void has_right_child(bool h)
    {
        node_has_right_child = h;
    }


    void has_sequence(bool s)
    {
        node_has_sequence = s;
    }
    bool has_sequence()
    {
        return node_has_sequence;
    }

    void prune_up();
    void prune_down();

    string name;
    void set_name(string n)
    {
        name = n;
    }
    string get_name()
    {
        return name;
    }

    int namehash;

    double dist_to_parent;
    void set_distance_to_parent(double d)
    {
        dist_to_parent = d;
    }
    double get_distance_to_parent()
    {
        return dist_to_parent;
    }

    void add_left_child(Node *child)
    {
        child0 = child;
        is_leaf(false);
        this->has_left_child(true);
    }
    void add_right_child(Node *child)
    {
        child1 = child;
        is_leaf(false);
        this->has_right_child(true);
    }

    void delete_left_child()
    {
        delete child0;
    }
    void delete_right_child()
    {
        delete child1;
    }

    string print_subtree()
    {
        if (!is_leaf())
        {
            stringstream ss;
            char num[10];
            sprintf(num,"%.5f",dist_to_parent);
            ss<<"("<<child0->print_subtree()<<","<<child1->print_subtree()<<"):"<<num; //dist_to_parent;
            return ss.str();
        }
        else
        {
            stringstream ss;
            char num[10];
            sprintf(num,"%.5f",dist_to_parent);
            ss<<tree<<":"<<num; //dist_to_parent;
            return ss.str();
        }
    }

//    int hash(string s)
//    {

//        int hash = 0;
//        int offset = 'a' - 1;
//        for(string::const_iterator it=s.begin(); it!=s.end(); ++it) {
//          hash = hash << 1 | (*it - offset);
//        }
//        return hash;
//    }

    unsigned int hash(const char* s,unsigned int seed = 0)
    {
        unsigned int hash = seed;
        while (*s)
        {
            hash = hash * 101  +  *s++;
        }
        return hash;
    }

    float halfLength;
    float maxSpan;

    static bool warnings;

    static std::string mpTree;
    static int count;
    static int warned;
public:
    Node();
    ~Node();

    Node(std::string t, bool root = true);
    std::string rootedTree();


    void mark_sequences(vector<string> *names);

    string markNodes();

    int findMarkedNode(int h0, int h1) {
        if (is_leaf())
            return 0;

//        cout<<name<<" "<<namehash<<endl;
        if(this->namehash == h0 || this->namehash == h1)
        {
//            cout<<name<<" "<<namehash<<endl;
            return this->namehash;
        }
        else
        {
            int f0 = child0->findMarkedNode(h0,h1);
            if(f0 == h0 || f0 == h1)
                return f0;
            else
                return child1->findMarkedNode(h0,h1);
        }
    }

    bool rootAtMarkedNode(int h0, int h1) {
        if (is_leaf())
            return false;

//        cout<<"root "<<name<<" "<<namehash<<endl;

        if(this->namehash == h1 || this->namehash == h0)
        {
            string rt = this->reverseTree;
            int pos = rt.find_last_of(':');
            if(pos != string::npos)
            {
//                cout<<"\nb1 "<<dist_to_parent<<endl;
//                cout<<"\nf "<<this->tree<<endl;
//                cout<<"r "<<this->reverseTree<<endl<<endl;

                string num = rt.substr(pos+1);
                rt = rt.substr(0,pos);
                stringstream ss(num);
                float f;
                ss >> f;
                f/=2;
                stringstream fs;
                fs << f;
                mpTree = "("+this->tree+":"+fs.str()+","+rt+":"+fs.str()+");";
            }
            else
            {
                cout<<"\nb "<<dist_to_parent<<endl;
                cout<<"\nf "<<this->tree<<endl;
                cout<<"r "<<this->reverseTree<<endl<<endl;

                mpTree = "("+this->tree+":0,"+this->reverseTree+");";
            }
            return true;
        }
//        else if(this->namehash == h0)
//        {
//            string rt = this->reverseTree;
//            int pos = rt.find_last_of(':');
//            if(pos != string::npos)
//            {
//                cout<<"\nb2 "<<dist_to_parent<<endl;
//                cout<<"\nf "<<this->tree<<endl;
//                cout<<"r "<<this->reverseTree<<endl<<endl;

//                string num = rt.substr(pos+1);
//                rt = rt.substr(0,pos);
//                stringstream ss(num);
//                float f;
//                ss >> f;
//                f/=2;
//                stringstream fs;
//                fs << f;
//                mpTree = "("+rt+":"+fs.str()+","+this->tree+":"+fs.str()+");";
//            }
//            else
//            {
//                cout<<"\nb "<<dist_to_parent<<endl;
//                cout<<"\nf "<<this->tree<<endl;
//                cout<<"r "<<this->reverseTree<<endl<<endl;

//                mpTree = "("+this->reverseTree+","+this->tree+":0);";
//            }
//            return true;
//        }


/*
        if(child0->namehash == h0 || child0->namehash == h1)
        {
            halfLength = maxSpan/2;

//            cout<<"root "<<child0->name<<" "<<child0->namehash<<endl;

            float b0 = halfLength-child0->maxLength;
            if(b0<0)
                b0=0.001;
            float b1 = subDistances[0]-b0;

            char num0[10];
            sprintf(num0,"%.5f",b0);
            char num1[10];
            sprintf(num1,"%.5f",b1);

            char num[10];
            sprintf(num,"%.5f",subDistances[1]);

            mpTree = "("+child0->tree+":"+num0+",("+subTrees[1]+":"+num+","+parent->revTrees[0]+"):"+num1+");";

            cout<<"\n\n1a: "<<child0->tree<<":"<<num0<<endl
               <<"2a: "
              <<subTrees[1]<<":"<<num<<endl
             <<"3a: "
            <<parent->revTrees[0]<<endl
            <<endl;

            cout<<mpTree<<endl<<endl;

            return true;
        }
        else if(child1->namehash == h0 || child1->namehash == h1)
        {
            halfLength = maxSpan/2;

//            cout<<"root "<<child1->name<<" "<<child1->namehash<<endl;

            float b0 = halfLength-child1->maxLength;
            if(b0<0)
                b0=0.001;
            float b1 = subDistances[1]-b0;

            char num0[10];
            sprintf(num0,"%.5f",b0);
            char num1[10];
            sprintf(num1,"%.5f",b1);

            char num[10];
            sprintf(num,"%.5f",subDistances[0]);

            int b = 0;
            if(parent->child1->name == this->name)
                b = 1;
            mpTree = "(("+parent->revTrees[b]+","+subTrees[0]+":"+num+"):"+num1+","+child1->tree+":"+num0+");";
//            cout<<"\n\n1b: "<<child1->tree<<":"<<num0<<endl<<"2b: "<<parent->revTrees[b]<<endl<<"3b: "<<subTrees[0]<<":"<<num<<endl<<endl;

//            cout<<mpTree<<endl;
            return true;

        }
*/
        else
        {
            if(child0->rootAtMarkedNode(h0,h1))
                return true;
            else
                return child1->rootAtMarkedNode(h0,h1);
        }
    }

    void prune_tree()
    {
        this->prune_down();
        this->prune_up();
    }

    void printTerminal()
    {
        if (!isLast)
        {
            child0->printTerminal();
            child1->printTerminal();
        }
        else
        {
            cout<<this->get_name()<<" "<<this->has_sequence()<<endl;
        }
    }

    void countMatchingLeaves(int *leaves,int *matches)
    {
        if (this->has_left_child())
            child0->countMatchingLeaves(leaves,matches);
        if (this->has_right_child())
            child1->countMatchingLeaves(leaves,matches);

        if (is_leaf())
        {
            (*leaves)++;
            if (this->has_sequence())
                (*matches)++;
        }
    }

    void collectUnmatchingLeaves(vector<string> *unmatching)
    {
        if (this->has_left_child())
            child0->collectUnmatchingLeaves(unmatching);
        if (this->has_right_child())
            child1->collectUnmatchingLeaves(unmatching);

        if (is_leaf())
        {
            if (!this->has_sequence())
                unmatching->push_back(name);
        }
    }


    string print_tree()
    {
        if (!is_leaf())
        {
            stringstream ss;
            ss<<"(";
            bool hasleft = false;
            if (this->has_left_child())
            {
                ss<<child0->print_subtree();
                hasleft = true;
            }
            if (this->has_right_child())
            {
                if (hasleft)
                    ss<<",";

                ss<<child1->print_subtree();
            }
            ss<<");";
            return ss.str();
        }
        else
        {
            return "";
        }
    }

    void getRootNamehashes(int *hash0,int *hash1)
    {
        *hash0 = child0->namehash;
        *hash1 = child1->namehash;
    }

    void getNamehashes(vector<int> *h)
    {
        if(!is_leaf())
        {
            child0->getNamehashes(h);
            child1->getNamehashes(h);
        }
        h->push_back(namehash);
        return;
    }

    void getNameNamehashes(map<string,int> *h)
    {
        if(!is_leaf())
        {
            child0->getNameNamehashes(h);
            child1->getNameNamehashes(h);
        }
        h->insert(make_pair(name,namehash));
        return;
    }

    void getNamehashNames(map<int,string> *h)
    {
        if(!is_leaf())
        {
            child0->getNamehashNames(h);
            child1->getNamehashNames(h);
        }
        h->insert(make_pair(namehash,name));
        return;
    }

    void getNamehashesDistance(map<int,float> *h)
    {
        if(!is_leaf())
        {
            child0->getNamehashesDistance(h);
            child1->getNamehashesDistance(h);
        }

        float dist;
        char num[10];
        sprintf(num,"%.5f",dist_to_parent);
        dist = strtof(num, NULL);
        h->insert(make_pair(namehash,dist));
        return;
    }

    int getNodeNumber()
    {
        if(!is_leaf())
            return child0->getNodeNumber()+child1->getNodeNumber()+1;
        return 1;
    }
};


#endif
