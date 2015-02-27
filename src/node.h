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
#include <iostream>
#include <sstream>

using namespace std;
class Node
{
    std::string tree;
    std::string subTrees[2];
    std::string revTrees[2];

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
            ss<<"("<<child0->print_subtree()<<","<<child1->print_subtree()<<"):"<<dist_to_parent;
            return ss.str();
        }
        else
        {
            stringstream ss;
            ss<<tree<<":"<<dist_to_parent;
            return ss.str();
        }
    }

    static int count;
public:
    Node();
    ~Node();

    Node(std::string t);
    std::string rootedTree();

    void mark_sequences(vector<string> *names);

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
};


#endif
