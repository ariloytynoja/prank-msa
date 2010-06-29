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
#ifndef READNEWICK_H
#define READNEWICK_H

/**
 * Reader for newick-format treefiles
 */

#include <string>
#include <map>
#include "treenode.h"

class ReadNewick{
    std::string s;
    std::string root;
    std::map<std::string,TreeNode*> nodes;
public:
    ReadNewick();
    ~ReadNewick();

    std::string readFile(const char* filename);
    void buildTree(std::string s,std::map<std::string,TreeNode*>* nodes);
    std::string getRoot(){ return root; }
};

#endif
