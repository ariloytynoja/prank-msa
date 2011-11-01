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
#ifndef GUIDETREE_H
#define GUIDETREE_H

#include <string>
#include <vector>
#include "flmatrix.h"
#include "intmatrix.h"

class GuideTree
{
    std::string tree;
public:
    GuideTree(std::vector<std::string>* seqs,std::vector<std::string>* names,IntMatrix* substScores);
    GuideTree(std::vector<std::string>* seqs,std::vector<std::string>* names,bool idDna);

    ~GuideTree();
    std::string getTree()
    {
        return tree;
    }

    void makeTree(FlMatrix* distance, std::vector<std::string>* names);
    void joinNeighbors(FlMatrix* distance, std::string* names,FlMatrix* newDistance, std::string* newNames,FlMatrix* rDist,int* no);

};

#endif
