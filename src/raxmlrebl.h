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

#ifndef RAXMLREBL_H
#define RAXMLREBL_H

#include "ancestralnode.h"
#include "config.h"
#include <sys/stat.h>
#include "node.h"

class RaxmlRebl
{
    std::string raxmlpath;

public:
    RaxmlRebl();
    bool testExecutable();
    bool inferBranchLengths(AncestralNode *root,vector<string> *names,vector<string> *sequences,bool isDna);
};

#endif // RAXMLREBL_H
