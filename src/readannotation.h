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
#ifndef READANNOTATION_H
#define READANNOTATION_H

#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include "flmatrix.h"

class ReadAnnotation
{
    int end,k;
public:

    ReadAnnotation(std::string annofile,std::map<std::string,FlMatrix*>*);

    ~ReadAnnotation();

    std::string nextNotComment(std::ifstream* in);
    int nextInt(std::string row);
    double nextDouble(std::string row);
};

#endif
