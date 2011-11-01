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
#include "readannotation.h"
#include "config.h"

using namespace std;

ReadAnnotation::ReadAnnotation(string annofile,map<string,FlMatrix*>* annotation)
{
    ifstream in(annofile.c_str());

    int nState = hmm->getNStates();

    if (!in)
    {
        cout<<"Could not read annotation file "<<annofile<<"!"<<endl;
        exit(-1);
    }

    string row = nextNotComment(&in);
    while (row.length()>0)
    {
        string name;
        if (row.at(0)=='>')
            name = row.substr(1);
        else
        {
            cout<<"Annotation file "<<annofile<<": erroneous name '"<<row<<"'"<<endl;
            exit(-1);
        }
        row = nextNotComment(&in);
        int nstate = nextInt(row);
        if (nstate != nState)
        {
            cout<<"Annotation "<<annofile<<": erroneous nstate '"<<nstate<<"'"<<endl;
            exit(-1);
        }
        int nsite = nextInt(row);
        FlMatrix* anno = new FlMatrix(nsite,nstate,"annotation");
        if (LOGVALUES)
            anno->initialise(0);
        else
            anno->initialise(1);
        row = nextNotComment(&in);
        int i;
        FOR(i,nsite)
        {
// 			cout<<i<<" "<<row<<endl;
            FOR(k,nstate)
            if (LOGVALUES)
                anno->s(log(nextDouble(row)),i,k);
            else
                anno->s(nextDouble(row),i,k);
            getline(in,row);
            end= -1;
        }
        annotation->insert(make_pair(name,anno));
        row = nextNotComment(&in);
    }
}


ReadAnnotation::~ReadAnnotation()
{
}

double ReadAnnotation::nextDouble(string row)
{
    int tmp=row.find_first_of("0123456789.",end+1);
    end=row.find_first_not_of("0123456789.",tmp);
    return atof(row.substr(tmp,end).c_str());
}

int ReadAnnotation::nextInt(string row)
{
    int tmp=row.find_first_of("0123456789",end+1);
    end=row.find_first_not_of("0123456789",tmp);
    return atoi(row.substr(tmp,end).c_str());
}

std::string ReadAnnotation::nextNotComment(ifstream *in)
{
    end = -1;

    string row = "";
    getline(*in,row);
    while (row[0]=='#')
    {
        getline(*in,row);
    }
    return row;
}

