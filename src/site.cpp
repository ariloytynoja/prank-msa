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

#include "site.h"
#include "hmmodel.h"
#include <iostream>

extern HMModel* hmm;
extern float initialMatrixSize;

Site::Site()
{
}
Site::Site(int i)
{
    in = i;
}

Site::~Site()
{
}

void Site::setMatrices(int longest,int )
{
    int s = (int)(initialMatrixSize*(float)longest);

    anc = new BoolMatrix(s,"site_anc");
    nus = new BoolMatrix(s,"site_nus");
    anc->allowResize(true);
    nus->allowResize(true);

    lSite = new IntMatrix(s,"site_lSite");
    rSite = new IntMatrix(s,"site_rSite");
    lSite->allowResize(true);
    rSite->allowResize(true);

    cIndex1 = new IntMatrix(s,"site_cIndex1");
    nIndex1 = new IntMatrix(s,"site_nIndex1");
    rIndex1 = new IntMatrix(s,"site_lIndex1");
    lIndex1 = new IntMatrix(s,"site_lIndex1");
    cIndex2 = new IntMatrix(s,"site_cIndex2");
    nIndex2 = new IntMatrix(s,"site_nIndex2");
    rIndex2 = new IntMatrix(s,"site_rIndex2");
    lIndex2 = new IntMatrix(s,"site_lIndex2");
    cIndex1->allowResize(true);
    nIndex1->allowResize(true);
    rIndex1->allowResize(true);
    lIndex1->allowResize(true);
    cIndex2->allowResize(true);
    nIndex2->allowResize(true);
    rIndex2->allowResize(true);
    lIndex2->allowResize(true);

    currMS = new IntMatrix(s,"site_currMS");
    currSS = new IntMatrix(s,"site_currSS");
    currMS->allowResize(true);
    currSS->allowResize(true);

    permIns = new IntMatrix(s,"site_permIns");
    permIns->allowResize(true);
    permIns->initialise(0);

    vf = new FlMatrix(s,"site_vf");
    vfM = new IntMatrix(s,"site_vfM");
    vfS = new IntMatrix(s,"site_vfS");
    vb = new FlMatrix(s,"site_vb");
    vbM = new IntMatrix(s,"site_vbM");
    vbS = new IntMatrix(s,"site_vbS");
    vf->allowResize(true);
    vfM->allowResize(true);
    vfS->allowResize(true);
    vb->allowResize(true);
    vbM->allowResize(true);
    vbS->allowResize(true);

    ffX = new FlMatrix(nState,s,"site_ffX");
    ffY = new FlMatrix(nState,s,"site_ffY");
    ffM = new FlMatrix(nState,s,"site_ffM");
    fbX = new FlMatrix(nState,s,"site_fbX");
    fbY = new FlMatrix(nState,s,"site_fbY");
    fbM = new FlMatrix(nState,s,"site_fbM");
    ffX->allowResize(false,true);
    ffY->allowResize(false,true);
    ffM->allowResize(false,true);
    fbX->allowResize(false,true);
    fbY->allowResize(false,true);
    fbM->allowResize(false,true);


    mcp = new DbMatrix(nState,aSize,s,"site_mcp");
    stp = new FlMatrix(nState,s,"site_stp");
    pop = new FlMatrix(s,"site_pop");
    mcp->allowResize(false,false,true);
    stp->allowResize(false,true);
    pop->allowResize(true);

}

void Site::deleteMatrices()
{
    delete anc;
    delete nus;

    delete lSite;
    delete rSite;

    delete cIndex1;
    delete nIndex1;
    delete rIndex1;
    delete lIndex1;
    delete cIndex2;
    delete nIndex2;
    delete rIndex2;
    delete lIndex2;

    delete currMS;
    delete currSS;

    delete permIns;

    delete vf;
    delete vfM;
    delete vfS;
    delete vb;
    delete vbM;
    delete vbS;

    delete ffX;
    delete ffY;
    delete ffM;
    delete fbX;
    delete fbY;
    delete fbM;

    delete mcp;
    delete stp;
    delete pop;

}

int Site::count = 2;

BoolMatrix *Site::anc;
BoolMatrix *Site::nus;

IntMatrix *Site::lSite;
IntMatrix *Site::rSite;

IntMatrix *Site::cIndex1;
IntMatrix *Site::nIndex1;
IntMatrix *Site::lIndex1;
IntMatrix *Site::rIndex1;
IntMatrix *Site::cIndex2;
IntMatrix *Site::nIndex2;
IntMatrix *Site::lIndex2;
IntMatrix *Site::rIndex2;

IntMatrix *Site::currMS;
IntMatrix *Site::currSS;

IntMatrix *Site::permIns;

FlMatrix *Site::vf;
IntMatrix *Site::vfM;
IntMatrix *Site::vfS;
FlMatrix *Site::vb;
IntMatrix *Site::vbM;
IntMatrix *Site::vbS;

FlMatrix *Site::ffX;
FlMatrix *Site::ffY;
FlMatrix *Site::ffM;
FlMatrix *Site::fbX;
FlMatrix *Site::fbY;
FlMatrix *Site::fbM;

DbMatrix *Site::mcp;
FlMatrix *Site::stp;
FlMatrix *Site::pop;

int Site::aSize=4;

int Site::nState=1;

