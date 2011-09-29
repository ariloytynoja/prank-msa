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
#include <iostream>
#include "pwsite.h"
extern bool PWMATRIXMAXSIZE;
extern float pwInitialMatrixSize;

PwSite::PwSite()
{
}
PwSite::PwSite(int i)
{
	in = i;
}
PwSite::~PwSite()
{
}

void PwSite::setMatrices(int longest,int slongest)
{
	int s;
	if(PWMATRIXMAXSIZE)
		s = longest+slongest+2;
	else
		s = (int)(pwInitialMatrixSize*(float)longest);

	lSite = new IntMatrix(s,"pwsite_lSite");
	rSite = new IntMatrix(s,"pwsite_rSite");
	lSite->allowResize(true);
	rSite->allowResize(true);

	cIndex1 = new IntMatrix(s,"pwsite_cIndex1");
	rIndex1 = new IntMatrix(s,"pwsite_rIndex1");
	lIndex1 = new IntMatrix(s,"pwsite_lIndex1");
	cIndex2 = new IntMatrix(s,"pwsite_cIndex2");
	rIndex2 = new IntMatrix(s,"pwsite_rIndex2");
	lIndex2 = new IntMatrix(s,"pwsite_lIndex2");
	cIndex1->allowResize(true);
	rIndex1->allowResize(true);
	lIndex1->allowResize(true);
	cIndex2->allowResize(true);
	rIndex2->allowResize(true);
	lIndex2->allowResize(true);

	currMS = new IntMatrix(s,"pwsite_currMS");
	currMS->allowResize(true);

	vfX = new IntMatrix(s,"pwsite_vfX");
	vfY = new IntMatrix(s,"pwsite_vfY");
	vfM = new IntMatrix(s,"pwsite_vfM");
	vbX = new IntMatrix(s,"pwsite_vbX");
	vbY = new IntMatrix(s,"pwsite_vbY");
	vbM = new IntMatrix(s,"pwsite_vbM");
	vfX->allowResize(true);
	vfY->allowResize(true);
	vfM->allowResize(true);
	vbX->allowResize(true);
	vbY->allowResize(true);
	vbM->allowResize(true);

}

void PwSite::deleteMatrices()
{
	delete lSite;
	delete rSite;

	delete cIndex1;
	delete rIndex1;
	delete lIndex1;
	delete cIndex2;
	delete rIndex2;
	delete lIndex2;

	delete currMS;

	delete vfX;
	delete vfY;
	delete vfM;
	delete vbX;
	delete vbY;
	delete vbM;

}

int PwSite::aSize=4;

int PwSite::count = 2;

IntMatrix *PwSite::lSite;
IntMatrix *PwSite::rSite;

IntMatrix *PwSite::cIndex1;
IntMatrix *PwSite::lIndex1;
IntMatrix *PwSite::rIndex1;
IntMatrix *PwSite::cIndex2;
IntMatrix *PwSite::lIndex2;
IntMatrix *PwSite::rIndex2;

IntMatrix *PwSite::currMS;

IntMatrix *PwSite::vfX;
IntMatrix *PwSite::vfY;
IntMatrix *PwSite::vfM;
IntMatrix *PwSite::vbX;
IntMatrix *PwSite::vbY;
IntMatrix *PwSite::vbM;

