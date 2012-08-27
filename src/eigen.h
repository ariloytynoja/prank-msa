/***************************************************************************
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

/*
  Copyright (C) by Ziheng Yang except where otherwise stated.
  The code is adapted from Ziheng Yang's software package PAML 3.14.
*/


#ifndef EIGEN_H
#define EIGEN_H

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <time.h>

class Eigen
{
public:
    Eigen();
    ~Eigen();

    int getpi_sqrt (double pi[], double pi_sqrt[], int n, int *npi0);
    int eigenQREV (double Q[], double pi[], double pi_sqrt[], int n, int npi0, double Root[], double U[], double V[]);
    int eigenRealSym(double A[], int n, double Root[], double work[]);

    void HouseholderRealSym(double a[], int n, double d[], double e[]);
    int EigenTridagQLImplicit(double d[], double e[], int n, double z[]);
    void EigenSort(double d[], double U[], int n);

/*
  Copyright (C) by Simon Whelan.
*/
    void computePMatrix(int n, double* pMat, double* U, double* V, double* Root, double time);

};

#endif
