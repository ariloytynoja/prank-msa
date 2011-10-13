/***************************************************************************
 *   Copyright (C) 2005 by Ari Loytynoja                                   *
 *   ari@sink                                                              *
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

#ifndef CONFIG_H
#define CONFIG_H

#ifndef RFOR
#define RFOR(i,n) for(i=n; i>=0; i--)
#endif
#ifndef FOR
#define FOR(i,n) for(i=0; i<n; i++)
#endif

#include "hmmodel.h"

using namespace std;

extern int NOISE;

extern string seqfile;
extern string treefile;
extern string outfile;

extern string hmmname;
extern HMModel *hmm;
extern bool HASHMM;

extern bool TREESTRING;
extern bool PARTLYALIGNED;
extern bool PREALIGNED;

extern string dnafile;

extern bool TREEFROMALIGNMENT;

extern int format;

extern bool SCREEN;
extern int reportLimit;

extern bool WRITEANC;
extern bool WRITEANCSEQ;
extern bool PRINTNODES;
extern bool PRINTTREE;
extern bool WRITEXML;

extern bool DOTS;

extern bool CONVERT;
extern bool SHORTNAMES;
extern bool BACKTRANSLATE;

extern std::string message;
extern std::string currentNode;

extern bool FOREVER;
extern bool FOREVER_FOR_PA;
extern bool FOREVER_OLD;

extern string dnaFreqs;
extern float gapRate;
extern float gapExt;
extern float kappa;
extern float rho;

extern float dnaGapRate;
extern float dnaGapExt;
extern float protGapRate;
extern float protGapExt;

extern bool CODON;
extern bool DNA;
extern bool PROTEIN;

extern bool DOPOST;
extern bool NOTGAP;
extern bool TERMF;
extern bool TWICE;
extern bool PRUNETREE;
extern bool LOGVALUES;

extern bool TRANSLATE;
extern bool MTTABLE;

extern bool NXis1;

extern string annofile;
extern bool PRIORS;

extern float pwDist;
extern float pwDnaDist;
extern float pwGapRate;
extern float pwGapExt;
extern float pwProtDist;
extern float pwProtGapRate;
extern float pwProtGapExt;

extern float minBrL;
extern float branchScalingFactor;
extern float dnaBranchScalingFactor;
extern float protBranchScalingFactor;

extern float fixedBranchLength;
extern bool FIXEDBRANCH;
extern bool MAXBRANCH;
extern bool REALBRANCHES;

extern float dnaMaxPairwiseLength;
extern float protMaxPairwiseLength;

extern bool ADJUSTMODEL;

extern bool CORRECTP;

extern bool FULLBAND;
extern bool FULLFULL;
extern int HBW;
extern int FBW;

extern bool SKIPINS;

extern bool ANCHORS;
extern bool EXONERATE;

extern int initialAnchDist;
extern int maxAnchDist;
extern int minAnchDist;
extern int anchSkipDist;
extern int anchDropDist;

extern bool PATCHMISSING;
extern int missingLimit;
extern bool SKIPGAPANCH;

extern float resizeFactor;
extern float initialMatrixSize;
extern bool PWMATRIXMAXSIZE;
extern float pwInitialMatrixSize;

extern bool PWGENOMIC;
extern float pwgendist;

extern int fOnNode;

extern std::string anchorfile;
extern bool HARDANCHORS;

extern double sumLogs(double a, double b);
extern std::string itos(int i);

#endif
