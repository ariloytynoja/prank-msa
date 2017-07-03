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
extern int version;

extern string seqfile;
extern string seqfile1;
extern string seqfile2;
extern string treefile;
extern string oldtreefile;
extern string treefile1;
extern string treefile2;
extern string outfile;
extern string tempdir;
extern string mafftpath;
extern string exoneratepath;

extern string hmmname;
extern HMModel *hmm;
extern bool HASHMM;

extern bool MERGE;
extern bool TREESTRING;
extern bool PARTLYALIGNED;
extern bool PREALIGNED;
extern bool PRINTSCOREONLY;
extern bool UPDATE;
extern bool UPDATESECOND;
extern float updateTolerance;

extern string dnafile;

extern bool MAFFTALIGNMENT;
extern bool TREEFROMALIGNMENT;
extern bool TREEONLY;
extern bool SCOREMAFFT;
extern bool BPPANCESTORS;

extern int format;

extern bool SCREEN;
extern int reportLimit;

extern bool WRITEANCSEQ;
extern bool LISTEVENTS;
extern bool WRITEITER;
extern bool PARSIMONYSCORE;
extern string INDELSCORE;
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

extern int rnd_seed;
extern bool REPRODUCIBLE;

extern bool CODON;
extern bool DNA;
extern bool PROTEIN;

extern bool DOPOST;
extern bool NOTGAP;
extern bool TERMF;
extern int iterations;

extern bool PRUNETREE;
extern bool PRUNEDATA;
extern bool LOGVALUES;

extern bool TRANSLATE;
extern bool MTTABLE;

extern bool NXis1;

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

extern float defaultBranchLength;
extern float fixedBranchLength;
extern float mergeBranchLength;
extern bool FIXEDBRANCH;
extern bool MAXBRANCH;
extern bool REALBRANCHES;

extern float dnaMaxPairwiseLength;
extern float protMaxPairwiseLength;

extern bool ADJUSTMODEL;

extern bool CORRECTP;
extern bool PENALISEGAPS;

extern bool FULLBAND;
extern bool FULLFULL;
extern int HBW;
extern int FBW;

extern bool SKIPINS;

extern bool EXONERATE;

extern int initialAnchDist;
extern int minAnchDist;
extern int anchSkipDist;
extern int anchDropDist;
extern bool dropRevAnch;

extern bool PATCHMISSING;
extern int missingLimit;
extern bool SKIPGAPANCH;

extern float resizeFactor;
extern float initialMatrixSize;
extern bool PWMATRIXMAXSIZE;
extern float pwInitialMatrixSize;

extern double sumLogs(double a, double b);
extern std::string itos(int i);

//global variable for temporary directory
extern char tmp_dir[1000];

//global variable for verbose flag
extern bool verbose;

#endif
