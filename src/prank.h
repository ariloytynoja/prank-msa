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

#ifndef PRANK_H
#define PRANK_H

#include <string>
#include <sstream>
#include "hmmodel.h"

using namespace std;

void readArguments(int argc, char *argv[]);
int parseFormat(string format);
void printHelp(bool complete);

int version;

int NOISE = 0;

/********* input/output: **********/

// sequence data file
string seqfile = "";
string seqfile1 = "";
string seqfile2 = "";

// guide tree file
string treefile = "";
string oldtreefile = "";
string treefile1 = "";
string treefile2 = "";

// alignment output file
string outfile = "output";

// alternative temp dir
string tempdir = "";

string mafftpath = "";
string exoneratepath = "";

// structure model file
string hmmname = "";
HMModel *hmm;
bool HASHMM = false;

// guide tree as a string
bool TREESTRING = false;

// mixture of existing and new alignments for Ziheng
bool PARTLYALIGNED = false;

bool PREALIGNED = false;
bool PRINTSCOREONLY = false;

bool UPDATE = false;
bool UPDATESECOND = true;
float updateTolerance = 0.1;

// DNA sequence data file for backtranslation
string dnafile = "";

/********* more input/output: **********/

// estimate guide tree from mafft alignment (before Prank realignment)
bool MAFFTALIGNMENT = true;

// estimate guide tree from input alignment (before realignment)
bool TREEFROMALIGNMENT = false;

// merge two alignments
bool MERGE = false;

// output alignment format
int format = 8;

// do output backtabs
bool SCREEN = true;

// reporting interval
int reportLimit = 100;

// write reconstructed ancestral seqs
bool WRITEANCSEQ = false;

// write iterations
bool WRITEITER = false;

// write evolutionary events
bool LISTEVENTS = false;

// compute parsimony score
bool PARSIMONYSCORE = true;
string INDELSCORE = "";

// write ancestral nodes as they are solved
bool PRINTNODES = false;

// don't print tree
bool PRINTTREE = false;

// don't write xml
bool WRITEXML = false;

// print dots for insertions
bool DOTS = false;

// no align, convert only
bool CONVERT = false;

// use short names (until first space)
bool SHORTNAMES = false;

// backtranslate existing alignment
bool BACKTRANSLATE = false;

// info on the screen
std::string message = "";
std::string currentNode = "";

/********* model: gaps and F **********/

// keep insertion forever
bool FOREVER = false;

// old implementation
bool FOREVER_OLD = false;

/********* model: substitutions, indels **********/

string dnaFreqs = "";
float gapRate = -1;
float gapExt = -1;
float kappa = 2;
float rho = 1;

float dnaGapRate = 0.025;
float dnaGapExt = 0.75;
float protGapRate = 0.005;
float protGapExt = 0.5;

/******* random number seed ********/

int rnd_seed = -1;

/********* model: other **********/

// codon alignment
bool CODON = false;
bool DNA = false;
bool PROTEIN = false;

// no posterior probabiliity calculation
bool DOPOST = false;

// penalise terminal gaps
bool NOTGAP = true;

// terminal gaps considered as normal gaps
bool TERMF = false;

// run iterations
int iterations = 5;

// prune the tree
bool PRUNETREE = false;

// prune the data
bool PRUNEDATA = false;

// use log values
bool LOGVALUES = false;


/********* more model: **********/

// translate DNA to protein, then backtranslate
bool TRANSLATE = false;

// translate mtDNA to protein, then backtranslate
bool MTTABLE = false;

// consider N or X identical to any
bool NXis1 = true;


/********* more model: pairwise alignment for guide tree **********/

// expected pairwise distance
float pwDist = -1;

// expected pairwise distance
float pwDnaDist = 0.25;

// defaults for pw alignment
float pwGapRate = -1;
float pwGapExt = -1;

float pwProtDist = 0.5;
float pwProtGapRate = 0.005;
float pwProtGapExt = 0.50;

/********* more model: branch lengths in guide tree **********/

// minimum length
float minBrL = 0.0001;

// scale branch lengths
float branchScalingFactor = -1;

float dnaBranchScalingFactor = 1.0;
float protBranchScalingFactor = 1.0;

// set branch lengths
float defaultBranchLength = 0.1; // used if nothing else provided
float fixedBranchLength = -1;
float mergeBranchLength = -1;
bool FIXEDBRANCH = false;
bool MAXBRANCH = false;

float dnaMaxPairwiseLength = 0.3;
float protMaxPairwiseLength = 0.5;

bool ADJUSTMODEL = false;

// use real guidetree distances
bool REALBRANCHES = false;

// correct guidetree distances
bool CORRECTP = false;


/********* technical: hirschberg, full probability **********/

// "band" full probability (less memory)
bool FULLBAND = false;

// complete full probability
bool FULLFULL = false;

// hirschberg band width (for hirschbergalignment)
int HBW = 50;

// full probability band width (for fullprobability)
int FBW = 50;

// skip insertions in postprobs
bool SKIPINS = false;

/********* technical: anchoring **********/

// use anchors
bool EXONERATE = true;

// initial anchor distance
int initialAnchDist = 500;

// minimum anchor distance
int minAnchDist = 30;

// anchor skip distance
int anchSkipDist = 200;

// anchor drop distance
int anchDropDist = 10;

// don't infer gaps caused by missing data
bool PATCHMISSING = true;

// length of gap deemed as missing data
int missingLimit = 1000;

// skip gaps in anchoring ancestral seqs (?)
bool SKIPGAPANCH = true;

/********* technical: memory & speed efficiency **********/

// matrix resize factor
float resizeFactor = 1.75;

// matrix initial factor
float initialMatrixSize = 1.5;

// use pwmatrix maximum size
bool PWMATRIXMAXSIZE = true;

float pwInitialMatrixSize = 1.5;


/************************************************/

//
double sumLogs(double a, double b)
{
    if (a==-HUGE_VAL && b==-HUGE_VAL)
    {
        return -HUGE_VAL;
    }
    else if (a==-HUGE_VAL)
    {
        return b;
    }
    else if (b==-HUGE_VAL)
    {
        return a;
    }
    if (b>a)
    {
        double c = a;
        a = b;
        b = c;
    }

    return (a+log(1+exp(b-a)));
}

std::string itos(int i)
{
    std::stringstream s;
    s << i;
    return s.str();
}

#endif

