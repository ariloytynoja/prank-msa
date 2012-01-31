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
#ifndef PROGRESSIVEALIGNMENT_H
#define PROGRESSIVEALIGNMENT_H

/**
 * Wrapper for sequence loading and multiple alignment.
 */

#include <string>
#include <map>
#include "config.h"
#include "ancestralnode.h"
#include "treenode.h"
#include "translatesequences.h"

class ProgressiveAlignment
{
    TranslateSequences *trseq;
public:
    ProgressiveAlignment(std::string treefile,std::string seqfile,std::string dnafile,int aMethod=0);
    ~ProgressiveAlignment();

    void printXml(AncestralNode *root,int iteration,bool translate);
    void getAlignmentMatrix(AncestralNode *root,char* alignment,bool translate);
    void printAncestral(AncestralNode *root,std::vector<std::string> *names,std::vector<std::string> *seqs,int iteration);
    void getAncestralAlignmentMatrix(AncestralNode *root,char* alignment);
    void getAncestralAlignmentSeqs(AncestralNode *root,map<string,string> *anc_seqs);
    void printAlignment(AncestralNode *root,std::vector<std::string> *nms,std::vector<std::string> *sqs,int iteration, bool isDna);
    std::string formatExtension(int format);
private:

    void makeSettings(bool isDna)
    {
        if (isDna)
        {
            if (gapRate<0)
                gapRate = dnaGapRate;

            if (gapExt<0)
                gapExt = dnaGapExt;

            if (pwDist<0)
                pwDist = pwDnaDist;

            if (pwGapRate<0)
                pwGapRate = dnaGapRate;

            if (pwGapExt<0)
                pwGapExt = dnaGapExt;

            if (branchScalingFactor<0)
                branchScalingFactor = dnaBranchScalingFactor;

        }
        else
        {
            if (gapRate<0)
                gapRate = protGapRate;

            if (gapExt<0)
                gapExt = protGapExt;

            if (pwDist<0)
                pwDist = pwProtDist;

            if (pwGapRate<0)
                pwGapRate = pwProtGapRate;

            if (pwGapExt<0)
                pwGapExt = pwProtGapExt;

            if (branchScalingFactor<0)
                branchScalingFactor = protBranchScalingFactor;

        }

    }
};

#endif
