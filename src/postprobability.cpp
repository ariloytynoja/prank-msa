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

#include <iostream>
#include <cmath>
#include "config.h"
#include "postprobability.h"

using namespace std;

PostProbability::~PostProbability(){
}

PostProbability::PostProbability(Sequence* ,Sequence* ,double fullScore,PhyloMatchScore *msr)
{
    int nState = hmm->getNStates();

	Site *pSite = new Site();
	pSite->index(0);
	Site *cSite = new Site();
	cSite->index(0);

    while (pSite->nullSite()){ // skip 'null' sites
        pSite->next();
    }
    while (cSite->nullSite()){
        cSite->next();
    }
    cSite->next();

    for (;cSite->index()!=1; cSite->next(),pSite->next()) {
        while (pSite->nullSite()){ // skip 'null' sites
            pSite->next();
        }

		while (cSite->nullSite()){
            cSite->postProb(-1);
            for (int k=0;k<nState;k++) {
                cSite->stateProb(-1,k);
            }

            cSite->next();

            if (cSite->index()==1) // stop if skipping brings to the end
                break;
        }

		if (cSite->index()==1) // stop if skipping brings to the end
            break;

        msr->computeFullFwd(cSite->nInd1(),cSite->nInd2());

        int moveTo = cSite->currMatchState();

        // Sum probabilities of all moves doing the same alignment
        // (i.e., the reliability of the solution)
        //
        double sumStates = -HUGE_VAL;

        for (int l=0;l<nState;l++) {
            double sumThis = -HUGE_VAL;
            for (int k=0;k<nState;k++) {

                if (moveTo==0 || moveTo==3 || moveTo==4 || moveTo==5) {

                    sumThis = sumLogs(sumThis,pSite->fullFwdX(k) + hmm->probXX(k,l) + msr->indelX(l) + cSite->fullBwdX(l));
                    sumThis = sumLogs(sumThis,pSite->fullFwdY(k) + hmm->probYX(k,l) + msr->indelX(l) + cSite->fullBwdX(l));
                    sumThis = sumLogs(sumThis,pSite->fullFwdM(k) + hmm->probMX(k,l) + msr->indelX(l) + cSite->fullBwdX(l));

                } else if (moveTo==1 || moveTo==6 || moveTo==7 || moveTo==8) {

                    sumThis = sumLogs(sumThis,pSite->fullFwdX(k) + hmm->probXY(k,l) + msr->indelY(l) + cSite->fullBwdY(l));
                    sumThis = sumLogs(sumThis,pSite->fullFwdY(k) + hmm->probYY(k,l) + msr->indelY(l) + cSite->fullBwdY(l));
                    sumThis = sumLogs(sumThis,pSite->fullFwdM(k) + hmm->probMY(k,l) + msr->indelY(l) + cSite->fullBwdY(l));

                } else if (moveTo==2) {

                    sumThis = sumLogs(sumThis,pSite->fullFwdX(k) + hmm->probXM(k,l) + msr->fullM(l) + cSite->fullBwdM(l));
                    sumThis = sumLogs(sumThis,pSite->fullFwdY(k) + hmm->probYM(k,l) + msr->fullM(l) + cSite->fullBwdM(l));
                    sumThis = sumLogs(sumThis,pSite->fullFwdM(k) + hmm->probMM(k,l) + msr->fullM(l) + cSite->fullBwdM(l));
                }
            }

            sumStates = sumLogs(sumStates,sumThis);
        }

        cSite->postProb(exp(sumStates-fullScore));
    }

	delete pSite;
	delete cSite;
}

