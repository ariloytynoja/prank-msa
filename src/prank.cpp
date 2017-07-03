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
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <unistd.h>
#include "progressivealignment.h"
#include "check_version.h"
#include "prank.h"

using namespace std;

// stores current temporary directory will be created and removed in main
char tmp_dir[1000] = "";
bool verbose = false;

//main
int main(int argc, char *argv[])
{
    version = 170427;

    readArguments(argc, argv);
    int time1 = time(0);

    // character array for assigning system tmp path
    char* tmpPath = NULL;  
    //mkdtemp template
    const char* mktemplate = "tmpdirprankmsaXXXXXX";

     //create temporary directory
     //get tmpPath from environment or assign tmpPath to NULL (if last getenv fails)
     if( (tmpPath = getenv("TMPDIR")) == NULL)
       {
       if( (tmpPath = getenv("TMP")) == NULL)
         {     
           if( (tmpPath = getenv("TEMPDIR")) == NULL)
             { 
               tmpPath = getenv("TEMP");
             } 
         }
       }
   
     // define temp dir based on whether tmpPath was found   
     if(tmpPath == NULL)
       {
       sprintf(tmp_dir, "%s", mktemplate);
       }
     else
       {
       sprintf(tmp_dir, "%s/%s", tmpPath, mktemplate); 
       }
   
     // call mkdtemp
     if( (mkdtemp(tmp_dir) == NULL) )
       {
       perror("'mkdtemp' failed to generate temporary directory while creating ProgressiveAlignbment object.\n");
       exit(EXIT_FAILURE);
       }


    ProgressiveAlignment* pa = new ProgressiveAlignment(treefile,seqfile,dnafile);
    if (NOISE>=0)
        cout<<endl<<"Analysis done. Total time "<<(time(0)-time1)<<"s"<<endl;

    delete pa;
    delete hmm;

    //unlink temporary directory
    if( rmdir(tmp_dir) != 0)
    {
        perror("ERROR! failed removing temporary directory.");
    }

    cout<<endl;
    exit(0);
}


void readArguments(int argc, char *argv[])
{

    // first see if 'noise' is defined.
    for (int i=1; i<argc; i++)
    {

        string s = string(argv[i]);

        // no debugging info
        if (s=="-quiet")
        {
            NOISE = -1;
            SCREEN = false;
        }
        // debugging info level
        else if (s.substr(0,7)=="-noise=")
        {
            NOISE = atoi(string(argv[i]).substr(7).c_str());
            SCREEN = false;
        }
    }

    // one argument only; should be sequence file!
    if (argc==2 && string(argv[1]).substr(0,1)!="-")
    {
        seqfile = string(argv[1]);
    }

    // otheriwise run in the normal manner
    else
    {
        for (int i=1; i<argc; i++)
        {

            string s = string(argv[i]);

            if (s=="-quiet")
            {
                // handled earlier
            }

            // debugging info level
            else if (s.substr(0,7)=="-noise=")
            {
                NOISE = atoi(string(argv[i]).substr(7).c_str());
                SCREEN = false;
            }

            // print help
            else if (s=="-help")
            {
                printHelp(true);
                exit(0);
            }

            // check version
            else if (s=="-version")
            {
                Check_version ch(version);
                exit(0);
            }


	    else if (s=="-verbose")
            {
	      verbose = true;
            }


            /********* input/output: **********/

            // sequence data file
            else if (s.substr(0,3)=="-d=")
            {
                seqfile = string(argv[i]).substr(3);
            }
            else if (s.substr(0,4)=="-d1=")
            {
                seqfile1 = string(argv[i]).substr(4);
                MERGE = true;
                PARTLYALIGNED = true;
            }
            else if (s.substr(0,4)=="-d2=")
            {
                seqfile2 = string(argv[i]).substr(4);
                MERGE = true;
                PARTLYALIGNED = true;
            }

            // guide tree file
            else if (s.substr(0,3)=="-t=")
            {
                treefile = string(argv[i]).substr(3);
            }
            else if (s.substr(0,4)=="-ot=")
            {
                oldtreefile = string(argv[i]).substr(4);
            }
            else if (s.substr(0,4)=="-t1=")
            {
                treefile1 = string(argv[i]).substr(4);
            }
            else if (s.substr(0,4)=="-t2=")
            {
                treefile2 = string(argv[i]).substr(4);
            }

            // alignment output file
            else if (s.substr(0,3)=="-o=")
            {
                outfile = string(argv[i]).substr(3);
            }

            // alignment output file
            else if (s.substr(0,5)=="-tmp=")
            {
                tempdir = string(argv[i]).substr(5);
            }

            // structure model file
            else if (s.substr(0,3)=="-m=")
            {
                hmmname = string(argv[i]).substr(3);
                hmm = new HMModel();
                hmm->readModel(hmmname.c_str());
                HASHMM = true;
            }

            // guide tree as a string
            else if (s.substr(0,6)=="-tree=")
            {
                treefile = string(argv[i]).substr(6).c_str();
                TREESTRING = true;
            }

            // mixture of existing and new alignments for Ziheng
            else if (s=="-partaligned")
            {
                PARTLYALIGNED = true;
            }

            // pre-aligned data, just reconstruct it (using a model if specified)
            else if (s=="-e")
            {
                PREALIGNED = true;
            }
            else if (s=="-keep")
            {
                PREALIGNED = true;
                PRINTSCOREONLY = false;
            }

            else if (s=="-score")
            {
                PREALIGNED = true;
                PRINTSCOREONLY = true;
            }

            else if (s=="-update")
            {
                UPDATE = true;
            }

            else if (s=="-realign")
            {
                UPDATESECOND = false;
            }

            else if (s.substr(0,13)=="-updatelimit=")
            {
                updateTolerance = atof(s.substr(13).c_str());
            }

            // backtranslate existing protein alignment to DNA
            else if (s.substr(0,5)=="-dna=")
            {
                dnafile = s.substr(5);
                BACKTRANSLATE = true;
            }

            /********* more input/output: **********/

            // do not estimate guide tree from mafft alignment
            else if (s=="-nomafft")
            {
                MAFFTALIGNMENT = false;
            }


            // compute score for mafft alignment
            else if (s=="-scoremafft")
            {
                SCOREMAFFT = true;
            }

            // estimate guide tree from input alignment (before realignment)
            else if (s=="-njtree")
            {
                TREEFROMALIGNMENT = true;
            }

            else if (s=="-treeonly")
            {
                TREEONLY = true;
            }

            // output alignment format
            else if (s.substr(0,11)=="-outformat=")
            {
                string tmp = string(argv[i]).substr(11);
                format = parseFormat(tmp);
            }

            // do not estimate ancestors with bppancestor
            else if (s=="-nobppa")
            {
                BPPANCESTORS = false;
            }

            // output alignment format
            else if (s.substr(0,3)=="-f=")
            {
                string tmp = string(argv[i]).substr(3);
                format = parseFormat(tmp);
            }

            // do output backtabs
            else if (s=="-cute")
            {
                SCREEN = true;
            }

            // reporting interval
            else if (s.substr(0,4)=="-rl=")
            {
                reportLimit = atoi(string(argv[i]).substr(4).c_str());
            }

            // write reconstructed ancestral seqs
            else if (s=="-showanc")
            {
                WRITEANCSEQ = true;
            }

            // write evolutionary events
            else if (s=="-showevents")
            {
                LISTEVENTS = true;
            }

            // write all iteration
            else if (s=="-showiter")
            {
                WRITEITER = true;
            }

            // write everything
            else if (s=="-showall")
            {
                WRITEANCSEQ = true;
                LISTEVENTS = true;
                PRINTTREE = true;
                WRITEXML = true;
//                WRITEITER = true;
            }

            // compute parsimony score
            else if (s=="-noscore")
            {
                PARSIMONYSCORE = false;
            }

            // parsimony score for indels
            else if (s.substr(0,12)=="-indelscore=")
            {
                INDELSCORE = string(argv[i]).substr(12);
            }

            // write ancestral nodes as they are solved
            else if (s=="-printnodes")
            {
                PRINTNODES = true;
            }

            // don't print tree
            else if (s=="-showtree")
            {
                PRINTTREE = true;
            }

            // don't write xml
            else if (s=="-showxml")
            {
                WRITEXML = true;
            }

            // print dots for insertions
            else if (s=="-dots" || s=="-esko")
            {
                DOTS = true;
            }

            // no align, convert only
            else if (s=="-convert")
            {
                CONVERT = true;
            }

            // use short names (until first space)
            else if (s=="-shortnames")
            {
                SHORTNAMES = true;
            }

            /********* model: gaps and F **********/

            // keep insertion forever
            else if (s=="-F" || s=="+F")
            {
                FOREVER = true;
                SKIPGAPANCH = true;
            }

            // not needed but still allowed option
            else if (s=="-no-F")
            {
                FOREVER = false;
            }

            // old implementation
            else if (s=="-F_old")
            {
                FOREVER_OLD = true;
                SKIPGAPANCH = true;
            }

            /********* model: substitutions, indels **********/

            else if (s.substr(0,10)=="-dnafreqs=")
            {
                dnaFreqs = string(argv[i]).substr(10);
            }

            else if (s=="-jc")
            {
                dnaFreqs = "1,1,1,1";
            }

            else if (s.substr(0,9)=="-gaprate=")
            {
                gapRate = atof(string(argv[i]).substr(9).c_str());
            }

            else if (s.substr(0,8)=="-gapext=")
            {
                gapExt = atof(string(argv[i]).substr(8).c_str());
            }

            else if (s.substr(0,7)=="-kappa=")
            {
                kappa = atof(string(argv[i]).substr(7).c_str());
            }

            else if (s.substr(0,5)=="-rho=")
            {
                rho = atof(string(argv[i]).substr(5).c_str());
            }

            /********* model: other **********/

            // codon alignment
            else if (s=="-codon")
            {
                CODON = true;
            }

            // force dna alignment
            else if (s=="-DNA")
            {
                DNA = true;
            }

            // force protein alignment
            else if (s=="-protein")
            {
                PROTEIN = true;
            }

            // no posterior probabiliity calculation
            else if (s=="-support")
            {
                DOPOST = true;
            }

            // penalise terminal gaps
            else if (s=="-termgap")
            {
                NOTGAP = false;
            }

            else if (s=="-nomissing")
            {
                TERMF = true;
            }

            // run once
            else if (s=="-once")
            {
                iterations = 1;
            }

            // run twice
            else if (s=="-twice")
            {
                iterations = 2;
            }

            // run many times
            else if (s.substr(0,9)=="-iterate=")
            {
                iterations = atoi(s.substr(9).c_str());
            }

            // prune the tree
            else if (s=="-prunetree")
            {
                PRUNETREE = true;
            }

            // prune the data
            else if (s=="-prunedata")
            {
                PRUNEDATA = true;
            }

            // use log values (slightly slower)
            else if (s=="-uselogs")
            {
                LOGVALUES = true;
            }

            // use log values (slightly slower)
            else if (s=="-nologs")
            {
                LOGVALUES = false;
            }

            // seed for random number generator
            else if (s.substr(0,6)=="-seed=")
            {
                rnd_seed = atoi(string(argv[i]).substr(6).c_str());
            }

            else if (s=="-reproducible")
            {
                REPRODUCIBLE = true;
            }


            /********* more model: **********/

            // translate DNA to protein, then backtranslate
            else if (s=="-translate")
            {
                TRANSLATE = true;
            }

            // translate mtDNA to protein, then backtranslate
            else if (s=="-mttranslate")
            {
                TRANSLATE = true;
                MTTABLE = true;
            }

            // consider N or X identical to any
            else if (s=="-NX")
            {
                NXis1 = true;
            }

            // split probbailities for N and X
            else if (s=="-splitNX")
            {
                NXis1 = false;
            }


            /********* more model: pairwise alignment for guide tree **********/

            // expected pairwise distance
            else if (s.substr(0,8)=="-pwdist=")
            {
                pwDist = atof(string(argv[i]).substr(8).c_str());;
            }

            // expected pairwise distance
            else if (s.substr(0,11)=="-pwdnadist=")
            {
                pwDnaDist = atof(string(argv[i]).substr(11).c_str());;
            }

            /********* more model: branch lengths in guide tree **********/

            // scale branch lengths
            else if (s.substr(0,15)=="-scalebranches=")
            {
                branchScalingFactor = atof(string(argv[i]).substr(15).c_str());;
            }

            // set branch lengths
            else if (s.substr(0,15)=="-fixedbranches=")
            {
                fixedBranchLength = atof(string(argv[i]).substr(15).c_str());;
                FIXEDBRANCH = true;
            }

            // set merge branch length
            else if (s.substr(0,11)=="-mergedist=")
            {
                mergeBranchLength = atof(string(argv[i]).substr(11).c_str());;
            }

            // set branch lengths
            else if (s.substr(0,13)=="-maxbranches=")
            {
                fixedBranchLength = atof(string(argv[i]).substr(13).c_str());;
                MAXBRANCH =true;
            }

            // set branch lengths
            else if (s.substr(0,13)=="-maxpairdist=")
            {
                dnaMaxPairwiseLength = atof(string(argv[i]).substr(13).c_str());
                protMaxPairwiseLength = atof(string(argv[i]).substr(13).c_str());
            }

            // use real guidetree distances
            else if (s=="-adjustmodel")
            {
                ADJUSTMODEL = true;
            }

            // use real guidetree distances
            else if (s=="-noadjustmodel")
            {
                ADJUSTMODEL = false;
            }

            // use real guidetree distances
            else if (s=="-realbranches")
            {
                REALBRANCHES = true;
            }

            // correct guidetree distances
            else if (s=="-correctp")
            {
                CORRECTP = true;
            }

            // penalise gaps in NJ distances
            else if (s=="-njgaps")
            {
                PENALISEGAPS = true;
            }

            /********* technical: hirschberg, full probability **********/

            // "band" full probability (less memory)
            else if (s=="-fb")
            {
                FULLBAND = true;
            }

            // complete full probability
            else if (s=="-ff")
            {
                FULLFULL = true;
            }

            // hirschberg band width (for hirschbergalignment)
            else if (s.substr(0,5)=="-hbw=")
            {
                HBW = atoi(string(argv[i]).substr(5).c_str());
            }

            // full probability band width (for fullprobability)
            else if (s.substr(0,5)=="-fbw=")
            {
                FBW = atoi(string(argv[i]).substr(5).c_str());
            }

            // skip insertions in postprobs
            else if (s=="-skipins")
            {
                SKIPINS = true;
            }

            /********* technical: anchoring **********/

            // use anchors
            else if (s=="-noanchors" || s=="-noa")
            {
                EXONERATE = false;
            }

            // anchor skip distance
            else if (s.substr(0,12)=="-anchorskip=")
            {
                anchSkipDist = atoi(string(argv[i]).substr(12).c_str());
            }

            // minimum anchor distance
            else if (s.substr(0,6)=="-mind=")
            {
                minAnchDist = atoi(string(argv[i]).substr(6).c_str());
            }

            // anchor skip distance
            else if (s.substr(0,7)=="-skipd=")
            {
                anchSkipDist = atoi(string(argv[i]).substr(7).c_str());
            }

            // anchor drop distance
            else if (s.substr(0,7)=="-dropd=")
            {
                anchDropDist = atoi(string(argv[i]).substr(7).c_str());
            }

            // ignore reverse anchors
            else if (s=="-droprevanch")
            {
                dropRevAnch = true;
            }

            // don't infer gaps caused by missing data
            else if (s=="-nopatchdata")
            {
                PATCHMISSING = false;
            }

            // length of gap deemed as missing data
            else if (s.substr(0,11)=="-misslimit=")
            {
                missingLimit = atoi(string(argv[i]).substr(11).c_str());
            }

            // skip gaps in anchoring ancestral seqs (?)
            else if (s=="-gapanch")
            {
                if (!FOREVER)
                    SKIPGAPANCH = false;
                else
                    cout<<endl<<"unaccepted combination: -F -gapanch; disabling -gapanch"<<endl;
            }

            /********* technical: memory & speed efficiency **********/

            // matrix resize factor
            else if (s.substr(0,11)=="-matresize=")
            {
                resizeFactor = atof(string(argv[i]).substr(11).c_str());
            }

            // matrix initial factor
            else if (s.substr(0,13)=="-matinitsize=")
            {
                initialMatrixSize = atof(string(argv[i]).substr(13).c_str());
            }

            // use pwmatrix maximum size
            else if (s=="-pwmatmax")
            {
                PWMATRIXMAXSIZE = true;
            }

            // use pwmatrix maximum size
            else if (s=="-longseq")
            {
                PWMATRIXMAXSIZE = false;
            }


            /************************************************/

            else
            {
                cout<<"Unknown option: "<<s<<endl<<endl;
                printHelp(false);
                exit(0);
            }
        }
    }

    if (seqfile=="" && (seqfile1=="" || seqfile2==""))
    {
        printHelp(false);
        exit(0);
    }

    // define a seed for random numbers
    if (rnd_seed>0)
        srand(rnd_seed);
    else
        srand(time(0));


    if (format!=8 && format!=11 && format!=12 && format!=17 && format!=18 && format!=19)
        format = 8;

}


void printHelp(bool complete)
{
    cout<<endl<<"prank v."<<version<<". ";
    cout<<"Minimal usage: 'prank sequence_file'"<<endl<<endl;;
    cout<<"Advanced usage: 'prank [optional parameters] -d=sequence_file [optional parameters]'"<<endl;;
    cout<<"\n input/output parameters:"<<endl;
    cout<<"  -d=sequence_file (in FASTA format)"<<endl;
    cout<<"  -t=tree_file [default: no tree, generate approximate NJ tree]"<<endl;
    if (complete)
        cout<<"  -tree=\"tree_string\" [tree in newick format; in double quotes]"<<endl;
    cout<<"  -o=output_file [default: 'output']"<<endl;
    cout<<"  -f=output_format ['fasta' (default), 'phylipi', 'phylips', 'paml', 'nexus']"<<endl;

    cout<<"  -showxml [output xml-files]"<<endl;
    cout<<"  -showtree [output dnd-files]"<<endl;
    cout<<"  -showanc [output ancestral sequences]"<<endl;
    cout<<"  -showevents [output evolutioanry events]"<<endl;
    cout<<"  -showall [output all of these]"<<endl;
    if (complete)
    {
        cout<<"  -noanchors [no Exonerate anchoring]"<<endl;
        cout<<"  -nomafft [no MAFFT guide tree]"<<endl;
        cout<<"  -scoremafft [score also MAFFT alignment]"<<endl;
    }
    cout<<"  -support [compute posterior support]"<<endl;
    cout<<"  -njtree [estimate tree from input alignment (and realign)]"<<endl;
    cout<<"  -treeonly [estimate tree only]"<<endl;
    cout<<"  -quiet"<<endl;
    if (complete)
    {
        cout<<"\n alignment merge parameters:"<<endl;
        cout<<"  -d1=sequence_file_1 (in FASTA format)"<<endl;
        cout<<"  -d2=sequence_file-2 (in FASTA format)"<<endl;
        cout<<"  -t1=tree_file_1 [if not provided, generate NJ tree]"<<endl;
        cout<<"  -t2=tree_file_2 [if not provided, generate NJ tree]"<<endl;
        cout<<"  -mergedist=# [if no tree provided; default: "<<defaultBranchLength<<"]"<<endl;
    }
    cout<<"\n model parameters:"<<endl;
    cout<<"  +F or -F [force insertions to be always skipped]"<<endl;
    if (complete)
        cout<<"  -dots [show insertion gaps as dots]"<<endl;
    if (complete)
        cout<<"  -m=model_file [default: HKY2/WAG]"<<endl;
    cout<<"  -gaprate=# [gap opening rate; default: dna "<<dnaGapRate<<" / prot "<<protGapRate<<"]"<<endl;
    cout<<"  -gapext=# [gap extension probability; default: dna "<<dnaGapExt<<" / prot "<<protGapExt<<"]"<<endl;
    if (complete)
    {
        cout<<"  -dnafreqs=#,#,#,# [ACGT; default: empirical]"<<endl;
        cout<<"  -kappa=# [ts/tv rate ratio; default:"<<kappa<<"]"<<endl;
        cout<<"  -rho=# [pur/pyr rate ratio; default:"<<rho<<"]"<<endl;
        cout<<"  -indelscore=#,#,#,# [1,2,3,>3; indel penalties for alignment score]"<<endl;
    }
    cout<<"  -codon [for coding DNA: use empirical codon model]"<<endl;
    cout<<"  -DNA / -protein [no autodetection: use dna or protein model]"<<endl;
    cout<<"  -termgap [penalise terminal gaps normally]"<<endl;
    cout<<"  -nomissing [no missing data, use -F for terminal gaps ]"<<endl;
    cout<<"\n other parameters:"<<endl;
    cout<<"  -keep [keep alignment \"as is\" (e.g. for ancestor inference)]"<<endl;
    if (complete)
        cout<<"  -pwdist=# [expected pairwise distance for computing guide tree; default: dna "<<pwDnaDist<<" / prot "<<pwProtDist<<"]"<<endl;
    cout<<"  -iterate=# [rounds of re-alignment iteration]"<<endl;
    cout<<"  -once [run only once; same as -iterate=1]"<<endl;
    cout<<"  -prunetree [prune guide tree branches with no sequence data]"<<endl;
    cout<<"  -prunedata [prune sequence data with no guide tree leaves]"<<endl;

    if (complete)
        cout<<"  -skipins [skip insertions in posterior support]"<<endl;
    cout<<"  -uselogs [slower but should work for a greater number of sequences]"<<endl;
    if (complete)
    {
        cout<<"  -shortnames [truncate names at first space]"<<endl;
        cout<<"  -anchorskip=# [min. sequence length for anchoring; default "<<anchSkipDist<<"]"<<endl;
        cout<<"  -scalebranches=# [scale branch lengths; default: dna "<<dnaBranchScalingFactor<<" / prot "<<protBranchScalingFactor<<"]"<<endl;
        cout<<"  -fixedbranches=# [use fixed branch lengths]"<<endl;
        cout<<"  -maxbranches=# [set maximum branch length]"<<endl;
        cout<<"  -realbranches [disable branch length truncation]"<<endl;
        cout<<"  -seed=# [set random number seed]"<<endl;
    }
    cout<<"  -translate [translate to protein]"<<endl;
    cout<<"  -mttranslate [translate to protein using mt table]"<<endl;
//    if (complete)
//        cout<<"  -maxpairdist=# [maximum pairwise distance for matrix computation]"<<endl;

    cout<<"\n other:"<<endl;
    cout<<"  -convert [no alignment, just convert to another format]"<<endl;
    if (complete)
        cout<<"  -dna=dna_sequence_file [DNA sequence file for backtranslation of protein alignment]"<<endl;
    cout<<"  -version [check for updates]"<<endl;
    cout<<"  -verbose [print progress etc. during runtime]"<<endl;
    cout<<"\n  -help [show more options]"<<endl;

    cout<<""<<endl;
}

int parseFormat(string Format)
{
    string format = "";
    for (unsigned int i = 0; i < Format.size(); i++)
    {
        format += tolower(Format[i]);
    }

    if (format == "fasta" || format == "8")
        return 8;
    else if (format == "phylipi" || format == "phylip" || format == "12")
        return 12;
    else if (format == "phylips" || format == "11")
        return 11;
    else if (format == "nexus" || format == "paup" || format == "17")
        return 17;
    else if (format == "paml")
        return 18;
    else if (format == "raxml")
        return 19;
    else
    {
        cout<<"Warning: output format not recognized, using FASTA.\n";
        return 8;
    }
}
