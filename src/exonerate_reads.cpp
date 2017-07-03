#include "exonerate_reads.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include "config.h"
#include "translatesequences.h"
#include <unistd.h>

#if defined (__APPLE__)
#include <mach-o/dyld.h>
#endif

using namespace std;

Exonerate_reads::Exonerate_reads()
{

}

bool Exonerate_reads::test_executable()
{

    int status = -1;

    #if defined (__CYGWIN__)
    char path[200] = "";
    int length = readlink("/proc/self/exe",path,200-1);
	
    string epath = string(path).substr(0,length);
    if (epath.find("/")!=std::string::npos)
        epath = epath.substr(0,epath.rfind("/")+1);
    exoneratepath = epath;
    epath = epath+"exonerate.exe > /dev/null 2>/dev/null";
    status = system(epath.c_str());

    #else

    if(WEXITSTATUS(status) != 1)
    {
        char path[200] = "";
        string epath;

        #if defined (__APPLE__)
        uint32_t size = sizeof(path);
        _NSGetExecutablePath(path, &size);
        epath = string(path);
        if (epath.find("/")!=std::string::npos)
            epath = epath.substr(0,epath.rfind("/")+1);
        //epath = "DYLD_LIBRARY_PATH="+epath+" "+epath;

        #else
        int length = readlink("/proc/self/exe",path,200-1);
        epath = string(path).substr(0,length);
        if (epath.find("/")!=std::string::npos)
            epath = epath.substr(0,epath.rfind("/")+1);

        #endif
        
        exoneratepath = epath;
        epath = epath+"exonerate >/dev/null 2>/dev/null";

        status = system(epath.c_str());
    }
    #endif

    if(WEXITSTATUS(status) == 1) {
        if(NOISE>0)
            cout<<"Using Exonerate to anchor alignments. Use option '-noanchors' to disable.\n";
        return true;
    }

    exoneratepath = "";
    status = system("`exonerate  >/dev/null 2>/dev/null`");

    if(WEXITSTATUS(status) == 1 && NOISE>0)
        cout<<"Using Exonerate to anchor alignments. Use option '-noanchors' to disable.\n";

    return WEXITSTATUS(status) == 1;
}

bool Exonerate_reads::split_sugar_string(const string& row,hit *h)
{
    string method; string lname; int lst; int len; char lstr; string rname; int rst; int ren; char rstr; int score;

    istringstream str(row);
    str >> method >> lname >> lst >> len >> lstr >> rname >> rst >> ren >> rstr >> score;

    if(method != "sugar:" || score < 0)
        return false;

    h->query    = lname;
    h->q_start  = lst;
    h->q_end    = len;
    h->q_strand = lstr;

    h->node     = rname;
    h->t_start  = rst;
    h->t_end    = ren;
    h->t_strand = rstr;

    h->score    = score;

    return true;

}



void Exonerate_reads::local_alignment(string* ls,string* rs, vector<hit> *hits, bool is_local)
{

    ofstream q_output;
    ofstream t_output;

    stringstream q_name;
    stringstream t_name;

    int r = rand();
    while(true)
    {
        
        q_name <<tmp_dir<<"/q"<<r<<".fas";
        t_name <<tmp_dir<<"/t"<<r<<".fas";

        ifstream q_file(q_name.str().c_str());
        ifstream t_file(t_name.str().c_str());

        if(!q_file && !t_file)
        {
            q_output.open( q_name.str().c_str(), (ios::out) );
            t_output.open( t_name.str().c_str(), (ios::out) );

            break;
        }
        r = rand();
    }

    char ic = 'X';
    if(DNA)
        ic = 'N';

    string left;
    for (int i=0; i<ls->length(); i++)
        if (ls->at(i)!='-')
            left+=(ls->at(i));
        else
            left+= ic;
    string right;
    for (int i=0; i<rs->length(); i++)
        if (rs->at(i)!='-')
            right+=(rs->at(i));
        else
            right+= ic;

    if(CODON)
    {
        TranslateSequences ts;
        vector<string> names;
        names.push_back("left");
        names.push_back("right");
        vector<string> sequences;
        sequences.push_back(left);
        sequences.push_back(right);
        std::map<std::string,std::string> dnaSeqs;
        ts.translateProtein(&names,&sequences,&dnaSeqs);

        left = sequences.at(0);
        right = sequences.at(1);

    }


    q_output<<">left"<<endl<<left<<endl;
    q_output.close();


    t_output<<">right"<<endl<<right<<endl;
    t_output.close();


    // exonerate command for local alignment

    stringstream command;
    command<<exoneratepath<<"exonerate  ";
//    #if defined (__CYGWIN__)
//    char path[200];
//    int length = readlink("/proc/self/exe",path,200-1);
	
//    string epath = string(path).substr(0,length);
//    epath.replace(epath.rfind("prank"),string("prank").size(),string("exonerate.exe "));
//    command<<epath;
 
//    # else
//    command<<"exonerate  ";
//    #endif
    command << " -q " <<tmp_dir<<"/q"<<r<<".fas -t "<<tmp_dir<<"/t"<<r<<".fas --showalignment no --showsugar yes --showvulgar no 2>&1";
    if(NOISE>0)
        cout<<"cmd: "<<command.str()<<endl;

    FILE *fpipe;
    if ( !(fpipe = (FILE*)popen(command.str().c_str(),"r")) )
    {
        perror("Problems with exonerate pipe.\nExiting.\n");
        exit(1);
    }

    // read exonerate output, summing the multiple hit scores

    char line[256];

    while ( fgets( line, sizeof line, fpipe))
    {
        hit h;
        bool valid = split_sugar_string(string(line),&h);

        if (valid)
        {
            hits->push_back(h);
        }
    }
    pclose(fpipe);

    if(dropRevAnch)
    {
        for(vector<hit>::iterator iter1 = hits->begin();iter1!=hits->end();)
        {
            if(iter1->t_strand == '-' || iter1->q_strand == '-')
                hits->erase(iter1);
            else
                iter1++;
        }
    }

    sort (hits->begin(), hits->end(), Exonerate_reads::q_earlier);

    vector<hit>::iterator iter1 = hits->begin();
    vector<hit>::iterator iter2 = hits->begin();
    if ( iter2 != hits->end() )
        iter2++;

    while ( iter2 != hits->end() )
    {
        if (iter1->t_start > iter2->t_start)
        {
            if (iter1->score > iter2->score)
            {
                hits->erase(iter2);

                iter1 = hits->begin();
                iter2 = hits->begin();
                if ( iter2 != hits->end() )
                    iter2++;

                continue;
            }
            else
            {
                if (iter1 == hits->begin())
                {
                    hits->erase(iter1);

                    iter1 = hits->begin();
                    iter2 = hits->begin();
                    if ( iter2 != hits->end() )
                        iter2++;

                    continue;
                }
                else
                {
                    hits->erase(iter1);

                    iter1 = hits->begin();
                    iter2 = hits->begin();
                    if ( iter2 != hits->end() )
                        iter2++;

                    continue;
                }
            }
        }
        else if (iter1->q_end > iter2->q_end)
        {
            hits->erase(iter2);

            iter1 = hits->begin();
            iter2 = hits->begin();
            if ( iter2 != hits->end() )
                iter2++;

            continue;
        }
        else
        {
            iter1++;
            iter2++;
        }
    }

    iter1 = hits->begin();
    iter2 = hits->begin();
    if ( iter2 != hits->end() )
        iter2++;

    while ( iter2 != hits->end() )
    {
        if (iter1->t_end > iter2->t_start)
        {
            int overlap = iter1->t_end - iter2->t_start;

            if(float(iter1->score/(iter1->t_end - iter1->t_start)) <
                    float(iter2->score/(iter2->t_end - iter2->t_start)))
            {
                iter1->score *= int ( float((iter1->t_end - iter1->t_start -overlap))/float((iter1->t_end - iter1->t_start)) );

                iter1->t_end-=overlap;
                iter1->q_end-=overlap;

                if(iter1->t_end <= iter1->t_start)
                {
                    hits->erase(iter1);

                    iter1 = hits->begin();
                    iter2 = hits->begin();
                    if ( iter2 != hits->end() )
                        iter2++;

                    continue;
                }
            }
            else
            {
                iter2->score *= int( float((iter2->t_end - iter2->t_start -overlap))/float((iter2->t_end - iter2->t_start)) );

                iter2->t_start+=overlap;
                iter2->q_start+=overlap;

                if(iter2->t_end <= iter2->t_start)
                {
                    hits->erase(iter2);

                    iter1 = hits->begin();
                    iter2 = hits->begin();
                    if ( iter2 != hits->end() )
                        iter2++;

                    continue;
                }
            }

        }
        if (iter1->q_end > iter2->q_start)
        {

            int overlap = iter1->q_end - iter2->q_start;
            if(float(iter1->score/(iter1->q_end - iter1->q_start)) <
                    float(iter2->score/(iter2->q_end - iter2->q_start)))
            {
                iter1->score *= int( float((iter1->t_end - iter1->t_start -overlap))/float((iter1->t_end - iter1->t_start)) );

                iter1->t_end-=overlap;
                iter1->q_end-=overlap;

                if(iter1->t_end <= iter1->t_start)
                {
                    hits->erase(iter1);

                    iter1 = hits->begin();
                    iter2 = hits->begin();
                    if ( iter2 != hits->end() )
                        iter2++;

                    continue;
                }
            }
            else
            {
                iter2->score *= int( float((iter2->t_end - iter2->t_start -overlap))/float((iter2->t_end - iter2->t_start)) );

                iter2->t_start+=overlap;
                iter2->q_start+=overlap;

                if(iter2->t_end <= iter2->t_start)
                {
                    hits->erase(iter2);

                    iter1 = hits->begin();
                    iter2 = hits->begin();
                    if ( iter2 != hits->end() )
                        iter2++;

                    continue;
                }
            }
        }


        iter1++;
        iter2++;
    }

    //delete files 
    if( remove( q_name.str().c_str() ) != 0 )
      perror("Error deleting temporary file in Exonerate_reads::local_alignment");
    if( remove( t_name.str().c_str() ) != 0 )
	perror("Error deleting temporary file in Exonerate_reads::local_alignment");

}
