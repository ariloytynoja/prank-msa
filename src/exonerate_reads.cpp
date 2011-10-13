#include "exonerate_reads.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <boost/regex.hpp>

using namespace std;

extern bool NOISE;

Exonerate_reads::Exonerate_reads()
{

}

bool Exonerate_reads::test_executable()
{
    int status = system("exonerate  >/dev/null");
    return WEXITSTATUS(status) == 1;
}

bool Exonerate_reads::split_sugar_string(const string& row,hit *h)
{

    const boost::regex pattern("sugar:\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+(\\S)\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+(\\S)\\s+(\\d+)\n");
    boost::match_results<string::const_iterator> result;
    bool valid = boost::regex_match(row, result, pattern);

    if(valid)
    {
        h->query    = result[1];
        h->q_start  = atoi( string(result[2]).c_str() );
        h->q_end    = atoi( string(result[3]).c_str() );
        h->q_strand = string(result[4]).at(0);

        h->node     = result[5];
        h->t_start  = atoi( string(result[6]).c_str() );
        h->t_end    = atoi( string(result[7]).c_str() );
        h->t_strand = string(result[8]).at(0);

        h->score    = atoi( string(result[9]).c_str() );
    }

    return valid;
}

bool Exonerate_reads::split_vulgar_string(const string& row,hit *h)
{

    const boost::regex pattern("vulgar:\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+(\\S)\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+(\\S)\\s+(\\d+)\\s+(.+)\n");
    boost::match_results<string::const_iterator> result;
    bool valid = boost::regex_match(row, result, pattern);

    if(valid)
    {
        h->query    = result[1];
        h->q_start  = atoi( string(result[2]).c_str() );
        h->q_end    = atoi( string(result[3]).c_str() );
        h->q_strand = string(result[4]).at(0);

        h->node     = result[5];
        h->t_start  = atoi( string(result[6]).c_str() );
        h->t_end    = atoi( string(result[7]).c_str() );
        h->t_strand = string(result[8]).at(0);

        h->score    = atoi( string(result[9]).c_str() );

    }



    return valid;
}


void Exonerate_reads::local_alignment(string* ls,string* rs, vector<hit> *hits, bool is_local)
{

    int r = rand();

    string left;
    for(int i=0;i<ls->length();i++)
        if(ls->at(i)!='-')
            left+=(ls->at(i));
    string right;
    for(int i=0;i<rs->length();i++)
        if(rs->at(i)!='-')
            right+=(rs->at(i));

    stringstream q_name;
    q_name <<"q"<<r<<".fas";

    ofstream q_output( q_name.str().c_str(), (ios::out));
    q_output<<">left"<<endl<<left<<endl;
    q_output.close();

    stringstream t_name;
    t_name <<"t"<<r<<".fas";

    ofstream t_output( t_name.str().c_str(), (ios::out));
    t_output<<">right"<<endl<<right<<endl;
    t_output.close();


    // exonerate command for local alignment

    stringstream command;
    if(is_local)
        command << "exonerate -q q"<<r<<".fas -t t"<<r<<".fas --showalignment no --showsugar yes --showvulgar no 2>&1";
    else
        command << "exonerate -q q"<<r<<".fas -t t"<<r<<".fas --showalignment no --showsugar yes --showvulgar no -m affine:local -E 2>&1";

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

        if(valid)
        {
            hits->push_back(h);
        }
    }
    pclose(fpipe);

    sort (hits->begin(), hits->end(), Exonerate_reads::q_earlier);

    vector<hit>::iterator iter1 = hits->begin();
    vector<hit>::iterator iter2 = hits->begin();
    if( iter2 != hits->end() )
        iter2++;

    while( iter2 != hits->end() )
    {
        if(iter1->t_start > iter2->t_start)
        {
            if(iter1->score > iter2->score)
            {
                hits->erase(iter2);
                iter2 = iter1;
                iter2++;
            } else {
                if(iter1 == hits->begin())
                {
                    hits->erase(iter1);
                    iter1 = iter2;
                    iter2++;
                }
                else
                {
                    hits->erase(iter1);
                    iter1 = iter2;
                    iter1--;
                }
            }
        }
        else
        {
            iter1++;
            iter2++;
        }
    }

    iter1 = hits->begin();
    iter2 = hits->begin();
    if( iter2 != hits->end() )
        iter2++;

    while( iter2 != hits->end() )
    {
        if(iter1->t_end > iter2->t_start)
        {
            int overlap = iter1->t_end - iter2->t_start;
            iter1->t_end-=overlap;
            iter1->q_end-=overlap;
            iter2->t_start+=overlap;
            iter2->q_start+=overlap;
        }
        if(iter1->q_end > iter2->q_start)
        {
            int overlap = iter1->q_end - iter2->q_start;
            iter1->t_end-=overlap;
            iter1->q_end-=overlap;
            iter2->t_start+=overlap;
            iter2->q_start+=overlap;
        }

        iter1++;
        iter2++;
    }

    if(NOISE==0)
        this->delete_files(r);
}

void Exonerate_reads::delete_files(int r)
{
    stringstream q_name;
    q_name <<"q"<<r<<".fas";

    stringstream t_name;
    t_name <<"t"<<r<<".fas";

    if( remove( q_name.str().c_str() ) != 0 )
       perror( "Error deleting file" );
    if( remove( t_name.str().c_str() ) != 0 )
       perror( "Error deleting file" );
}