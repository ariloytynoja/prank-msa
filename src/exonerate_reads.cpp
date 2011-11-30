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

using namespace std;


Exonerate_reads::Exonerate_reads()
{

}

bool Exonerate_reads::test_executable()
{
    int status = system("exonerate  >/dev/null 2>/dev/null");
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

    int r = rand();

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

        ts.translateProtein(&names,&sequences);

        left = sequences.at(0);
        right = sequences.at(1);

    }

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
    command << "exonerate -q q"<<r<<".fas -t t"<<r<<".fas --showalignment no --showsugar yes --showvulgar no 2>&1";

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
                iter1->score *= float((iter1->t_end - iter1->t_start -overlap))/float((iter1->t_end - iter1->t_start));

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
                iter2->score *= float((iter2->t_end - iter2->t_start -overlap))/float((iter2->t_end - iter2->t_start));

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
                iter1->score *= float((iter1->t_end - iter1->t_start -overlap))/float((iter1->t_end - iter1->t_start));

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
                iter2->score *= float((iter2->t_end - iter2->t_start -overlap))/float((iter2->t_end - iter2->t_start));

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

    this->delete_files(r);
}

void Exonerate_reads::delete_files(int r)
{
    stringstream q_name;
    q_name <<"q"<<r<<".fas";

    stringstream t_name;
    t_name <<"t"<<r<<".fas";

    if ( remove( q_name.str().c_str() ) != 0 )
        perror( "Error deleting file" );
    if ( remove( t_name.str().c_str() ) != 0 )
        perror( "Error deleting file" );
}
