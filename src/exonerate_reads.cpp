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
using namespace ppa;

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

//        string match = result[10];

//        std::string::const_iterator start, end;
//        start = match.begin();
//        end = match.end();

//        const boost::regex triple_pattern("([MG])\\s+(\\d+)\\s+(\\d+)");
//        boost::match_results<std::string::const_iterator> what;
//        boost::match_flag_type flags = boost::match_default;
//        while(regex_search(start, end, what, triple_pattern, flags))
//        {
//            cout<<what[0]<<endl;
//            start = what[0].second;
//            // update flags:
//            flags |= boost::match_prev_avail;
//            flags |= boost::match_not_bob;
//        }
    }



    return valid;
}

void Exonerate_reads::write_exonerate_input(Node *root, Fasta_entry *read, map<string,string> *names, int r)
{
    vector<Fasta_entry> aligned_sequences;
    root->get_alignment(&aligned_sequences,true);

    // create exonerate input
    stringstream q_name;
    q_name <<"q"<<r<<".fas";

    ofstream q_output( q_name.str().c_str(), (ios::out));
    q_output<<">"<<read->name<<endl<<read->sequence<<endl;
    q_output.close();

    stringstream t_name;
    t_name <<"t"<<r<<".fas";

    ofstream t_output( t_name.str().c_str(), (ios::out));
    vector<Fasta_entry>::iterator it = aligned_sequences.begin();
    for(;it!=aligned_sequences.end();it++)
    {
        if(names->find(it->name) != names->end())
        {
            string seq = it->sequence;
            for (string::iterator si = seq.begin();si != seq.end();)
                if(*si == '-')
                    seq.erase(si);
                else
                    si++;
            t_output<<">"<<it->name<<endl<<seq<<endl;
        }
    }
    t_output.close();

    return;
}

void Exonerate_reads::local_alignment(Node *root, Fasta_entry *read, multimap<string,string> *tid_nodes, map<string,hit> *hits, bool is_local)
{

    int r = rand();

//    set<string> names;
//    multimap<string,string>::iterator it = tid_nodes->begin();
//    for(;it!=tid_nodes->end();it++)
//        names.insert(it->second);

    map<string,string> names;
    multimap<string,string>::iterator it = tid_nodes->begin();
    for(;it!=tid_nodes->end();it++)
    {
        if(it->first==read->tid)
            names.insert(pair<string,string>(it->second,it->first));
    }

    this->write_exonerate_input(root,read,&names,r);

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
    map<string,hit> all_hits;
    vector<string> hit_names;

    while ( fgets( line, sizeof line, fpipe))
    {
        hit h;
        bool valid = split_sugar_string(string(line),&h);

        if(valid)
        {
            map<string,hit>::iterator iter = all_hits.find(h.node);
            if( iter != all_hits.end() )
            {
                if(iter->second.t_strand == h.t_strand && iter->second.q_strand == h.q_strand)
                {
                    iter->second.score += h.score;

                    if(iter->second.q_start > h.q_start)
                        iter->second.q_start = h.q_start;
                    if(iter->second.q_end < h.q_end)
                        iter->second.q_end = h.q_end;
                    if(iter->second.t_start > h.t_start)
                        iter->second.t_start = h.t_start;
                    if(iter->second.t_end < h.t_end)
                        iter->second.t_end = h.t_end;
//                    cout<<"i "<<h.query<<" "<<h.node<<" "<<h.score<<" "<<h.q_start<<" "<<h.q_end<<" "<<h.q_strand<<" "<<h.t_start<<" "<<h.t_end<<" "<<h.t_strand<<"\n";
                }
                else if(iter->second.score < h.score)
                {
                    iter->second = h;
//                    cout<<"b "<<h.query<<" "<<h.node<<" "<<h.score<<" "<<h.q_start<<" "<<h.q_end<<" "<<h.q_strand<<" "<<h.t_start<<" "<<h.t_end<<" "<<h.t_strand<<"\n";
                }
            }
            else
            {
//                cout<<"n "<<h.query<<" "<<h.node<<" "<<h.score<<" "<<h.q_start<<" "<<h.q_end<<" "<<h.q_strand<<" "<<h.t_start<<" "<<h.t_end<<" "<<h.t_strand<<"\n";

                all_hits.insert( make_pair(h.node, h) );
                hit_names.push_back(h.node);
            }
        }
    }
    pclose(fpipe);


    if(Settings::noise>1)
        cout<<"\nExonerate_reads: "<<read->name<<" has "<<hit_names.size()<<" hits\n";

    if(hit_names.size()>0)
    {
        tid_nodes->clear();
        hits->clear();

        vector<hit> best_hits;
        vector<string>::iterator iter = hit_names.begin();

        for(;iter!=hit_names.end();iter++)
        {
             map<string,hit>::iterator iter2 = all_hits.find(*iter);
             if( iter2 != all_hits.end() )
             {
               best_hits.push_back(iter2->second);
             }
        }

        sort (best_hits.begin(), best_hits.end(), Exonerate_reads::better);


        // keep hits that are above a relative threshold

        if( ( is_local && Settings_handle::st.is("exonerate-local-keep-above") &&
                Settings_handle::st.get("exonerate-local-keep-above").as<float>()>0 ) ||

            ( !is_local && Settings_handle::st.is("exonerate-gapped-keep-above") &&
                Settings_handle::st.get("exonerate-gapped-keep-above").as<float>()>0 ) )
        {
            int lim = best_hits.at(0).score;
            if(is_local)
                lim = int (lim * Settings_handle::st.get("exonerate-local-keep-above").as<float>() );
            else
                lim = int (lim * Settings_handle::st.get("exonerate-gapped-keep-above").as<float>() );


            for(int i=0; i<(int)hit_names.size(); i++)
            {
                if(best_hits.at(i).score > lim)
                {
                    string tid = names.find(best_hits.at(i).node)->second;
                    tid_nodes->insert(pair<string,string>(tid,best_hits.at(i).node));
                    hits->insert(pair<string,hit>(best_hits.at(i).node,best_hits.at(i)));

                    if(Settings::noise>2)
                        cout<<"adding "<<best_hits.at(i).node<<" "<<best_hits.at(i).score<<endl;
                }
            }
        }


        // keep a fixed number of hits

        else
        {
            int lim = hit_names.size();
            if( is_local && Settings_handle::st.is("exonerate-local-keep-best") &&
                    Settings_handle::st.get("exonerate-local-keep-best").as<int>()>0 )
                lim = Settings_handle::st.get("exonerate-local-keep-best").as<int>();

            if( !is_local && Settings_handle::st.is("exonerate-gapped-keep-best") &&
                    Settings_handle::st.get("exonerate-gapped-keep-best").as<int>()>0 )
                lim = Settings_handle::st.get("exonerate-gapped-keep-best").as<int>();

            for(int i=0; i<lim && i<(int)hit_names.size(); i++)
            {
                string tid = names.find(best_hits.at(i).node)->second;
                tid_nodes->insert(pair<string,string>(tid,best_hits.at(i).node));
                hits->insert(pair<string,hit>(best_hits.at(i).node,best_hits.at(i)));

                if(Settings::noise>2)
                    cout<<"adding "<<best_hits.at(i).node<<" "<<best_hits.at(i).score<<endl;
            }
        }
    }
    else if(!Settings_handle::st.is("keep-despite-exonerate-fails"))
    {
        tid_nodes->clear();
        read->node_to_align = "discarded_read";
    }

    if(!Settings_handle::st.is("keep-exonerate-files"))
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
