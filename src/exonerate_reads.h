#ifndef EXONERATE_READS_H
#define EXONERATE_READS_H

#include <fstream>
#include <string>
#include <vector>
#include <map>

struct hit {
    std::string query;
    std::string node;
    int score;
    int q_start;
    int q_end;
    char q_strand;
    int t_start;
    int t_end;
    char t_strand;
};

class Exonerate_reads
{
    static bool better (hit i,hit j) { return (i.score>j.score); }

    bool split_sugar_string(const std::string& row,hit *h);
    bool split_vulgar_string(const std::string& row,hit *h);
    void write_exonerate_input(Node *root, Fasta_entry *read, map<string,string> *names, int r);
    void delete_files(int r);

public:
    Exonerate_reads();
    bool test_executable();

    void local_alignment(Node *root, Fasta_entry *read, std::multimap<std::string,std::string> *good_hits, std::map<std::string,hit> *hits, bool is_local);

};

#endif // EXONERATE_READS_H
