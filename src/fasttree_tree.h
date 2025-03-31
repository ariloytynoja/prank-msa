#ifndef FASTTREE_TREE_H
#define FASTTREE_TREE_H

#include <fstream>
#include <string>
#include <vector>
#include <sys/stat.h>
#include "config.h"

using namespace std;


class FastTree_tree
{
    string progpath;

public:
    FastTree_tree();
    bool test_executable();
    string infer_phylogeny(std::vector<string> *names, std::vector<string> *sequences, bool is_protein);
};

#endif // FASTTREE_TREE_H
