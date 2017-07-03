#ifndef MAFFT_ALIGNMENT_H
#define MAFFT_ALIGNMENT_H

#include <fstream>
#include <string>
#include <vector>
#include <sys/stat.h>

extern std::string tempdir;


class Mafft_alignment
{

    std::string remove_last_whitespaces(const std::string & s)
    {
        // Copy sequence
        std::string st (s);

        while (st.size() > 0 && this->is_whitespace_character(st[st.size() - 1]))
        {
            st.erase(st.end() - 1);
        }

        // Send result
        return st;
    }

    std::string remove_whitespaces(const std::string & s)
    {
        std::string st="";

        for (unsigned int i = 0; i < s.size(); i++)
        {
            if (!this->is_whitespace_character(s[i]))
            {
                st+=s[i];
            }
        }
        return st;
    }

    bool is_whitespace_character(char c)
    {
        return (c == ' ')
               || (c == '\t')
               || (c == '\n')
               || (c == '\r')
               || (c == '\f');
    }


public:
    Mafft_alignment();
    bool test_executable();
    void align_sequences(std::vector<std::string> *names,std::vector<std::string> *sequences);
};

#endif // MAFFT_ALIGNMENT_H
