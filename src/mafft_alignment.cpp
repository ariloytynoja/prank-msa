#include "mafft_alignment.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include "config.h"
#include <algorithm>
#include <unistd.h>
#include <sys/wait.h>

#if defined (__APPLE__)
#include <mach-o/dyld.h>
#endif

using namespace std;

Mafft_alignment::Mafft_alignment()
{
}

bool Mafft_alignment::test_executable()
{
    #if defined (__CYGWIN__)
    char path[200] = "";
    int length = readlink("/proc/self/exe",path,200-1);

    string epath = string(path).substr(0,length);
    if (epath.find("/")!=std::string::npos)
        epath = epath.substr(0,epath.rfind("/")+1);
    mafftpath = epath;
    epath = epath+"sh.exe "+epath+"mafft -h >/dev/null 2>/dev/null";
    int status = system(epath.c_str());
    return WEXITSTATUS(status) == 1;

    # else

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

    mafftpath = epath;
    epath = epath+"mafft -h >/dev/null 2>/dev/null";
    int status = system(epath.c_str());

    if(WEXITSTATUS(status) == 1)
        return true;

    mafftpath = "";
    status = system("mafft -h >/dev/null 2>/dev/null");

    return WEXITSTATUS(status) == 1;

    #endif
}

void Mafft_alignment::align_sequences(vector<string> *names,vector<string> *sequences)
{
    ofstream m_output;

    stringstream m_name;

    int r = rand();
    while(true)
    {
        m_name <<tmp_dir<<"/m"<<r<<".fas";
        ifstream m_file(m_name.str().c_str());

        if(!m_file)
        {
            m_output.open( m_name.str().c_str(), (ios::out) );
            break;
        }
        r = rand();
    }

    vector<string>::iterator ni = names->begin();
    vector<string>::iterator si = sequences->begin();

    for(;ni!=names->end();ni++,si++)
    {
        m_output<<">"<<*ni<<endl<<*si<<endl;
    }
    m_output.close();

    // mafft writes its progress to stderr, so stderr cannot be folded into
    // the stdout this function parses -- but discarding it means a failure
    // arrives with no reason attached. Keep it in a file beside the input
    // instead: it is read only if the alignment comes back empty, and it is
    // then the stderr OF THE RUN THAT FAILED rather than of a second run.
    stringstream e_name;
    e_name <<tmp_dir<<"/m"<<r<<".err";

    stringstream command;
    if (PROTEIN) {
        command << mafftpath<<"mafft --amino "<<tmp_dir<<"/m"<<r<<".fas 2> "<<e_name.str();
    }else{
        command << mafftpath<<"mafft "<<tmp_dir<<"/m"<<r<<".fas 2> "<<e_name.str();
    }
    if(NOISE>0)
        cout<<"cmd: "<<command.str()<<endl;

//    #if defined (__CYGWIN__)
//    char path[200];
//    int length = readlink("/proc/self/exe",path,200-1);

//    string epath = string(path).substr(0,length);
//    epath.replace(epath.rfind("prank"),string("prank").size(),string(""));
//    command << epath<<"sh.exe "<<epath<<"mafft "<<tmp_dir<<"m"<<r<<".fas 2>/dev/null";
//    # else
//    command << "mafft "+tmp_dir+"m"<<r<<".fas  2>/dev/null";
//    #endif

    FILE *fpipe;
    if ( !(fpipe = (FILE*)popen(command.str().c_str(),"r")) )
    {
        perror("Problems with mafft pipe.\nExiting.\n");
        exit(1);
    }

    names->clear();
    sequences->clear();

    // read mafft output
    string name, sequence = "";  // Initialization
    char temp[256];

    while ( fgets( temp, sizeof temp, fpipe))
    {
       string line(temp);

        if (line[0] == '>')
        {
            line = this->remove_last_whitespaces(line);

            // If a name and a sequence were found
            if ((name != "") && (sequence != ""))
            {
                names->push_back(name);
                sequence = this->remove_whitespaces(sequence);
                transform( sequence.begin(), sequence.end(), sequence.begin(), (int(*)(int))toupper );
                sequences->push_back(sequence);
                name = "";
                sequence = "";
            }
            name = line;
            name.erase(name.begin());  // Character > deletion
        }
        else
        {
            sequence += temp;  // Sequence isolation
        }

    }

    // Addition of the last sequence in file
    if ((name != "") && (sequence != ""))
    {
        names->push_back(name);
        sequence = this->remove_whitespaces(sequence);
        transform( sequence.begin(), sequence.end(), sequence.begin(), (int(*)(int))toupper );
        sequences->push_back(sequence);
    }

    // A non-zero exit or a signal is a failure even if something was printed:
    // a partial alignment must not be used as if it were the whole one.
    int mafft_status = pclose(fpipe);
    bool mafft_ok = mafft_status != -1 && WIFEXITED(mafft_status)
                    && WEXITSTATUS(mafft_status) == 0;

    if(sequences->size()==0 || !mafft_ok)
    {

        cerr<<"\nError: Initial alignment with Mafft failed";
        if(mafft_status == -1)
            cerr<<" (its exit status could not be read)";
        else if(WIFSIGNALED(mafft_status))
            cerr<<" (killed by signal "<<WTERMSIG(mafft_status)<<")";
        else if(WEXITSTATUS(mafft_status) != 0)
            cerr<<" (exit status "<<WEXITSTATUS(mafft_status)<<")";
        else
            cerr<<" (it returned no sequences)";
        cerr<<". mafft's error output was:\n";

        // Report the stderr already captured above. Running mafft a SECOND
        // time to find out why the first failed costs a whole alignment and
        // can print a different failure than the one that happened.
        // The shell creates the file even when mafft writes nothing to it.
        ifstream e_file(e_name.str().c_str());
        if(e_file && e_file.peek() != EOF)
        {
            string line;
            while(getline(e_file,line))
                cerr<<line<<endl;
            e_file.close();
        }
        else
        {
            cerr<<"(mafft produced no diagnostic output)"<<endl;
        }

        cerr<<"\nNow exiting.\n";

        remove( m_name.str().c_str() );
        remove( e_name.str().c_str() );
        // exit() below skips main()'s own rmdir() of this run's directory.
        rmdir( tmp_dir );

        // Exit NON-ZERO. This is a failure, and prank is invoked per
        // sequence pair by callers that branch on the exit status; exiting 0
        // reported a misconfigured mafft as a successful run that happened
        // to produce nothing.
        exit(1);
    }

    //remove file
    if( remove( m_name.str().c_str() ) != 0)
      perror("Error deleting temporary file in Mafft_alignment::align_sequences");
    remove( e_name.str().c_str() );

}
