#ifndef TOOL_PROBE_H
#define TOOL_PROBE_H

/*
 * Ask whether a helper tool is available by LOOKING FOR IT, not by running it.
 *
 * Every `test_executable()` in this tree answered that question the same way:
 * build `<dir>/<tool> <flag> </dev/null >/dev/null 2>/dev/null`, hand it to
 * system(3), and test the exit status -- then, when the binary is not beside
 * prank's own, retry with a bare name. Six tools are probed per invocation and
 * most probes run twice, so a prank process paid ~10 of those before it
 * aligned anything.
 *
 * That is expensive in a way the source does not show, because `mafft` is a
 * SHELL SCRIPT. Measured against an unpatched build of this same commit,
 * 20 real 3822 nt vs 3810-3816 nt pairs, `-once -codon -F -reproducible`,
 * glibc/Linux, local SSD:
 *
 *     execve  217 -> 137   (-80, -36.9%)
 *     clone   189 -> 111   (-78, -41.3%)
 *     wall    0.220 s -> 0.155 s median   (-29.7%)
 *
 * Of the 80 execs that go, 42 grep, 37 awk, 30 mktemp, 16 cat, 15 sh, 12 tr
 * and 11 file are mafft's wrapper unpacking itself -- to print a help text
 * nobody reads. The probe is issued twice per process (the banner in
 * progressivealignment.h, then again before the initial alignment), and each
 * issue runs mafft up to twice. prank is invoked once per sequence pair by
 * some callers, so that fraction is paid per pair, and on a network
 * filesystem each of those execs is a remote lookup as well.
 *
 * HOW IT LOOKS. Like the shell's own command search: the first REGULAR file
 * with execute permission, beside prank's binary first, then along $PATH (an
 * empty element is the current directory; an unset $PATH is the system
 * default path, as the shell uses). A directory of the same name does not
 * count. access(2) checks the real rather than the effective user, which only
 * differs for a set-uid prank.
 *
 * WHAT IS GIVEN UP, STATED PLAINLY. What is no longer established is "...and
 * it exits with the status this particular probe expects". A tool that is present and
 * executable but broken now reports AVAILABLE and fails at the real invocation,
 * where the failure is visible, instead of reporting absent and having prank
 * silently continue without it. That direction was chosen deliberately: a
 * silent fallback is the harder failure to diagnose, and this tree already
 * carries a fix for a probe whose expected status was simply wrong (FastTree's,
 * which read as absent on every machine but one named host).
 *
 * The answer is memoised per (dir, tool). The filesystem is not expected to
 * change under a running alignment, and the probes are re-issued from several
 * call sites within one process.
 */

#include <string>
#include <map>
#include <cstdlib>
#include <sys/stat.h>
#include <unistd.h>

namespace prank_tool_probe
{

inline bool is_executable_file(const std::string &path)
{
    struct stat st;
    return stat(path.c_str(), &st) == 0 && S_ISREG(st.st_mode)
           && access(path.c_str(), X_OK) == 0;
}

/*
 * Is <tool> executable in <dir>, or on $PATH?
 *
 * On success *found_dir is set the way the callers' own variables expect: the
 * directory (with its trailing '/') when the tool sits beside prank's binary,
 * and "" when it was found on $PATH and must therefore be invoked by bare name.
 * That is exactly what the system(3) probes assigned, so the callers' later
 * command building is unchanged.
 */
inline bool find_tool(const std::string &dir, const std::string &tool,
                      std::string *found_dir)
{
    static std::map<std::string, std::pair<bool, std::string> > cache;

    std::string key = dir;
    key += '\x1f';                       // not a path character
    key += tool;

    std::map<std::string, std::pair<bool, std::string> >::const_iterator it
        = cache.find(key);
    if (it != cache.end())
    {
        if (found_dir != 0)
            *found_dir = it->second.second;
        return it->second.first;
    }

    bool ok = false;
    std::string where = "";

    if (!dir.empty())
    {
        if (is_executable_file(dir + tool))
        {
            ok = true;
            where = dir;
        }
    }

    if (!ok)
    {
        std::string path;
        const char *path_env = getenv("PATH");
        if (path_env != 0)
            path = path_env;
        else
        {
            size_t n = confstr(_CS_PATH, NULL, 0);
            if (n > 0)
            {
                std::string def(n, '\0');
                confstr(_CS_PATH, &def[0], n);
                path = def.c_str();
            }
            else
                path = "/bin:/usr/bin";
        }
        {
            std::string::size_type start = 0;
            while (start <= path.size())
            {
                std::string::size_type end = path.find(':', start);
                if (end == std::string::npos)
                    end = path.size();

                std::string element = path.substr(start, end - start);
                // POSIX: an empty $PATH element means the current directory.
                if (element.empty())
                    element = ".";
                if (element[element.size() - 1] != '/')
                    element += '/';

                if (is_executable_file(element + tool))
                {
                    ok = true;
                    where = "";          // invoked by bare name, as before
                    break;
                }

                if (end == path.size())
                    break;
                start = end + 1;
            }
        }
    }

    cache[key] = std::make_pair(ok, where);
    if (found_dir != 0)
        *found_dir = where;
    return ok;
}

} // namespace prank_tool_probe

#endif // TOOL_PROBE_H
