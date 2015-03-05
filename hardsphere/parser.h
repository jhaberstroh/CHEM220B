#ifndef __PARSER_H__
#define __PARSER_H__
#include <sstream>
#include <string.h>
#include <vector>
#include <sstream>
#include <iostream>

std::string default_line = "OPTIONS:";
std::vector<std::string> PARSER_HELP(1, default_line);
int HELP_PADDING=15;

bool log_in_help(std::string log)
{
    for (int i = 0 ; i < PARSER_HELP.size() ; i++)
    {
        if (log.compare(PARSER_HELP[i]) == 0)
        {
            return true;
        }
    }
    return false;
}

bool help(int argc, char * argv[], int i)
{
    if (!strcmp(argv[i], "-h"))
    {
        for (int help_line = 0 ; help_line < PARSER_HELP.size(); help_line++)
        {
            std::cout << PARSER_HELP[help_line] << std::endl;
        }
        exit(1);
    }
}


bool parse_int(int argc, char * argv[], int i, const char * argname, int * arg,
    const std::string& help="", bool bVerbose = false)
{
    std::ostringstream help_line;
    int numspace = HELP_PADDING - strlen(argname);
    help_line << argname;
    for (int sp = 0 ; sp < numspace ; sp++)
    {
        help_line << " ";
    }  
    help_line << ": " << help;
    std::string help_line_str = help_line.str();
    if (!log_in_help(help_line_str))
    {
        PARSER_HELP.push_back(help_line_str);
    }

    int output;
    if (!strcmp(argv[i], argname))
    {
        if (i + 1 >= argc)
        {
            std::cerr << "Invalid format for option " << argname << std::endl;
            exit(100); //My exit code for arg position failure 
        }
        std::istringstream iss(argv[i+1]);
        iss >> output;
        if (iss.fail())
        {
            std::cerr << "Non-integer argument for " << argname << std::endl;
            exit(101); //My exit code for arg malformed failure
        }
        *arg = output;
        if (bVerbose)
        {
            std::cout << "Found option " << argname << ": " << *arg <<std::endl;
        }
        return true;
    }
    return false;
}

bool parse_double(int argc, char * argv[], int i, const char * argname, double * arg,
    const std::string& help="", bool bVerbose=false)
{
    std::ostringstream help_line;
    int numspace = HELP_PADDING - strlen(argname);
    help_line << argname;
    for (int sp = 0 ; sp < numspace ; sp++)
    {
        help_line << " ";
    }  
    help_line << ": " << help;
    std::string help_line_str = help_line.str();
    if (!log_in_help(help_line_str))
    {
        PARSER_HELP.push_back(help_line_str);
    }

    double output;
    if (!strcmp(argv[i], argname))
    {
        if (i + 1 >= argc)
        {
            std::cerr << "Invalid format for option " << argname << std::endl;
            exit(100); //My exit code for arg position failure 
        }
        std::istringstream iss(argv[i+1]);
        iss >> output;
        if (iss.fail())
        {
            std::cerr << "Non-integer argument for " << argname << std::endl;
            exit(101); //My exit code for arg malformed failure
        }
        *arg = output;
        if (bVerbose)
        {
            std::cout << "Found option " << argname << ": " << *arg <<std::endl;
        }
        return true;
    }
    return false;
}
#endif //__PARSER_H__
