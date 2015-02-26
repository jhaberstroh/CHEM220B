#include <sstream>

bool parse_int(int argc, char * argv[], int i, const char * argname, int * arg,
    bool bVerbose = false)
{
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
    bool bVerbose=false)
{
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
