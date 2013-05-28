#include "cqca.hpp"
#include "version.hpp"

const char PROGRAM_NAME[] = "runQca";

void printVersion ()
{
    std::cout << PROGRAM_NAME << ", version " << PROGRAM_VERSION << std::endl;
}

void printUsage ()
{
    printVersion();
    std::cout 
        << "Usage: " << PROGRAM_NAME 
        << " -h | --help | " << std::endl
        << "              -v | --version | " << std::endl 
        << "              [configfile] |" << std::endl
        << "              [json-configuration]" << std::endl;
}

int main(int argc, const char** argv)
{
    std::string config;
    if (argc == 2 && (std::string(argv[1]) == "-h" || 
                      std::string(argv[1]) == "--help"))
    {
        printUsage();
        std::exit(EXIT_SUCCESS);
    }
    else if (argc == 2 && (std::string(argv[1]) == "-v" || 
                           std::string(argv[1]) == "--version"))
    {
        printVersion();
        std::exit(EXIT_SUCCESS);
    }
    else if (argc == 2 && argv[1][0] != '-')
        config = argv[1];
    else
    {
        printUsage();
        std::exit(EXIT_FAILURE);
    }
    
    Runner::run(config);
    return EXIT_SUCCESS;
}
