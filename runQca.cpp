#include "qca.hpp"
#include "utilities.hpp"
#include <map>
#include "version.hpp"
#include <ctime>

/*
 * TODO
 * ----
 * 
 * Output for the following command does not look correct to me:
 * ./runQca -m bondDP -p 1 -P 0 -Pext 0,1,2,3 -a 1,2 -E 
 */

class EpsilonLess
{
public:
    EpsilonLess (double e_ = 1E-10)
    : e(e_)
    {}

    bool operator() (double a, double b)
    {
        return b-a > e;
    }
private:
    double e;
};

CommandLineOptions setupCLOptions ()
{
    CommandLineOptions o;
    o.add("help", "h", "Print this help message.")
     .add("version", "Print program version.")
     .add("model", "m", "Which QCA model to use. Options are: bond, bondDP, fixed2, fixed6, fixed2DP, fixed6DP, grand, grand2DP, grand6DP.")
     .add("", "p", "Number of plaquets.")
     .add("", "t", "Hopping parameter.")
     .add("", "td", "Diagonal hopping parameter.")
     .add("", "ti", "Inter-plaquet hopping parameter.")
     .add("", "V0", "On-site coulomb repulsion (Hubbard U).")
     .add("", "mu", "Chemical potential.")
     .add("", "Vext", "External potential.")
     .add("", "Pext", "External, 'driving' polarisation of a 'dead plaquet' to the left of the linear chain system.")
     .add("", "a", "Intra-plaquet spacing.")
     .add("", "b", "Inter-plaquet spacing.")
     .add("", "beta", "Inverse temperature.")
     .add("", "epsilonr", "Relative permattiviy.")
     .add("", "lambdaD", "Debye screening length. Set to zero to disable screening.")
     .add("energy-spectrum", "E", "Calculate the energy spectrum.")
     .add("polarisation", "P", "Calculate the polarisation for the specified plaquet(s).")
     .add("particle-number", "N", "Calculate the particle number for the specified plaquet(s).");
    
    o["p"].setDefault(1);
    o["t"].setDefault(1);
    o["td"].setDefault(0);
    o["ti"].setDefault(0);
    o["V0"].setDefault(0);
    o["mu"].setDefault(0);
    o["Vext"].setDefault(0);
    o["Pext"].setDefault(0);
    o["a"].setDefault(1);
    o["b"].setDefault(1.75);
    o["beta"].setDefault(1);
    o["epsilonr"].setDefault(QCA_EPSILON_R_DEFAULT_VALUE);
    o["lambdaD"].setDefault(0);

    return o;
}

void printUsage (CommandLineOptions& o)
{
    std::cerr << "Usage: " << std::endl;
    std::cerr << o.usage();
}

enum OutputMode {Header, Line, None};

template<class System>
class Measurement
{
public:
    Measurement (OptionSection o_)
    : o(o_), s(o_), needHeader(true), globalFirst(true)
    {
        for (OptionSection::OptionsType::iterator i=o.getOptions().begin(); i!=o.getOptions().end(); i++)
            outputConfig[i->getName()] = Header;
        outputConfig["energy-spectrum"] = None;
        outputConfig["polarisation"] = None;
        outputConfig["particle-number"] = None;
        outputConfig["help"] = None;
        outputConfig["version"] = None;
    }

    void setParam (const std::string& param, const OptionValue& value)
    {
        o[param] = value;
        if (outputConfig[param] == Header)
            needHeader = true;
    }

    void measure (const std::string& param, const OptionValue& value)
    {
        setParam(param, value);
        measure();
    }

    void measure () 
    {
        s.setParameters(o);
        s.update();

        if (o["energy-spectrum"].isSet())
            measureEnergies();
        if (o["polarisation"].isSet())
            measurePolarisation();
        if (o["particle-number"].isSet())
            measureParticleNumber();

        store();
    }

    void measureEnergies ()
    {
        E.resize(0);
        typedef std::map<double, int, EpsilonLess> myMap;
        typedef myMap::const_iterator mapIt;
        myMap evs;
        for (int i=0; i<s.energies().size(); i++)
            evs[s.energies()(i)]++;
        for (mapIt i=evs.begin(); i!=evs.end(); i++)
        {
            std::vector<double> v(3);
            v[0] = i->first;
            v[1] = i->first - s.Emin();
            v[2] = i->second;
            E.push_back(v);
        }
    }

    void measurePolarisation ()
    {
        P.resize(0);
        std::vector<size_t> ps = o["polarisation"];
        for (size_t i=0; i<ps.size(); i++)
        {
            if (ps[i] >= s.N_p) 
            {
                std::cerr << std::endl << "Polarisation: There is no plaquet " 
                    << ps[i] << " in this system." << std::endl;
                std::exit(EXIT_FAILURE);
            }
            P.push_back(s.measure(o["beta"], s.P(ps[i])));
        }
    }

    void measureParticleNumber ()
    {
        N.resize(0);
        std::vector<size_t> ns = o["particle-number"];
        for (size_t i=0; i<ns.size(); i++)
        {
            if (ns[i] >= s.N_p) 
            {
                std::cerr << std::endl << "Particle number: There is no plaquet " 
                    << ns[i] << " in this system." << std::endl;
                std::exit(EXIT_FAILURE);
            }
            N.push_back(s.measure(o["beta"], s.N(ns[i])));
        }
    }

    void store ()
    {
        // read parameters back from system, so we are sure we are printing the
        // parameters that are actually used
        OptionSection oUsed = s.getParameters();
        for (OptionSection::OptionsType::iterator i=oUsed.getOptions().begin(); 
                                                  i!=oUsed.getOptions().end(); i++)
            o[i->getName()] = i->getValue();

        // always print energy spectrum as a separate block
        if (o["energy-spectrum"].isSet())
        {
            printHeaderAll();
            std::cout << "# " << std::endl << "# E\tE-Emin\tDegeneracy" << std::endl;
            for (size_t i=0; i<E.size(); i++)
                std::cout << E[i][0] << "\t" << E[i][1] << "\t" << E[i][2] << std::endl;
            std::cout << std::endl;
        }

        if (o["polarisation"].isSet() || o["particle-number"].isSet())
        {
            if (needHeader)
            {
                printHeader();

                std::cout << "# " << std::endl << "# ";
                bool first = true;
                for (typename OutputMap::const_iterator i=outputConfig.begin(); i!=outputConfig.end(); i++)
                    if (i->second == Line)
                    {
                        if (!first) std::cout << "\t";
                        else first = false;
                        std::cout << i->first;
                    }
                if (o["polarisation"].isSet())
                {
                    std::vector<size_t> ps = o["polarisation"];
                    for (size_t i=0; i<ps.size(); i++)
                    {
                        if (!first) std::cout << "\t";
                        else first = false;
                        std::cout << "P" << ps[i];
                    }
                }
                if (o["particle-number"].isSet())
                {
                    std::vector<size_t> ns = o["particle-number"];
                    for (size_t i=0; i<ns.size(); i++)
                    {
                        if (!first) std::cout << "\t";
                        else first = false;
                        std::cout << "N" << ns[i];
                    }
                }
                std::cout << std::endl;
                
                needHeader = false;
            }
            
            bool first = true;
            for (typename OutputMap::const_iterator i=outputConfig.begin(); i!=outputConfig.end(); i++)
                if (i->second == Line)
                {
                    if (!first) std::cout << "\t";
                    else first = false;
                    std::cout << o[i->first];
                }
            if (o["polarisation"].isSet())
                for (size_t j=0; j<P.size(); j++)
                {
                    if (!first) std::cout << "\t";
                    else first = false;
                    std::cout << P[j];
                }
            if (o["particle-number"].isSet())
                for (size_t j=0; j<N.size(); j++)
                {
                    if (!first) std::cout << "\t";
                    else first = false;
                    std::cout << N[j];
                }
            std::cout << std::endl;
        }
    }

    void printHeader ()
    {
        if (!globalFirst) std::cout << std::endl;
        else globalFirst = false;
        std::cout << "# program version: " << GIT_PROGRAM_VERSION << std::endl 
                  << "# date: " << getDate() << std::endl
                  << "# " << std::endl;
        for (OptionSection::OptionsType::iterator i=o.getOptions().begin(); i!=o.getOptions().end(); i++)
            if (outputConfig[i->getName()] == Header)
                std::cout << "# " << i->getName() << " = " << i->getValue() << std::endl;
    }

    void printHeaderAll ()
    {
        if (!globalFirst) std::cout << std::endl;
        else globalFirst = false;
        std::cout << "# program version: " << GIT_PROGRAM_VERSION << std::endl
                  << "# date: " << getDate() << std::endl
                  << "# " << std::endl;
        for (OptionSection::OptionsType::iterator i=o.getOptions().begin(); i!=o.getOptions().end(); i++)
            if (outputConfig[i->getName()] != None)
                std::cout << "# " << i->getName() << " = " << i->getValue() << std::endl;
    }

    void setOutputMode (std::string param, OutputMode mode)
    {
        outputConfig[param] = mode;
    }

    std::string getDate () const
    {
        time_t t;
        tm localTime;
        char cstr[100];
        time(&t);
        tzset();
        localtime_r(&t, &localTime);
        strftime(cstr, 100, "%a %b %d %T %Y %z", &localTime);
        return cstr;
    }
        
private:
    OptionSection o;
    System s;

    std::vector<std::vector<double> > E;
    std::vector<double> P;
    std::vector<double> N;
    
    bool needHeader, globalFirst;
    typedef std::map<std::string, OutputMode> OutputMap;
    OutputMap outputConfig;
    std::vector<std::string> changed;
};

bool comparePair (const std::pair<std::string, size_t>& p1, 
                  const std::pair<std::string, size_t>& p2)
{
    return p1.second < p2.second;
}

template<class Measurement>
void runMeasurement (Measurement& M, CommandLineOptions& opts, const std::vector<std::pair<std::string, size_t> >& params, size_t pos)
{
    const std::string pName = params[pos].first;
    if (pos+1 == params.size() && params[pos].second > 1)
        M.setOutputMode(pName, Line);

    std::vector<double> v = opts[pName];
    for (size_t i=0; i<v.size(); i++)
    {
        M.setParam(pName, v[i]);
        if (pos+1<params.size())
            runMeasurement(M, opts, params, pos+1);
        else
            M.measure();
    }
}

template<class System>
void run (CommandLineOptions& opts)
{
    const char* pnamesv[] = {"t","td","ti","V0","mu","Vext","Pext","a","b",
                             "beta","epsilonr","lambdaD"};
    size_t pnamesc = 12;
    std::vector<std::pair<std::string, size_t> > params;
    for (size_t i=0; i<pnamesc; i++)
    {
        std::pair<std::string, size_t> p;
        p.first = pnamesv[i];
        std::vector<double> v = opts[pnamesv[i]];
        p.second = v.size();
        params.push_back(p);
    }

    std::sort(params.begin(), params.end(), comparePair);

    OptionSection cOpts(opts);
    cOpts["t"] = 1;
    cOpts["td"] = 0;
    cOpts["ti"] = 0;
    cOpts["V0"] = 0;
    cOpts["mu"] = 0;
    cOpts["Vext"] = 0;
    cOpts["Pext"] = 0;
    cOpts["a"] = 1;
    cOpts["b"] = 1.75;
    cOpts["beta"] = 1;
    cOpts["epsilonr"] = QCA_EPSILON_R_DEFAULT_VALUE;
    cOpts["lambdaD"] = 0;
    Measurement<System> M(cOpts);
    runMeasurement(M, opts, params, 0);
}

int main(int argc, const char** argv)
{
    CommandLineOptions opts = setupCLOptions();
    try {
        opts.parse(argc, argv);
    }
    catch (CommandLineOptionsException e)
    {
        std::cerr << e.what() << std::endl << std::endl;
        printUsage(opts);
        std::exit(EXIT_FAILURE);
    }

    if (opts["help"].isSet())
    {
        printUsage(opts);
        std::exit(EXIT_SUCCESS);
    }
    if (opts["version"].isSet())
    {
        std::cout << GIT_PROGRAM_VERSION << std::endl;
        std::exit(EXIT_SUCCESS);
    }

    try 
    {
        if (opts["model"] == "bond")
            run<DQcaBondPlain>(opts);
        else if (opts["model"] == "fixedcharge2" || opts["model"] == "fixed2")
            run<DQcaFixedCharge2Plain>(opts);
        else if (opts["model"] == "fixedcharge6" || opts["model"] == "fixed6")
            run<DQcaFixedCharge6Plain>(opts);
        else if (opts["model"] == "grandcanonical" || opts["model"] == "grand")
            run<DQcaGrandCanonicalPlain>(opts);
        else if (opts["model"] == "bondDP")
            run<DQcaBondDeadPlaquet>(opts);
        else if (opts["model"] == "fixedcharge2DP" || opts["model"] == "fixed2DP")
            run<DQcaFixedCharge2DeadPlaquet>(opts);
        else if (opts["model"] == "fixedcharge6DP" || opts["model"] == "fixed6DP")
            run<DQcaFixedCharge6DeadPlaquet>(opts);
        else if (opts["model"] == "grandcanonical2DP" || opts["model"] == "grand2DP")
            run<DQcaGrandCanonical2DeadPlaquet>(opts);
        else if (opts["model"] == "grandcanonical6DP" || opts["model"] == "grand6DP")
            run<DQcaGrandCanonical6DeadPlaquet>(opts);
        else
        {
            printUsage(opts);
            std::cerr << std::endl << "Please specify a model." << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    catch (ConversionException e)
    {
        std::cerr << "One or more of the specified command line arguments are "
                  << "of the wrong type." << std::endl;
        std::cerr << "ConversionException: " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}
