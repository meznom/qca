#include "qca.hpp"
#include "utilities.hpp"
#include <map>
#include "version.hpp"
#include <ctime>
#include <iomanip>

class EpsilonLess
{
public:
    EpsilonLess (double e_ = 1E-20)
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
     .add("model", "m", "Which QCA model to use. Options are: bond, fixed, grand.")
     .add("natural", "n", "Use natural units. Effectively uses a Coulomb potential 1/r by setting epsilon_r = e^2/(4pi epsilon_0) (in eV).")
     .add("", "p", "Number of plaquets.")
     .add("", "t", "Hopping parameter.")
     .add("", "td", "Diagonal hopping parameter.")
     .add("", "V0", "On-site coulomb repulsion (Hubbard U).")
     .add("", "mu", "Chemical potential.")
     .add("", "Vext", "External potential.")
     .add("", "beta", "Inverse temperature.")
     .add("", "epc", "Number of electrons per cell, either 2 or 6.")
     .add("", "q", "Compensation charge. Defaults to 0.")
     .add("", "epsilonr", "Relative permattiviy.")
     .add("", "lambdaD", "Debye screening length. Set to zero to disable screening.")
     .add("", "layout", "Which cell layout to use: wire or nonuniformwire")
     .add("", "a", "Intra-cell spacing.")
     .add("", "b", "Inter-cell spacing. For the (uniform) wire.")
     .add("", "bs", "List of inter-cell spacings. Has to equal the number of cells. For the non-uniform wire.")
     .add("", "Pext", "Polarization of the driver cell (left-most cell in the wire).")
     .add("energy-spectrum", "E", "Calculate the energy spectrum.")
     .add("polarization", "P", "Calculate the polarization for the specified plaquet(s).")
     .add("polarization2", "P2", "Calculate the polarization for the specified plaquet(s). Uses a different definition for the polarization.")
     .add("particle-number", "N", "Calculate the particle number for the specified plaquet(s).");
    
    o["p"].setDefault(1);
    o["t"].setDefault(1);
    o["td"].setDefault(0);
    o["V0"].setDefault(0);
    o["mu"].setDefault(0);
    o["Vext"].setDefault(0);
    o["beta"].setDefault(1);
    o["q"].setDefault(0);
    o["epsilonr"].setDefault(1);
    o["lambdaD"].setDefault(0);
    o["layout"].setDefault("wire");
    o["a"].setDefault(1);
    o["b"].setDefault(1.75);
    std::vector<double> bs;
    bs.push_back(1.75);
    //o["bs"].setDefault(bs); TODO
    o["Pext"].setDefault(0);

    return o;
}

void printUsage (CommandLineOptions& o)
{
    std::cerr << "Usage: " << std::endl;
    std::cerr << o.usage();
}

enum OutputMode {Header, Line, None};

/**
 * Run and store measurements.
 *
 * This simple implementation does not really store results of measurements, but
 * prints them straight to the standard output.
 */
template<class System>
class Measurement
{
public:
    Measurement (OptionSection o_)
    : o(o_), s(o_), needHeader(true), globalFirst(true), lineCount(0)
    {
        for (OptionSection::OptionsType::iterator i=o.getOptions().begin(); i!=o.getOptions().end(); i++)
            outputConfig[i->getName()] = Header;
        outputConfig["energy-spectrum"] = None;
        outputConfig["polarization"] = None;
        outputConfig["polarization2"] = None;
        outputConfig["particle-number"] = None;
        outputConfig["help"] = None;
        outputConfig["version"] = None;
        outputConfig["natural"] = None;
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
        if (o["polarization"].isSet())
            measurePolarization();
        if (o["polarization2"].isSet())
            measurePolarization2();
        if (o["particle-number"].isSet())
            measureParticleNumber();

        store();
    }

    void setOutputMode (std::string param, OutputMode mode)
    {
        outputConfig[param] = mode;
    }

private:
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

    void measurePolarization ()
    {
        P.resize(0);
        std::vector<int> ps = o["polarization"];
        for (size_t i=0; i<ps.size(); i++)
        {
            if (ps[i] >= s.N_p()) 
            {
                std::cerr << std::endl << "Polarization: There is no plaquet " 
                    << ps[i] << " in this system." << std::endl;
                std::exit(EXIT_FAILURE);
            }
            P.push_back(s.measure(o["beta"], s.P(ps[i])));
        }
    }

    void measurePolarization2 ()
    {
        P2.resize(0);
        std::vector<int> ps = o["polarization2"];
        for (size_t i=0; i<ps.size(); i++)
        {
            if (ps[i] >= s.N_p()) 
            {
                std::cerr << std::endl << "Polarization2: There is no plaquet " 
                    << ps[i] << " in this system." << std::endl;
                std::exit(EXIT_FAILURE);
            }
            P2.push_back(s.measurePolarization2(o["beta"], ps[i]));
        }
    }

    void measureParticleNumber ()
    {
        N.resize(0);
        std::vector<int> ns = o["particle-number"];
        for (size_t i=0; i<ns.size(); i++)
        {
            if (ns[i] >= s.N_p()) 
            {
                std::cerr << std::endl << "Particle number: There is no plaquet " 
                    << ns[i] << " in this system." << std::endl;
                std::exit(EXIT_FAILURE);
            }
            N.push_back(s.measureParticleNumber(o["beta"], ns[i]));
        }
    }

    /**
     * "Store" measurements.
     *
     * This simple implementation just prints the results to stdout.
     */
    void store ()
    {
        // configure output formatting
        std::cout << std::setprecision(6);
        std::cout << std::left;

        // read parameters back from system, so we are sure we are printing the
        // parameters that are actually used
        updateParametersFromSystem();

        if (needHeader)
        {
            //begin a new block in output
            printHeader();
            lineCount = 0;
            needHeader = false;
        }

        // print table heading every twenty lines for easier readability
        if (lineCount%20 == 0)
            printTableHeading();
        
        // print all results (each parameter in a separate column, all in one
        // line)
        std::cout << "  ";
        for (typename OutputMap::const_iterator i=outputConfig.begin(); i!=outputConfig.end(); i++)
            if (i->second == Line)
            {
                std::cout << std::setw(width) << o[i->first];
            }
        if (o["polarization"].isSet())
            for (size_t j=0; j<P.size(); j++)
            {
                std::cout << std::setw(width) << P[j];
            }
        if (o["polarization2"].isSet())
            for (size_t j=0; j<P2.size(); j++)
            {
                std::cout << std::setw(width) << P2[j];
            }
        if (o["particle-number"].isSet())
            for (size_t j=0; j<N.size(); j++)
            {
                // by convention the last entry is the total number of
                // particles for each plaquet
                const size_t n = N[j].size();
                std::cout << std::setw(width) << N[j][n-1];
                for (size_t k=0; k<n-1; k++)
                    std::cout << std::setw(width) << N[j][k];
            }
        if (o["energy-spectrum"].isSet())
            for (size_t i=0; i<E.size(); i++)
                std::cout << std::setw(width) << E[i][0]
                          << std::setw(width) << E[i][1]
                          << std::setw(width) << E[i][2];
        std::cout << std::endl;
        lineCount++;
    }

    void updateParametersFromSystem ()
    {
        OptionSection oUsed = s.getParameters();
        for (OptionSection::OptionsType::iterator i=oUsed.getOptions().begin(); 
                                                  i!=oUsed.getOptions().end(); i++)
            o[i->getName()] = i->getValue();
    }

    /**
     * Prints the heading for the table of results.
     *
     * Here table means the columns of all measured observables (e.g.
     * polarization, particle number).
     */
    void printTableHeading ()
    {
        std::cout << "# ";
        for (typename OutputMap::const_iterator i=outputConfig.begin(); i!=outputConfig.end(); i++)
            if (i->second == Line)
            {
                std::cout << std::setw(width) << i->first;
            }
        if (o["polarization"].isSet())
        {
            std::vector<size_t> ps = o["polarization"];
            for (size_t i=0; i<ps.size(); i++)
            {
                std::stringstream ss;
                ss << "P_" << ps[i];
                std::string s(ss.str());
                std::cout << std::setw(width) << s;
            }
        }
        if (o["polarization2"].isSet())
        {
            std::vector<size_t> ps = o["polarization2"];
            for (size_t i=0; i<ps.size(); i++)
            {
                std::stringstream ss;
                ss << "P2_" << ps[i];
                std::string s(ss.str());
                std::cout << std::setw(width) << s;
            }
        }
        if (o["particle-number"].isSet())
        {
            std::vector<size_t> ns = o["particle-number"];
            for (size_t i=0; i<ns.size(); i++)
            {
                std::stringstream ss;
                ss << "N_" << ns[i];
                std::string s(ss.str());
                std::string s0 = s + "_0";
                std::string s1 = s + "_1";
                std::string s2 = s + "_2";
                std::string s3 = s + "_3";
                std::cout << std::setw(width) << s
                          << std::setw(width) << s0
                          << std::setw(width) << s1
                          << std::setw(width) << s2
                          << std::setw(width) << s3;
            }
        }
        if (o["energy-spectrum"].isSet())
        {
            for (size_t i=0; i<E.size(); i++)
            {
                std::stringstream ss;
                ss << "E_" << i;
                std::string s1(ss.str());
                ss.str("");
                ss << "E_" << i << "-E_0";
                std::string s2(ss.str());
                ss.str("");
                ss << "n_deg_" << i; //degeneracy
                std::string s3(ss.str());
                std::cout << std::setw(width) << s1
                          << std::setw(width) << s2
                          << std::setw(width) << s3;
            }
        }
        std::cout << std::endl;
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
        std::cout << "# " << std::endl;
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
    std::vector<double> P2;
    std::vector<std::vector<double> > N;
    
    bool needHeader, globalFirst;
    typedef std::map<std::string, OutputMode> OutputMap;
    OutputMap outputConfig;
    unsigned int lineCount;

    enum {width=14}; //column width for output
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
    if (opts["natural"].isSet())
        opts["epsilonr"] = QCA_NATURAL_EPSILON_R;

    //TODO: include bs here as well
    const char* pnamesv[] = {"t","td","V0","mu","Vext","Pext","a","b",
                             "beta","q","epsilonr","lambdaD"};
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
    cOpts["V0"] = 0;
    cOpts["mu"] = 0;
    cOpts["Vext"] = 0;
    cOpts["beta"] = 1;
    cOpts["q"] = 0;
    cOpts["epsilonr"] = 1;
    cOpts["lambdaD"] = 0;
    cOpts["a"] = 1;
    cOpts["b"] = 1.75;
    std::vector<double> bs;
    bs.push_back(1.75);
    //cOpts["bs"] = bs; TODO
    cOpts["Pext"] = 0;
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
            run<DQcaGeneric<QcaBond> >(opts);
        else if (opts["model"] == "fixedcharge" || opts["model"] == "fixed")
            run<DQcaGeneric<QcaFixedCharge> >(opts);
        else if (opts["model"] == "grandcanonical" || opts["model"] == "grand")
            run<DQcaGeneric<QcaGrandCanonical> >(opts);
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
