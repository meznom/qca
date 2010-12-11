#include "qca.hpp"
#include "utilities.hpp"
#include <map>

/*
 * TODO
 * ----
 * 
 * * Dead Cell Systems
 * * grand canonical system
 * * thoroughly test dead plaquet and grand canonical / new hopping etc; compare with Mathematica
 * * P-P for three plaquets (bond system)
 * * P-P with dead cell (bond and qf)
 * * energies with dead cell? -- shouldn't make much difference
 * * grand canonical
 *    * number of electrons per cell for different a and b
 *    * number of electrons, considering Vext or dead cell
 *    * P-P
 * * compare energy, P-P for bond, qf, gc
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

template<class System>
void printEnergySpectrum (const System& s) 
{
    typedef std::map<double, int, EpsilonLess> myMap;
    typedef myMap::const_iterator mapIt;
    myMap evs;
    for (int i=0; i<s.H.eigenvalues.size(); i++)
        evs[s.H.eigenvalues(i)]++;
    for (mapIt i=evs.begin(); i!=evs.end(); i++)
        std::cout << i->first << "\t" << i->first - s.H.Emin << "\t" << i->second << std::endl;
}

template<class System>
void printPolarisationPolarisation (System& s, size_t i, size_t j, 
                                    double beta, const std::vector<double>& Vexts)
{
    for (size_t k=0; k<Vexts.size(); k++)
    {
        s.H.Vext = Vexts[k];
        s.H.construct();
        s.H.diagonalise();
        std::cout
            << s.ensembleAverage(beta, s.P(i)) << "\t"
            << s.ensembleAverage(beta, s.P(j)) << std::endl;
    }
}

CommandLineOptions setupCLOptions ()
{
    CommandLineOptions o;
    o.add("help", "h", "Print this help message.")
     .add("model", "m", "Which QCA model to use.")
     .add("N_p", "p", "Number of plaquets.")
     .add("t", "t", "Hopping parameter.")
     .add("td", "td", "Diagonal hopping parameter.")
     .add("ti", "ti", "Inter-plaquet hopping parameter.")
     .add("V0", "V0", "On-site coulomb repulsion (Hubbard U).")
     .add("mu", "mu", "Chemical potential.")
     .add("Vext", "Vext", "External potential.")
     .add("Pext", "Pext", "External, 'driving' polarisation of a 'dead plaquet' to the left of the linear chain system.")
     .add("a", "a", "Intra-plaquet spacing.")
     .add("b", "b", "Inter-plaquet spacing.")
     .add("beta", "beta", "Inverse temperature.")
     .add("energy-spectrum", "E", "Calculate the energy spectrum.")
     .add("polarisation", "P", "Calculate the polarisation for specified plaquet(s).")
     .add("particle-number", "N", "Calculate the particle number for the specified plaquet or the total particle number.");
    
    o["N_p"].setDefault(1);
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

    return o;
}

void printUsage (CommandLineOptions& o)
{
    std::cerr << "Usage: " << std::endl;
    std::cerr << o.optionsDescription();
}

enum OutputMode {Header, Line, None};

template<class System>
class Measurement
{
public:
    Measurement (Description o_)
    : o(o_), s(o_), needHeader(true), globalFirst(true)
    {
        for (Description::ItemsType::const_iterator i=o.items.begin(); i!=o.items.end(); i++)
            outputConfig[i->first] = Header;
        outputConfig["energy-spectrum"] = None;
        outputConfig["polarisation"] = None;
        outputConfig["particle-number"] = None;
        outputConfig["help"] = None;
    }

    void setParam (const std::string& param, const DescriptionItem& value)
    {
        o[param] = value;
        if (outputConfig[param] == Header)
            needHeader = true;
    }

    void measure (const std::string& param, const DescriptionItem& value)
    {
        setParam(param, value);
        measure();
    }

    void measure () 
    {
        s.setParameters(o);
        s.H.construct();
        s.H.diagonalise();

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
        for (int i=0; i<s.H.eigenvalues.size(); i++)
            evs[s.H.eigenvalues(i)]++;
        for (mapIt i=evs.begin(); i!=evs.end(); i++)
        {
            std::vector<double> v(3);
            v[0] = i->first;
            v[1] = i->first - s.H.Emin;
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
            P.push_back(s.ensembleAverage(o["beta"], s.P(ps[i])));
        }
    }

    void measureParticleNumber ()
    {
        //TODO o["particle-number"].isSize_t() or isBool() would be uselful
        size_t p;
        try
        {
            p = o["particle-number"];
            if (p >= s.N_p) 
            {
                std::cerr << std::endl << "Particle number: There is no plaquet " 
                    << p << " in this system." << std::endl;
                std::exit(EXIT_FAILURE);
            }
            N = s.ensembleAverage(o["beta"], s.N(p));
        }
        catch (ConversionException e)
        {
            N = s.ensembleAverage(o["beta"], s.N());
        }
    }

    void store ()
    {
        // always print energy spectrum as a separate block
        if (o["energy-spectrum"].isSet())
        {
            printHeader();
            std::cout << "# " << std::endl << "# E\tE-Emin\tDegeneracy" << std::endl;
            for (size_t i=0; i<E.size(); i++)
                std::cout << E[i][0] << "\t" << E[i][1] << "\t" << E[i][2] << std::endl;
            std::cout << std::endl;
        }

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
                if (!first) std::cout << "\t";
                else first = false;
                
                try
                {
                    size_t p;
                    p = o["particle-number"];
                    std::cout << "N" << N;
                }
                catch (ConversionException e)
                {
                    std::cout << "N";
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
        {
            if (!first) std::cout << "\t";
            else first = false;
            std::cout << N;
        }
        std::cout << std::endl;
    }

    void printHeader ()
    {
        if (!globalFirst) std::cout << std::endl;
        else globalFirst = false;
        for (Description::ItemsType::const_iterator i=o.items.begin(); i!=o.items.end(); i++)
            if (outputConfig[i->first] == Header)
                std::cout << "# " << i->first << " = " << i->second << std::endl;
    }

    void setOutputMode (std::string param, OutputMode mode)
    {
        outputConfig[param] = mode;
    }

private:
    Description o;
    System s;

    std::vector<std::vector<double> > E;
    std::vector<double> P;
    double N;
    
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
    if (pos+1 == params.size())
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
    const char* pnamesv[] = {"t","td","ti","V0","mu","Vext","Pext","a","b","beta"};
    size_t pnamesc = 10;
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

    Description cOpts(opts);
    cOpts["N_p"] = 1;
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
    Measurement<System> M(cOpts);
    runMeasurement(M, opts, params, 0);


    /*
    std::vector<double> Vexts = opts["Vext"];
    if (Vexts.size() > 1) opts["Vext"] = 0; 

    System qca(opts);

    for (size_t j=0; j<Vexts.size(); j++)
    {
        qca.H.Vext = Vexts[j];
    
        qca.H.construct();
        qca.H.diagonalise();

        if (opts["energy-spectrum"].isSet())
        {
            printEnergySpectrum(qca);
            std::cout << std::endl;
        }

        if (opts["polarisation"].isSet())
        {
            std::cout << qca.H.Vext;
            std::vector<size_t> ps = opts["polarisation"];
            for (size_t i=0; i<ps.size(); i++)
            {
                if (ps[i] >= qca.N_p) 
                {
                    std::cerr << std::endl << "Polarisation: There is no plaquet " 
                        << ps[i] << " in this system." << std::endl;
                    std::exit(EXIT_FAILURE);
                }
                std::cout << "\t" << qca.ensembleAverage(opts["beta"], qca.P(ps[i]));
            }
            std::cout << std::endl;
        }

        if (opts["particle-number"].isSet())
        {
            //TODO opts["particle-number"].isSize_t() or isBool() would be uselful
            size_t p;
            try
            {
                p = opts["particle-number"];
                if (p >= qca.N_p) 
                {
                    std::cerr << std::endl << "Particle number: There is no plaquet " 
                        << p << " in this system." << std::endl;
                    std::exit(EXIT_FAILURE);
                }
                std::cout << qca.ensembleAverage(opts["beta"], qca.N(p)) << std::endl;
            }
            catch (ConversionException e)
            {
                std::cout << qca.ensembleAverage(opts["beta"], qca.N()) << std::endl;
            }
        }
    }
    */
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

    if (opts["model"] == "bond")
        run<DQcaBondPlain>(opts);
    else if (opts["model"] == "quarterfilling" || opts["model"] == "quarter" || 
             opts["model"] == "qf")
        run<DQcaQuarterFillingPlain>(opts);
    else if (opts["model"] == "grandcanonical" || opts["model"] == "grand")
        run<DQcaGrandCanonicalPlain>(opts);
    else if (opts["model"] == "bondDP")
        run<DQcaBondDeadPlaquet>(opts);
    else if (opts["model"] == "quarterfillingDP" || opts["model"] == "quarterDP" || 
             opts["model"] == "qfDP")
        run<DQcaQuarterFillingDeadPlaquet>(opts);
    else if (opts["model"] == "grandcanonicalDP" || opts["model"] == "grandDP")
        run<DQcaGrandCanonicalDeadPlaquet>(opts);
    else
    {
        std::cerr 
            << "Please specify a model. Options are: " << std::endl 
            << "\t" << "'bond'" << std::endl
            << "\t" << "'bondDP'" << std::endl
            << "\t" << "'quarterfilling' or 'quarter'" << std::endl
            << "\t" << "'quarterfillingDP' or 'quarterDP'" << std::endl
            << "\t" << "'grandcanonical' or 'grand'" << std::endl
            << "\t" << "'grandcanonicalDP' or 'grandDP'" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}
