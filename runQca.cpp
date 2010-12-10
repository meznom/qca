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
    o["beta"].setDefault(1);
    o["energy-spectrum"].setDefault(false);
    o["particle-number"].setDefault(false);
    o["help"].setDefault(false);

    return o;
}

void printUsage (CommandLineOptions& o)
{
    std::cerr << "Usage: " << std::endl;
    std::cerr << o.optionsDescription();
}

template<class System>
void run (CommandLineOptions& opts)
{
    std::vector<double> Vexts = opts["Vext"];
    if (Vexts.size() > 1) opts["Vext"] = 0; 

    System qca(opts);

    for (size_t j=0; j<Vexts.size(); j++)
    {
        qca.H.Vext = Vexts[j];
    
        qca.H.construct();
        qca.H.diagonalise();

        if (opts["energy-spectrum"])
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

    if (opts["help"])
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
