#include "qca.hpp"
#include "utilities.hpp"
#include <map>

/*
 * TODO
 * ----
 * 
 * * Dead Cell Systems
 * * functions to measure / output energy spectra (over various parameters)
 * * function to output polarisation (over Vext or for dead cell's P)
 * * polarisation-polarisation response function
 * * set system, parameters from command line or configuration file
 * * nice / meaningful output functions
 * * grand canonical system
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
        std::cout << i->first << "\t" << i->second << std::endl;
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
    /*
    o.add("help", "h", "Print this help message.")
     .add("model", "m", "Which QCA model to use.")
     .add("N_p", "p", "Number of plaquets.")
     .add("hopping", "t", "Hopping parameter.")
     .add("hopping-diagonal", "t_d", "Diagonal hopping parameter.")
     .add("on-site", "V_0", "On-site coulomb repulsion (Hubbard U).")
     .add("external", "V_ext", "External potential.")
     .add("intra-plaquet-spacing", "a", "Intra-plaquet spacing.")
     .add("inter-plaquet-spacing", "b", "Inter-plaquet spacing.");
     */
    o.add("help", "h", "Print this help message.")
     .add("model", "m", "Which QCA model to use.")
     .add("N_p", "p", "Number of plaquets.")
     .add("t", "t", "Hopping parameter.")
     .add("td", "td", "Diagonal hopping parameter.")
     .add("V0", "V0", "On-site coulomb repulsion (Hubbard U).")
     .add("Vext", "Vext", "External potential.")
     .add("a", "a", "Intra-plaquet spacing.")
     .add("b", "b", "Inter-plaquet spacing.")
     .add("beta", "beta", "Inverse temperature.")
     .add("test-list", "Test list parsing.")
     .add("energy-spectrum", "E", "Calculate the energy spectrum.")
     .add("polarisation", "P", "Calculate the polarisation for specified plaquet(s).");
    
    o["N_p"] = 1;
    o["beta"] = 1;

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
    //TODO: opts["Vext"] = 0 does not work
    if (Vexts.size() > 1) opts["Vext"] = "0"; 

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

        /*
         * TODO: -P 0 is a valid specification for the polarisation, but does not
         * work (yields false below) -> fix this in DescriptionItem
         */
        if (opts["polarisation"])
        {
            std::vector<size_t> ps = opts["polarisation"];
            for (size_t i=0; i<ps.size(); i++)
            {
                if (ps[i] >= qca.N_p) 
                {
                    std::cerr << std::endl << "Polarisation: There is no plaquet " 
                        << ps[i] << " in this system." << std::endl;
                    std::exit(EXIT_FAILURE);
                }
                if (i>0) std::cout << "\t";
                std::cout << qca.ensembleAverage(opts["beta"], qca.P(ps[i]));
            }
            std::cout << std::endl;
        }
    }
            
    /*
    std::cout << qca.ensembleAverage(1, qca.P(0)) << std::endl;
    if (qca.N_p > 1)
        std::cout << qca.ensembleAverage(1, qca.P(1)) << std::endl;
    
    qca.H.Vext = -1;
    qca.H.construct();
    qca.H.diagonalise();
    std::cerr << qca.ensembleAverage(1, qca.P(0)) << std::endl;
    if (qca.N_p > 1)
        std::cerr << qca.ensembleAverage(1, qca.P(1)) << std::endl;
    
    if (qca.N_p > 1)
    {
        double Vexts[] = {-1,0,1};
        std::vector<double> vVexts(Vexts, Vexts + sizeof(Vexts)/sizeof(double));
        printPolarisationPolarisation(qca, 0, 1, 1.0, vVexts);
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

    if (opts["help"])
    {
        printUsage(opts);
        std::exit(EXIT_SUCCESS);
    }

    //TODO --test-list '' does not work as expected
    if (opts["test-list"])
    {
        std::vector<double> tv = opts["test-list"];
        std::cerr << "Testing list parsing: " << std::endl;
        for (size_t i=0; i<tv.size(); i++)
            std::cerr << tv[i] << std::endl;
        std::exit(EXIT_SUCCESS);
    }

    if (opts["model"] == "bond")
        run<DQcaBond>(opts);
    else if (opts["model"] == "quarterfilling" || opts["model"] == "qf")
        run<DQcaQuarterFilling>(opts);
    else
    {
        std::cerr 
            << "Please specify a model. Options are: " << std::endl 
            << "\t" << "'bond'" << std::endl
            << "\t" << "'quarterfilling' or 'qf'" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}
