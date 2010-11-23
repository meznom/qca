#include "qca.hpp"
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

int main(int argc, const char *argv[])
{

    //QCABond qca(2);
    QCAQuarterFilling qca(2);
    //Filter::And<Filter::NElectronsPerPlaquet, Filter::Spin> f(Filter::NElectronsPerPlaquet(2,8), Filter::Spin(0));
    //Filter::NElectronsPerPlaquet f(2,8);
    //Filter::NElectronsPerPlaquet f(2,4);
    //qca.basis.mask(f);
    
    qca.H.t = 0.2;
    qca.H.td= 0.2 * qca.H.t;
    qca.H.Vext = 1;
    qca.H.V0 = 10;
    qca.H.a = 1;
    qca.H.b = 3 * qca.H.a;

    qca.H.construct();
    qca.H.diagonalise();
    //printEnergySpectrum(qca);
    std::cout << qca.ensembleAverage(1, qca.P(0)) << std::endl;
    std::cout << qca.ensembleAverage(1, qca.P(1)) << std::endl;
    
    qca.H.Vext = -1;
    qca.H.construct();
    qca.H.diagonalise();
    std::cerr << qca.ensembleAverage(1, qca.P(0)) << std::endl;
    std::cerr << qca.ensembleAverage(1, qca.P(1)) << std::endl;
    
    double Vexts[] = {-1,0,1};
    std::vector<double> vVexts(Vexts, Vexts + sizeof(Vexts)/sizeof(double));
    printPolarisationPolarisation(qca, 0, 1, 1.0, vVexts);
    
    return 0;
}
