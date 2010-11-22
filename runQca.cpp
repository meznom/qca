#include "qca.hpp"

int main(int argc, const char *argv[])
{

    //QCABond qca(1);
    QCAQuarterFilling qca(1);
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
    //qca.selectBlock();
    //std::cerr << qca.sparseBlock(qca.H, 2,2,4,4) << std::endl;
    //std::cerr << qca.H << std::endl;
    //qca.selectBlock(qca.H, FilterNElectronsPerPlaquet(2));
    //std::cerr << qca.H.H << std::endl;
    qca.H.diagonalise();
    std::cerr << qca.H.eigenvalues.size() << std::endl;
    std::cerr << qca.ensembleAverage(1, qca.P(0)) << std::endl;
    
    return 0;
}
