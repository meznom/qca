//#include "diagonaliser.hpp"

//#include "system.hpp"
#include "anderson.hpp"
#include <cstdlib>

int main (int argc, char** argv)
{
    //System s(4);
    //s.dumpBasis();
    //FermionicState<> s2(40);
    //std::cerr << s.getIndexFromState(s.basis[7]) << std::endl;
    //std::cerr << s.cs[2] << std::endl;


    /*
     * Test Anderson model
     */
    /*
    AndersonModel am(1);
    am.t = 1;
    am.V = 1;
    am.mu = 0;
    am.Ed = 0;
    am.U = 2;
    am.constructH();
    am.diagonalise();
    //std::cerr << am.eigenvalues << std::endl;

    //DoubleOccupancy<AndersonModel> d(am);
    //std::cerr << am.eigenvectors.col(2).dot( d(0) * am.eigenvectors.col(2) ) << std::endl;

    std::cerr << am.ensembleAverage(1, am.doubleOccupancy(0)) << std::endl;
    std::cerr << am.ensembleAverage(1, am.particleNumber()) << std::endl;
    //std::cerr << am.H << std::endl;
    */

    /*
     * QCA bond
     */
    /*
    QCABond qca(1);
    qca.param.t = 0.2;
    qca.param.td= 0.2 * qca.param.t;
    qca.param.Vext = 1;
    qca.param.V0 = 10;
    qca.param.a = 1;
    qca.param.b = 3 * qca.param.a;
    qca.constructH();
    //qca.selectBlock();
    //std::cerr << qca.sparseBlock(qca.H, 2,2,4,4) << std::endl;
    //std::cerr << qca.H << std::endl;
    qca.selectBlock(qca.H, FilterNElectronsPerPlaquet(2));
    std::cerr << qca.H << std::endl;
    qca.diagonalise();
    std::cerr << qca.eigenvalues.size() << std::endl;
    //this is correct! compared with Mathematica

    QCABond::SMatrix P = qca.polarisation(0);
    qca.selectBlock(P, FilterNElectronsPerPlaquet(2));
    //std::cerr << P << std::endl;
    std::cerr << qca.ensembleAverage(1.0, P) << std::endl;
    */

    /*
     * QCA quarter filling
     */
    /*
    QCAQuarterFilling qca(2);
    std::cerr << "QCA constructed." << std::endl;
    qca.param.t = 1.0;
    qca.param.td= 0.2 * qca.param.t;
    qca.param.Vext = 10;
    qca.param.a = 1.0 / 50.0;
    qca.param.b = 1.75 * qca.param.a;
    qca.param.V0 = 10.0 * 1.0 / qca.param.a;
    qca.constructH();
    std::cerr << "H constructed." << std::endl;
    //FilterAnd<FilterNElectronsPerPlaquet, FilterSpin> f(FilterNElectronsPerPlaquet(2,8), FilterSpin(0));
    FilterNElectronsPerPlaquet f(2,8);
    qca.selectBlock(qca.H, f);
    std::cerr << qca.H.cols() << std::endl;
    qca.diagonalise();
    std::cerr << "H diagonalised." << std::endl;
    std::cerr << qca.eigenvalues << std::endl;
    QCAQuarterFilling::SMatrix P1 = qca.polarisation(0);
    //QCAQuarterFilling::SMatrix P2 = qca.polarisation(1);
    qca.selectBlock(P1, f);
    //qca.selectBlock(P2, f);
    std::cerr << qca.ensembleAverage(1000, P1) << std::endl;
    //std::cerr << qca.ensembleAverage(1000, P2) << std::endl;
    */

    //System<Filter::SelectAll, Sorter::DontSort> s(4, Filter::SelectAll(), Sorter::DontSort());
    //s.basis.dump();
    //std::cerr << s.creator(0) << std::endl;
    //s.testPolarisation();

    //System2 s2;
    //s2.testPolarisation();

    // TODO: I get a segmentation fault with this:
    //AndersonModel am(7);
    // whereas the following works
    AndersonModel am(3);
    //std::cerr << am.c(0,UP).cols() << std::endl;
    am.H.t = 1;
    am.H.V = 1;
    am.H.mu = 0;
    am.H.Ed = 0;
    am.H.U = 2;
    am.H.construct();
    am.H.diagonalise();
    std::cerr << am.H.eigenvalues << std::endl;

    //DoubleOccupancy<AndersonModel> d(am);
    //std::cerr << am.eigenvectors.col(2).dot( d(0) * am.eigenvectors.col(2) ) << std::endl;

    std::cerr << am.ensembleAverage(1, am.doubleOccupancy(0)) << std::endl;
    std::cerr << am.ensembleAverage(1, am.particleNumber()) << std::endl;


    std::exit(0);
}
