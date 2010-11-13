#include "diagonaliser.hpp"

int main (int argc, char** argv)
{
    //System s(4);
    //s.dumpBasis();
    //FermionicState<> s2(40);
    //std::cerr << s.getIndexFromState(s.basis[7]) << std::endl;
    //std::cerr << s.cs[2] << std::endl;


    AndersonModel am(3);
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


    std::exit(0);
}
