#include "diagonaliser.hpp"

int main (int argc, char** argv)
{
    System s(4);
    s.dumpBasis();
    FermionicState<> s2(40);
    std::cerr << s.getIndexFromState(s.basis[7]) << std::endl;
    std::cerr << s.cs[2] << std::endl;


    AndersonModel a(2);
    a.H.getMatrix();
    a.H.diagonalise();
    std::cerr << a.H.eigenvalues << std::endl;


    std::exit(0);
}
