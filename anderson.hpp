#ifndef __ANDERSON_HPP__
#define __ANDERSON_HPP__

#include "system.hpp"

template<class System>
class DoubleOccupancy
{
public:
    DoubleOccupancy (const System& s_)
    : s(s_)
    {}

    SMatrix operator() (size_t i) const
    {
        //return s.c(i,UP)*s.a(i,UP) * s.c(i,DOWN)*s.a(i,DOWN);
        return s.n(i,UP) * s.n(i,DOWN);
    }
private:
    const System& s;
};

template<class System>
class ParticleNumber
{
public:
    ParticleNumber (const System& s_)
    : s(s_)
    {}

    SMatrix operator() () const
    {
        SMatrix m(s.basis.size(), s.basis.size());
        for (size_t i=0; i<s.N_sites; i++)
            m += s.n(i);
        return m;
    }
private:
    const System& s;
};

template<class System>
class AndersonHamiltonian : public Hamiltonian<System>
{
public:
    AndersonHamiltonian (const System& s_)
    : Hamiltonian<System>(s_), H(Hamiltonian<System>::H), s(Hamiltonian<System>::s)
    {}

    void construct() 
    {
        H.setZero();
        for (size_t i = 1; i<s.N_sites; i++)
        {
            H += - mu * ( s.c(i,UP) * s.a(i,UP) + s.c(i,DOWN) * s.a(i,DOWN) );
            
            // this is physically correct, I think
            //size_t j = (i+1==N_sites)?1:i+1;
            //if (i==j) continue;
            //H += -t * ( 
            //    s.c(i,UP) * s.a(j,UP) + s.c(j,UP) * s.a(i,UP) +
            //    s.c(i,DOWN) * s.a(j,DOWN) + s.c(j,DOWN) * s.a(i,DOWN)
            //);

            // this is what my python diagonaliser does (and probably my ctqmc
            // code as well)
            for (size_t j=i-1; j<=i+1; j+=2)
            {
                size_t k = j;
                if (k==0) k = s.N_sites-1;
                else if (k==s.N_sites) k = 1;
                H += -t * ( 
                    s.c(i,UP) * s.a(k,UP) + s.c(i,DOWN) * s.a(k,DOWN)
                );
            }
        }
        H += Ed * (s.c(0,UP)*s.a(0,UP) + s.c(0,DOWN)*s.a(0,DOWN));
        H += V * (
            s.c(0,UP)*s.a(1,UP) + s.c(1,UP)*s.a(0,UP) + 
            s.c(0,DOWN)*s.a(1,DOWN) + s.c(1,DOWN)*s.a(0,DOWN)
        );
        //SMatrix I(DMatrix::Identity(basis.size(),basis.size()));
        //TODO
        SMatrix Id(s.basis.size(), s.basis.size());
        Id.reserve(s.basis.size());
        for (size_t i=0; i<s.basis.size(); i++)
        {
            Id.startVec(i);
            Id.insertBack(i,i) = 1;
        }
        Id.finalize();
        H += U * (s.c(0,UP)*s.a(0,UP) - 0.5*Id) * (s.c(0,DOWN)*s.a(0,DOWN) - 0.5*Id);
    }

    double t, mu, U, V, Ed;
    //TODO: it is too annoying having to reintroduce H and s!!!
    SMatrix& H;
    const System& s;
};

class AndersonModel : public BasicSystem<Filter::SelectAll, Sorter::DontSort>
{
public:
    AndersonModel (size_t N_bath_)
    : BasicSystem<Filter::SelectAll, Sorter::DontSort>(2 + 2*N_bath_, Filter::SelectAll(), Sorter::DontSort()),
      N_bath(N_bath_), N_sites(N_bath+1), H(*this), ensembleAverage(*this), 
      doubleOccupancy(*this), particleNumber(*this)
    {}

    size_t I (size_t i, Spin s) const
    {
        return 2*i + s;
    }

    const SMatrix& c (size_t i, Spin s) const
    {
        return creator(I(i,s));
    }

    const SMatrix& a (size_t i, Spin s) const
    {
        return annihilator(I(i,s));
    }

    SMatrix n (size_t i, Spin s) const
    {
        return c(i,s)*a(i,s);
    }

    SMatrix n (size_t i) const
    {
        return n(i,UP) + n(i,DOWN);
    }

    size_t N_bath, N_sites;
    AndersonHamiltonian<AndersonModel> H;
    EnsembleAverage<AndersonModel> ensembleAverage;
    DoubleOccupancy<AndersonModel> doubleOccupancy;
    ParticleNumber<AndersonModel> particleNumber;
};

#endif // __ANDERSON_HPP__
