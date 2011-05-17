#ifndef __BASIS_HPP__
#define __BASIS_HPP__

#include <stdint.h>
#include <vector>
#include <map>
#include <bitset>
#include <algorithm>
#include <stdexcept>
#include <cassert>
#include <cmath>
#include <Eigen/Sparse>
#include <Eigen/Dense>
using namespace Eigen;


typedef SparseMatrix<double> SMatrix;
typedef MatrixXd DMatrix;
typedef VectorXd DVector;

enum Spin {UP=0, DOWN=1};

template<typename T=uint64_t>
class FermionicState
{
public:
    typedef T num_t; 
    typedef std::bitset<sizeof(T)*8> bitset;

    FermionicState (size_t N_) : N(N_)
    {
        assert(N<=sizeof(T)*8);
    }

    void operator= (const T& v)
    {
        s = v;
    }

    bool operator== (const FermionicState<T>& fs) const
    {
        return s == fs.s && N == fs.N;
    }

    bool operator[] (size_t pos) const
    {
        return s[pos];
    }

    typename bitset::reference operator[] (size_t pos)
    {
        //return s[pos];
        return typename bitset::reference(s, pos);
    }

    size_t size() const
    {
        return N;
    }

    size_t count(size_t fromPos, size_t toPos) const
    {
        if (toPos<fromPos) std::swap(toPos,fromPos);
        assert (fromPos <= toPos);
        size_t sum = 0;
        for (size_t i=fromPos; i<toPos; i++)
            sum += s[i];
        return sum;
    }

    size_t count(size_t toPos) const
    {
        size_t sum = 0;
        for (size_t i=0; i<toPos; i++)
            sum += s[i];
        return sum;
    }

    size_t count() const
    {
        return s.count();
    }

    bool operator< (const FermionicState<T>& state2) const
    {
        return s.to_ulong() < state2.s.to_ulong();
    }

private:
    size_t N;
    bitset s;
};

template<typename T>
std::ostream& operator<< (std::ostream& o, const FermionicState<T>& fs)
{
    for (size_t i=0; i<fs.size(); i++)
        o << fs[i];
    return o;
}

typedef FermionicState<uint16_t> State;

class SymmetryOperator
{
public:
    virtual ~SymmetryOperator() {}
    virtual double operator() (const State&) const = 0;
};

class ParticleNumberSymmetryOperator : public SymmetryOperator
{
public:
    virtual ~ParticleNumberSymmetryOperator() {}

    virtual double operator() (const State& s) const
    {
        return s.count();
    }
};

class SpinSymmetryOperator : public SymmetryOperator
{
public:
    virtual ~SpinSymmetryOperator() {}

    virtual double operator() (const State& s) const
    {
        // N_down = N - N_up with N the total particle number
        // spin = N_up - N_down = 2 N_up - N
        const int N = s.count();
        int N_up=0;
        for (size_t i=0; i<s.size(); i+=2)
            N_up += s[i];
        return 2*N_up - N;
    }
};

class BasisException : public std::logic_error
{
public:
    BasisException (const std::string& what)
    : std::logic_error(what)
    {}

    virtual ~BasisException() throw() {}
};

template<class Filter, class Sorter>
class Basis
{
public:
    typedef State::num_t num_t;

    typedef std::pair<State, std::vector<double> > SortableState;

    class SymmetrySorter
    {
    public:
        bool operator() (const SortableState& s1, const SortableState& s2) const
        {
            const std::vector<double>& ev1 = s1.second;
            const std::vector<double>& ev2 = s2.second;
            const size_t size = (ev1.size()<ev2.size())?ev1.size():ev2.size();
            for (size_t i=0; i<size; i++)
            {
                //TODO: is this dangerous with doubles? 
                // -- it should be okay, because doubles that come from 
                // equal integers should equal each other (5 == 5 => 5.0 == 5.0)
                if (ev1[i] < ev2[i]) 
                    return true;
                else if (ev1[i] > ev2[i])
                    return false;
            }
            return false;
        }
    };

    Basis (size_t N_orbital_, const Filter& filter_, const Sorter& sorter_)
    : N_orbital(N_orbital_), maskStart(-1), maskEnd(-1), filter(filter_), sorter(sorter_)
    {}

    void construct ()
    {
        /*
         * Note: sizeof(N_basis) must be bigger than sizeof(num_t), because a
         * FermionicState<num_t> can hold N = 2^(sizeof(num_t)) different
         * states. However this number N is not representable in num_t (num_t's
         * biggest number is N-1).
         */
        assert(sizeof(N_basis) > sizeof(num_t));
        N_basis = std::pow(2.0, static_cast<int>(N_orbital));
        State s(N_orbital);
        
        //TODO: disabled filter and sorter for now
        
        std::vector<SortableState> sStates(N_basis, SortableState(s, std::vector<double>(sos.size())));

        for (size_t num=0; num<N_basis; num++)
        {
            s = static_cast<num_t>(num);
            sStates[num] = SortableState(s, applySymmetryOperators(s));
            //if (filter(s))
            //    states.push_back(s);
        }

        std::sort(sStates.begin(), sStates.end(), SymmetrySorter());

        for (size_t i=0; i<sStates.size(); i++)
            states.push_back(sStates[i].first);
            

        //N_basis = states.size();
        //std::sort(states.begin(), states.end(), sorter);
        for (size_t i=0; i<states.size(); i++)
            indices[states[i]] = i;
    }

    Basis& addSymmetryOperator (const SymmetryOperator* so)
    {
        sos.push_back(so);
        return *this;
    }

    std::vector<double> applySymmetryOperators (const State& s) const
    {
        std::vector<double> evs(sos.size());
        for (size_t i=0; i<sos.size(); i++)
            evs[i] = sos[i]->operator()(s);
        return evs;
    }

    const State& operator() (num_t index) const
    {
        return states[index];
    }

    num_t operator() (const State& state) const
    {
        typedef std::map<State, num_t>::const_iterator CIT;
        CIT i = indices.find(state);
        if (i != indices.end())
            return i->second;
        throw BasisException("No such state in basis.");
    }

    size_t size() const
    {
        return states.size();
    }

    void dump () const
    {
        for (size_t i=0; i<states.size(); i++)
            std::cerr << states[i] << std::endl;
    }

    template<class MaskFilter>
    void mask (const MaskFilter& filter)
    {
        maskStart = -1;
        maskEnd = -1;
        for (size_t i=0; i<states.size(); i++)
        {
            if (maskStart == -1 && maskEnd == -1)
                if (filter(states[i])) maskStart = i;
                else continue;
            else if (maskStart != -1 && maskEnd == -1)
                if (!filter(states[i])) maskEnd = i;
                else continue;
            else if (filter(states[i]))
            {
                assert(maskStart != -1 && maskEnd != -1);
                throw BasisException("Non-contiguous block selected.");
            }
        }
        if (maskStart == -1) throw BasisException("Zero-size block selected.");
        if (maskEnd == -1) maskEnd = states.size();
    }

    void unmask () 
    {
        maskStart = -1;
        maskEnd = -1;
    }

    void applyMask (SMatrix& m) const
    {
        if (maskStart == -1 || maskEnd == -1)
            return;
        m = sparseBlock(m, maskStart, maskStart, maskEnd-maskStart, maskEnd-maskStart);
    }

    SMatrix applyMask (const SMatrix& m) const
    {
        if (maskStart == -1 || maskEnd == -1)
            return m;
        return sparseBlock(m, maskStart, maskStart, maskEnd-maskStart, maskEnd-maskStart);
    }

    /**
     * Block method for sparse matrices. 
     * 
     * Eigen doesn't implement block() for sparse matrices, only for dense
     * matrices. Hence this small helper function.
     */
    SMatrix sparseBlock (const SMatrix& om, int i, int j, int p, int q) const
    {
        SMatrix nm(p,q);
        for (int l=0; l<q; l++)
        {
            nm.startVec(l);
            for (SMatrix::InnerIterator it(om,l+j); it; ++it)
                if (it.row()>=i && it.row()<i+p)
                    nm.insertBack(it.row()-i,it.col()-j) = it.value();
        }
        nm.finalize();
        return nm;
    }

private:
    size_t N_orbital, N_basis;
    std::vector<State> states;
    std::map<State, num_t> indices;
    int maskStart, maskEnd;
    std::vector<const SymmetryOperator*> sos;
    const Filter& filter;
    const Sorter& sorter;
};

#endif // __BASIS_HPP__
