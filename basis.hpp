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

    Basis (size_t N_orbital_, const Filter& filter, const Sorter& sorter)
    : N_orbital(N_orbital_)
    {
        /*
         * Note: N_basis must be bigger than num_t, because a
         * FermionicState<num_t> can hold N = 2^(sizeof(num_t)) different
         * states. However this number N is not representable in num_t (num_t's
         * biggest number is N-1).
         */
        assert(sizeof(N_basis) > sizeof(num_t));
        N_basis = std::pow(2.0, static_cast<int>(N_orbital));
        states = std::vector<State>(N_basis, State(N_orbital));
        State s(N_orbital);
        size_t i=0;
        for (size_t num=0; num<N_basis; num++)
        {
            s = static_cast<num_t>(num);
            if (filter(s))
            {
                states[i] = static_cast<num_t>(num);
                i++;
            }
        }
        std::sort(states.begin(), states.end(), sorter);
        for (size_t i=0; i<states.size(); i++)
            indices[states[i]] = i;
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

private:
    size_t N_orbital, N_basis;
    std::vector<State> states;
    std::map<State, num_t> indices;
};

#endif // __BASIS_HPP__
