/*
 * Copyright (c) 2011-2014 Burkhard Ritter
 * This code is distributed under the two-clause BSD License.
 */
#ifndef __BASIS_HPP__
#define __BASIS_HPP__

#include <stdint.h>
#include <vector>
#include <map>
#include <bitset>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <cassert>
#include <cmath>

class FermionicStateException : public std::logic_error
{
public:
    FermionicStateException (const std::string& what)
    : std::logic_error(what)
    {}
};

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

    FermionicState (const std::string& sString)
    {
        N = sString.size();
        for (size_t i=0; i<N; i++)
            if (sString[i] == '1')
                s[i] = 1;
            else if (sString[i] == '0')
                s[i] = 0;
            else
                throw FermionicStateException("Can't construct FermionicState from string '" + 
                                              sString + "'.");
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
        return s[pos];
        //return typename bitset::reference(s, pos);
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

    std::string toString() const
    {
        std::stringstream ss;
        for (size_t i=0; i<N; i++)
            ss << s[i];
        return ss.str();
    }

private:
    size_t N;
    bitset s;
};

template<typename T>
std::ostream& operator<< (std::ostream& o, const FermionicState<T>& fs)
{
    o << fs.toString();
    return o;
}

#ifndef STORAGE_TYPE_OF_FERMIONIC_STATE
#define STORAGE_TYPE_OF_FERMIONIC_STATE uint32_t
#endif

typedef FermionicState<STORAGE_TYPE_OF_FERMIONIC_STATE> State;

class SymmetryOperator
{
public:
    virtual ~SymmetryOperator() {}
    virtual int operator() (const State&) const = 0;
};

class ParticleNumberSymmetryOperator : public SymmetryOperator
{
public:
    virtual ~ParticleNumberSymmetryOperator() {}

    virtual int operator() (const State& s) const
    {
        return s.count();
    }
};

class SpinSymmetryOperator : public SymmetryOperator
{
public:
    virtual ~SpinSymmetryOperator() {}

    virtual int operator() (const State& s) const
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

struct Range
{
    Range () : a(0), b(0) {}
    Range (int a_, int b_) : a(a_), b(b_) {}

    bool operator== (const Range& r)
    {
        return a==r.a && b==r.b;
    }

    int a;
    int b;
};

typedef std::vector<int> Sector;

/*
 * Some small helper functions for convenient construction of Sector.
 */
Sector constructSector (int a)
{
    Sector s;
    s.push_back(a);
    return s;
}

Sector constructSector (int a, int b)
{
    Sector s;
    s.push_back(a);
    s.push_back(b);
    return s;
}

Sector constructSector (int a, int b, int c)
{
    Sector s;
    s.push_back(a);
    s.push_back(b);
    s.push_back(c);
    return s;
}

Sector constructSector (int a, int b, int c, int d)
{
    Sector s;
    s.push_back(a);
    s.push_back(b);
    s.push_back(c);
    s.push_back(d);
    return s;
}

Sector constructSector (int a, int b, int c, int d, int e)
{
    Sector s;
    s.push_back(a);
    s.push_back(b);
    s.push_back(c);
    s.push_back(d);
    s.push_back(e);
    return s;
}

Sector constructSector (int a, int b, int c, int d, int e, int f)
{
    Sector s;
    s.push_back(a);
    s.push_back(b);
    s.push_back(c);
    s.push_back(d);
    s.push_back(e);
    s.push_back(f);
    return s;
}

class BasisException : public std::logic_error
{
public:
    BasisException (const std::string& what)
    : std::logic_error(what)
    {}
};

class Basis
{
private:
    struct SectorAndState 
    {
        SectorAndState (Sector sector_, State state_)
            : sector(sector_), state(state_) {}
        
        Sector sector;
        State state;
    };

    class SymmetrySorter
    {
    public:
        bool operator() (const SectorAndState& sas1, const SectorAndState& sas2) const
        {
            const Sector& s1 = sas1.sector;
            const Sector& s2 = sas2.sector;
            assert(s1.size() == s2.size());
            for (size_t i=0; i<s1.size(); i++)
            {
                if (s1[i] < s2[i]) 
                    return true;
                else if (s1[i] > s2[i])
                    return false;
            }
            //They are equal. Has to be false, otherwise sort keeps on
            //reordering indefinitely.
            return false;
        }
    };

public:
    typedef State::num_t num_t;

    Basis () 
    : N_orbital(0), filterSet(false)
    {}

    void construct (size_t N_orbital_)
    {
        N_orbital = N_orbital_;

        /*
         * N_basis should be of type size_t. If the size of the basis cannot be
         * represented by size_t then we potentially have a problem anyway,
         * because we are using std::vector which uses size_t. 
         * The assertion checks that the number of states we are going to use
         * is actually representable by N_basis (and thus size_t).
         */
        size_t N_basis;
        assert(std::pow(2.0, static_cast<int>(N_orbital)) <= 
               std::pow(2.0, static_cast<int>(sizeof(N_basis)*8)) - 1);
        N_basis = std::pow(2.0, static_cast<int>(N_orbital));
        State state(N_orbital);
        
        std::vector<SectorAndState> sas;
        for (size_t num=0; num<N_basis; num++)
        {
            state = static_cast<num_t>(num);
            Sector sector = computeSectorForState(state);
            if (!filterSet || applyFilter(sector))
                sas.push_back(SectorAndState(sector, state));
        }
        std::sort(sas.begin(), sas.end(), SymmetrySorter());
        constructSectorsAndRanges(sas);

        states.resize(sas.size(), State(N_orbital));
        indices.clear();
        for (size_t i=0; i<sas.size(); i++)
        {
            states[i] = sas[i].state;
            indices[states[i]] = i;
        }
    }

    Basis& addSymmetryOperator (const SymmetryOperator* so)
    {
        symmetryOperators.push_back(so);
        return *this;
    }

    void setFilter (Sector filter_)
    {
        if (filter_.size() > symmetryOperators.size())
            throw BasisException("Invalid filter specified.");
        filter = filter_;
        if (filter.size() == 0) filterSet = false;
        else filterSet = true;
    }

    Range getRangeOfSector (const Sector& s) const
    {
        assert(sectors.size() == ranges.size());
        for (size_t i=0; i<sectors.size(); i++)
            if (sectors[i] == s)
                return ranges[i];
        return Range(0,0);
    }

    Range getRangeOfSector (int a, int b=0, int c=0, int d=0, int e=0) const
    {
        size_t size = symmetryOperators.size();
        Sector s;
        if (size > 0) s.push_back(a);
        if (size > 1) s.push_back(b);
        if (size > 2) s.push_back(c);
        if (size > 3) s.push_back(d);
        if (size > 4) s.push_back(e);
        return getRangeOfSector(s);
    }

    const std::vector<Range>& getRanges () const
    {
        return ranges;
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

    size_t numberOfOrbitals () const
    {
        return N_orbital;
    }

    void dump () const
    {
        for (size_t i=0; i<states.size(); i++)
            std::cerr << states[i] << std::endl;
    }

private:
    Sector computeSectorForState (const State& s) const
    {
        Sector sector(symmetryOperators.size());
        for (size_t i=0; i<symmetryOperators.size(); i++)
            sector[i] = symmetryOperators[i]->operator()(s);
        return sector;
    }

    void constructSectorsAndRanges(const std::vector<SectorAndState>& sas)
    {
        if (sas.size() == 0)
            return;
        Sector const * currentSector = &(sas[0].sector);
        sectors.clear();
        ranges.clear();
        int a = 0;
        int b = 0;
        for (size_t i=1; i<sas.size(); i++)
        {
            Sector const * newSector = &(sas[i].sector);
            if (*newSector != *currentSector)
            {
                b = i;
                sectors.push_back(*currentSector);
                ranges.push_back(Range(a,b));
                a = i;
                currentSector = newSector;
            }
        }
        sectors.push_back(*currentSector);
        ranges.push_back(Range(b, sas.size()));
    }

    bool applyFilter (const Sector& sector) const
    {
        assert(filter.size() <= sector.size());
        for (size_t i=0; i<filter.size(); i++)
            if (filter[i] != sector[i])
                return false;
        return true;
    }

    size_t N_orbital;
    std::vector<State> states;
    std::map<State, num_t> indices;
    std::vector<Sector> sectors;
    std::vector<Range> ranges;
    std::vector<const SymmetryOperator*> symmetryOperators;
    Sector filter;
    bool filterSet;
};

#endif // __BASIS_HPP__
