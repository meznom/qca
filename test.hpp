#ifndef __TEST_HPP__
#define __TEST_HPP__

#include <iostream>

template<class System>
class Polarisation
{
public:
    Polarisation (System s_)
    : s(s_)
    {}

    int operator() (int i) const
    {
        return s.c(i+5);
    }

private:
    const System& s;
};

class System
{
public:
    System ()
    : polarisation(*this) 
    {}

    int c (int i) const
    {
        return 4+i;
    }

    void testPolarisation() const
    {
        std::cerr << polarisation(3) << std::endl;
    }

private:
    Polarisation<System> polarisation;
};

class System2 : public Polarisation<System2>
{
public:
    System2 ()
    : Polarisation<System2>(*this)
    {}

    int c (int i) const
    {
        return 4+i;
    }

    void testPolarisation() const
    {
        std::cerr << operator()(3) << std::endl;
    }

};

#endif // __TEST_HPP__
