#ifndef __CQCA_HPP__
#define __CQCA_HPP__

#include <boost/property_tree/ptree.hpp>
#include "qca.hpp"

class ConfigurationException: public std::runtime_error
{
public:
    ConfigurationException(const std::string& msg)
        : std::runtime_error(msg) {}
};

namespace ptreeHelpers
{
    using boost::property_tree::ptree;

    template<typename T>
    std::vector<T> getArray (const ptree& c)
    {
        std::vector<T> v;
        for (ptree::const_iterator i=c.begin(); i!=c.end(); i++)
            v.push_back(i->second.get_value<T>());
        return v;
    }

    template <typename T>
    ptree constructArray (const std::vector<T>& v)
    {
        ptree c;
        for (int i=0; i<v.size(); i++)
        {
            ptree p;
            p.put_value(v[i]);
            c.push_back(ptree::value_type("", p));
        }
        return c;
    }

    // rtree = result tree
    // TODO: rtree really should have double as the value type, but there are a
    // couple of issues that need to be resolved before this will work/compile
    //typedef boost::property_tree::basic_ptree<std::string, double> rtree;
    typedef ptree rtree;
};

using boost::property_tree::ptree;
using namespace ptreeHelpers;

class CLayout : public Layout
{
private:
    typedef Layout Base;
    enum LayoutType {wire, nonuniformWire};

    LayoutType type;
    int cells;
    double a, b, Pext;
    std::vector<double> bs;
    ElectronsPerCell epc;

public:
    CLayout ()
    {}

    void setConfig (const ptree& c)
    {
        // read in parameters
        type = getType(c);
        cells = c.get("cells", 1);
        a = c.get("a", 1.0);
        b = c.get("b", 3.0);
        bs = getBs(c);
        Pext = c.get("Pext", 0.0);
        epc = getEpc(c);

        // validate parameter values
        if (cells<1)
            throw ConfigurationException("A minimum of 1 cell is required");
        if (Pext<-1 || Pext>1)
            throw ConfigurationException("Pext must be in the range [-1:1]");
        bool bsPos=true;
        for (int i=0; i<bs.size(); i++)
            if (bs[i]<0) bsPos=false;
        if (a<0 || b<0 || !bsPos)
            throw ConfigurationException("Wire dimensions (a, b, bs) must be positive");
        if (cells != bs.size())
            throw ConfigurationException("The number of b-values in 'bs' must be"
                        " equal to the number of cells for the non-uniform wire.");

        // construct the wire
        if (type == wire)
            Base::wire(cells, a, b, Pext, epc);
        else if (type == nonuniformWire)
            Base::nonuniformWire(cells, a, bs, Pext, epc);
    }

    ptree getConfig () const
    {
        ptree c;

        if (type == wire)
        {
            c.put("type", "wire");
            c.put("cells", cells);
            c.put("a", a);
            c.put("b", b);
            c.put("Pext", Pext);
            c.put("epc", epc);
        }
        else if (type == nonuniformWire)
        {
            c.put("type", "nonuniformwire");
            c.put("cells", cells);
            c.put("a", a);
            c.put_child("bs", ptreeArray(bs));
            c.put("Pext", Pext);
            c.put("epc", epc);
        }
        return c;
    }

private:
    LayoutType getType (const ptree& c) const
    {
        std::string typeString = c.get("type", "wire");
        if (typeString == "wire")
            return wire;
        else if(typeString == "nonuniformwire")
            return nonuniformWire;
        else
            throw ConfigurationException("Unknown layout type: '" + typeString + "'");
    }

    ElectronsPerCell getEpc (const ptree& c) const
    {
        int epcInt = c.get("epc", 2);
        if (epcInt == 2)
            return epc2;
        else if(epcInt == 6)
            return epc6;
        else
            throw ConfigurationException("Electrons per Cell (epc) must be either 2 or 6");
    }

    std::vector<double> getBs (const ptree& c) const
    {
        std::vector<double> myBs;
        ptree bsc = c.get_child("bs", ptree());
        for (ptree::const_iterator i=bsc.begin(); i!=bsc.end(); i++)
            myBs.push_back(i->second.get_value<double>());
        return myBs;
    }

    template <typename T>
    ptree ptreeArray (const std::vector<T>& v) const
    {
        ptree c;
        for (int i=0; i<v.size(); i++)
        {
            ptree p;
            p.put_value(v[i]);
            c.push_back(ptree::value_type("", p));
        }
        return c;
    }
};

template<class QcaSystem>
class CQca : public QcaSystem
{
private:
    typedef QcaSystem Base;
    typedef CQca<QcaSystem> Self;

    CLayout l;
    ptree os;
public:
    CQca ()
    {}

    void setConfig (const ptree& c)
    {
        //TODO: model, natural
        Base::t = c.get("t", 1.0);
        Base::td = c.get("td", 0.0); 
        Base::Vext = c.get("Vext", 0.0);
        Base::V0 = c.get("V0", 1000.0); 
        Base::mu = c.get("mu", 0.0);
        Base::epsilonr = c.get("epsilonr", 1.0);
        Base::lambdaD = c.get("lambdaD", 0.0);
        Base::q = c.get("q", 0);
        Base::beta = c.get("beta", 1);
        os = c.get_child("observables", ptree());
        l.setConfig(c.get_child("layout", ptree()));
        Base::l = l;
    }

    ptree getConfig () const
    {
        /*
         * Ideally we would read the layout configuration back from Base::l. However, the
         * the implementation of the class Layout currently does not know about its
         * high-level layout (e.g. is it a wire or a nonuniform wire). That's
         * why we store an instance of CLayout in CQca and use CLayout to get
         * the layout configuration.
         */
        ptree c;
        c.put("t", Base::t);
        c.put("td", Base::td);
        c.put("Vext", Base::Vext);
        c.put("V0", Base::V0);
        c.put("mu", Base::mu);
        c.put("epsilonr", Base::epsilonr);
        c.put("lambdaD", Base::lambdaD);
        c.put("q", Base::q);
        c.put("beta", Base::beta);
        c.put_child("layout", l.getConfig());
        c.put_child("observables", os);
        return c;
    }

    rtree measure ()
    {
        //TODO: only call update when necessary -- e.g. not necessary for
        //changed beta
        Base::update();

        rtree r;
        for (ptree::const_iterator i=os.begin(); i!=os.end(); i++)
        {
            if (i->first == "P")
                r.put_child("P", measureP(i->second));
            else if (i->first == "P2")
                r.put_child("P2", measureP2(i->second));
            else if (i->first == "N")
                r.put_child("N", measureN(i->second));
            else if (i->first == "E" && i->second.get_value<std::string>() == "yes")
                r.put_child("E", measureE(i->second));
            else
                throw ConfigurationException("Unknown observable: '" + i->first + "'");
        }
        return r;
    }

private:
    std::vector<int> getCells (const ptree& c) const
    {
        std::vector<int> allCells;
        for (int i=0; i<Base::N_p(); i++)
            allCells.push_back(i);

        std::vector<int> cells;
        if (c.get_value<std::string>() == "all")
            cells = allCells;
        else if (c.get_value<std::string>() != "")
            cells.push_back(c.get_value<int>());
        else
            cells = getArray<int>(c);

        for (int i=0; i<cells.size(); i++)
            if (cells[i] < 0)
                throw ConfigurationException("In observable specification: "
                                             "Cell numbers must be positive");
        return cells;
    }

    rtree measureP (const ptree& c) const
    {
        std::vector<int> cells = getCells(c);
        rtree r;
        for (int i=0; i<cells.size(); i++)
            r.put(toString(cells[i]), Base::measurePolarization(cells[i]));
        return r;
    }

    rtree measureP2 (const ptree& c) const
    {
        std::vector<int> cells = getCells(c);
        rtree r;
        for (int i=0; i<cells.size(); i++)
            r.put(toString(cells[i]), Base::measurePolarization2(cells[i]));
        return r;
    }

    rtree measureN (const ptree& c) const
    {
        std::vector<int> cells = getCells(c);
        rtree r;
        for (int i=0; i<cells.size(); i++)
        {
            std::vector<double> ns = Base::measureParticleNumber(cells[i]);
            r.put(toString(cells[i]) + ".total", ns[4]);
            r.put(toString(cells[i]) + ".0", ns[0]);
            r.put(toString(cells[i]) + ".1", ns[1]);
            r.put(toString(cells[i]) + ".2", ns[2]);
            r.put(toString(cells[i]) + ".3", ns[3]);
        }
        return r;
    }

    class EpsilonLess
    {
    public:
        EpsilonLess (double e_ = 1E-20)
        : e(e_)
        {}

        bool operator() (double a, double b)
        {
            return b-a > e;
        }
    private:
        double e;
    };

    rtree measureE (const ptree& c) 
    {
        rtree r;
        typedef std::map<double, int, EpsilonLess> myMap;
        typedef typename myMap::const_iterator mapIt;
        myMap evs;
        for (int i=0; i<Base::energies().size(); i++)
            evs[Base::energies()(i)]++;
        for (mapIt i=evs.begin(); i!=evs.end(); i++)
        {
            std::vector<double> v(3);
            v[0] = i->first;
            v[1] = i->first - Base::Emin();
            v[2] = i->second;
            r.push_back(ptree::value_type("", constructArray(v)));
        }
        return r;
    }

    template<typename T>
    std::string toString (T v) const
    {
        std::stringstream ss;
        ss << v;
        return ss.str();
    }
};

#endif // __CQCA_HPP__
